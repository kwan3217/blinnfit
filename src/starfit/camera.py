"""
Functions to calculate the camera matrix

Created: 1/31/25
"""
from dataclasses import dataclass
from enum import Enum

import numpy as np
from sqlite3 import Connection

from kwanmath.geodesy import llr2xyz
from kwanmath.interp import linterp
from kwanmath.vector import vlength, vnormalize, vcross


class Source(Enum):
    # Sources are, in order of increasing confidence:
    # unknown source
    UNKNOWN=0
    # Interpolated from higher-confidence sources
    INTERPOLATED=1
    # Fit manually
    MANUAL=2
    # Fit via a least-squares model
    FIT=3
    # 4 - Modeled. Certain parameters like right_denom and probably clock are actually constant
    #     over the video. Also, the intent is that there is eventually a spline model for each of
    #     the viewpoint variables. "Modeled" means either constant or splined.
    MODELED=4


@dataclass
class Camera:
    lat:float         # Camera axis latitude in B1950ECLIP, deg [0,360)
    lon:float         # camera axis longitude, deg [-90,90]
    angle:float       # Horizontal field of view in degrees
    clock:float       # Camera axis barrel twist, with 0deg such that up points as close to north pole as possible,
                      # positive twists camera to right, negative to left, deg [-180,180]
    right_denom:float # Denominator of right vector length. Right vector length is 4/right_denom, so right_denom=3 gives
                      # a 4:3 aspect ratio.
    width:int         # Width of image in pixels
    height:int        # Height of image in pixels
    lat_sig:float=None
    lon_sig:float=None
    angle_sig:float=None
    clock_sig:float=None
    right_denom_sig:float=None
    lat_source:Source=None
    lon_source:Source=None
    angle_source:Source=None
    clock_source:Source=None
    right_denom_source:Source=None
    def __post_init__(self):
        """
        Calculate internal parameters from the values collected from __init__()

        self.dir - direction unit vector that the camera is looking in
        self.sky - sky unit vector for camera
        self.right - actual aspect ratio
        self.C - Camera matrix that transforms from world to camera coordinates
        """
        self.dir = llr2xyz(lat=self.lat, lon=self.lon,deg=True)
        self.sky = self._make_sky()
        self.right=4/self.right_denom
        self._cmatrix()
    def _cmatrix(self)->np.ndarray:
        """

        :param loc: Location of camera, equivalent to camera{location...}
        :param look: Look-at point of camera, equivalent to camera{look_at...}
        :param sky: Sky vector of camera, equivalent to camera{sky...}
        :return: Camera matrix which transforms a global vector to a vector in camera space
        We will return a 4x4 matrix which will transform a vector in homogeneous coordinates into
        the camera frame. This matrix has the form characteristic of homogeneous transformations:
         [a b c xt]
         [d e f yt]
         [g h i zt]
         [0 0 0  1]
        Multiply this matrix by a 4-element column vector, with the last element being 1
        if the point is a finite distance from the camera (and thus affected by translation)
        zero if the point is infinitely far away, like a star, and thus unaffected by translation.
        The camera frame is set up such that the look direction is +z, right is +x, and up is +y.

        The coefficients a through i could in theory describe a matrix with arbitrary scaling,
        shearing, mirroring, etc. In practice this code will only ever return a matrix which describes
        a pure rotation.
        :raises: ValueError if location and look are the same vector, or direction is parallel or anti-parallel
                 to sky.
        """
        loc = np.zeros((3, 1))
        look=self.dir
        #Calculate the relative look direction
        look_rel = look - loc
        if vlength(look_rel)==0.0:
            raise ValueError("Camera look_at same as location")
        look_rel=vnormalize(look_rel)

        #Calculate the right vector as the cross product of the relative look and sky vector
        right=vcross(look_rel,self.sky)
        if vlength(right)==0:
            raise ValueError("Camera looking at sky")
        right=vnormalize(right)

        down=vcross(look_rel,right) #guaranteed to be unit-length since product of two perpendicular unit-length vectors

        # Build the pure rotation part of the matrix. This one transforms
        # from camera to world, so the first column is the direction of
        # the first camera basis vector (right) in world space, the second
        # is the second camera basis vector (down) in world space, and the
        # third is the third camera basis vector (look_rel) in world space.
        # r_bw=[right.x down.x look_rel.x 0]
        #      [right.y down.y look_rel.y 0]
        #      [right.z down.z look_rel.z 0]
        #      [   0       0       0      1]
        r_wc=np.identity(4)
        r_wc[0:3,0,None]=right
        r_wc[0:3,1,None]=down
        r_wc[0:3,2,None]=look_rel
        # Build the translation part of the matrix. This one translate (IE shears in 4D space)
        # and can be thought of as where the fourth (w) basis vector=[0,0,0,1]^T in camera space lands
        # in world space.
        # t_bw=[1 0 0 loc.x]
        #      [0 1 0 loc.y]
        #      [0 0 1 loc.z]
        #      [0 0 0   1  ]
        t_wc=np.identity(4)
        t_wc[0:3,3,None]=loc
        self.M_wc=t_wc@r_wc
        self.M_cw=np.linalg.inv(self.M_wc)
    def project(self,vs_w:np.ndarray,*,out_nan:bool=True,w:float=0):
        """
        Project the target into the camera field of view
        :param v: Column vector (or stack of column vectors of shape (3,N)) in world coordinates
                       to project into camera coordinates.
        :return: 2D position on camera, in the form of a shape (2,N) numpy array. Row 0 is
                 horizontal coordinate, row 1 is vertical
        """
        if vs_w.shape[0]==3:
            # Add homogeneous coordinate if needed
            vs_w=np.vstack((vs_w,np.zeros(vs_w.shape[1])+w))
        vs_c = self.M_cw @ vs_w
        #Convert to normalized screen coordinates. In this frame, the screen is on a plane perpendicular and
        #out along the z axis The edges of the screen are at +-0.5*up and +-0.5*right. Angle determines the distance
        #between the camera and the plane of the screen. Using the image at http://www.povray.org/documentation/view/3.7.0/246/
        #as a reference, tan(angle/2)=0.5*right/direction. We can solve this for direction:
        # tan(angle/2)*direction=0.5*right
        # direction=0.5*right/tan(angle/2)
        direction=0.5*self.right/np.tan(np.radians(self.angle)/2)
        #if the z component is negative, we don't want to plot. Do this by setting the z component to NaN if it was negative.
        vs_c[2,vs_c[2,:]<0]=float('NaN')
        #If the target is at the screen, then the x and y coordinates are already what we want. If it is twice as far,
        #then we need to divide x and y by 2. If half as far, then they need to multiply by two. In general, multiply
        # by direction/z. If we do this right, the z coordinate will become equal to direction, which indicates the other
        # components are normalized screen coordinates
        target_scl=vs_c[0:2,:]*direction/vs_c[2,:]
        result=np.zeros(target_scl.shape)
        cx=self.width/2
        cy=self.height/2
        up=1
        rx=linterp(-0.5*self.right,-self.width /2,0.5*self.right,self.width /2,target_scl[0,:])
        ry=linterp(-0.5*up        ,-self.height/2,0.5*up        ,self.height/2,target_scl[1,:])
        result[0,:]=rx+cx
        result[1,:]=ry+cy
        if out_nan:
            result[:,result[0,:]<0]=float('NaN')
            result[:,result[1,:]<0]=float('NaN')
            result[:,result[0,:]>self.width]=float('NaN')
            result[:,result[1,:]>self.height]=float('NaN')
        return result,np.isfinite(result[0,:])
    def _make_sky(self):
        ssky=vnormalize(vcross(self.dir,np.array([[0.0],[0.0],[1.0]])))
        csky=vnormalize(vcross(ssky,self.dir))
        sky=np.cos(np.deg2rad(self.clock))*csky+np.sin(np.deg2rad(self.clock))*ssky
        return sky
    def to_params(self)->np.ndarray:
        """
        Convert the Camera's parameters to a numpy array.

        :return: Numpy array with [lat, lon, clock, angle, right_denom]
        """
        return np.array([self.lat, self.lon, self.angle, self.clock, self.right_denom])
    @classmethod
    def from_params(cls, params: np.ndarray,width:int,height:int)->'Camera':
        """
        Factory method to create a Camera instance from a numpy array of 5 elements.

        :param params: Array containing [lat, lon, angle, clock, right_denom]
        :type params: np.ndarray
        :return: A new Camera instance
        """
        return cls(lat=params[0], lon=params[1], angle=params[2], clock=params[3], right_denom=params[4],width=width,height=height)
    @classmethod
    def from_frame_db(cls, conn:Connection, i_frame:int, width:int,height:int)->'Camera':
        """
        Factory method to create a Camera instance from a numpy array of 5 elements.

        :param params: Array containing [lat, lon, clock, angle, right_denom]
        :type params: np.ndarray
        :return: A new Camera instance
        :rtype: Camera
        :raises ValueError: If the array does not have exactly 5 elements
        """
        # Order of first five fields must match the canonical param order
        fields=("lat_c","lon_c","angle","clock","right_denom",
                "lat_c_sig","lon_c_sig","angle_sig","clock_sig","right_denom_sig",
                "lat_c_source","lon_c_source","angle_source","clock_source","right_denom_source")
        members=("lat","lon","angle","clock","right_denom",
                "lat_sig","lon_sig","angle_sig","clock_sig","right_denom_sig",
                "lat_source","lon_source","angle_source","clock_source","right_denom_source")
        sql=f"select {','.join(fields)} from frames where framenum=?;"
        params=conn.execute(sql,(i_frame,)).fetchone()
        result=cls.from_params(params[:5],width=width,height=height)
        # Get the optional parameters
        result.__data__.update({k:v for k,v in zip(members,params) if v is not None})
    def write_db(self,conn:Connection,i_frame:int):
        fields={"lat_c":self.lat,
                "lon_c":self.lon,
                "angle":self.angle,
                "clock":self.clock,
                "right_denom":self.right_denom}
        optional_fields = {
            "lat_c_sig": self.lat_sig,
            "lon_c_sig": self.lon_sig,
            "angle_sig": self.angle_sig,
            "clock_sig": self.clock_sig,
            "right_denom_sig": self.right_denom_sig,
            "lat_c_source": int(self.lat_source),
            "lon_c_source": int(self.lon_source),
            "angle_source": int(self.angle_source),
            "clock_source": int(self.clock_source),
            "right_denom_source": int(self.right_denom_source)
        }
        for k,v in optional_fields.items():
            if v is not None:
                fields[k]=v
        values=tuple([v for k,v in fields.items()]+[i_frame])
        set_clause=",".join([f"{k}=?" for k,v in fields.items()])
        sql=f"update frames set {set_clause} where framenum=?"
        with conn:
            conn.execute(sql,values)



