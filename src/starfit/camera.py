"""
Functions to calculate the camera matrix

Created: 1/31/25
"""
import numpy as np
from kwanmath.interp import linterp
from kwanmath.vector import vlength, vnormalize, vcross


def cmatrix(loc=None,look=None,sky=None):
    """

    :param loc: Location of camera, equivalent to camera{location...}
    :param look: Look-at point of camera, equivalent to camera{look_at...}
    :param sky: Sky vector of camera, equivalent to camera{sky...}
    :return: Camera matrix which transforms a global vector to a vector in camera space
    We will return a 4x4 matrix which will transform a vector in homogeneous coordinates into
    the camera frame. This matrix is:
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
    """

    #Calculate the relative look direction
    look_rel = look - loc
    assert vlength(look_rel)>0,"Camera look_at same as location"
    look_rel=vnormalize(look_rel)

    #Calculate the right vector as the cross product of the relative look and sky vector
    right=vcross(look_rel,sky)
    assert vlength(right)>0,"Camera looking at sky"
    right=vnormalize(right)

    down=vcross(look_rel,right) #guaranteed to be unit-length since product of two perpendicular unit-length vectors

    r=np.zeros((4,4))
    r[0:3,0,None]=right
    r[0:3,1,None]=down
    r[0:3,2,None]=look_rel
    #result[0:3,3]=-loc
    r[3,3]=1
    t=np.zeros((4,4))
    t[0,0]=1
    t[1,1]=1
    t[2,2]=1
    t[3,3]=1
    t[0:3,3,None]=loc
    result=t@r
    result=np.linalg.inv(result)
    return result


def project(right=None,angle=None,width=None,height=None,target_c=None,out_nan:bool=True):
    """
    Project the target into the camera field of view
    :param right: Length of Right vector, equivalent to camera{right -x*...}.
                  This controls the aspect ratio. Note that a right-handed
                  camera should use a positive value for right. The up vector is
                  implicitly camera{up y*1 ...} .
    :param angle: Field-of-view angle in degrees, equivalent to camera{angle...}
                  Length of direction vector is calculated from this and right
    :param width: Width of image in pixels, equivalent to image_width
    :param height: Height of image in pixels, equivalent to image_height
    :param target: 3D position of point to project in camera coordinates.
                   If your point is in world coordinates, transform it first
                   with cmatrix(...)@target
    :return: 2D position on camera, in the form of a 2xN numpy array. Row 0 is
             horizontal coordinate, row 1 is vertical
    """
    #Convert to normalized screen coordinates. In this frame, the screen is on a plane perpendicular and
    #out along the z axis The edges of the screen are at +-0.5*up and +-0.5*right. Angle determines the distance
    #between the camera and the plane of the screen. Using the image at http://www.povray.org/documentation/view/3.7.0/246/
    #as a reference, tan(angle/2)=0.5*right/direction. We can solve this for direction:
    # tan(angle/2)*direction=0.5*right
    # direction=0.5*right/tan(angle/2)
    direction=0.5*right/np.tan(np.radians(angle)/2)
    #if the z component is negative, we don't want to plot. Do this by setting the z component to NaN if it was negative.
    target_c[2,target_c[2,...]<0]=float('NaN')
    #If the target is at the screen, then the x and y coordinates are already what we want. If it is twice as far,
    #then we need to divide x and y by 2. If half as far, then they need to multiply by two. In general, multiply
    # by direction/z. If we do this right, the z coordinate will become equal to direction, which indicates the other
    # components are normalized screen coordinates
    target_scl=target_c[0:2,...]*direction/target_c[2,...]
    result=np.zeros(target_scl.shape)
    cx=width/2
    cy=height/2
    up=1
    rx=linterp(-0.5*right,-width /2,0.5*right,width /2,target_scl[0,...])
    ry=linterp(-0.5*up,   -height/2,0.5*up   ,height/2,target_scl[1,...])
    result[0,...]=rx+cx
    result[1,...]=ry+cy
    if out_nan:
        result[:,result[0,...]<0]=float('NaN')
        result[:,result[1,...]<0]=float('NaN')
        result[:,result[0,...]>width]=float('NaN')
        result[:,result[1,...]>height]=float('NaN')
    return result


def make_sky(clock,dir):
    ssky=vnormalize(vcross(dir,np.array([[0.0],[0.0],[1.0]])))
    csky=vnormalize(vcross(ssky,dir))
    sky=np.cos(np.deg2rad(clock))*csky+np.sin(np.deg2rad(clock))*ssky
    return sky
