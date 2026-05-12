"""
Describe purpose of this script here

Created: 2/9/25
"""
import sqlite3
from contextlib import closing
from dataclasses import dataclass

import numpy as np
from kwanmath.geodesy import xyz2llr
from kwanmath.interp import linterp
from kwanmath.matrix import point_toward, m_to_aa
from kwanmath.vector import vlength, vnormalize, vangle
from matplotlib import pyplot as plt
from spiceypy import spkezr, furnsh, scs2e, etcal

from starfit.camera import make_sky, Camera


names={
    701:"Ariel",
    702:"Umbriel",
    703:"Titania",
    704:"Oberon",
    705:"Miranda"
}

@dataclass
class Zoom:
    """
    f0 - Last frame at normal angle, therefore first frame of zoom is f0+1
    f1 - first frame at maximum zoom, therefore last frame before maximum zoom is f1-1
    f2 - last frame at maximum zoom, therefore first frame of zoom out is f2+1
    f3 - first frame at normal angle, therefore last frame of motion is f3-1
    """
    f0:int
    f1:int
    f2:int
    f3:int
    flaststars:int=None
    fease01:int=None
    zoomangle1:float=None


zoom_VOBEST=Zoom(f0=1293,f1=1382,f2=1533,f3=1622,zoomangle1=2.30407199274772)
zoom_VUBEST=Zoom(f0=2240,f1=2344,f2=2480,f3=2570,zoomangle1=2.3)
zoom_VMBEST=Zoom(f0=4190,f1=4250,f2=4580,f3=4680,zoomangle1=2.30653115847091)


def patch_vmbest_actual(*,conn:sqlite3.Connection):
    """
    Make sure there is an actual time for each FOV in the VMBEST sequence
    :param conn:
    :return:
    """
    patch=[(4310,scs2e(-32,"3/26846:07:768")), # C2684608, clear filter
           (4340,scs2e(-32,"3/26846:10:768")), # C2684611, clear filter
           (4370,scs2e(-32,"3/26846:13:768")), # C2684614, clear filter
           (4400,scs2e(-32,"3/26846:16:768")), # C2684617, clear filter
           (4430,scs2e(-32,"3/26846:19:768")), # C2684620, clear filter
           (4460,scs2e(-32,"3/26846:22:768")), # C2684623, clear filter
           (4490,scs2e(-32,"3/26846:25:768")), # C2684626, clear filter
           (4520,scs2e(-32,"3/26846:28:768"))] # C2684629, clear filter
    with conn, closing(conn.cursor()) as cur:
        sql = "insert or replace into frames (framenum,et,et_source,angle,angle_source,width,height,right_num) values (?,?,?,?,?,?,?,?)"
        for framenum,et in patch:
            cur.execute(sql,(framenum,et,2,zoom_VMBEST.zoomangle1,4,1280,720,16.0))
        frames=[row[0] for row in patch]
        ets=[row[1] for row in patch]
        et1,et2=linterp(frames[0],ets[0],frames[-1],ets[-1],np.array((zoom_VMBEST.f1,zoom_VMBEST.f2)))
        cur.execute(sql, (zoom_VMBEST.f1, et1, 2, zoom_VMBEST.zoomangle1, 4,1280,720,16.0))
        cur.execute(sql, (zoom_VMBEST.f2, et2, 2, zoom_VMBEST.zoomangle1, 4,1280,720,16.0))


def interp_ellipse(*, conn, ells):
    sparse_framenums=np.array([ell[0] for ell in ells])
    sparse_set=set(sparse_framenums)
    dense_framenums=np.arange(sparse_framenums[0],sparse_framenums[-1]+1,dtype=np.uint16)
    spice_ids=np.array([ell[1] for ell in ells])
    if not np.all(spice_ids==spice_ids[0]):
        raise ValueError("Different Spice ids got mixed in")
    # Note that in all the interps, we are going to get dense results, not just
    # missing results. This is fine, we just won't write results that are in
    # the sparse results.
    cxs=[float(x) for x in np.interp(dense_framenums,sparse_framenums,np.array([ell[2] for ell in ells]))]
    cys=[float(x) for x in np.interp(dense_framenums,sparse_framenums,np.array([ell[3] for ell in ells]))]

    # Special handling for major and minor axes. Since many ellipses are near-circular,
    # they can vary wildly in direction from frame to frame but still describe nearly
    # the same near-circle. We will interpolate them on a polar coordinate so that
    # the interpolated results are sensible length even while their azimuths vary
    # wildly.
    axs=np.array([ell[4] for ell in ells])
    ays=np.array([ell[5] for ell in ells])
    bxs=np.array([ell[6] for ell in ells])
    bys=np.array([ell[7] for ell in ells])
    # If the major axis was on the left side, flip it. Don't worry about
    # the minor axis because we only use its length, which doesn't depend
    # on beign flipped.
    # Index w identifies those cases.
    w=axs<0
    axs[w]=-axs[w]
    ays[w]=-ays[w]
    # All axes now flipped, don't refer to w from now on
    avs=np.vstack((axs,ays))
    bvs=np.vstack((bxs,bys))
    # Axis lengths. Use 'aa' instead of 'as' because
    # that's a reserved word. Use 'bb' instead of 'bs'
    # for consistency.
    bb=vlength(bvs)
    aa=vlength(avs)
    # Interpolate vectors in polar coordinates
    thetas=np.arctan(ays/axs)
    aa=np.interp(dense_framenums,sparse_framenums,aa)
    bb=np.interp(dense_framenums,sparse_framenums,bb)
    thetas=np.interp(dense_framenums,sparse_framenums,thetas)
    # Reconstruct a vector from length and direction
    axs=[float(x) for x in aa*np.cos(thetas)]
    ays=[float(x) for x in aa*np.sin(thetas)]
    # Since a and b axes are perpendicular, reconstruct b vector
    # from b length and a direction +90deg.
    bxs=[float(x) for x in bb*np.cos(thetas+np.pi/2)]
    bys=[float(x) for x in bb*np.sin(thetas+np.pi/2)]
    if conn is not None:
        with conn, closing(conn.cursor()) as cur:
            sql1="insert or replace into ellipses (framenum,spice_id,cx,cy,ax,ay,bx,by,source) values (?,?,?,?,?,?,?,?,?);"
            sql2="insert or ignore into frames (framenum,width,height,right_num) values (?,?,?,?);"
            for framenum,cx,cy,ax,ay,bx,by in zip(dense_framenums,cxs,cys,axs,ays,bxs,bys):
                if framenum in sparse_set:
                    continue
                row1=(int(framenum),int(spice_ids[0]),cx,cy,ax,ay,bx,by,1)
                row2=(int(framenum),1280,720,16.0)
                print(row1,row2)
                cur.execute(sql1,row1)
                cur.execute(sql2,row2)
    new_ells=[(int(framenum),int(spice_ids[0]),cx,cy,ax,ay,bx,by)
              for framenum,cx,cy,ax,ay,bx,by in zip(dense_framenums,cxs,cys,axs,ays,bxs,bys)]
    return new_ells


def cubic_bridge(*,x0:float,y0:float,x1:float,y1:float,x2:float,y2:float,x3:float,y3:float,
                   x:np.ndarray=None)->tuple[float,float,float,float]|np.ndarray:
    """
    Bridge between two linear segments
    :param x0: independent variable at point 0
    :param y0:   dependent variable at point 0
    :param x1: independent variable at point 1
    :param y1:   dependent variable at point 1
    :param x2: independent variable at point 2
    :param y2:   dependent variable at point 2
    :param x3: independent variable at point 3
    :param y3:   dependent variable at point 3
    :param x:  Optional independent variable values at which to evaluate
    :return: Either:
      * Coefficients a,b,c,d of cubic ax**3+bx**2+cx+d that bridges
        the two segments with a continuous value and derivative at each side
        OR
      * cubic evaluated at points passed in x
    """
    y1p=(y1-y0)/(x1-x0)
    y2p=(y3-y2)/(x3-x2)
    bv=np.array([[y1 ],
                 [y1p],
                 [y2 ],
                 [y2p]])
    A=np.array([[x1**3,x1**2,x1,1],
                [3*x1**2,2*x1,1,0],
                [x2**3,x2**2,x2,1],
                [3*x2**2,2*x2,1,0]])
    xv=np.linalg.solve(A,bv)
    (a,),(b,),(c,),(d,)=xv
    assert np.allclose(A@xv,bv)



def fill_out_zoom(*,conn:sqlite3.Connection,zoom:Zoom):
    """
    Figure out and fill in a zoom from everything else we have. We can do:
    * Maximum zoom from FOV analysis
    * Angular size of object from FOV and ellipse analysis
    * Handover from star angles to ellipse angles
    * Stage vector across zoom

    Before this is done: Make sure that the zoom has been measured as well as possible:
    * Align the predicted and observed position of the moon as well as possible in time.
      In starfit, it's not always possible to perfectly overlay the two positions,
      so prefer the one where the projected orbit bisects the moon equally.
    * Make sure that the ellipse is plotted for the first FOV image
    * Make sure there is a frame record for each FOV, with a good estimated ET.
      If necessary, use something like patch_vmbest_actual above. Go ahead and
      set the angle, width, height, and right_num too. Set the angle based on
      FOVs, width, heigh, and right_num from known prototype video properties
    * A patch will only be generated for frames which have an entry in the frames table
      so
    :param conn: Database
    :param zoom: Zoom parameters
    :return:
    """
    # Maximum zoom from FOV analysis of frames within [f1-f2].
    # Consider just using the first one - in principle since the FOV
    # is attached to the body, it moves during the rest of the zoom.

    # Look up the FOV - Since the best imaging sequences don't overlap,
    # we don't need to know which sequence.
    with closing(conn.cursor()) as cur:
        sql="select framenum,fov_size,angle from fovs where framenum>=? and framenum<=? order by framenum asc;"
        fovs=cur.execute(sql,(zoom.f0,zoom.f3)).fetchall()
        fov_framenums=np.array([fov[0] for fov in fovs])
        sql = "select framenum,spice_id,cx,cy,ax,ay,bx,by from ellipses where framenum>=? and framenum<=? and source=2 and spice_id>0"
        ells = cur.execute(sql, (zoom.f0, zoom.f3)).fetchall()
        sql = "select framenum,et, lat,lon,angle,clock,right_denom,width, height,right_num from frames where framenum>=? and framenum<=?"
        frames = cur.execute(sql, (zoom.f0,zoom.f3)).fetchall()
    ells=interp_ellipse(conn=conn, ells=ells)
    ell_framenums = np.array([ell[0] for ell in ells])
    print("FOV but no ellipse for frames ", set([int(x) for x in fov_framenums]) - set([int(x) for x in ell_framenums]))
    print("Ellipse but no FOV for frames ", set([int(x) for x in ell_framenums]) - set([int(x) for x in fov_framenums]))
    et_framenums=np.array([frame[0] for frame in frames])
    r_moons=[]
    for (fov_framenum, fov_size,angle) in fovs:
        i_et=np.where(et_framenums==fov_framenum)[0][0]
        i_ell=np.where(ell_framenums==fov_framenum)[0][0]
        ell_framenum, spice_id, cx, cy, ax, ay, bx, by=ells[i_ell]
        et_framenum,et, lat,lon,angle,clock,right_denom,width, height,right_num=frames[i_et]
        print(f"{spice_id=}")
        print(f"{fov_framenum=}, {fov_size=}, {angle=}")
        print(et)

        tan_fov_over_2=np.tan(np.deg2rad(angle/2))
        # Tangent subtended by one pixel
        tan_per_pix=tan_fov_over_2/(width/2)
        # Get distance to target moon
        moon_state,_=spkezr(str(spice_id),et,"ECLIPB1950","NONE","-32")
        moon_pos=moon_state[:3]
        moon_dist=vlength(moon_pos)
        print(f"{moon_dist=} km")
        #angle subtended by moon
        r_pix=(np.sqrt(ax**2+ay**2)+np.sqrt(bx**2+by**2))/2
        print(f"{r_pix=} pix")
        r_tan=tan_per_pix*r_pix
        #Photogrammetric radius of moon
        r_moon=moon_dist*np.arcsin(r_tan) # Yes, arcsin, because we are measuring to the horizon not to the plane through the center
        print(r_moon)
        r_moons.append(r_moon)
    r_moons=np.array(r_moons)
    r_moon=np.mean(r_moons)
    print(f"Mean {r_moon=}")
    # Now for the zoom in and zoom out. Use the moon size
    w_frames_f23=et_framenums>zoom.f2
    w_ells_f23=ell_framenums>zoom.f2
    print("Frame but no ellipse ", set(et_framenums[w_frames_f23]) - set(ell_framenums[w_ells_f23]))
    print("Ellipse but no frame ", set(ell_framenums[w_ells_f23]) - set(et_framenums[w_frames_f23]))
    # Now build a fresh zoom
    framenums=np.arange(np.max(np.hstack((fov_framenums, ell_framenums, et_framenums))))
    ets=np.zeros_like(framenums,dtype=np.float64)*np.nan
    lats=np.zeros_like(framenums,dtype=np.float64)*np.nan
    lons=np.zeros_like(framenums,dtype=np.float64)*np.nan
    angles=np.zeros_like(framenums,dtype=np.float64)*np.nan
    clocks=np.zeros_like(framenums,dtype=np.float64)*np.nan
    right_denoms=np.zeros_like(framenums,dtype=np.float64)*np.nan
    # For segment 12, use the first angle from the first image
    clock_stars0=frames[ 0][5]
    clock_stars3=frames[-1][5]
    clock_stars=(clock_stars0+clock_stars3)/2
    right_num=frames[0][9]
    right_denom_stars0=frames[0][6]
    right_denom_stars3=frames[-1][6]
    right_denom_stars=(right_denom_stars0+right_denom_stars3)/2
    # Angle from the first FOV
    angles[np.logical_and(zoom.f1<framenums,framenums<zoom.f2)]=fovs[0][2]
    # Clock from the mean of the first and last star images
    clocks[np.logical_and(zoom.f1<framenums,framenums<zoom.f2)]=clock_stars
    # Right denominator from first and last star images
    right_denoms[np.logical_and(zoom.f1<framenums,framenums<zoom.f2)]=right_denom_stars
    # Find the camera boresight at over segment 12 from the moon ellipse position
    for framenum_moon in et_framenums[np.logical_and(zoom.f1<et_framenums,et_framenums<zoom.f2)]:
        framenums[framenum_moon]=framenum_moon
        i_ell_moon=list(ell_framenums).index(framenum_moon )
        i_et_moon =list( et_framenums).index(framenum_moon )
        angle_moon=angles[framenum_moon]
        _,et_moon,*_=frames[i_et_moon]
        ets[framenum_moon]=et_moon
        _, spice_id, cx, cy, ax, ay, bx, by=ells[i_ell_moon]
        # Calculate position of moon relative to spacecraft in universe space
        x_moon,*_=spkezr(str(spice_id),et_moon,"ECLIPB1950","NONE","-32")
        moon_pos=x_moon[:3].reshape(-1,1)
        r_moon=vlength(moon_pos)
        # Pixel scale of an image at the same distance as the moon.
        tan_per_pix=np.arctan(np.deg2rad(angle_moon/2))/(width/2)
        km_per_pix=tan_per_pix*r_moon
        # Position of moon in camera space
        moon_pos_cam=np.array([[km_per_pix*(cx-width/2)],    # right component
                               [km_per_pix*(cy-height/2)],   # down component
                               [r_moon]])                    # boresight component
        # Now we use the power of point_toward(). The point body
        # vector is moon_pos_cam, while the toward body vector is
        # moon_pos. The toward body vector is Up (or maybe Down)
        # and the toward vector is the sky vector calculated from
        # Clock from the stars
        dir=vnormalize(moon_pos)
        for i in range(3):
            # The sky vector is weakly dependent on the direction, which in turn is
            # dependent on the sky, so iterate a few times. Once iteration is written,
            # it's just computer time to iterate once, 3 times, or 100 times.
            sky=make_sky(dir=dir,clock=clock_stars)
            M_uc=point_toward(p_b=moon_pos_cam,p_r=moon_pos,t_b=np.array([[0.0],[-1.0],[0.0]]),t_r=sky)
            # Now where does the boresight point in universe space?
            new_dir=M_uc @ np.array([[0.0],[0.0],[1.0]])
            print(f"{vangle(dir,new_dir,deg=True)=}")
            dir=new_dir
        lons[framenum_moon],lats[framenum_moon],_=xyz2llr(dir,deg=True)
        lons[framenum_moon]%=360.0
        # Test the round-trip. The position of the moon should fall right on the
        # center of the ellipse. Also M_uc should match the corresponding matrix
        # in the camera.
        camera=Camera(et=et_moon,
                      lat=lats[framenum_moon],lon=lons[framenum_moon],
                      angle=angles[framenum_moon],clock=clocks[framenum_moon],
                      right_denom=right_denoms[framenum_moon],width=width,height=height)
        camera._update()
        # Difference between matrices -- Suppose we have frames A and B, encoded as matrices M_Ab and M_Bb.
        # We can rotate from body to A, and then from A to B, by rotating with a difference matrix delta_BA
        # M_Bb=delta_BA @ M_Ab
        # We can find the difference by solving for it:
        # M_Bb @ M_Ab.inv = delta_BA @ M_Ab @ M_Ab.inv=delta_BA*I=delta_BA
        # Then we can show the difference in human-readable form by showing the magnitude of the axis-angle
        # representation of the difference
        delta=M_uc @ np.linalg.inv(camera.M_wc[0:3,0:3])
        print(vlength(m_to_aa(delta,deg=True)))
        print(camera)
        # Now round-trip the moon vector and the ellipse coordinates
        print(camera.project(vs_w=dir)) # This should be the exact center
        print(camera.project(vs_w=moon_pos))
        print(cx,cy)

    # Write stuff to a scratch file that overwrites the pointing from pointing_{casename}.inc
    with open(f"data/scratch/patch_pointing_VoyagerUranusHD_V{names[spice_id][0].upper()}BEST.inc","wt") as ouf:
        for i_frame,et,lon,lat,angle,clock,right_denom in zip(framenums,ets,lons,lats,angles,clocks,right_denoms):
            if not np.all((np.isfinite(et),np.isfinite(lat))):
                continue
            print(f"#declare FET [{i_frame:4d}]={et:14.3f};"
                  f"#declare FLat[{i_frame:4d}]={lat:9.6f};"
                  f"#declare FLon[{i_frame:4d}]={lon:9.6f};"
                  f"#declare FAngle[{i_frame:4d}]={angle:9.6f};"
                  f"#declare FTwist[{i_frame:4d}]={clock:9.6f};"
                  f"#declare FRight[{i_frame:4d}]={16.0/right_denom:9.6f};"
                  ,file=ouf)
    plt.figure("lon,lat")
    plt.plot(framenums,lons)
    plt.plot(framenums,lats)
    plt.figure("angle")
    plt.plot(framenums,angles,'.')
    plt.figure("clock")
    plt.plot(framenums,clocks,'.')
    plt.figure("right_denom")
    plt.plot(framenums,right_denoms,'.')
    plt.show()


def main():
    furnsh("data/spice/vgr2.tm")
    with closing(sqlite3.connect("data/db/frame_index_VoyagerUranusHD.sqlite")) as conn:
        patch_vmbest_actual(conn=conn)
        fill_out_zoom(conn=conn,zoom=zoom_VMBEST)


if __name__ == "__main__":
    main()


