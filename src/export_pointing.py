"""
Export the pointing from one of the fit databases. Ultimately
this will do fitting, interpolation, and shifting from the
star fit to the moon and FOV fit

Created: 2/7/25
"""
import argparse
import sqlite3
from contextlib import closing
from dataclasses import dataclass
from datetime import datetime

import numpy as np
from kwanmath.interp import linterp

from starfit.camera import Source, Camera


@dataclass
class frame_rowtype:
    # framenum:int # Use this as a dictionary key instead
    timestamp:datetime
    nstars:Source
    rmsdiff:float
    et:float
    et_sig:float
    et_source:Source
    lat:float
    lat_sig:float
    lat_source:Source
    lon:float
    lon_sig:float
    lon_source:Source
    angle:float
    angle_sig:float
    angle_source:Source
    clock:float
    clock_sig:float
    clock_source:Source
    right_denom:float
    right_denom_sig:float
    right_denom_source:Source
    width:int
    height:int
    right_num:float

@dataclass
class ellipse_rowtype:
    framenum:int
    cx:float
    cy:float
    r:float


@dataclass
class ZoomRawData:
    width:int
    height:int
    frame_number:int
    spice_id:int
    x0:int
    y0:int
    x1:int
    y1:int
    x2:int
    y2:int
    x3:int
    y3:int
    fov_size:float
    cx:float=None
    cy:float=None
    ax:float=None
    ay:float=None
    bx:float=None
    by:float=None


def figure_zoom(*,conn:sqlite3.Connection,zoom_raw_data:ZoomRawData):
    """
    Figure out the zoom level of a zoomed-in image with no stars based on the
    square marked field of view
    :param conn:
    :return:
    """
    # For Oberon, the moon is very near the center of the field of view.
    # It's also zoomed in a long way so the small angle approximation
    # is king.
    dx10=zoom_raw_data.x1-zoom_raw_data.x0
    dy10=zoom_raw_data.y1-zoom_raw_data.y0
    dx21=zoom_raw_data.x2-zoom_raw_data.x1
    dy21=zoom_raw_data.y2-zoom_raw_data.y1
    dx32=zoom_raw_data.x3-zoom_raw_data.x2
    dy32=zoom_raw_data.y3-zoom_raw_data.y2
    dx03=zoom_raw_data.x0-zoom_raw_data.x3
    dy03=zoom_raw_data.y0-zoom_raw_data.y3
    dr10=np.sqrt(dx10**2+dy10**2)
    dr21=np.sqrt(dx21**2+dy21**2)
    dr32=np.sqrt(dx32**2+dy32**2)
    dr03=np.sqrt(dx03**2+dy03**2)
    # All of these will be in degrees per pixel
    pix_scale_10=zoom_raw_data.fov_size/dr10
    pix_scale_21=zoom_raw_data.fov_size/dr21
    pix_scale_32=zoom_raw_data.fov_size/dr32
    pix_scale_03=zoom_raw_data.fov_size/dr03
    pix_scale=(pix_scale_10+pix_scale_21+pix_scale_32+pix_scale_03)/4
    print(f"{pix_scale=}")
    # amount of tangent at the center over one pixel
    tan_pix_scale=np.tan(np.deg2rad(pix_scale))
    # Tangent from center to horizontal edge
    tan_img_over_2=zoom_raw_data.width/2*tan_pix_scale
    # Horizontal render FOV is then 2*the angle with the
    # above tangent. Express it in degrees.
    angle=2*np.rad2deg(tan_img_over_2)
    print(f"{angle=}")


def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Export camera data for the derived case")
    parser.add_argument("case_name", default="VoyagerUranusHD",type=str, help="The name of the case")

    args = parser.parse_args()

    print(f"Exporting pointing for the following case: {args.case_name}")

    dbname = f"data/db/frame_index_{args.case_name}.sqlite"
    oufn=f"data/scratch/pointing_{args.case_name}.inc"
    # Suck in all rows
    frames:dict[int,frame_rowtype]={}
    with closing(sqlite3.connect(dbname)) as conn, open(oufn,"wt") as ouf:
        with closing(conn.cursor()) as cur:
            frames={row[0]:frame_rowtype(*row[1:]) for row in cur.execute("select * from frames;")}
        n=np.max(np.array(list(frames.keys())))+1
        print(f"#declare NominalImageWidth={frames[1].width:4d};",file=ouf)
        print(f"#declare NominalImageHeight={frames[1].height:4d};",file=ouf)
        print(f"#declare FET=array[{n}]",file=ouf)
        print(f"#declare FLat=array[{n}]",file=ouf)
        print(f"#declare FLon=array[{n}]",file=ouf)
        print(f"#declare FAngle=array[{n}]",file=ouf)
        print(f"#declare FTwist=array[{n}]",file=ouf)
        print(f"#declare FRight=array[{n}]",file=ouf)
        print(f"#declare Fpixc_sc=array[{n}] // Pixel vector of center of Voyager, and size of dish",file=ouf)
        for i_frame,frame in frames.items():
            if frame.lat is None:
                print(f"// Some elements of frame {i_frame:4d} are none: {frame}",file=ouf)
            else:
                camera=Camera.from_db(conn=conn,framenum=i_frame)
                pixs_c=np.vstack((np.array([0,1,0,1])*camera.width,
                                  np.array([0,0,1,1])*camera.height))
                vs_u=camera.project_inv(pixs_c)
                print(f"#declare FET [{i_frame:4d}]={frame.et:14.3f};"
                      f"#declare FLat[{i_frame:4d}]={frame.lat:9.6f};"
                      f"#declare FLon[{i_frame:4d}]={frame.lon:9.6f};"
                      f"#declare FAngle[{i_frame:4d}]={frame.angle:9.6f};"
                      f"#declare FTwist[{i_frame:4d}]={frame.clock:9.6f};"
                      f"#declare FRight[{i_frame:4d}]={frame.right_num/frame.right_denom:9.6f};"
                      ,file=ouf)
            with closing(conn.cursor()) as cur:
                sc_pos=[ellipse_rowtype(*row) for row in cur.execute("select framenum,cx,cy,(sqrt(ax*ax+ay*ay)+sqrt(bx*bx+by*by))/2 as r from ellipses where spice_id=? order by abs(framenum-?) asc",(-32,i_frame)).fetchmany(2)]
                if sc_pos[0].framenum==i_frame:
                    # Exact match, use it
                    cx=sc_pos[0].cx
                    cy=sc_pos[0].cy
                    r=sc_pos[0].r
                elif len(sc_pos)==2:
                    # No exact match, interpolate
                    cx=linterp(sc_pos[0].framenum,sc_pos[0].cx,sc_pos[1].framenum,sc_pos[1].cx,i_frame)
                    cy=linterp(sc_pos[0].framenum,sc_pos[0].cy,sc_pos[1].framenum,sc_pos[1].cy,i_frame)
                    r =linterp(sc_pos[0].framenum,sc_pos[0].r ,sc_pos[1].framenum,sc_pos[1].r ,i_frame)
                else:
                    # No rows at all, shouldn't happen
                    raise ValueError("No ellipse rows at all in database")
                print(f"#declare Fpixc_sc[{i_frame:4d}]=<{cx:10.4f},{cy:10.4f},{r:10.4f}>;",file=ouf)



if __name__ == "__main__":
    main()


