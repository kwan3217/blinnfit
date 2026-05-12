"""
Test and debug entry point for Camera class

Created: 2/3/25
"""
import sqlite3
from contextlib import closing
from itertools import product

import numpy as np
import pytest
from kwanmath.geodesy import llr2xyz
from matplotlib import image as mpimg

from bsc import load_catalog, parse_stars
from starfit.camera import Camera


def test_camera():
    vs_w,names,mags,colors=parse_stars(load_catalog(),frame='ECLIPB1950')
    #vs_w=llr2xyz(lat=np.arange(0,10,1),lon=0)
    camera=Camera(lat=0,lon=0,clock=0,angle=45,right_denom=3,width=1280,height=720)
    pixs,on_screen=camera.project(vs_w)
    pixs=pixs[:,on_screen]
    names=names[on_screen]
    mags=mags[on_screen]
    colors=colors[:,on_screen]
    for i_star in range(pixs.shape[1]):
        print(f"{names[i_star]}: x={pixs[0,i_star]}, y={pixs[1,i_star]}, mag={mags[i_star]}, color={colors[:,i_star]}")


def test_camera_from_db(framenum:int=3500,casename:str="VoyagerUranusHD"):
    with closing(sqlite3.connect(f"data/db/frame_index_{casename}.sqlite")) as conn:
        camera=Camera.from_db(conn=conn, framenum=framenum, width=1280, height=720)
    vs_w=llr2xyz(lat=camera.lat,lon=camera.lon,r=1,deg=True)
    names=np.array(("Boresight star",))
    #v_w,names,mags,colors=parse_stars(load_catalog(limit_mag=6,count=4000))
    pixs,on_screen=camera.project(vs_w)
    print(pixs,on_screen)


def test_camera_sgr(framenum:int=3500,casename:str="VoyagerUranusHD"):
    vs_w,names,mags,colors=parse_stars(load_catalog(),frame='ECLIPB1950')
    infn = f"data/frames/{casename}/frame{framenum:04d}.png"
    img = mpimg.imread(infn)[:,:,0]
    with closing(sqlite3.connect(f"data/db/frame_index_{casename}.sqlite")) as conn:
        camera=Camera.from_db(conn=conn, framenum=framenum, width=1280, height=720)
    vs_w=llr2xyz(lat=camera.lat,lon=camera.lon,r=1,deg=True)
    names=np.array(("Boresight star",))
    #v_w,names,mags,colors=parse_stars(load_catalog(limit_mag=6,count=4000))
    pixs,on_screen=camera.project(vs_w)
    print(pixs,on_screen)

@pytest.mark.parametrize(
    "lat,lon,angle,clock,right_num,right_denom,width,height",
    [( 0,  0,   45,    0,       16,          9, 1280,720),
     ( 0,  0,   45,    0,        4,          3,  640,480),
     (90, 20,   45,    0,       16,          9, 1280,720),
     ]
)
def test_camera_round_trip(lat,lon,angle,clock,right_num,right_denom,width,height):
    camera=Camera(lat=lat,lon=lon,angle=angle,clock=clock,right_num=right_num,right_denom=right_denom,width=width,height=height)
    pixs_x=np.array([0,1,2,3,4])*width/4
    pixs_y=np.array([0,2,4])*height/4
    pixs_c=np.column_stack(list(product(pixs_x,pixs_y)))
    print(pixs_c.shape)
    print(pixs_c)
    vs_w=camera.project_inv(pixs_c)
    pixs_c_round_trip,*_=camera.project(vs_w,out_nan=False)
    assert np.allclose(pixs_c,pixs_c_round_trip)
