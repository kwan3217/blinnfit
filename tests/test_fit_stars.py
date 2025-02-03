"""
Test and debug entry point for fit_stars

Created: 2/3/25
"""
import sqlite3
from contextlib import closing

import pytest
from matplotlib import image as mpimg
from matplotlib.pyplot import figure

from bsc import load_catalog, parse_stars
from starfit.camera import Camera
from starfit.fit_stars import fit_stars

def test_fit_stars(framenum:int=3500,casename:str="VoyagerUranusHD"):
    fig=figure("test_fit_stars box")
    ax_box=fig.gca()
    fig=figure("test_fit_stars img")
    ax_img=fig.gca()
    dbname = f"data/db/frame_index_{casename}.sqlite"
    infn = f"data/frames/{casename}/frame{framenum:04d}.png"
    img = mpimg.imread(infn)[:,:,0]
    height,width=img.shape
    with closing(sqlite3.connect(dbname)) as conn:
        camera=Camera.from_frame_db(conn,framenum,width=width,height=height)
    v_w,names,mags,colors=parse_stars(load_catalog(limit_mag=6,count=4000),frame='ECLIPB1950')
    camera_opt=fit_stars(img=img,star_vs=v_w,star_names=names,camera0=camera,ax_box=ax_box,ax_img=ax_img)
    print(camera_opt)


if __name__=="__main__":
    test_fit_stars()