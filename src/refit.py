import sqlite3
from contextlib import closing
from multiprocessing import Pool

import numpy as np
from matplotlib import image as mpimg

from bsc import parse_stars, load_catalog
from starfit.camera import Camera
from starfit.fit_stars import fit_stars

casename:str="VoyagerUranusHD"
dbname:str=f"data/db/frame_index_{casename}.sqlite"
v_w:np.ndarray=None
names:np.ndarray[str]=None
mags:np.ndarray=None
colors:np.ndarray=None


def process_frame(framenum: int):
    print(f"Frame {framenum:04d}")
    infn = f"data/frames/{casename}/frame{framenum:04d}.png"
    img = mpimg.imread(infn)[:, :, 0]
    height, width = img.shape
    with closing(sqlite3.connect(dbname)) as conn:
        old_camera = Camera.from_db(conn=conn, framenum=framenum, width=width, height=height)
        try:
            new_camera = fit_stars(img=img, star_vs=v_w, star_names=names, camera0=old_camera)
            new_camera.write_db(conn=conn, framenum=framenum)
        except Exception:
            import traceback
            traceback.print_exc()


def refit():
    """
    For all frames which have an automatic fit solution, run the fit solution again

    """
    # Get a list of all the images to re-fit. These will be those with a finite sigma on lat.
    global v_w,names,mags,colors
    v_w,names,mags,colors=parse_stars(load_catalog(limit_mag=6,count=4000),frame='ECLIPB1950')
    with closing(sqlite3.connect(dbname)) as conn:
        sql = "select framenum from frames where lat_source<=3;"
        cur = conn.cursor()
        frames = [x[0] for x in cur.execute(sql).fetchall()]
        print(len(frames))
    with Pool(24) as p:
        p.map(process_frame,frames)


def main():
    refit()


if __name__=="__main__":
    main()