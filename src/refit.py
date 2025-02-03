import sqlite3
from contextlib import closing

from matplotlib import image as mpimg

from starfit.camera import Camera


def refit(casename:str):
    """
    For all frames which have an automatic fit solution, run the fit solution again

    """
    # Get a list of all the images to re-fit. These will be those with a finite sigma on lat.
    dbname = f"data/db/frame_index_{casename}.sqlite"
    with closing(sqlite3.connect(dbname)) as conn:
        sql = "select framenum from frames where lat_c_source<=3;"
        cur = conn.cursor()
        frames = [x[0] for x in cur.execute(sql).fetchall()]
        for i_frame in frames:
            infn = f"data/frames/{casename}/frame{i_frame:04d}.png"
            img = mpimg.imread(infn)[:, :, 0]
            height, width = img.shape
            old_camera=Camera.from_frame_db(i_frame,)
            new_camera=self.fit()