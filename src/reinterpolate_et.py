"""
Some frames have a carefully (manually) marked ET

Created: 1/30/25
"""
import sqlite3
from contextlib import closing

import numpy as np
from matplotlib import pyplot as plt


def main():
    with closing(sqlite3.connect("data/db/frame_index_VoyagerUranusHD.sqlite")) as conn:
        # Load all the keyframes
        keyframes=[]
        sql = f"select framenum,et from frames where et_source>1 order by framenum asc"
        print(sql)
        with closing(conn.cursor()) as cur:
            for this_row in cur.execute(sql):
                keyframes.append(this_row)
        keyframes=np.array(keyframes)
        print(keyframes)
        # Build the interpolator
        frames=keyframes[:,0]
        ets=keyframes[:,1]
        kf_interp=lambda frame:np.interp(frame,frames,ets)
        # Plot the interpolator
        plt.plot(frames,ets,'+')
        # Load the original data
        sql = f"select framenum,et from frames where et_source=1 order by framenum asc"
        print(sql)
        interpframes=[]
        with closing(conn.cursor()) as cur:
            for this_row in cur.execute(sql):
                interpframes.append(this_row)
        interpframes=np.array(interpframes)
        iframes=interpframes[:,0]
        iets=interpframes[:,1]
        plt.plot(iframes,iets,'-')
        plt.plot(iframes,kf_interp(iframes),'-')
        # Write back the interpolated data
        with closing(conn.cursor()) as cur:
            for framenum,et in zip(iframes,kf_interp(iframes)):
                sql=f"update frames set et=? where framenum=?"
                print(sql,(int(framenum),et))
                cur.execute(sql,(et,int(framenum)))
            conn.commit()

        plt.show()

if __name__ == "__main__":
    main()
