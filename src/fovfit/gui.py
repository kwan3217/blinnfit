import sqlite3
from contextlib import closing

import cv2
from kwanmath.conic import fit_conic, eval2_conic, identify_conic
from matplotlib.axes import Axes
from matplotlib.backend_tools import Cursors
import spiceypy as cspice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.widgets as widgets

cspice.furnsh('data/spice/vgr1.tm')
cspice.furnsh('data/spice/vgr2.tm')
cspice.furnsh('data/spice/lsk/naif0012.tls')
cspice.furnsh('data/spice/spk/de430.bsp')


class FovFit(object):
    def __init__(self,*,casename:int,initframe:int,seqid:str,imgid:int,fov_tan:float=0.003700098):
        #initialize from parameters
        self.casename=casename
        self.initframe=initframe
        self.framenum=self.initframe
        self.seqid=seqid
        self.imgid=imgid
        self.framepat=f"data/frames/{casename}/frame%04d.png"
        # Data from vg2_issna_v02.ti, field of view vector. Field of view is a square
        # with a boresight of 1 and an x and y of +-0.0037... . Convert this to
        # a full-width field of view in degrees
        self.has_img=False
        self.fov_tan=fov_tan
        self.fov_size = 2 * np.degrees(np.arctan(fov_tan))

        # Open the database
        self.open_fov_index()

        #set up graphics
        self.fig = plt.figure("Main fitting")
        self.ax = self.fig.add_subplot(111)
        self.xclicks=[]
        self.yclicks=[]
        self.load_image()
        self.step_size=10

        #set up buttons
        self.control_axs:list[Axes]=[]
        self.btns:list[widgets.Button]=[]

        # Function callbacks should have BTNxxx for buttons, CHKxxx fot checkboxes,
        # like from my old VB5 days. Then minimize the amount of code in the callbacks.
        self.fig_controls=plt.figure("Controls")
        self.makebtn(0.00, 0.05,'<frame', self.BTNframem)
        self.makebtn(0.10, 0.05,'>frame', self.BTNframep)
        self.makebtn(0.00, 0.15,'<10frame', self.BTNframem10)
        self.makebtn(0.10, 0.15,'>10frame', self.BTNframep10)
        self.makebtn(0.00, 0.25,'<100frame', self.BTNframem100)
        self.makebtn(0.10, 0.25,'>100frame', self.BTNframep100)
        self.fig.canvas.mpl_connect('button_press_event',self.PICclick)
        self.use_par={}
        #Once everything is loaded, do a replot to make sure it's visible
        self.ax.cursor_to_use = Cursors.SELECT_REGION

        def hover(event):
            if self.fig.canvas.widgetlock.locked():
                # Don't do anything if the zoom/pan tools have been enabled.
                return
            self.fig.canvas.set_cursor(
                event.inaxes.cursor_to_use if event.inaxes else Cursors.POINTER)

        self.fig.canvas.mpl_connect('motion_notify_event', hover)
        self.imshow()
    def open_fov_index(self):
        """
        Open an SQLite database and make sure that the appropriate table(s) are present
        in the database
        """
        dbname=f"data/db/frame_index_{self.casename}.sqlite"
        self.conn = sqlite3.connect(dbname)
        sql = ("create table if not exists fovs (" +
               "framenum    integer not null," +
               "timestamp datetime default CURRENT_TIMESTAMP,"
               "seqid string,"
               "imgid int,"
               "fov_size real,"
               "x0 real,y0 real,"
               "x1 real,y1 real,"
               "x2 real,y2 real,"
               "x3 real,y3 real,"
               "angle real,"
               "primary key (framenum,seqid,imgid))")
        with self.conn:
            self.conn.execute(sql)
    def makebtn(self,x,y,name,f):
        ax = self.fig_controls.add_axes((x, y, 0.05, 0.05))
        self.control_axs.append(ax)
        bx = widgets.Button(ax, name)
        self.btns.append(bx)
        bx.on_clicked(f)
    def load_image(self):
        infn = self.framepat % self.framenum
        self.backimg = mpimg.imread(infn)
        self.read_record()
    def imshow(self):
        if self.has_img:
            x0,x1,y0,y1=self.ax.axis()
        self.ax.clear()
        self.ax.imshow(self.backimg)
        if self.has_img:
            self.ax.axis([x0,x1,y0,y1])
        self.clicks,*_=self.ax.plot([],[],'g+',markersize=20,mew=3)
        self.y0plot,*_=self.ax.plot([],[],'b-')
        self.y1plot,*_=self.ax.plot([],[],'r-')
        self.center,*_=self.ax.plot([],[],'k*')
        self.ax.set_title(f"Frame {self.framenum}")
        self.plot_fov()
        self.has_img=True
        plt.pause(0.001)
    def set_frame(self,framenum):
        self.framenum=framenum
        self.load_image()
        self.read_record()
        self.imshow()
        self.plot_fov()
    def d_frame(self,d_frame):
        self.set_frame(self.framenum+d_frame)
    def figure_zoom(self):
        """
        Figure out the zoom level of a zoomed-in image with no stars based on the
        square marked field of view
        :param conn:
        :return:
        """
        # For Oberon, the moon is very near the center of the field of view.
        # It's also zoomed in a long way so the small angle approximation
        # is king.
        dx10 = self.xclicks[1]-self.xclicks[0]
        dy10 = self.yclicks[1]-self.yclicks[0]
        dx21 = self.xclicks[2]-self.xclicks[1]
        dy21 = self.yclicks[2]-self.yclicks[1]
        dx32 = self.xclicks[3]-self.xclicks[2]
        dy32 = self.yclicks[3]-self.yclicks[2]
        dx03 = self.xclicks[0]-self.xclicks[3]
        dy03 = self.yclicks[0]-self.yclicks[3]
        dr10 = np.sqrt(dx10 ** 2 + dy10 ** 2)
        dr21 = np.sqrt(dx21 ** 2 + dy21 ** 2)
        dr32 = np.sqrt(dx32 ** 2 + dy32 ** 2)
        dr03 = np.sqrt(dx03 ** 2 + dy03 ** 2)
        # All of these will be in degrees per pixel
        pix_scale_10 = self.fov_size / dr10
        pix_scale_21 = self.fov_size / dr21
        pix_scale_32 = self.fov_size / dr32
        pix_scale_03 = self.fov_size / dr03
        pix_scale = (pix_scale_10 + pix_scale_21 + pix_scale_32 + pix_scale_03) / 4
        print(f"{pix_scale=}")
        # amount of tangent at the center over one pixel
        tan_pix_scale = np.tan(np.deg2rad(pix_scale))
        # Tangent from center to horizontal edge
        tan_img_over_2 = self.backimg.shape[1] / 2 * tan_pix_scale
        # Horizontal render FOV is then 2*the angle with the
        # above tangent. Express it in degrees.
        angle = 2 * np.rad2deg(tan_img_over_2)
        print(f"{angle=}")
        return angle
    def update_click(self):
        this_xclicks=np.array(self.xclicks)
        this_yclicks=np.array(self.yclicks)
        self.clicks.set_data(this_xclicks,this_yclicks)
        if len(this_xclicks)==4:
            fields = {"seqid":self.seqid,
                      "imgid":self.imgid,
                      "fov_size":self.fov_size,
                      "x0":self.xclicks[0],
                      "y0":self.yclicks[0],
                      "x1":self.xclicks[1],
                      "y1":self.yclicks[1],
                      "x2":self.xclicks[2],
                      "y2":self.yclicks[2],
                      "x3":self.xclicks[3],
                      "y3":self.yclicks[3],
                      "angle":self.figure_zoom()
                     }
            values = tuple([v for k, v in fields.items()] + [self.framenum])
            sql = (f"insert or replace into fovs ({','.join([k for k, v in fields.items()])},timestamp,framenum) "
                   f"values ({','.join(['?' for k, v in fields.items()])},datetime('now'),?)")
            with self.conn:
                self.conn.execute(sql,values)
        self.plot_fov()
        plt.pause(0.001)
    def plot_fov(self):
        this_xclicks=np.array(self.xclicks)
        this_yclicks=np.array(self.yclicks)
        self.clicks.set_data(this_xclicks,this_yclicks)
        if len(this_xclicks)==4:
            self.y0plot.set_data(np.hstack((this_xclicks,this_xclicks[0])),np.hstack((this_yclicks,this_yclicks[0])))
    def read_record(self):
        sql=("select x0,y0,x1,y1,x2,y2,x3,y3,et,lat,lon,frames.angle as frame_angle,clock,right_num,right_denom "
             "from fovs left join frames on fovs.framenum=frames.framenum where fovs.framenum=? and fovs.seqid=? and fovs.imgid=?;")
        with self.conn:
            with closing(self.conn.cursor()) as cur:
                row=cur.execute(sql,(self.framenum,self.seqid,self.imgid)).fetchone()
        if row is not None:
            x0,y0,x1,y1,x2,y2,x3,y3,et,lat,lon,angle,clock,right_num,right_denom=row
            this_xclicks=np.array((x0,x1,x2,x3))
            this_yclicks=np.array((y0,y1,y2,y3))
            self.xclicks=list(this_xclicks)
            self.yclicks=list(this_yclicks)
        else:
            self.xclicks=[]
            self.yclicks=[]
    # Callbacks only below this point. Callbacks should be trivial.
    # If they are more than about 1 line, or if they are called from
    # somewhere else, break out the functionality into its own function
    # and call it from the callback.
    def BTNframem(self, event):
        self.d_frame(-1)
    def BTNframep(self, event):
        self.d_frame(+1)
    def BTNframem10(self, event):
        self.d_frame(-10)
    def BTNframep10(self, event):
        self.d_frame(+10)
    def BTNframem100(self, event):
        self.d_frame(-100)
    def BTNframep100(self, event):
        self.d_frame(+100)
    def PICclick(self,event):
        if self.fig.canvas.widgetlock.locked():
            # Don't do anything if the zoom/pan tools have been enabled.
            return
        if event.button==3:
            # Right mouse button - clear the ellipse
            self.xclicks=[]
            self.yclicks=[]
        else:
            self.xclicks.append(event.xdata)
            self.yclicks.append(event.ydata)
        self.update_click()
        print(event)


