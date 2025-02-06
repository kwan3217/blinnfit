from contextlib import closing
from functools import partial

from matplotlib.axes import Axes
from spiceypy import furnsh, spkezr, timout, edlimb, pxform, SpiceFRAMEDATANOTFOUND
import spiceypy as cspice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.widgets as widgets
import sqlite3

from starfit.camera import Camera
from starfit.fit_stars import fit_stars
from starfit.videos import projects, ProjectBody
from bsc import load_catalog, parse_stars

cspice.furnsh('data/spice/vgr1.tm')
cspice.furnsh('data/spice/vgr2.tm')
cspice.furnsh('data/spice/lsk/naif0012.tls')
cspice.furnsh('data/spice/spk/de430.bsp')


class CameraMount(object):
    q:np.ndarray=np.arange(0,np.pi*2,0.01)
    c:np.ndarray=np.cos(q)
    s:np.ndarray=np.sin(q)
    # Resistor color code
    colors = {0: '#404040', 1: '#804000', 2: '#ff0000', 3: '#ff8000', 4: '#ffff00', 5: '#00ff00', 6: '#0000ff',
             7: '#8000ff', 8: '#c0c0c0', 9: '#ffffff'}

    def __init__(self,*,casename:int,initframe:int):
        #initialize from parameters
        self.casename=casename
        self.initframe=initframe

        #initialize by table lookup
        self.project=projects[self.casename]
        self.framepat:str=self.project.framepat
        self.goodstars:set[int]=self.project.goodstars
        self.badstars:set[int]=self.project.badstars
        self.spiceobjs:dict[int,ProjectBody]=self.project.bodies
        self.ringcenter:int=list(self.spiceobjs.keys())[0]
        self.rings:dict[float]=self.project.rings
        self.framenum:int=self.initframe

        #set up graphics
        self.fig = plt.figure("Main fitting")
        self.ax = self.fig.add_subplot(111)

        #initialize by reading and calculation
        furnsh("data/spice/vgr2.tm")
        self.open_frame_index()
        self.read_record()
        self.figimg=None
        self.load_image()
        self.figimg = self.ax.imshow(self.backimg)
        self.step_size=10

        # Load stars
        self.star_vec,self.starnames,*_=parse_stars(load_catalog(limit_mag=6,count=6000),frame="ECLIPB1950")

        #set up buttons
        self.control_axs:list[Axes]=[]
        self.btns:list[widgets.Button]=[]

        # Function callbacks should have BTNxxx for buttons, CHKxxx fot checkboxes,
        # like from my old VB5 days. Then minimize the amount of code in the callbacks.
        self.fig_controls=plt.figure("Controls")
        self.makebtn(0.55, 0.05,'-camlon', self.BTNcamlonm)
        #self.makebtn(0.75, 0.00,'-var', self.BTNvarm)
        #self.makebtn(0.75, 0.10,'+var', self.BTNvarp)
        self.makebtn(0.80, 0.00,'-clock', self.BTNclockm)
        self.makebtn(0.80, 0.10,'+clock', self.BTNclockp)
        self.makebtn(0.45, 0.05,'+camlon', self.BTNcamlonp)
        self.makebtn(0.50, 0.00,'-camlat', self.BTNcamlatm)
        self.makebtn(0.50, 0.10,'+camlat', self.BTNcamlatp)
        self.makebtn(0.95, 0.00,'-angle', self.BTNanglem)
        self.makebtn(0.95, 0.10,'+angle', self.BTNanglep)
        self.makebtn(0.90, 0.00,'-time', self.BTNtimem)
        self.makebtn(0.90, 0.05,'$time', self.BTNtimeconf)
        self.makebtn(0.90, 0.10,'+time', self.BTNtimep)
        self.makebtn(0.85, 0.00,'-right', self.BTNrightm)
        self.makebtn(0.85, 0.10,'+right', self.BTNrightp)
        self.makebtn(0.00, 0.00,'/step', self.BTNsm)
        self.makebtn(0.10, 0.00,'*step', self.BTNbig)
        self.makebtn(0.00, 0.05,'<frame', self.BTNframem)
        self.makebtn(0.10, 0.05,'>frame', self.BTNframep)
        self.makebtn(0.05, 0.00,'fit', self.BTNfit)
        self.makebtn(0.00, 0.10,'<auto', self.BTNautom)
        self.makebtn(0.10, 0.10,'auto>', self.BTNautop)
        self.use_par={}
        self.makechk(0.50,0.15,'lat')
        self.makechk(0.45,0.15,'lon')
        self.makechk(0.95,0.15,'angle')
        self.makechk(0.80,0.15,'clock')
        self.makechk(0.85,0.15,'right')
        self.makebtn(0.55, 0.35,'-dlon', self.BTNcamlonm)
        self.makebtn(0.75, 0.30,'-size', self.BTNvarm)
        self.makebtn(0.75, 0.40,'+size', self.BTNvarp)
        self.makebtn(0.80, 0.30,'-twist', self.BTNclockm)
        self.makebtn(0.80, 0.40,'+twist', self.BTNclockp)
        self.makebtn(0.45, 0.35,'+dlon', self.BTNcamlonp)
        self.makebtn(0.50, 0.30,'-dlat', self.BTNcamlatm)
        self.makebtn(0.50, 0.40,'+dlat', self.BTNcamlatp)
        #Once everything is loaded, do a replot to make sure it's visible
        self.replot()
    def makebtn(self,x,y,name,f):
        ax = self.fig_controls.add_axes((x, y, 0.05, 0.05))
        self.control_axs.append(ax)
        bx = widgets.Button(ax, name)
        self.btns.append(bx)
        bx.on_clicked(f)
    def makechk(self,x,y,name):
        ax = self.fig_controls.add_axes((x, y, 0.05, 0.05))
        self.control_axs.append(ax)
        bx = widgets.CheckButtons(ax, [name],[True])
        self.use_par[name]=bx
    def open_frame_index(self):
        """
        Open an SQLite database and make sure that the appropriate table(s) are present
        in the database
        """
        dbname=f"data/db/frame_index_{self.casename}.sqlite"
        self.conn = sqlite3.connect(dbname)
        sql = ("create table if not exists frames (" +
               "framenum    integer not null," +
               "timestamp datetime default CURRENT_TIMESTAMP,"
               "nstars          integer,"+
               "rmsdiff         real," +
               "et              real," +
               "et_sig          real," +
               "et_source          real," +
               "lat           real," +
               "lat_sig       real," +
               "lat_source    integer," +
               "lon           real," +
               "lon_sig       real," +
               "lon_source    integer," +
               "angle           real," +
               "angle_sig       real," +
               "angle_source    integer," +
               "clock           real," +
               "clock_sig       real," +
               "clock_source    integer," +
               "right_denom     real," +
               "right_denom_sig real," +
               "right_denom_source integer," +
               "primary key (framenum))")
        cur = self.conn.cursor()
        cur.execute(sql)
        self.conn.commit()
    def read_record(self):
        #Check if this frame is already recorded
        self.camera=Camera.from_db(conn=self.conn, framenum=self.framenum)
    def write(self):
        self.camera.write_db(conn=self.conn,framenum=self.framenum)
    def load_image(self):
        infn = self.framepat % self.framenum
        self.backimg = mpimg.imread(infn)
        if self.figimg is not None:
            self.figimg.set_data(self.backimg)
        self.ax.set_title(f"Frame {self.framenum}")
    def autofit(self,dframe,limit=8000):
        for i in range(limit):
            self.d_frame(dframe)
            self.fit()
    def plotspice(self,this_camera:Camera=None):
        if this_camera is None:
            this_camera=self.camera
        # print(f"Voyager 2 in kernel {which_kernel('SPK',-32,self.camera.et)}")
        for i_spice,body in self.spiceobjs.items():
            try:
                if body.dt is not None and body.dt>0:
                    rofs=body.dt*spkezr(str(i_spice),this_camera.et,"ECLIPB1950","NONE",body.parent)[0][3:].reshape(-1,1)
                else:
                    rofs=np.zeros((3,1))
                bodyframe=f"IAU_{body.name.upper()}"
                #Vector labels have two letters:
                # * First is center, one of:
                #    - v: Voyager Spacecraft
                #    - b: center of body in question
                # * Second is coordinate frame
                #    - b: body-fixed frame of body in question
                #    - i: Global inertial frame (Ecliptic B1950)
                # Observe the target natural spice object from Voyager
                # at the given time, in the body's own body-fixed frame
                try:
                    xbody_vb,lt=spkezr(str(i_spice),this_camera.et,bodyframe,"NONE","-32")
                except SpiceFRAMEDATANOTFOUND:
                    # Frame not found -- no frame is included for Nereid
                    bodyframe="ECLIPB1950"
                    xbody_vb,lt=spkezr(str(i_spice),this_camera.et,bodyframe,"NONE","-32")
                # position and velocity of Voyager relative to body is
                # reverse of position of body relative to Voyager. Likewise
                # for velocity.
                xvoy_bb=-xbody_vb
                rvoy_bb=xvoy_bb[:3]
                if True:
                    # Use pre-encounter spherical radius for all ellipsoid radii
                    a,b,c=body.r,body.r,body.r
                else:
                    a,b,c=body.a,body.b,body.c
                if a>0:
                    ell_bb=edlimb(a,b,c,rvoy_bb)
                    # limb points in body centered body frame
                    rlimb_bb=self.c*ell_bb.semi_major.reshape(-1,1)+self.s*ell_bb.semi_minor.reshape(-1,1)+ell_bb.center.reshape(-1,1)
                else:
                    # If body size is zero, then put all the limb points at the center
                    rlimb_bb=np.zeros((3,len(self.c)))
                rvoy_bb = rvoy_bb.reshape(-1, 1)
                # limb points in Voyager-centered body frame
                rlimb_vb=rlimb_bb-rvoy_bb
                M_ib=pxform(bodyframe,"ECLIPB1950",this_camera.et)
                # limb points in Voyager-centered inertial frame
                rlimb_vi=M_ib @ rlimb_vb+rofs
                # Stack up into homogeneous vectors
                rlimb_vi=np.vstack((rlimb_vi,np.zeros((1,rlimb_vi.shape[1]))))
                pixlimb,_= this_camera.project(vs_w=rlimb_vi)

                # print(f"{i_spice} in kernel {which_kernel('SPK',i_spice,self.camera.et)}")
                self.ax.plot(pixlimb[0,:],pixlimb[1,:],'-',color=self.colors[i_spice%10])
                if body.parent is not None:
                    et_orbit=(np.arange(100)-50)*60+this_camera.et
                    v_orbit=np.ones((4,100))
                    rmoon_vi=spkezr(str(i_spice),this_camera.et,"ECLIPB1950","NONE","-32")[0][0:3].reshape(-1,1)
                    for i_et,this_et in enumerate(et_orbit):
                        v_orbit[:3,i_et]=spkezr(str(i_spice),this_et,"ECLIPB1950","NONE",str(body.parent))[0][0:3]
                    v_orbit[0:3,:]-=v_orbit[0:3,None,50]
                    v_orbit[0:3,:]+=rmoon_vi
                    pixorbit,*_= this_camera.project(vs_w=v_orbit)
                    self.ax.plot(pixorbit[0,:],pixorbit[1,:],color=self.colors[i_spice%10])
                    self.ax.text(pixorbit[0,50],pixorbit[1,50],body.name,color=self.colors[i_spice%10])
            except Exception:
                import traceback
                traceback.print_exc()
                print(f"{i_spice} not in kernels")
        for i_ring,ring_r in enumerate(self.rings):
            name="Uranus"
            i_spice=self.ringcenter
            bodyframe = f"IAU_{name.upper()}"
            xbody_vb, lt = spkezr(str(i_spice), this_camera.et, bodyframe, "NONE", "-32")
            xvoy_bb = -xbody_vb
            rvoy_bb = xvoy_bb[:3]
            # ring points in body centered body frame
            rring_bb = self.c*np.array([[ring_r],[0.0],[0.0]]) + self.s*np.array([[0.0],[ring_r],[0]])
            rvoy_bb = rvoy_bb.reshape(-1, 1)
            # ring points in Voyager-centered body frame
            rring_vb = rring_bb - rvoy_bb
            M_ib = pxform(bodyframe, "ECLIPB1950", this_camera.et)
            # ring points in Voyager-centered inertial frame
            rring_vi = M_ib @ rring_vb
            # Stack up into homogeneous vectors
            rring_vi = np.vstack((rring_vi, np.zeros((1, rring_vi.shape[1]))))
            pixring,_ = this_camera.project(vs_w=rring_vi)
            self.ax.plot(pixring[0, :],pixring[1, :],'c-')
    def imshow(self):
        self.ax.clear()
        self.ax.imshow(self.backimg)
        self.ax.axis('scaled')
        self.ax.set_xlabel(timout(self.camera.et, "YYYY-MM-DD HR:MN:SC.###::UTC"))
        self.ax.set_title(f"Frame {self.framenum}")
        self.ax.set_ylabel(f"lat={self.camera.lat:.1f}, lon={self.camera.lon:.1f}, "
                           f"angle={self.camera.angle:.1f}, right={self.camera.right_num:.0f}/{self.camera.right_denom:.3f}")
    def plot_predicted_stars(self,this_camera:Camera=None):
        if this_camera is None:
            this_camera=self.camera
        star_pix,w = this_camera.project(vs_w=self.star_vec)
        self.ax.plot(star_pix[0,w],star_pix[1,w],'.',markerfacecolor='none',color='#ffff00')
        for x,y,name in zip(star_pix[0,w],star_pix[1,w],self.starnames[w]):
            self.ax.text(x,y,name,color='#ffff00')
    def replot(self):
        self.imshow()
        self.plot_predicted_stars()
        def set_btn(name,sig):
            if (self.use_par[name].get_status()[0] ^ np.isfinite(sig)):
                self.use_par["lat"].set_active(0)
        set_btn("lat",self.camera.lat_sig)
        set_btn("lon",self.camera.lon_sig)
        set_btn("angle",self.camera.angle_sig)
        set_btn("clock",self.camera.clock_sig)
        set_btn("right",self.camera.right_denom_sig)
        self.plotspice()
        try:
            plt.pause(0.001)
        except Exception:
            pass
    def f_findstars(self,*,starpixo:np.ndarray,starpixc:np.ndarray,w:np.ndarray[bool]=None,names:np.ndarray[str]=None):
        self.imshow()
        self.plotspice()
        if w is None:
            w=np.array([True]*starpixo.shape[1])
        self.ax.plot(starpixo[0,w],starpixo[1,w],'gx')
        self.ax.plot(starpixc[0,w],starpixc[1,w],'r+')
        if names is not None:
            for x,y,name in zip(starpixo[0,w],starpixc[1,w],names[w]):
                self.ax.text(x,y,name,color='#ffff00')
        plt.pause(1)
    def fit(self):
        self.camera=fit_stars(img=self.backimg[:,:,0],star_vs=self.star_vec,star_names=self.starnames,
                              camera0=self.camera,
                              f_findstars=self.f_findstars)
    def set_frame(self,framenum):
        self.write()
        self.framenum=framenum
        self.read_record()
        self.load_image()
        self.replot()
    def d_frame(self,d_frame):
        self.set_frame(self.framenum+d_frame)
    # Callbacks only below this point. Callbacks should be trivial.
    # If they are more than about 1 line, or if they are called from
    # somewhere else, break out the functionality into its own function
    # and call it from the callback.
    def BTNbig(self, event):
        self.step_size*=10
    def BTNsm(self, event):
        self.step_size/=10
    def BTNcamlonp(self, event):
        self.camera.lon_source=2
        self.camera.lon+=self.step_size
        self.replot()
    def BTNcamlonm(self, event):
        self.camera.lon_source=2
        self.camera.lon-=self.step_size
        self.replot()
    def BTNcamlatp(self, event):
        self.camera.lat_source=2
        self.camera.lat+=self.step_size
        self.replot()
    def BTNcamlatm(self, event):
        self.camera.lat_source=2
        self.camera.lat-=self.step_size
        self.replot()
    def BTNanglep(self, event):
        self.camera.angle_source=2
        self.camera.angle += self.step_size
        self.replot()
    def BTNanglem(self, event):
        self.camera.angle_source=2
        self.camera.angle -= self.step_size
        self.replot()
    def BTNrightp(self, event):
        self.camera.right_denom_source=2
        self.camera.right_denom += self.step_size
        self.replot()
    def BTNrightm(self, event):
        self.camera.right_denom_source=2
        self.camera.right_denom -= self.step_size
        self.replot()
    def BTNvarp(self, event):
        self.spiceobjs[801][5] += self.step_size
        print(self.spiceobjs[801])
        self.replot()
    def BTNvarm(self, event):
        self.spiceobjs[801][5] -= self.step_size
        print(self.spiceobjs[801])
        self.replot()
    def BTNclockp(self, event):
        self.camera.clock_source=2
        self.camera.clock += self.step_size
        self.replot()
    def BTNclockm(self, event):
        self.camera.clock_source=2
        self.camera.clock -= self.step_size
        self.replot()
    def BTNtimep(self, event):
        self.camera.et_source=2
        self.camera.et += self.step_size*60
        self.replot()
    def BTNtimem(self, event):
        self.camera.et_source=2
        self.camera.et -= self.step_size*60
        self.replot()
    def BTNtimeconf(self, event):
        self.camera.et_source=2
        self.replot()
        self.write()
    def BTNframem(self, event):
        self.d_frame(-1)
    def BTNframep(self, event):
        self.d_frame(+1)
    def BTNfit(self, event):
        self.fit()
    def BTNautop(self, event):
        self.autofit(+1)
    def BTNautom(self, event):
        self.autofit(-1)


