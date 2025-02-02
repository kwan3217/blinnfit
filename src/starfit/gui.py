from contextlib import closing

from spiceypy import furnsh, spkezr, str2et, timout, edlimb, pxform, SpiceFRAMEDATANOTFOUND
import spiceypy as cspice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.widgets as widgets
import sqlite3

from kwanmath.gaussian import correlation_matrix, infamily
from kwanmath.geodesy import llr2xyz
from kwanmath.optimize import curve_fit, bounded, positive, rbounded

from starfit.camera import cmatrix, project, make_sky
from starfit.videos import projects, ProjectBody
from bsc import load_catalog, GetMag, GetDec, GetRA, GetName

from kwanmath.interp import linterp

from find_stars import find_stars
from which_kernel import which_kernel

cspice.furnsh('data/spice/vgr1.tm')
cspice.furnsh('data/spice/vgr2.tm')
cspice.furnsh('data/spice/lsk/naif0012.tls')
cspice.furnsh('data/spice/spk/de430.bsp')


def curve_fitsky_interface(starvec,lat_c,lon_c,angle,clock,right_denom,*,width,height):
    """
    Calculate the pixel positions of the given stars, given these camera parameters
    :param starvec: List of star vectors with homogeneous coordinates of shape (4,M//2)
    :param camlat: Scalar camera latitude in degrees
    :param camlon: Scalar camera longitude in degrees
    :param dist:   Scalar camera distance in AU
    :param angle:  Scalar camera FOV angle in degrees
    :param right_denom: Scalar camera aspect ratio constant
    :param cx:     Scalar image distortion center horizontal coordinate
    :param cy:     Scalar image distortion center vertical coordinate
    :param x_sky: x coordinate of unit sky vector in degrees
    :param ym_sky: modified y coordinate of unit sky vector in degrees
      So that the curve fitter may freely choose any x_sky,ym_sky pair without the bounds
      of y dependent on x, we use ym_sky=
    :return: 1D array of shape M, representing a 2D array of pixel coordinates [x_or_y,star] shape (2,M//2)
             raveled so as to work with scipy.optimize.curve_fit. This will be all the x coordinates first,
             then all the y coordinates
    """
    # The curve fitter scipy.optimize.curve_fit takes a function to fit f,
    # independent xdata (can be any object, but f(xdata,*p) must return an
    # array of shape M), dependent ydata (shape M), and an initial guess at
    # a set of parameters p0 (shape N). It returns a set of parameters popt
    # which best fits the data. In our case, the ydata is pixel positions of
    # the stars, and therefore M is twice the number of stars we are trying
    # to fit. The p is camera parameters, and therefore by process of elimination
    # the xdata must be the positions of the stars. In our case, it's easiest to take
    # the vectors of the stars as inputs, so xdata will be an array of shape
    # (4,M//2) and we will return a 1D array of raveled x and y pixel coordinates
    # of each star
    starvec=starvec.reshape(4,-1)
    dir = llr2xyz(lat=lat_c, lon=lon_c)
    sky = make_sky(clock, dir)
    C = cmatrix(loc=np.zeros((3, 1)), look=dir, sky=sky)
    Cv = C @ starvec
    prj = project(right=4 / right_denom, angle=angle, width=width, height=height, target_c=Cv, out_nan=False)
    if not np.all(np.isfinite(prj)):
        print("Dir: \n",dir)
        print("Sky: \n",sky)
        print("C:   \n",C)
        print("Cv:  \n",Cv)
        print("prj: \n",prj)
        print("Nonfinite: \n",np.where(np.isnan(prj[0,:])))
        prj = project(right=4 / right_denom, angle=angle, width=width, height=height, target_c=Cv, out_nan=False)
        raise AssertionError("Not all projected star locations are finite")
    return prj.ravel()


class CameraMount(object):
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
        self.rings:dict[float]=self.project.rings
        self.framenum:int=self.initframe

        #set up graphics
        self.fig = plt.figure("Main fitting")
        self.ax = self.fig.add_subplot(111)

        #initialize by reading and calculation
        furnsh("data/spice/vgr2.tm")
        self.et=None
        self.open_frame_index()
        self.read_record()
        self.figimg=None
        self.load_image()
        self.figimg = self.ax.imshow(self.backimg)
        self.width = self.backimg.shape[1]
        self.height = self.backimg.shape[0]
        self.step_size=10
        self.rmsdiff=None
        self.loadstars()

        #set up buttons
        self.axs=[]
        self.btns=[]

        # Function callbacks should have BTNxxx for buttons, CHKxxx fot checkboxes,
        # like from my old VB5 days. Then minimize the amount of code in the callbacks.
        self.fig_controls=plt.figure("Controls")
        self.makebtn(0.55, 0.05,'-camlon', self.BTNcamlonm)
        self.makebtn(0.75, 0.00,'-var', self.BTNvarm)
        self.makebtn(0.75, 0.10,'+var', self.BTNvarp)
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
        self.makebtn(0.00,0.10,'<auto',self.autom)
        self.makebtn(0.10,0.10,'auto>',self.autop)
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
        self.makebtn(0.00, 0.75,'refit', self.BTNrefit)
        #Once everything is loaded, do a replot to make sure it's visible
        self.replot()
    def makebtn(self,x,y,name,f):
        ax = self.fig_controls.add_axes((x, y, 0.05, 0.05))
        self.axs.append(ax)
        bx = widgets.Button(ax, name)
        self.btns.append(bx)
        bx.on_clicked(f)
    def makechk(self,x,y,name):
        ax = self.fig_controls.add_axes((x, y, 0.05, 0.05))
        self.axs.append(ax)
        bx = widgets.CheckButtons(ax, [name],[True])
        self.use_par[name]=bx
    def loadstars(self):
        # Load stars
        LimitMag = 6
        self.catalog=load_catalog()
        self.starnames=[]
        ras=[]
        decs=[]
        for i, this_star in enumerate(self.catalog):
            if i>4000:
                break
            star = " " + this_star
            if GetMag(star) > LimitMag:
                break
            decs.append(GetDec(star))
            ras.append(GetRA(star))
            self.starnames.append(f"{i:4d} "+GetName(star))
        #vector above is in J2000. Natural frame of Voyager animations is Ecliptic.
        #Since we are transforming anyway, rotate to B1950 as well.
        star_vecs=cspice.pxform("J2000","ECLIPB1950",0) @ llr2xyz(lat=np.array(decs),lon=np.array(ras),deg=True)
        hom_comps=np.zeros((1, len(self.starnames)))
        self.star_vec = np.vstack((star_vecs,hom_comps))
        self.nameobj=[None]*self.star_vec.shape[1]
        if True:
            for i in range(self.star_vec.shape[1]):
                if self.starnames[i]!="":
                    print(i, self.starnames[i])
                    self.nameobj[i]=self.ax.text(0,0,self.starnames[i],
                                             color=('yellow' if i in self.goodstars else
                                                    ('#808080' if i in self.badstars else 'blue')),
                                             visible=False)

        self.starplot, = self.ax.plot(np.zeros((self.star_vec.shape[1],)),
                                      np.zeros((self.star_vec.shape[1],)), 'r+')
        self.starplot.set_visible(False)
        q=np.arange(0,2*np.pi,0.01)
        self.c=np.cos(q)
        self.s=np.sin(q)
        self.diskplot={}
        self.orbplot={}
        self.ringplot=[]
        self.nametext={}
        self.ringcenter=list(self.spiceobjs.keys())[0]
        for k,body in self.spiceobjs.items():
            color={0:'#404040',1:'#804000',2:'#ff0000',3:'#ff8000',4:'#ffff00',5:'#00ff00',6:'#0000ff',7:'#8000ff',8:'#c0c0c0',9:'#ffffff'}
            self.diskplot[k],=self.ax.plot(self.c*0,self.s*0,'-',color=color[k%10])
            self.diskplot[k].set_visible(False)
            self.orbplot[k], = self.ax.plot(np.zeros(100), np.zeros(100),'-', color=color[k%10])
            self.orbplot[k].set_visible(False)
            self.nametext[k]=self.ax.text(0,0,body.name,color=color[k%10])
        for ring_r in self.rings:
            self.ringplot.append(self.ax.plot(self.c*0,self.s*0,'c-')[0])
            self.ringplot[len(self.ringplot)-1].set_visible(False)
        self.fitplot,  = self.ax.plot(np.zeros((self.star_vec.shape[1],)),
                                      np.zeros((self.star_vec.shape[1],)), 'g*')
        self.fitplot.set_visible(False)
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
               "lat_c           real," +
               "lat_c_sig       real," +
               "lat_c_source    integer," +
               "lon_c           real," +
               "lon_c_sig       real," +
               "lon_c_source    integer," +
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
        def fill_in_value(fieldname):
            # Sources are, in order of decreasing confidence:
            # 4 - constrained. Certain parameters like right_denom and probably clock are actually constant
            #     over the video. Also, the intent is that there is eventually a spline model for each of
            #     the viewpoint variables. "Constrained" means either constant or splined.
            # 3 - Fit via a least-squares model
            # 2 - Fit manually
            # 1 - Interpolated from higher-confidence sources
            # 0 - unknown source
            if self.__dict__[fieldname] is None or (fieldname=="et" and self.__dict__[fieldname+"_source"]<2):
                sql=f"select framenum,{fieldname} from frames where {fieldname}_source>1 order by abs(framenum-?) asc"
                print(sql)
                with closing(self.conn.cursor()) as cur:
                    this_has_row = False
                    for this_row in cur.execute(sql,(self.framenum,)):
                        if not this_has_row:
                            fn0,val0=this_row
                            this_has_row=True
                        else:
                            fn1,val1=this_row
                            break
                self.__dict__[fieldname]=linterp(fn0,val0,fn1,val1,self.framenum)
                self.__dict__[fieldname+"_source"]=1
        has_row=False
        #Check if this frame is already recorded
        sql=("select framenum,nstars,et,et_source,"+
                                       "lat_c       ,lon_c       ,angle       ,clock       ,right_denom,"+
                                       "lat_c_sig   ,lon_c_sig   ,angle_sig   ,clock_sig   ,right_denom_sig,"+
                                       "lat_c_source,lon_c_source,angle_source,clock_source,right_denom_source "+
                                       "from frames order by abs(framenum-?) asc")
        cur = self.conn.cursor()
        old_et = self.et
        self.lat_c=None
        self.lon_c=None
        self.angle=None
        self.clock=None
        self.right_denom=None
        self.et=None
        has_row=False
        for row in cur.execute(sql, (self.framenum,)):
            has_row=True
            if row[0]==self.framenum:
                (self.nstars, self.et,self.et_source,
                 self.lat_c       ,self.lon_c       ,self.angle       ,self.clock       ,self.right_denom       ,
                 self.lat_c_sig   ,self.lon_c_sig   ,self.angle_sig   ,self.clock_sig   ,self.right_denom_sig   ,
                 self.lat_c_source,self.lon_c_source,self.angle_source,self.clock_source,self.right_denom_source,)=row[1:]
            break
        if has_row:
            #There was at least one row, so we can interpolate
            fill_in_value("lat_c")
            fill_in_value("lon_c")
            fill_in_value("angle")
            fill_in_value("clock")
            fill_in_value("right_denom")
            fill_in_value("et")
            self.lat_c_sig = float('inf')
            self.lon_c_sig = float('inf')
            self.angle_sig = float('inf')
            self.clock_sig = float('inf')
            self.right_denom_sig = float('inf')
            self.nstars = None
        else:
            # No rows at all -- use initial conditions
            #initial conditions valid for frame 675
            if False:
                self.lat_c = 25.186          #Spherical coordinate latitude of camera position relative to its look point, in degrees
                self.lon_c = -35.169         #Spherical coordinate longitude of camera position, in degrees
                self.angle = 45.987          #Angle parameter of camera, equivalent to angle keyword in POV-Ray perspective camera
            else:
                self.lat_c = -2.4          #Spherical coordinate latitude of camera position relative to its look point, in degrees
                self.lon_c = 98.6-180         #Spherical coordinate longitude of camera position, in degrees
                self.angle = 45          #Angle parameter of camera, equivalent to angle keyword in POV-Ray perspective camera
            self.right_denom = 2.897004  #Camera right vector denominator -- right vector in POV-Ray perspective camera is right -x*4/right_denom
            self.clock=0.0
            self.lat_c_sig = float('inf')
            self.lon_c_sig = float('inf')
            self.angle_sig = float('inf')
            self.clock_sig = float('inf')
            self.right_denom_sig = float('inf')
            self.nstars=None
            if self.et is None:
                self.et=old_et
            if self.et is None:
                self.et=str2et("1986-01-25 04:00:00 UTC")-11*3600-20*60
    def write(self):
        sql=("insert or replace into frames (framenum,nstars,rmsdiff,et,et_source,"+
             "lat_c    ,lon_c    ,angle    ,clock    ,right_denom,   "+
             "lat_c_sig,lon_c_sig,angle_sig,clock_sig,right_denom_sig,"
             "lat_c_source,lon_c_source,angle_source,clock_source,right_denom_source) "+
             "values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)")
        cur = self.conn.cursor()
        cur.execute(sql, (self.framenum,
                          self.nstars,
                          self.rmsdiff,
                          self.et,
                          self.et_source,
                          self.lat_c,
                          self.lon_c,
                          self.angle,
                          self.clock,
                          self.right_denom,
                          self.lat_c_sig,
                          self.lon_c_sig,
                          self.angle_sig,
                          self.clock_sig,
                          self.right_denom_sig,
                          self.lat_c_source,
                          self.lon_c_source,
                          self.angle_source,
                          self.clock_source,
                          self.right_denom_source,
                          ))
        self.conn.commit()
    def project_stuff(self):
        dir=llr2xyz(lat=self.lat_c,lon=self.lon_c)
        sky= make_sky(self.clock, dir)
        self.C = cmatrix(loc=np.zeros((3, 1)), look=dir, sky=sky)
        self.star_pix = project(right=4 / self.right_denom, angle=self.angle, width=self.width,
                                height=self.height, target_c=self.C @ self.star_vec)

        print("image_width  ", self.width)
        print("image_height ", self.height)
        print("angle        ", self.angle)
        print("lat_c(deg)   ", self.lat_c)
        print("lon_c(deg)   ", self.lon_c)
        print("clock(deg)   ", self.clock)
        print("right    x*4/", self.right_denom)
        print("framenum     ", self.framenum)
        print("time         ", timout(self.et,"YYYY-MM-DD HR:MN:SC.###::UTC"))
    def load_image(self):
        infn = self.framepat % self.framenum
        self.backimg = mpimg.imread(infn)
        if self.figimg is not None:
            self.figimg.set_data(self.backimg)
        self.ax.set_title(f"Frame {self.framenum}")
    def refit(self):
        """
        For all frames which have an automatic fit solution, run the fit solution again

        """
        # Get a list of all the images to re-fit. These will be those with a finite sigma on lat.
        sql = "select framenum from frames where lat_c_sig<1;"
        # frames=cur.execute(sql).
        # for row in cur.execute(sql):
    def autofit(self,dframe,limit=8000):
        for i in range(limit):
            self.d_frame(dframe)
            self.fit()
    def plotspice(self):
        print(f"Voyager 2 in kernel {which_kernel('SPK',-32,self.et)}")
        for i_spice,body in self.spiceobjs.items():
            try:
                if body.dt is not None and body.dt>0:
                    rofs=body.dt*spkezr(str(i_spice),self.et,"ECLIPB1950","NONE",body.parent)[0][3:].reshape(-1,1)
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
                    xbody_vb,lt=spkezr(str(i_spice),self.et,bodyframe,"NONE","-32")
                except SpiceFRAMEDATANOTFOUND:
                    # Frame not found -- no frame is included for Nereid
                    bodyframe="ECLIPB1950"
                    xbody_vb,lt=spkezr(str(i_spice),self.et,bodyframe,"NONE","-32")
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
                M_ib=pxform(bodyframe,"ECLIPB1950",self.et)
                # limb points in Voyager-centered inertial frame
                rlimb_vi=M_ib @ rlimb_vb+rofs
                # Stack up into homogeneous vectors
                rlimb_vi=np.vstack((rlimb_vi,np.zeros((1,rlimb_vi.shape[1]))))
                pixlimb= project(right=4 / self.right_denom, angle=self.angle, width=self.width,
                                 height=self.height, target_c=self.C @ rlimb_vi)

                print(f"{i_spice} in kernel {which_kernel('SPK',i_spice,self.et)}")
                self.diskplot[i_spice].set_xdata(pixlimb[0,:])
                self.diskplot[i_spice].set_ydata(pixlimb[1,:])
                self.diskplot[i_spice].set_visible(True)
                if body.parent is not None:
                    et_orbit=(np.arange(100)-50)*60+self.et
                    v_orbit=np.ones((4,100))
                    rmoon_vi=spkezr(str(i_spice),self.et,"ECLIPB1950","NONE","-32")[0][0:3].reshape(-1,1)
                    for i_et,this_et in enumerate(et_orbit):
                        v_orbit[:3,i_et]=spkezr(str(i_spice),this_et,"ECLIPB1950","NONE",str(body.parent))[0][0:3]
                    v_orbit[0:3,:]-=v_orbit[0:3,None,50]
                    v_orbit[0:3,:]+=rmoon_vi
                    pixorbit= project(right=4 / self.right_denom, angle=self.angle, width=self.width,
                                      height=self.height, target_c=self.C @ v_orbit)
                    self.orbplot[i_spice].set_xdata(pixorbit[0,:])
                    self.orbplot[i_spice].set_ydata(pixorbit[1,:])
                    self.orbplot[i_spice].set_visible(True)
                    self.nametext[i_spice].set_x(pixorbit[0,50])
                    self.nametext[i_spice].set_y(pixorbit[1,50])
            except Exception:
                import traceback
                traceback.print_exc()
                print(f"{i_spice} not in kernels")
        for i_ring,ring_r in enumerate(self.rings):
            name="Uranus"
            i_spice=self.ringcenter
            bodyframe = f"IAU_{name.upper()}"
            xbody_vb, lt = spkezr(str(i_spice), self.et, bodyframe, "NONE", "-32")
            xvoy_bb = -xbody_vb
            rvoy_bb = xvoy_bb[:3]
            # ring points in body centered body frame
            rring_bb = self.c*np.array([[ring_r],[0.0],[0.0]]) + self.s*np.array([[0.0],[ring_r],[0]])
            rvoy_bb = rvoy_bb.reshape(-1, 1)
            # ring points in Voyager-centered body frame
            rring_vb = rring_bb - rvoy_bb
            M_ib = pxform(bodyframe, "ECLIPB1950", self.et)
            # ring points in Voyager-centered inertial frame
            rring_vi = M_ib @ rring_vb
            # Stack up into homogeneous vectors
            rring_vi = np.vstack((rring_vi, np.zeros((1, rring_vi.shape[1]))))
            pixring = project(right=4 / self.right_denom, angle=self.angle, width=self.width,
                              height=self.height, target_c=self.C @ rring_vi)
            self.ringplot[i_ring].set_xdata(pixring[0, :])
            self.ringplot[i_ring].set_ydata(pixring[1, :])
            self.ringplot[i_ring].set_visible(True)
    def replot(self):
        self.project_stuff()
        for i,n_o in enumerate(self.nameobj):
            if n_o is not None:
                n_o.set_visible(np.isfinite(self.star_pix[0,i]))
                if np.isfinite(self.star_pix[0,i]):
                    n_o.set_position((self.star_pix[0,i],self.star_pix[1,i]))
        self.starplot.set_xdata(self.star_pix[0,...])
        self.starplot.set_ydata(self.star_pix[1,...])
        self.starplot.set_visible(True)
        def set_btn(name,sig):
            if (self.use_par[name].get_status()[0] ^ np.isfinite(sig)):
                self.use_par["lat"].set_active(0)
        set_btn("lat",self.lat_c_sig)
        set_btn("lon",self.lon_c_sig)
        set_btn("angle",self.angle_sig)
        set_btn("clock",self.clock_sig)
        set_btn("right",self.right_denom_sig)
        self.plotspice()
        self.ax.set_xlabel(timout(self.et,"YYYY-MM-DD HR:MN:SC.###::UTC"))
        try:
            plt.pause(0.001)
        except Exception:
            pass
    def fit(self):
        """
        Given the current position as an initial guess, find the optimum
        camera parameters and position to fit the stars.
        """
        self.ax.set_ylabel("Fitting...")
        plt.pause(0.001)
        done=False
        #All of the following arrays will be edited down as we edit the data
        #Array of star indices for stars under consideration
        this_goodstars=np.array([False]*self.star_pix.shape[1])
        for i_star in self.goodstars:
            if i_star<len(this_goodstars):
                this_goodstars[i_star]=True
        while not done:
            print("Number of stars on-screen:      ",np.sum(np.isfinite(self.star_pix[0,:])))
            goodstar_pix=self.star_pix.copy()
            goodstar_pix[:,np.logical_not(this_goodstars)]=np.array([[float('nan')],[float('nan')]])
            print("Number of good stars on-screen: ",np.sum(np.logical_and(this_goodstars,np.isfinite(goodstar_pix[0,:]))))
            (findx,findy),(sigx,sigy),rho=find_stars(self.backimg,goodstar_pix,names=self.starnames,ax=None,boxr=10)#self.ax_controls)
            print("Number of good stars found:     ",np.sum(np.isfinite(findx)))
            self.fitplot.set_visible(True)
            self.fitplot.set_xdata(findx)
            self.fitplot.set_ydata(findy)
            plt.pause(0.001)
            w=np.where(np.isfinite(findx))
            names=[name for i_name,name in enumerate(self.starnames) if np.isfinite(findx[i_name])]
            indices=[i_name for i_name,name in enumerate(self.starnames) if np.isfinite(findx[i_name])]
            findx=findx[w]
            findy=findy[w]
            sigx=sigx[w]
            sigy=sigy[w]
            rho=rho[w]
            cov=np.zeros((len(findx)*2,len(findx)*2))
            for i in range(len(findx)):
                cov[i*2  ,i*2  ]=sigx[i]**2
                cov[i*2+1,i*2+1]=sigy[i]**2
                cov[i*2+1,i*2  ]=sigx[i]*sigy[i]*rho[i]
                cov[i*2  ,i*2+1]=sigx[i]*sigy[i]*rho[i]
            weights=1.0/np.sqrt(sigx**2+sigy**2)
            # make an array [findx
            #                findy] then ravel it. The result is [findx|findy]
            pixdata=np.vstack((findx,findy)).ravel()
            fitv=self.star_vec[:,w[0]] #If we do fitv[:,w] we get a shape (4,1,nstars) instead of the (4,nstars) we want
            #fiti=np.array(goodstars)[w]

            p0 = np.array((self.lat_c                    ,self.lon_c                       ,self.angle                  ,self.clock                      ,self.right_denom))
            vary=[         bounded(-90.0,90.0,self.lat_c),rbounded(-180.0,180.0,self.lon_c),bounded(0.0,120.0,self.angle),bounded(-180.0,180.0,self.clock),positive()       ]
            if len(findx)<2:
                p0[2]=False
                p0[3]=False
                p0[4]=False
                print("Few usable stars, only fitting pointing")
            elif len(findx)<5:
                p0[3]=False
                p0[4]=False
                print("Few usable stars, only fitting pointing and angle")
            for v,fieldname in zip(vary,["lat_c_source","lon_c_source","angle_source","clock_source","right_denom_source"]):
                if v:
                    self.__dict__[fieldname]=3
            (popt,pcov)=curve_fit(curve_fitsky_interface,fitv,pixdata,p0=p0,vary=vary,sigma=cov,absolute_sigma=True,f_kwargs={'width':self.width,'height':self.height})
            fit_pixdata=curve_fitsky_interface(fitv,*popt,width=self.width,height=self.height)
            fitx=fit_pixdata[:len(fit_pixdata)//2]
            fity=fit_pixdata[len(fit_pixdata)//2:]
            total_lensq=0
            for name, ox, oy, cx, cy in zip(names, findx, findy,fitx,fity):
                lensq=(ox-cx)**2+(oy-cy)**2
                total_lensq+=lensq
                print(f"{name},{ox:8.3f},{oy:8.3f},{cx:8.3f},{cy:8.3f},{np.sqrt(lensq):8.3f}")
            self.rmsdiff=np.sqrt(total_lensq/len(name))
            print(f"RMS diff: {self.rmsdiff:8.5f}")
            while popt[1]>360:
                popt[1]-=360
            while popt[1]<0:
                popt[1]+=360
            cor = correlation_matrix(pcov)
            (self.lat_c, self.lon_c, self.angle,self.clock,self.right_denom) = popt
            self.lat_c_sig,self.lon_c_sig,self.angle_sig,self.clock_sig,self.right_denom_sig=tuple([cor[i,i] for i in range(len(popt))])
            infam=infamily(np.vstack((fitx,fity)),np.vstack((findx,findy)),weights=weights)
            if len(infam)>5:
                for i_infam in range(len(infam)):
                    if not infam[i_infam]:
                        print(f"Star {names[i_infam]} not in family")
                        this_goodstars[indices[i_infam]]=False
                done=np.all(infam)
            else:
                done=True
            self.replot()
        self.nstars = len(infam)
        self.write()
        self.ax.set_ylabel("")
        plt.pause(0.001)
    def set_frame(self,i_frame):
        self.write()
        self.framenum=i_frame
        self.read_record()
        self.load_image()
        if self.fitplot is not None:
            self.fitplot.set_visible(False)
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
    def BTNrefit(self, event):
        self.refit()
    def BTNcamlonp(self, event):
        self.lon_c_source=2
        self.lon_c+=self.step_size
        self.replot()
    def BTNcamlonm(self, event):
        self.lon_c_source=2
        self.lon_c-=self.step_size
        self.replot()
    def BTNcamlatp(self, event):
        self.lat_c_source=2
        self.lat_c+=self.step_size
        self.replot()
    def BTNcamlatm(self, event):
        self.lat_c_source=2
        self.lat_c-=self.step_size
        self.replot()
    def BTNanglep(self, event):
        self.angle_source=2
        self.angle += self.step_size
        self.replot()
    def BTNanglem(self, event):
        self.angle_source=2
        self.angle -= self.step_size
        self.replot()
    def BTNrightp(self, event):
        self.right_denom_source=2
        self.right_denom += self.step_size
        self.replot()
    def BTNrightm(self, event):
        self.right_denom_source=2
        self.right_denom -= self.step_size
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
        self.clock_source=2
        self.clock += self.step_size
        self.replot()
    def BTNclockm(self, event):
        self.clock_source=2
        self.clock -= self.step_size
        self.replot()
    def BTNtimep(self, event):
        self.et_source=2
        self.et += self.step_size*60
        self.replot()
    def BTNtimem(self, event):
        self.et_source=2
        self.et -= self.step_size*60
        self.replot()
    def BTNtimeconf(self, event):
        self.et_source=2
        self.replot()
        self.write()
    def BTNframem(self, event):
        self.d_frame(-1)
    def BTNframep(self, event):
        self.d_frame(+1)
    def BTNfit(self, event):
        self.fit()
    def autop(self, event):
        self.autofit(+1)
    def autom(self, event):
        self.autofit(-1)


