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


class EllFit(object):
    q:np.ndarray=np.arange(0,np.pi*2,0.01)
    c:np.ndarray=np.cos(q)
    s:np.ndarray=np.sin(q)
    # Resistor color code
    colors = {0: '#404040', 1: '#804000', 2: '#ff0000', 3: '#ff8000', 4: '#ffff00', 5: '#00ff00', 6: '#0000ff',
             7: '#8000ff', 8: '#c0c0c0', 9: '#ffffff'}

    def __init__(self,*,casename:int,initframe:int,spice_id:int):
        #initialize from parameters
        self.casename=casename
        self.initframe=initframe
        self.framenum=self.initframe
        self.spice_id=spice_id
        self.framepat=f"data/frames/{casename}/frame%04d.png"

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
    def makebtn(self,x,y,name,f):
        ax = self.fig_controls.add_axes((x, y, 0.05, 0.05))
        self.control_axs.append(ax)
        bx = widgets.Button(ax, name)
        self.btns.append(bx)
        bx.on_clicked(f)
    def load_image(self):
        infn = self.framepat % self.framenum
        self.backimg = mpimg.imread(infn)
        self.edge=cv2.Canny((self.backimg[:,:,2]*255.0).astype(np.uint8), 50, 150)
        self.backimg+=self.edge[:,:,np.newaxis]/255.0
    def imshow(self):
        self.ax.clear()
        self.ax.imshow(self.backimg)
        self.ax.axis('scaled')
        self.clicks,*_=self.ax.plot([],[],'g+')
        self.y0plot,*_=self.ax.plot([],[],'b-')
        self.y1plot,*_=self.ax.plot([],[],'r-')
        self.center,*_=self.ax.plot([],[],'k*')
        self.ax.set_title(f"Frame {self.framenum}")
        plt.pause(0.001)
    def set_frame(self,framenum):
        self.framenum=framenum
        self.load_image()
        self.imshow()
    def d_frame(self,d_frame):
        self.set_frame(self.framenum+d_frame)
    def update_click(self):
        this_xclicks=np.array(self.xclicks)
        this_yclicks=np.array(self.yclicks)
        self.clicks.set_data(this_xclicks,this_yclicks)
        if len(this_xclicks)>=5:
            Ap,Bp,Cp,Dp,Ep=fit_conic(np.row_stack((this_xclicks,this_yclicks)),scale=1000.0)
            r=eval2_conic(Ap,Bp,Cp,Dp,Ep,scale=1000.0)
            self.y0plot.set_data(r[0,:],r[1,:])
            cv,av,bv=identify_conic(Ap,Bp,Cp,Dp,Ep,scale=1000.0)
            cv=cv.reshape(-1)
            axes=np.column_stack((cv+av,cv,cv+bv))
            self.y1plot.set_data(axes[0,:],axes[1,:])
        plt.pause(0.001)
    # Callbacks only below this point. Callbacks should be trivial.
    # If they are more than about 1 line, or if they are called from
    # somewhere else, break out the functionality into its own function
    # and call it from the callback.
    def BTNframem(self, event):
        self.d_frame(-1)
    def BTNframep(self, event):
        self.d_frame(+1)
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


