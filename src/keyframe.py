"""
Generate a table of function values given a descrption of the
function in a piecewise fashion. Effectively each knot between
the pieces describes a keyframe, and the transition between keyframes
is attached to the left knot

Created: 2/13/25
"""
import sqlite3
from contextlib import closing
from typing import Iterable

import numpy as np
from kwanmath.interp import linterp
from matplotlib import pyplot as plt


class Knot:
    t0:float
    y0:float
    def __init__(self,*,t0:float,y0:float|np.ndarray):
        """

        :param t0: Independent variable at this knot
        :param y0: Dependent variable at this knot,
                   may be either a scalar (float) or column vector (N,1 array)
        """
        self.t0=t0
        self.y0=y0
    def eval(self,t:float|np.ndarray,next_knot:'Knot')->float|np.ndarray:
        """
        Evaluate the segment between the given knot and the next knot
        :param t: Independent variable at which to evaluate. Should
                  be either a scalar (float), a 1D (M,) array or a (1,M) row vector
        :param next_knot: Next knot in the sequence, used to define the
                          segment end point
        :return: Segment evaluated at each independent variable. Can be:
          * If t is scalar, result is scalar if y0 is scalar, or Nx1 column vector if y0 is a vector
          * if t is a row vector or 1D vector, result is a row vector if y0 is scalar or an N,M set of
            column vectors.
        """
        raise NotImplementedError


class ConstantSegment(Knot):
    def eval(self,t:float|np.ndarray,next_knot:Knot)->float|np.ndarray:
        # Do it this way so that it works whether:
        #   * t is a scalar, we just get y0*1=y0, either column vector or scalar
        #   * t is an (M,) 1D array. We get an (N,1) y0 broadcast with a (M,) t
        #     so the result is (N,M)
        #   * t is an (1,M) 2D row vector. We still get a result of (N,M)
        return self.y0*(0*t+1)


class LinearSegment(Knot):
    def eval(self,t:float|np.ndarray,next_knot:Knot)->float|np.ndarray:
        u=(t-self.t0)/(next_knot.t0-self.t0)
        return self.y0*(1-u)+next_knot.y0*u


class QuadraticEaseLeft(Knot):
    def __init__(self,*,t0:float,y0:float|np.ndarray,dt:float):
        """

        :param t0: Independent knot coordinate
        :param y0: Dependent knot coordinate (vector or scalar)
        :param dt: Distance of ease
        """
        super().__init__(t0=t0,y0=y0)
        self.dt=dt
    def eval(self,t:float|np.ndarray,next_knot:Knot)->float|np.ndarray:
        y0=self.y0
        t0=self.t0
        y1=next_knot.y0
        t1=next_knot.t0
        yp1=y1-y0
        tp1=t1-t0
        a=yp1/(2*tp1*self.dt-self.dt**2)
        tp=t-t0
        # Coordinates of the green point
        tgp=self.dt
        ygp=a*self.dt**2
        # Do the curve to the left of the green point and a line
        # to the right of the green point
        yp=np.where(tp<=tgp,a*tp**2,linterp(tgp,ygp,tp1,yp1,tp))
        y=yp+y0
        return y


class QuadraticEaseRight(Knot):
    def __init__(self,*,t0:float,y0:float|np.ndarray,dt:float):
        """

        :param t0: Independent knot coordinate
        :param y0: Dependent knot coordinate (vector or scalar)
        :param dt: Distance of ease
        """
        super().__init__(t0=t0,y0=y0)
        self.dt=dt
    def eval(self,t:float|np.ndarray,next_knot:Knot)->float|np.ndarray:
        # There's lots of derivation and lots of signs, so we will make
        # it super-explicit to help minimize that.
        y1=next_knot.y0
        t1=next_knot.t0
        y0=self.y0
        t0=self.t0
        yp0=y0-y1
        tp0=t0-t1
        a=-yp0/(self.dt**2+2*tp0*self.dt)
        tp=t-t1
        # Coordinates of the green point
        tgp=-self.dt
        ygp=a*self.dt**2
        # Do the curve to the right of the green point and a line
        # to the left of the green point
        yp=np.where(tp>=tgp,a*tp**2,linterp(tp0,yp0,tgp,ygp,tp))
        y=yp+y1
        return y

class HermiteEase(Knot):
    """
    Do easing with a
    """
    def __init__(self,*,t0:float,y0:float|np.ndarray,yd0:float|np.ndarray,yd1:float|np.ndarray):
        super().__init__(t0=t0,y0=y0)
        self.yd0=yd0
        self.yd1=yd1
    def eval(self,t:float|np.ndarray,next_knot:Knot)->float|np.ndarray:
        t0=self.t0
        y0=self.y0
        t1=next_knot.t0
        y1=next_knot.y0
        yd0=self.yd0
        yd1=self.yd1
        b=np.array([[y0],
                    [yd0],
                    [y1],
                    [yd1]])
        A=np.array([[  t0**3,  t0**2,t0,1],
                    [3*t0**2,2*t0   , 1,0],
                    [  t1**3,  t1**2,t1,1],
                    [3*t1**2,2*t1   , 1,0]])
        x=np.linalg.solve(A,b)
        a,b,c,d=x
        # Check that these coefficients actually work
        assert np.isclose(  a*t0**3+  b*t0**2+c*t0+d,y0)
        assert np.isclose(3*a*t0**2+2*b*t0   +c     ,yd0)
        assert np.isclose(  a*t1**3+  b*t1**2+c*t1+d,y1)
        assert np.isclose(3*a*t1**2+2*b*t1   +c     ,yd1)
        # Now we can just evaluate it
        return a*t**3+b*t**2+c*t+d


class CubicFitEase(Knot):
    """
    Do easing with a
    """
    def __init__(self,*,ts:np.ndarray,ys:np.ndarray):
        # Do the fit now
        ts=ts.reshape(-1,1)
        A=np.hstack((ts**3,ts**2,ts,np.ones_like(ts)))
        b=ys.reshape(-1,1)
        x,*_=np.linalg.lstsq(A,b)
        self.a,self.b,self.c,self.d=x
        ycs=self.a*ts**3+self.b*ts**2+self.c*ts+self.d
        plt.plot(ts,ycs,'-')
        plt.plot(ts,ys,'.')
        plt.show()
        self.t0=ts[0]
        self.y0=ycs[0]
    def eval(self,t:float|np.ndarray,next_knot:Knot)->float|np.ndarray:
        t0=self.t0
        y0=self.y0
        t1=next_knot.t0
        y1=next_knot.y0
        yd0=self.yd0
        yd1=self.yd1
        b=np.array([[y0],
                    [yd0],
                    [y1],
                    [yd1]])
        A=np.array([[  t0**3,  t0**2,t0,1],
                    [3*t0**2,2*t0   , 1,0],
                    [  t1**3,  t1**2,t1,1],
                    [3*t1**2,2*t1   , 1,0]])
        x=np.linalg.solve(A,b)
        a,b,c,d=x
        # Check that these coefficients actually work
        assert np.isclose(  a*t0**3+  b*t0**2+c*t0+d,y0)
        assert np.isclose(3*a*t0**2+2*b*t0   +c     ,yd0)
        assert np.isclose(  a*t1**3+  b*t1**2+c*t1+d,y1)
        assert np.isclose(3*a*t1**2+2*b*t1   +c     ,yd1)
        # Now we can just evaluate it
        return a*t**3+b*t**2+c*t+d


def hp00(tp):
    """
    Hermite basis function associated with the zero derivative of the point at t=t0
    :param tp: Normalized time coordinate t-prime
    """
    return 2*tp**3-3*tp**2+1
def hp10(tp):
    """
    Hermite basis function associated with the first derivative of the point at t=t0
    :param tp: Normalized time coordinate t-prime
    """
    return tp**3-2*tp**2+tp
def hp01(tp):
    """
    Hermite basis function associated with the zero derivative of the point at t=t1
    :param tp: Normalized time coordinate t-prime
    """
    return -2*tp**3+3*tp**2
def hp11(tp):
    """
    Hermite basis function associated with the first derivative of the point at t=t0
    :param tp: Normalized time coordinate t-prime
    """
    return tp**3-tp**2


def hermitep(y0:float,ydp0:float,y1:float,ydp1:float,tp:float|np.ndarray)->float|np.ndarray:
    """
    Run the Hermite interpolation using normalized time (t-prime) which runs from 0 to 1
    :param y0:
    :param ydp0:
    :param y1:
    :param ydp1:
    :param tp:
    :return:
    """
    hp00y0=hp00(tp)*y0
    hp10ydp0=hp10(tp)*ydp0
    hp01y1=hp01(tp)*y1
    hp11ydp1=hp11(tp)*ydp1
    return hp00y0+hp10ydp0+hp01y1+hp11ydp1



class HermiteFitEase(Knot):
    """
    Do easing with a Hermite polynomial fit to data. We take a definite
    t0,y0,t1,y1 and fit the slope yd0 and yd1.
    """
    def __init__(self,*,y0:float|np.ndarray,ts:np.ndarray,ys:np.ndarray):
        self.ts=ts.copy().reshape(-1)
        self.ys=ys.copy().reshape(-1)
        self.t0=ts[0]
        self.y0=y0
    def eval(self,t:float|np.ndarray,next_knot:Knot)->float|np.ndarray:
        t0=self.t0
        y0=self.y0
        t1=next_knot.t0
        y1=next_knot.y0
        ydp1=0.0
        tpdata=linterp(t0,0,t1,1,self.ts)
        # General linear least squares. Unknowns are slopes ydp0 and ydp1
        # equations are
        A=np.stack((hp10(tpdata),)).T
        b=(self.ys-hp00(tpdata)*y0-hp01(tpdata)*y1-hp11(tpdata)*ydp1).reshape(-1,1)
        x,*_=np.linalg.lstsq(A,b)
        ydp0,=x
        # Switch from
        w=np.logical_and(t0<=t,t<=t1)
        tp=linterp(t0,0,t1,1,t)
        ycs=hp00(tp)*y0+hp10(tp)*ydp0+hp01(tp)*y1+hp11(tp)*ydp1
        plt.plot(t[w],ycs[w],'-')
        plt.plot(self.ts,self.ys,'.')
        plt.show()
        return ycs


class Keyframe:
    def __init__(self,knots:list[Knot]):
        self.knots=knots
    def eval(self,t:float|np.ndarray):
        t=np.array(t)
        y=np.zeros_like(t)
        for this_knot,next_knot in zip(self.knots[:-1],self.knots[1:]):
            t0=this_knot.t0
            t1=next_knot.t0
            y=np.where(np.logical_and(t0<=t,t<=t1),this_knot.eval(t,next_knot),y)
        return y


def exercise_hermite():
    tp=np.arange(0.0,1.0,0.01)
    y=hermitep(y0=0,ydp0=0,y1=1,ydp1=0,tp=tp)
    plt.plot(tp,y)
    plt.show()


def main():
    exercise_hermite()
    with closing(sqlite3.connect("data/db/frame_index_VoyagerUranusHD.sqlite")) as conn, closing(conn.cursor()) as cur:
        sql="select framenum,angle from frames order by framenum asc;"
        result=cur.execute(sql).fetchall()
    result=np.array(result).T
    tdata=result[0,:]
    ydata=result[1,:]
    keyframe=Keyframe([
      # Frame Angle Ease Body cx cy scr scx scy
      ConstantSegment(t0=0,y0=97.8), # Show all orbits
      #HermiteEase(t0=100,y0=97.8,yd0=-0.1,yd1=0), # Start zoom in to standard zoom
      HermiteFitEase(y0=97.8,ts=tdata[np.logical_and(tdata>=100,tdata<=220)],
                   ys=ydata[np.logical_and(tdata>=100,tdata<=220)]),
      ConstantSegment(t0=220, y0=59.0), # Finish zoom to standard
      QuadraticEaseLeft(t0=401,y0=59.0,dt=5), # Start zoom into SigSag
      QuadraticEaseRight(t0=420,y0=40.0,dt=150), # Start zoom into SigSag
      ConstantSegment(t0=520,y0=6.8), # Finish zoom in
      QuadraticEaseLeft(t0=910,y0=6.8,dt=150), # Start zoom out
      ConstantSegment(t0=1033,y0=59.0), # Finish zoom out
      #{1292, 59.0,   0}, # Start zoom into 704, VOBEST
      #{1380,  2.3, -10}, # Finish zoom in
      #{1555,  2.3,  10}, # Start zoom out
      #{1622, 59.0,   0}, # Finish zoom out
      #{2240, 59.0,   0}, # Start zoom into 702, VUBEST
      #{2325,  2.3,   0}, # Finish zoom in
      #{2480,  2.3,  20}, # Start zoom out
      #{2570, 59.0,   0}, # Finish zoom out
      #{2988, 59.0,   0}, # Start zoom into 703, VTBEST
      #{3075,  2.3, -10}, # Finish zoom in
      #{3230,  2.3,   0}, # Start zoom out
      #{3318, 59.0,   0}, # Finish zoom out
      #{3740, 59.0,   0}, # Start zoom into 701, VABEST
      #{3830,  2.3, -15}, # Finish zoom in
      #{4040,  2.3,  10}, # Start zoom out
      #{4100, 59.0,   0}, # Finish zoom out
      #{4190, 59.0,   0}, # Start zoom into 705, VMBEST
      #{4240,  2.3, -15}, # Finish zoom in. Note that distance to Miranda changes significantly over sequence
      #{4578,  2.3,  15}, # Start zoom out
      #{4682, 59.0,   0}, # Finish zoom out */
      #{4790, 59.0,   0,799,686,426} #
      #{6906, 59.0,   0,799,}  # Constant to the end of the video
    ])
    dt=1
    t=np.arange(1,6906+1,dt)
    y=keyframe.eval(t)
    plt.figure("fov")
    plt.plot(t,y,'-')
    plt.plot(tdata,ydata,'.')
    dt=tdata-np.roll(tdata,1)
    dt[0]=dt[1]
    dy=ydata-np.roll(ydata,1)
    dy[0]=dy[1]
    dydt=dy/dt
    plt.figure("slope")
    plt.plot(tdata,dydt,'.')
    plt.show()


if __name__ == "__main__":
    main()
