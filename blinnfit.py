from spiceypy import furnsh, spkezr, str2et, timout, sxform, edlimb, pxform
import spiceypy as cspice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.widgets as widgets
import os
import sqlite3

from kwanmath.gaussian import correlation_matrix, infamily
from kwanmath.geodesy import llr2xyz
from kwanmath.vector import vlength, vnormalize, vcross
from kwanmath.optimize import curve_fit, bounded, positive, rbounded

from bsc import load_catalog, GetMag, LimitMag, GetDec, GetRA, GetName

from kwanmath.interp import linterp

from find_stars import find_stars
from which_kernel import which_kernel

cspice.furnsh('data/spice/vgr1.tm')
cspice.furnsh('data/spice/vgr2.tm')
cspice.furnsh('data/spice/lsk/naif0012.tls')
cspice.furnsh('data/spice/spk/de430.bsp')

goodstarsSuperTraj=[
           1553,1097, 191, 765, 613, 289, 456, 472, 968,1000,
           1292,1443, 411,1400,1243, 979, 801,1013,1492, 791,
            741,1149,1020,1496,1399, 724,1569,1148,1108, 467,
           1475, 887,  47, 666, 133,1414, 241, 909, 845,1109,
           1021, 501,1045, 436,  88,1323, 848,  65, 277,  22,
            980,1611, 379, 495, 755,1585, 931, 300,1585,1371,
           1291, 205, 278,  92, 197, 886, 565, 357, 159,1491,
            494,  37, 391, 515,   1, 175,1290,1320,1035,1199,
            237     , 790,1060, 908, 227, 510, 877, 753,1136,
            433,1269, 750,1147, 489, 516, 386,1474,1318, 226,

            580, 397, 661, 850,1146, 432, 220,1058, 308, 567,
            879,1497,1416,1324,1002,1124,1159, 301, 172,1444,
            530, 172,  90, 694,1125, 107, 821, 517,1557,1526,
            122, 458, 682,  94,1189, 811, 744, 263, 737,  98,
           1046,1446, 173,1139, 228,  57,1062, 727, 111,1229,
            216, 128, 318, 491, 552,1007, 563, 454, 252, 346,
            949, 309, 460,1064, 648,1039,1595, 871, 614, 649,
            387,  75, 200, 240, 142, 243, 334, 163,1387, 555,
           1217,1559,  24, 115, 706,  39, 992, 347,  83,1376,

            214, 193, 594, 972, 951, 924, 208, 239,  36,1510,
            584, 505     , 595, 194, 983, 784,1421,  13,1390,
           1378, 586, 603,1092,1308,1222,1434,1142, 177, 860,
            618, 633, 862,1360, 192, 566,1589, 591, 738, 100,
            154,1602, 955, 135, 392, 900, 564, 118, 376, 746,
            147, 711, 686,  16,  30,1433,  60,  43, 796, 815,
            785, 718,  52, 285, 320, 230, 976, 835,1513, 960,
             81,1267,1145,1181, 284,1118, 305,1194,1422,1350,
            513, 689, 461, 182, 161, 448, 543, 872, 470, 678,
            479, 371, 660,1183, 275, 557,  93, 700, 375, 732,

            722, 855,1071, 646,1604,1423, 907, 576, 307, 977,
           1184, 677, 805, 953, 306, 659, 843,1317, 106, 143,
            714, 840, 328,1001,1214,1495,1260,1614,1242,1173]


goodstarsVoyagerUranus={
    1015, 481,1545,1127, 675, 302, 739,  85, 244, 444,
    1154, 154, 392, 288,1333,1128, 924, 983,1008,1378,
    1359,  51, 208, 239, 100, 135,1067, 972, 118,1112,
    1280,1573,1206,1602,1222, 584, 585, 994, 505,1390,
     214,  24, 115, 706,  83,1279, 181,1404,1387     ,
    1392,  36, 193,1329, 168, 951, 229, 617,1278,1389,
     793,  39, 190,1530,1026,  75,1065,1207,1434, 784,
    1092, 586, 603, 136, 240, 595,1542, 271, 533,1233,
    1433,1508, 871,1040, 950, 224, 774,1420,1453,1308,
     992, 649, 648, 935, 992, 194,1142, 900,1348,1421,

    1391, 824, 730,1164, 395, 564,1376,1192,1544, 287,
     440,  96, 695, 834, 899, 851,  13, 991,1230, 460,
     221,1218, 912,1510, 200,1047,1506, 177,1687,1923,
         2155,1785,2834,2628,2157,2314,2212,2313,2977,
    2445,2950,2883,1852, 614,2623,1871,2574,1417, 264,
    2389,2776,1595,2439,2388,1922,1153, 128, 738,2149,
     252,2190,2745, 594, 347,1219,1874, 163,2128,1948,
    2189, 142, 468,1083,1824,2836,1688,2553,2080,2064,
    2412,2555,1879,2972,1924,1823,2552,2944,1099,2837,
    1664,1745,2782,2911,2881,2882,2261,2411,2410,1898,

    1822,1966,2652,2946,2945,2654,2409,2521,2971,1764,
    2441,2970,1684,2154,2522,1482,2678,2063,1290,2915,
     823, 522, 498,  79,1683, 152, 460, 823,2775,2813,
    1246,1661, 957, 153,1896,2519,1619, 346, 309, 188,
    2362,2393,2265, 656,1466,2835,2487,1052,2779,2744,
     333,2711,2101,2812,1432,2153,1997,1967,1726,2444,
    2391, 243, 555,1288,2312,2784,2193, 387, 684,1299,
     491,1191,1090, 216, 531,1597,1783,2711, 454,2151,
     948, 949,1963,1152, 263,1007, 199, 111, 552, 334,
     550, 822, 563, 910, 744,1051, 304, 756, 429,2594,

      78, 355,1804,2079, 743, 453,  94,  77,2594,1450,
    2672,2360, 517,1682, 821, 548, 811,   3, 439, 583,
    1375,  10,  76, 551,1462,  57,1062, 757, 335,1897,
    2233, 141, 767, 262, 183,1941, 735,1461,1189, 736,
    1459,1993,1161, 964, 622, 647, 683,1139,  63,1080,
    1140,1061,1004,  99, 442,  41,1231,1451, 151, 211,
    2941,1039,2062,  17, 737,1761,1962, 474,  17, 933,
     122, 458,1356,1803, 682,  23, 934,2058, 727,1216,
    1740, 420,

}

#Stars which are visible but shouldn't be used (for instance too close to other stars)
# For instance -- if two stars are within about 10 pixels of each other, they will inevitably
# be in each other's boxes. If one is much brighter than the other, the dimmer star will still
# get the position of the brighter one, with a good fit measurement. In this case, put the
# brighter one (with the lower index) in the list above, and the dimmer one in the list below
# so that we know that it was excluded on purpose, rather than neglected.
badstarsVoyagerUranus={
    729,313,315,1162,1298,965,20,870,2238
}

goodstarsVoyagerNeptune=goodstarsVoyagerUranus|{
   2105,2021,1826,2913,2488, 157, 914,1115,1807,2841,
   2266, 696, 147, 915, 376,2159,2043, 354,1166,2602,
   1102,1647, 746,1168, 901,2023,2656, 686,2749,1786,
    711,2657,1901,2364, 618,1084, 431, 164,2316,1117,
    825,1691, 615,2494,1309,2414, 526, 230,1118,2213,
   2066, 872,  16,2918,2629,1468, 377, 642,1747, 614,
    797,1576,1250,2820, 525,2660, 348,1119, 507, 798,
   1144,1310,2415,2162,1790,2067, 959,1791,1973,1436,
    926,1882,1974, 937, 258,1730, 698, 598,1512, 902,
   1709, 527, 535, 358,2578,1379,1010, 816,1857, 305,

   1513,1628, 891,2267,1408, 960,2025,2296, 835,1792,
   1252,2497,1731, 137,2887,2318,1748,2888,2498, 826,
    687,1211,2039,2633,2164,2633,1194,2027,2319,1884,
    272, 320, 626, 259, 297,1349,2495,2559,2108,1195,
   1054, 553,2751, 817,2448,1749,1069,2135,2499,2851,
    330,2134, 668,1335,1955, 836,1017,2634,1579,2397,
     91, 261,1565,1236,  48,1668, 105,2215,1808,2752,
    688, 778,1365,2166, 275,2245,2609,2846,1831, 170,
   1930, 518,1380,2244,2852,1487,2269,1669, 905,1070,
    337,2299,2635,2006, 676,1121,  93,1030,2849, 777,

   1954,2606, 904,2085,2822,2848,3497,3242,3610,3743,
   3519,3883,3344,3208,3380,1851,3421,3660,3988,3833,
   3449,3931,3562,3700,3378,3207,3162,3417,3559,3013,
   3487,3241,3048,3932,2466,3012,3608,3239,3785,3452,

}

badstarsVoyagerNeptune=badstarsVoyagerUranus|{
   1331,1598,1603,2492, 997,1031,3879,3933,1765,3276
}

goodstarsVoyagerUranus-=badstarsVoyagerUranus
goodstarsVoyagerUranus=list(goodstarsVoyagerUranus)

goodstarsVoyagerNeptune-=badstarsVoyagerNeptune
goodstarsVoyagerNeptune=list(goodstarsVoyagerNeptune)


projects={
    "SuperTrajectory":(675,800,"data/frames/SuperTrajectory/frame%04d.png",goodstarsSuperTraj),
    "VoyagerUranusQuick":(201,533,'data/frames/VoyagerUranus/frame%04d.png',None,),
    "VoyagerUranusDetailed": (891, 533, 'data/frames/VoyagerUranus/frame%04d.png',goodstarsVoyagerUranus,badstarsVoyagerUranus,{
        #Spice ID
        #     Name     Pre-encounter radius (km)
        #                    a        b        c (polar) from pck00010.tpc
        799:("Uranus",25550,25559  ,25559  , 24973),
        701:("Ariel",   665,  581.1,  577.9,   577.7),
        702:("Umbriel", 555,  584.7,  584.7,   584.7),
        703:("Titania", 800,  788.9,  788.9,   788.9),
        704:("Oberon",  815,  761.4,  761.4,   761.4),
        705:("Miranda", 250,  240.4,  234.2,   232.9)
    }),
    "VoyagerNeptune": (1758, 533, 'data/frames/VoyagerNeptune/frame%04d.png', goodstarsVoyagerNeptune, badstarsVoyagerNeptune,{
        #Spice ID
        #     Name     Pre-encounter radius (km)
        #                    a        b        c (polar) from pck00010.tpc
        899:("Neptune",24781,24764  , 24764  ,24341),
        801:("Triton" , 1500, 1352.6,  1352.6,  1352.6),
        802:("Neried" ,  400,  170  ,   170  ,   170), #Radius from Neptune Travel Guide p146
        803:("Naiad N6",   0,   29  ,    29  ,    29),
        804:("Thalassa N5",0,   40  ,    40  ,    40),
        805:("Despina N3" ,0,   74  ,    74  ,    74),
        806:("Galatea N4" ,0,   79  ,    79  ,    79),
        807:("Larissa N2" ,0,  104  ,   104  ,    89),
        808:("Proteus N1" ,0,  218  ,   208  ,   201),
    })
}

def cmatrix(loc=None,look=None,sky=None):
    """

    :param loc: Location of camera, equivalent to camera{location...}
    :param look: Look-at point of camera, equivalent to camera{look_at...}
    :param sky: Sky vector of camera, equivalent to camera{sky...}
    :return: Camera matrix which transforms a global vector to a vector in camera space
    We will return a 4x4 matrix which will transform a vector in homogeneous coordinates into
    the camera frame. This matrix is:
     [a b c xt]
     [d e f yt]
     [g h i zt]
     [0 0 0  1]
    Multiply this matrix by a 4-element column vector, with the last element being 1
    if the point is a finite distance from the camera (and thus affected by translation)
    zero if the point is infinitely far away, like a star, and thus unaffected by translation.
    The camera frame is set up such that the look direction is +z, right is +x, and up is +y.

    The coefficients a through i could in theory describe a matrix with arbitrary scaling,
    shearing, mirroring, etc. In practice this code will only ever return a matrix which describes
    a pure rotation.
    """

    #Calculate the relative look direction
    look_rel = look - loc
    assert vlength(look_rel)>0,"Camera look_at same as location"
    look_rel=vnormalize(look_rel)

    #Calculate the right vector as the cross product of the relative look and sky vector
    right=vcross(look_rel,sky)
    assert vlength(right)>0,"Camera looking at sky"
    right=vnormalize(right)

    down=vcross(look_rel,right) #guaranteed to be unit-length since product of two perpendicular unit-length vectors

    r=np.zeros((4,4))
    r[0:3,0,None]=right
    r[0:3,1,None]=down
    r[0:3,2,None]=look_rel
    #result[0:3,3]=-loc
    r[3,3]=1
    t=np.zeros((4,4))
    t[0,0]=1
    t[1,1]=1
    t[2,2]=1
    t[3,3]=1
    t[0:3,3,None]=loc
    result=t@r
    result=np.linalg.inv(result)
    return result


def project(right=None,angle=None,width=None,height=None,target_c=None,out_nan:bool=True):
    """
    Project the target into the camera field of view
    :param right: Length of Right vector, equivalent to camera{right -x*...}.
                  This controls the aspect ratio. Note that a right-handed
                  camera should use a positive value for right. The up vector is
                  implicitly camera{up y*1 ...} .
    :param angle: Field-of-view angle in degrees, equivalent to camera{angle...}
                  Length of direction vector is calculated from this and right
    :param width: Width of image in pixels, equivalent to image_width
    :param height: Height of image in pixels, equivalent to image_height
    :param target: 3D position of point to project in camera coordinates.
                   If your point is in world coordinates, transform it first
                   with cmatrix(...)@target
    :return: 2D position on camera, in the form of a 2xN numpy array. Row 0 is
             horizontal coordinate, row 1 is vertical
    """
    #Convert to normalized screen coordinates. In this frame, the screen is on a plane perpendicular and
    #out along the z axis The edges of the screen are at +-0.5*up and +-0.5*right. Angle determines the distance
    #between the camera and the plane of the screen. Using the image at http://www.povray.org/documentation/view/3.7.0/246/
    #as a reference, tan(angle/2)=0.5*right/direction. We can solve this for direction:
    # tan(angle/2)*direction=0.5*right
    # direction=0.5*right/tan(angle/2)
    direction=0.5*right/np.tan(np.radians(angle)/2)
    #if the z component is negative, we don't want to plot. Do this by setting the z component to NaN if it was negative.
    target_c[2,target_c[2,...]<0]=float('NaN')
    #If the target is at the screen, then the x and y coordinates are already what we want. If it is twice as far,
    #then we need to divide x and y by 2. If half as far, then they need to multiply by two. In general, multiply
    # by direction/z. If we do this right, the z coordinate will become equal to direction, which indicates the other
    # components are normalized screen coordinates
    target_scl=target_c[0:2,...]*direction/target_c[2,...]
    result=np.zeros(target_scl.shape)
    cx=width/2
    cy=height/2
    up=1
    rx=linterp(-0.5*right,-width /2,0.5*right,width /2,target_scl[0,...])
    ry=linterp(-0.5*up,   -height/2,0.5*up   ,height/2,target_scl[1,...])
    result[0,...]=rx+cx
    result[1,...]=ry+cy
    if out_nan:
        result[:,result[0,...]<0]=float('NaN')
        result[:,result[1,...]<0]=float('NaN')
        result[:,result[0,...]>width]=float('NaN')
        result[:,result[1,...]>height]=float('NaN')
    return result


def make_sky(clock,dir):
    ssky=vnormalize(vcross(dir,np.array([[0.0],[0.0],[1.0]])))
    csky=vnormalize(vcross(ssky,dir))
    sky=np.cos(np.deg2rad(clock))*csky+np.sin(np.deg2rad(clock))*ssky
    return sky


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
    def __init__(self,casename):
        #initialize from parameters
        self.casename=casename

        #initialize by table lookup
        self.framenum0,self.framenum1,self.framepat,self.goodstars,self.badstars,self.spiceobjs=projects[casename]
        self.framenum=self.framenum0

        #set up graphics
        self.fig = plt.figure("Main fitting")
        self.ax = self.fig.add_subplot(111)

        #initialize by reading and calculation
        furnsh("data/spice/vgr2.tm")
        self.et=None
        self.open_frame_index(casename)
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

        self.fig_controls=plt.figure("Controls")
        self.makebtn(0.55,0.05,'-camlon',self.camlonm)
        self.makebtn(0.80,0.00,'-clock',self.clockm)
        self.makebtn(0.80,0.10,'+clock',self.clockp)
        self.makebtn(0.45,0.05,'+camlon',self.camlonp)
        self.makebtn(0.50,0.00,'-camlat',self.camlatm)
        self.makebtn(0.50,0.10,'+camlat',self.camlatp)
        self.makebtn(0.95,0.00,'-angle',self.anglem)
        self.makebtn(0.95,0.10,'+angle',self.anglep)
        self.makebtn(0.90,0.00,'-time',self.timem)
        self.makebtn(0.90,0.10,'+time',self.timep)
        self.makebtn(0.85,0.00,'-right',self.rightm)
        self.makebtn(0.85,0.10,'+right',self.rightp)
        self.makebtn(0.00,0.00,'/step',self.sm)
        self.makebtn(0.10,0.00,'*step',self.big)
        self.makebtn(0.00,0.05,'<frame',self.framem)
        self.makebtn(0.10,0.05,'>frame',self.framep)
        self.makebtn(0.05,0.00,'fit',self.fit)
        self.makebtn(0.00,0.10,'<auto',self.autom)
        self.makebtn(0.10,0.10,'auto>',self.autop)
        self.use_par={}
        self.makechk(0.50,0.15,'lat')
        self.makechk(0.45,0.15,'lon')
        self.makechk(0.95,0.15,'angle')
        self.makechk(0.80,0.15,'clock')
        self.makechk(0.85,0.15,'right')

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
        LimitMag = 10
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
        for k,(name,pre,a,b,c) in self.spiceobjs.items():
            self.diskplot[k],=self.ax.plot(self.c*0,self.s*0,'m-')
            self.diskplot[k].set_visible(False)
        self.fitplot,  = self.ax.plot(np.zeros((self.star_vec.shape[1],)),
                                      np.zeros((self.star_vec.shape[1],)), 'g*')
        self.fitplot.set_visible(False)
    def open_frame_index(self,casename):
        """
        Open an SQLite database and make sure that the appropriate table(s) are present
        in the database
        """
        dbname=f"data/db/frame_index_{casename}.sqlite"
        self.conn = sqlite3.connect(dbname)
        sql = ("create table if not exists frames (" +
               "framenum    integer not null," +
               "timestamp datetime default CURRENT_TIMESTAMP,"
               "nstars          integer,"+
               "et              real," +
               "lat_c           real," +
               "lon_c           real," +
               "angle           real," +
               "clock           real," +
               "right_denom     real," +
               "lat_c_sig       real," +
               "lon_c_sig       real," +
               "angle_sig       real," +
               "clock_sig       real," +
               "right_denom_sig real," +
               "primary key (framenum))")
        cur = self.conn.cursor()
        cur.execute(sql)
        self.conn.commit()
    def read_record(self):
        has_row=False
        #Check if this frame is already recorded
        sql=("select framenum,nstars,et,lat_c    ,lon_c    ,angle    ,clock    ,right_denom,"+
                                       "lat_c_sig,lon_c_sig,angle_sig,clock_sig,right_denom_sig "+
                                       "from frames order by abs(framenum-?) asc")
        cur = self.conn.cursor()
        old_et = self.et
        for row in cur.execute(sql, (self.framenum,)):
            if row[0]==self.framenum:
                (self.nstars, self.et,
                 self.lat_c    ,self.lon_c    ,self.angle,    self.clock    ,self.right_denom,
                 self.lat_c_sig,self.lon_c_sig,self.angle_sig,self.clock_sig,self.right_denom_sig)=row[1:]
                has_row = True
                break
            elif not has_row:
                fn0=row[0]
                lat_c0       = row[3]
                lon_c0       = row[4]
                angle0       = row[5]
                clock0       = row[6]
                right0       = row[7]
                has_row=True
            else:
                fn1          = row[0]
                lat_c1       = row[3]
                lon_c1       = row[4]
                angle1       = row[5]
                clock1       = row[6]
                right1       = row[7]
                self.lat_c       = linterp(fn0,lat_c0,fn1,lat_c1,self.framenum)
                self.lon_c       = linterp(fn0,lon_c0,fn1,lon_c1,self.framenum)
                self.angle       = linterp(fn0,angle0,fn1,angle1,self.framenum)
                self.clock       = linterp(fn0,clock0,fn1,clock1,self.framenum)
                self.right_denom = linterp(fn0,right0,fn1,right1,self.framenum)
                self.lat_c_sig=float('inf')
                self.lon_c_sig=float('inf')
                self.angle_sig=float('inf')
                self.clock_sig=float('inf')
                self.right_denom_sig=float('inf')
                self.nstars=None
                break
        if not has_row:
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
            self.et=str2et("1989-08-25 04:00:00 UTC")-11*3600-20*60
    def write(self):
        sql=("insert or replace into frames (framenum,nstars,rmsdiff,et,"+
             "lat_c    ,lon_c    ,angle    ,clock    ,right_denom,   "+
             "lat_c_sig,lon_c_sig,angle_sig,clock_sig,right_denom_sig) "+
             "values (?,?,?,?,?,?,?,?,?,?,?,?,?,?)")
        cur = self.conn.cursor()
        cur.execute(sql, (self.framenum,
                          self.nstars,
                          self.rmsdiff,
                          self.et,
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
                          ))
        self.conn.commit()
    def project_stuff(self):
        dir=llr2xyz(lat=self.lat_c,lon=self.lon_c)
        sky=make_sky(self.clock,dir)
        self.C = cmatrix(loc=np.zeros((3,1)), look=dir, sky=sky)
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
    def big(self, event):
        self.step_size*=10
    def sm(self, event):
        self.step_size/=10
    def plotspice(self):
        print(f"Voyager 2 in kernel {which_kernel('SPK',-32,self.et)}")
        for i_spice,(name,r,a,b,c) in self.spiceobjs.items():
            try:
                bodyframe=f"IAU_{name.upper()}"
                #Vector labels have two letters:
                # * First is center, one of:
                #    - v: Voyager Spacecraft
                #    - b: center of body in question
                # * Second is coordinate frame
                #    - b: body-fixed frame of body in question
                #    - i: Global inertial frame (Ecliptic B1950)
                # Observe the target natural spice object from Voyager
                # at the given time, in the body's own body-fixed frame
                xbody_vb,lt=spkezr(str(i_spice),self.et,bodyframe,"LT+S","-32")
                # position and velocity of Voyager relative to body is
                # reverse of position of body relative to Voyager. Likewise
                # for velocity.
                xvoy_bb=-xbody_vb
                rvoy_bb=xvoy_bb[:3]
                if True:
                    # Use pre-encounter spherical radius for all ellipsoid radii
                    a,b,c=r,r,r
                if a>0:
                    ell_bb=edlimb(a,b,c,rvoy_bb)
                    rvoy_bb=rvoy_bb.reshape(-1,1)
                    # limb points in body centered body frame
                    rlimb_bb=self.c*ell_bb.semi_major.reshape(-1,1)+self.s*ell_bb.semi_minor.reshape(-1,1)+ell_bb.center.reshape(-1,1)
                else:
                    # If body size is zero, then put all the limb points at the center
                    rlimb_bb=np.zeros((3,len(self.c)))
                # limb points in Voyager-centered body frame
                rlimb_vb=rlimb_bb-rvoy_bb
                M_ib=pxform(bodyframe,"ECLIPB1950",self.et)
                # limb points in Voyager-centered inertial frame
                rlimb_vi=M_ib @ rlimb_vb
                # Stack up into homogeneous vectors
                rlimb_vi=np.vstack((rlimb_vi,np.zeros((1,rlimb_vi.shape[1]))))
                pixlimb= project(right=4 / self.right_denom, angle=self.angle, width=self.width,
                                 height=self.height, target_c=self.C @ rlimb_vi)

                print(f"{i_spice} in kernel {which_kernel('SPK',i_spice,self.et)}")
                self.diskplot[i_spice].set_xdata(pixlimb[0,:])
                self.diskplot[i_spice].set_ydata(pixlimb[1,:])
                self.diskplot[i_spice].set_visible(True)
            except Exception:
                print(f"{i_spice} not in kernels")
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
        plt.pause(0.001)
    def camlonp(self, event):
        self.lon_c+=self.step_size
        self.replot()
    def camlonm(self, event):
        self.lon_c-=self.step_size
        self.replot()
    def camlatp(self, event):
        self.lat_c+=self.step_size
        self.replot()
    def camlatm(self, event):
        self.lat_c-=self.step_size
        self.replot()
    def lookxp(self, event):
        self.look[0]+=self.step_size
        self.replot()
    def lookxm(self, event):
        self.look[0]-=self.step_size
        self.replot()
    def lookyp(self, event):
        self.look[1]+=self.step_size
        self.replot()
    def lookym(self, event):
        self.look[1]-=self.step_size
        self.replot()
    def anglep(self, event):
        self.angle += self.step_size
        self.replot()
    def anglem(self, event):
        self.angle -= self.step_size
        self.replot()
    def rightp(self, event):
        self.right_denom += self.step_size
        self.replot()
    def rightm(self, event):
        self.right_denom -= self.step_size
        self.replot()
    def clockp(self, event):
        self.clock += self.step_size
        self.replot()
    def clockm(self, event):
        self.clock -= self.step_size
        self.replot()
    def timep(self, event):
        self.et += self.step_size*3600
        self.replot()
    def timem(self, event):
        self.et -= self.step_size*3600
        self.replot()
    def framem(self, event):
        self.framenum-=1
        self.read_record()
        self.load_image()
        if self.fitplot is not None:
            self.fitplot.set_visible(False)
        self.replot()
    def framep(self, event):
        self.framenum+=1
        self.read_record()
        self.load_image()
        if self.fitplot is not None:
            self.fitplot.set_visible(False)
        self.replot()
    def fit(self, event):
        """
        Given the current position as an initial guess, find the optimum
        camera parameters and position to fit the stars.
        """
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
            (findx,findy),(sigx,sigy),rho=find_stars(self.backimg,goodstar_pix,names=self.starnames,ax=None)#self.ax_controls)
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
            vary=[         bounded(-90.0,90.0,self.lat_c),rbounded(-180.0,180.0,self.lon_c),bounded(0.0,90.0,self.angle),bounded(-180.0,180.0,self.clock),positive()       ]
            if len(findx)<2:
                p0[2]=False
                p0[3]=False
                p0[4]=False
                print("Few usable stars, only fitting pointing")
            elif len(findx)<5:
                p0[3]=False
                p0[4]=False
                print("Few usable stars, only fitting pointing and angle")
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
    def autop(self, event):
        for i in range(60):
            self.framep(event)
            self.fit(event)
    def autom(self, event):
        for i in range(100):
            self.framem(event)
            self.fit(event)

def main():
    boxfig = None
    boxax = None
    boximg = None
    callback=CameraMount(casename="VoyagerNeptune")

    plt.show()


if __name__=="__main__":
    main()

