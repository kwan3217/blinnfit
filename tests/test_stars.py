"""
Test and debug entry point for star catalog handler

Created: 2/3/25
"""
import numpy as np
from kwanmath.geodesy import xyz2llr
from matplotlib import pyplot as plt
from numpy import where

from bsc import load_catalog, parse_stars


def test_parse_stars():
    lines=load_catalog()
    vs_w,names,mags,colors=parse_stars(lines)
    w_ori=np.array(["Sgr" in name for name in names])
    lon,lat,r=xyz2llr(vs_w,deg=True)
    xs=vs_w[0,:]
    ys=vs_w[1,:]
    zs=vs_w[2,:]
    names=names[:]
    plt.figure("Map")
    ras=lon/15
    ras[ras<0]+=24
    decs=lat
    plt.plot(ras,decs,'*')
    plt.axis([24,0,-90,90])
    for ra,dec,name in zip(ras[w_ori],decs[w_ori],names[w_ori]):
        plt.text(ra,dec,name)
    plt.show()
