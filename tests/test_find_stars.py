"""

"""

import pytest
from kwanmath.gaussian import twoD_Gaussian,fit_twoD_Gaussian
import numpy as np
import matplotlib.pyplot as plt
from find_stars import find_star, find_stars


@pytest.mark.parametrize(
    "xs,ys,xc,yc,amp,sigx,sigy,rho,ofs",
    [(20,20,10,10,10,3,3,0,10)]
)
def test_twoD_Gaussian(xs,ys,xc,yc,amp,sigx,sigy,rho,ofs):
    x = np.linspace(0, xs-1, xs)
    y = np.linspace(0, ys-1, ys)
    x, y = np.meshgrid(x, y)
    img=twoD_Gaussian(x,y,amp,xc,yc,sigx,sigy,rho,ofs)
    plt.imshow(img)
    plt.show()


@pytest.mark.parametrize(
    "xs,ys,amp,a_xc,a_yc,sigx,sigy,rho,ofs",
    [(20,20,10,10,10,3,3,0,10)]
)
def test_twoD_Gaussian_fit(xs,ys,amp,a_xc,a_yc,sigx,sigy,rho,ofs):
    x = np.linspace(0, xs-1, xs)
    y = np.linspace(0, ys-1, ys)
    x, y = np.meshgrid(x, y)
    img=twoD_Gaussian(x,y,amp,a_xc,a_yc,sigx,sigy,rho,ofs)
    img+=np.random.default_rng().normal(0,1,(xs,ys))
    plt.figure("test_twoD_Gaussian_fit")
    plt.imshow(img)
    amp,xc,yc,sigx,sigy,rho,ofs,fit_img=fit_twoD_Gaussian(img)
    plt.contour(fit_img,color='w')
    print(amp,xc,yc,sigx,sigy,rho,ofs)
    plt.show()


def test_find_star():
    xs=20
    ys=20
    x = np.linspace(0, xs-1, xs)
    y = np.linspace(0, ys-1, ys)
    x, y = np.meshgrid(x, y)
    amp=100
    a_xc=8
    a_yc=12
    sigx=2.5
    sigy=1.5
    rho=0
    ofs=5
    img=twoD_Gaussian(x,y,amp,a_xc,a_yc,sigx,sigy,rho,ofs)
    img+=np.random.default_rng().normal(0,1,(xs,ys))
    img+=np.random.default_rng().normal(0,1,(xs,ys))
    fig,ax=plt.subplots(num="test_find_star")
    ax.imshow(img)
    print(find_star(img,ax=ax))
    plt.show()

def test_find_stars():
    xs=200
    ys=200
    x = np.linspace(0, xs-1, xs)
    y = np.linspace(0, ys-1, ys)
    x, y = np.meshgrid(x, y)
    amp=100
    rng=np.random.default_rng(seed=3217)
    a_xcs=rng.uniform(0,xs,size=20)
    a_ycs=rng.uniform(0,ys,size=20)
    sigx=2.5
    sigy=1.5
    rho=0
    ofs=0
    img=np.zeros((xs,ys))
    for a_xc,a_yc in zip(a_xcs,a_ycs):
        img+=twoD_Gaussian(x,y,amp,a_xc,a_yc,sigx,sigy,rho,ofs)
    img+=rng.normal(0,1,(xs,ys))
    fig,ax=plt.subplots(num="test_find_star")
    ax.imshow(img)
    plt.pause(0.001)
    print(find_stars(img,np.vstack((a_xcs+rng.normal(0,3,a_xcs.shape),a_ycs+rng.normal(0,3,a_xcs.shape))),ax=ax))
    plt.show()



if __name__=="__main__":
    test_twoD_Gaussian_fit(100,100,10,30,30,30,10,0.5,10)