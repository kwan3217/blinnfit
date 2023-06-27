"""
Given an image and a set of pixels of where there are *supposed* to be stars,
extract the positions of the *actual* stars
"""
from collections import namedtuple

import numpy as np
from kwanmath.gaussian import fit_twoD_Gaussian
from matplotlib import pyplot as plt
from matplotlib.axes import Axes


def subset_img(img:np.array,xc:int,yc:int,boxr:int=10)->np.array:
    """
    Given an image, get the square subimage centered on the given pixel
    with the given square "radius"

    :param img: image to subset. If the image is more than 2D (say RGB),
                the first two dimensions are y and x, and the result is summed
                across any other dimensions
    :param xc:  horizontal coordinate of center of subset
    :param yc:  vertical coordinate of center of subset
    :param boxr: radius of box
    :return: tuple of:
      * subset of image
      * x and y coordinate of requested center pixel relative to index 0,0 in this box

    If square hangs over the edge of the image, return just the
    part of the image that is in bounds (might not be square).

    If the square is completely off the image, return None.
    """
    def constrain(xx,xmin,xmax):
        if xx<xmin:
            return xmin
        if xx>=xmax:
            return xmax-1
        return xx
    if xc<-boxr:
        return None,None,None
    elif xc>=img.shape[1]+boxr:
        return None,None,None
    elif yc<-boxr:
        return None,None,None
    elif yc>=img.shape[0]+boxr:
        return None,None,None
    x0=constrain(xc-boxr,0,img.shape[1])
    x1=constrain(xc+boxr,0,img.shape[1])
    y0=constrain(yc-boxr,0,img.shape[0])
    y1=constrain(yc+boxr,0,img.shape[0])
    box = img[y0:y1, x0:x1,...]
    if len(box.shape)>2:
        box = np.sum(box, tuple(range(2,len(box.shape))))
    box = box - np.min(box)
    return box,xc-x0,yc-y0


def find_star(box:np.array,name:str=None,ax:Axes=None)->np.array:
    """
    Given a subset of an image, find the one star in it

    :param box: subset of image to search for a star
    :param ax: Axes to draw the results on, or None to not draw the results
    :return: Coordinates of centroid of image subset

    Raises an exception if we can't find the star
    """
    if ax is not None:
        ax.clear()
        ax.imshow(box)
        if name is not None:
            ax.set_title(name)
        plt.pause(0.001)

    # Try the 2D Gaussian fit on the subset
    reject=False
    data_fitted=None
    try:
        fit= fit_twoD_Gaussian(box,dn_size=1.0/255.0)
        data_fitted=fit.eval
        residual_stdev=np.std(box-fit.eval)
        print(f"{name},amp={fit.amp:8.3f},xc={fit.xc:8.3f}+-{np.sqrt(fit.cov[1,1]):5.3f},yc={fit.yc:8.3f}+-{np.sqrt(fit.cov[2,2]):5.3f},sigx={fit.sigx:8.3f},sigy={fit.sigy:8.3f},rho={fit.rho:8.5f},ofs={fit.ofs:8.3f},res_stdev={residual_stdev:8.3f}")
    except Exception:
        print(f"{name} Reject: fit failed")
        reject=True
    if not reject and np.abs(fit.sigx) > 5:
        print("Reject: bad sigx %f" % fit.sigx)
        reject = True
    elif not reject and np.abs(fit.sigy) > 5:
        print("Reject: bad sigy %f" % fit.sigy)
        reject = True
    elif not reject:
        if fit.amp<residual_stdev*3:
            print(f"Reject: low amplitude {fit.amp}<{3*residual_stdev}")
            reject=True
        elif fit.amp==0.00:
            print(f"Reject: low abs amplitude {fit.amp}")
            reject=True
        elif not np.isfinite(fit.cov[1,1]):
            print(f"Reject: bad uncertainty on xc {fit.cov[1,1]}")
            reject=True
    if reject:
        fit=namedtuple("fit_failed","xc yc cov")(xc = float('nan'),yc = float('nan'),cov=np.zeros((7,7))*float('nan'))
        residual_stdev=float('nan')
    if ax is not None:
        color = 'r' if reject else 'g'
        if data_fitted is not None:
            ax.contour(data_fitted, 8, colors=color)
        plt.pause(0.001)
        # plt.waitforbuttonpress()
    return fit.xc, fit.yc, fit.cov[1:3,1:3],residual_stdev


def find_stars(img:np.array,g_cs:np.array,boxr:int=10,names:list=None,ax:Axes=None):
    """
    Given an image with stars on it and the predicted location of
    a bunch of stars, find the actual position of each star

    :param img: Image to search
    :param g_cs: guess coordinates of stars, in the form of a 2xN numpy array
    :param boxr:
    :return:
    """
    result=np.zeros(g_cs.shape)
    sig_c=result*0
    rho_c=np.zeros(g_cs.shape[1])
    std_c=rho_c*0
    min_std=float('inf')
    for i,(g_xc,g_yc) in enumerate(zip(g_cs[0,:],g_cs[1,:])):
        name=names[i] if names is not None else None
        if np.isfinite(g_xc):
            g_xc = int(g_xc)
            g_yc = int(g_yc)
            box,g_bxc,g_byc=subset_img(img,g_xc,g_yc,boxr=boxr)
            bxc,byc,pxcyc,resid_stdev=find_star(box,name=name,ax=ax)
            result[0,i]=bxc+g_xc-g_bxc
            result[1,i]=byc+g_yc-g_byc
            sig_c[0,i]=np.sqrt(pxcyc[0,0])
            sig_c[1,i]=np.sqrt(pxcyc[1,1])
            rho_c[i]=pxcyc[0,1]/sig_c[0,i]/sig_c[1,i]
            std_c[i]=resid_stdev
            if resid_stdev<min_std:
                min_std=resid_stdev
        else:
            result[0,i]=float('nan')
            result[1,i]=float('nan')
    std_c/=min_std
    sig_c*=std_c
    return result,sig_c,rho_c




