"""
Refactor the star fitting stuff to here. Ultimately it is a function
which takes a list of stars, a star image, and an initial guess for
the camera paramters and returns an optimized set of camera parameters.


Created: 2/2/25
"""
from copy import copy
from typing import Callable

import numpy as np
from kwanmath.gaussian import correlation_matrix, infamily
from matplotlib.figure import Figure
from scipy.optimize import curve_fit, Bounds

from starfit.find_stars import find_stars
from starfit.camera import Camera, Source


def fit_stars(*,img:np.ndarray,
                star_vs:np.ndarray,
                star_names:np.ndarray[str],
                camera0:Camera,
                f_prefit:Callable=None,
                f_findstars:Callable=None,
                f_fitstars:Callable=None,
                ax_box:Figure=None,
                ax_img:Figure=None
              )->Camera:
    """
    Given the current position as an initial guess, find the optimum
    camera parameters and position to fit the stars.
    :param img: Image to fit
    :param star_vs: Unit vector towards each star to fit, should already be filtered by goodstars
    :param star_names: Names of stars
    :param camera0: Initial guess of camera pointing
    :param i_freeze: scipy.optimize.curve_fit() is advertised to allow Bounds such that
                     the upper and lower bounds of the parameter are the same,
                     which is how we describe freezing a parameter. Unfortunately that
                     doesn't work so we just truncate parameters that aren't used
                     and pass the initial guesses as arguments instead. This means
                     that only parameters from the end of the list lat,lon,angle,clock,right_denom
                     can be frozen. The order of the parameters is chosen specifically
                     with this constraint in mind. i_freeze is the index of the first
                     frozen parameter. This and ones after will not be solved for.
    """
    if f_prefit is not None:
        f_prefit()
    camera=copy(camera0)
    done = False
    dropset=set()
    while not done:
        # Trim stars in the dropset
        keep=[star_name not in dropset for star_name in star_names]
        this_star_pixs,on_screen=camera.project(star_vs)
        goodstars=np.logical_and(on_screen,keep)
        print("Number of stars on-screen:      ", np.sum(on_screen))
        print("Number of stars kept:      ", np.sum(goodstars))
        this_star_pixs=this_star_pixs[:,goodstars]
        this_star_vs=star_vs[:,goodstars]
        this_star_names=star_names[goodstars]
        findc, sigc, rho = find_stars(img, this_star_pixs, names=this_star_names,
                                      ax_box=ax_box,ax_img=ax_img, boxr=10)
        w = np.isfinite(findc[0,:])
        print("Number of good stars found:     ", np.sum(w))
        if f_findstars is not None:
            f_findstars(starpixo=findc,starpixc=this_star_pixs,w=w,names=this_star_names)
        # Trim down star list again. Those that weren't found get NaN for their pixel locations
        this_star_names = this_star_names[w]
        findc = findc[:,w]
        findx=findc[0,:]
        findy=findc[1,:]
        sigc = sigc[:,w]
        rho = rho[w]
        this_star_vs=this_star_vs[:,w]
        cov = np.zeros((findc.shape[1] * 2, findc.shape[1] * 2))
        sigx = sigc[0, :]
        sigy = sigc[1, :]
        for i in range(findc.shape[1]):
            cov[i * 2, i * 2] = sigx[i] ** 2
            cov[i * 2 + 1, i * 2 + 1] = sigy[i] ** 2
            cov[i * 2 + 1, i * 2] = sigx[i] * sigy[i] * rho[i]
            cov[i * 2, i * 2 + 1] = sigx[i] * sigy[i] * rho[i]
        weights = 1.0 / np.sqrt(sigx ** 2 + sigy ** 2)
        # make an array [findx
        #                findy] then ravel it. The result is [findx|findy]
        pixdata = findc.ravel()
        # fiti=np.array(goodstars)[w]
        cutoff = 5
        fix_right_denom = True
        mean_right_denom = 2.19884598665809*(camera0.right_num/4)  # Average of fit values from all frames where right_denom could be fit
        if fix_right_denom:
            cutoff = 4
            right_denom0 = mean_right_denom
        p0 = camera0.to_params()

        vary = [[-90.0, 90.0], [-180.0, 180.0], [0.0, 120.0], [-180.0, 180.0], [0.0, np.inf]]
        # We intend to fit all the parameters at first, mark them as such
        sources={k:Source.FIT for k in ("lat_source", "lon_source", "angle_source", "clock_source", "right_denom_source")}
        # Now check if we have enough stars. If not enough, turn off some parameters
        if findc.shape[1] < 2:
            cutoff = 2
            del sources["angle_source"]
            del sources["clock_source"]
            del sources["right_denom_source"]
            print("Very few usable stars, only fitting pointing")
        elif findc.shape[1] < 2:
            cutoff = 3
            del sources["clock_source"]
            del sources["right_denom_source"]
            print("Few usable stars, only fitting pointing and angle")
        if fix_right_denom:
            # If right_denom is fixed, mark it as constrained despite above
            camera.right_denom = mean_right_denom
            sources["right_denom_source"] = Source.MODELED
        # Now handle longitude. It wraps, so we hand the estimator
        # a zero guess and give the fit interface an offset to add to get the
        # actual longitude.
        lon0 = camera.lon
        p0[1] = 0.0
        # Package up the bounds
        lb = [a for a, b in vary]
        ub = [b for a, b in vary]
        bounds = Bounds(lb=lb[:cutoff], ub=ub[:cutoff])
        # Set up a partial to handle everything my old wrapper did
        count=0
        def curve_fitsky_interface(starvec, lat, lon, this_angle=None, this_clock=None, this_right_denom=None):
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
            nonlocal count
            this_camera = copy(camera)
            this_camera.lat=lat
            this_camera.lon=lon=lon+lon0
            this_camera.angle=camera.angle if this_angle is None else this_angle
            this_camera.clock=camera.clock if this_clock is None else this_clock
            this_camera.right_denom=camera.right_denom if this_right_denom is None else this_right_denom
            #print(f"{count=}")
            #print(f"old {camera=}")
            #print(f"{this_camera=}")
            # We set out_nan=False so stars off the edge will still have
            # finite values, just off the edge of the image.
            prj,_ = this_camera.project(starvec,out_nan=False)
            count+=1
            return prj.ravel()
        popt, pcov, *_ = curve_fit(curve_fitsky_interface, this_star_vs, pixdata, p0=p0[:cutoff], bounds=bounds, sigma=cov, absolute_sigma=True)
        # Calculated star coordinates from the best fit that curve_fit came up with
        fitc = curve_fitsky_interface(this_star_vs, *popt)
        fitx = fitc[:len(fitc) // 2]
        fity = fitc[len(fitc) // 2:]
        # Reravel the pixel coordinates
        fitc=np.vstack((fitx,fity))
        # Extend popt and pcov
        ext_popt = p0 * 1.0  # Use the pre-optimized values as default
        ext_popt[:cutoff] = popt  # Replace with optimized values
        ext_popt[1] += lon0
        ext_pcov = np.zeros((5, 5))  # Fixed values get covariance rows and cols of 0 (IE as if perfectly known)
        ext_pcov[:cutoff, :cutoff] = pcov  # Replace upper left corner with optimized values
        popt = ext_popt
        pcov = ext_pcov
        total_lensq = 0
        # Compare the observed (o*) and calculated (c*) positions of the pixels.
        #  * Observed is where find_stars() found them
        #  * Calculated is where they are supposed to be based on the camera fit
        for name, ox, oy, cx, cy in zip(this_star_names, findx, findy, fitx, fity):
            dx=ox-cx
            dy=oy-cy
            lensq = dx ** 2 + dy ** 2
            total_lensq += lensq
            #print(f"{name},{ox:8.3f},{oy:8.3f},{cx:8.3f},{cy:8.3f},{np.sqrt(lensq):8.3f}")
        rmsdiff = np.sqrt(total_lensq / len(this_star_names))
        print(f"RMS diff: {rmsdiff:8.5f}")
        while popt[1] > 360:
            popt[1] -= 360
        while popt[1] < 0:
            popt[1] += 360
        cor = correlation_matrix(pcov)
        camera.lat,camera.lon,camera.angle,camera.clock,camera.right_denom=popt
        camera.lat_sig, camera.lon_sig, camera.angle_sig, camera.clock_sig, camera.right_denom_sig = tuple(
            [cor[i, i] for i in range(len(popt))])
        camera.__dict__.update(sources)
        camera.nstars=len(this_star_names)
        camera.rmsdiff=rmsdiff
        infam = infamily(fitc, findc, weights=weights)
        if len(infam) > 5:
            for i_infam in range(len(infam)):
                if not infam[i_infam]:
                    print(f"Star {this_star_names[i_infam]} not in family")
                    dropset.add(this_star_names[i_infam])
            done = np.all(infam)
        else:
            done = True
        if f_findstars is not None:
            f_findstars(starpixo=findc, starpixc=fitc)
    return camera


