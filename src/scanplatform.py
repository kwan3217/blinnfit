"""
We have a scan platform CK file with discrete pointing instances at each observation.
Convert this into a CK which has continuous coverage and slews the scan platform
such that it hits each pointing instance.
"""
from datetime import datetime, timedelta, UTC

import numpy as np
from kwanmath.interp import linterp
from kwanmath.matrix import rot_y, rot_z
from spiceypy import furnsh, ckobj, ckcov, wncard, wnfetd, str2et, sct2e, scdecd, etcal, pxform, spkezr, timout, ckgp, \
    sce2t,scencd
from spiceypy.utils.support_types import SPICEDOUBLE_CELL
import matplotlib.pyplot as plt
import os
import subprocess

from which_kernel import ls_spice


def et_to_dt(et):
    cal=timout(et,"YYYY-MM-DD HR:MN:SC.###### ::UTC")
    dt = datetime(int(cal[0:4]), int(cal[5:7]), int(cal[8:10]),
                  int(cal[11:13]), int(cal[14:16]), int(cal[17:19]), int(cal[20:26],10),tzinfo=UTC)
    return dt


def furnish_kernels(vgr:int,ScanPlatformCK:str):
    furnsh("data/spice/lsk/naif0012.tls")
    if vgr==2:
        furnsh("data/spice/fk/vg2_v02.tf")
        furnsh("data/spice/sclk/vg200045.tsc")
        furnsh("data/spice/ck/vgr2_super.bc")
        furnsh("data/spice/spk/nep095.bsp")
        furnsh("data/spice/spk/nep097.bsp")
        furnsh("data/spice/spk/vgr2_nep097.bsp")
    else:
        furnsh("data/spice/vgr1.tm")
    furnsh(ScanPlatformCK)


def get_pointing_instances(ScanPlatformCK:str, sclkid:int=-32)->tuple[list[float],list[str],list[float]]:
    """
    Get times of pointing instances from a discrete pointing kernel

    :param ScanPlatformCK:
    :return: list of tuples, one for each pointing. Returned value is
      * encoded sclk of pointing
      * string sclk of pointing
      * translated ET using best available sclk kernel
      * todo: transformation matrix at the given time
    """
    ckids=ckobj(ScanPlatformCK)
    print(ckids[0])
    cover = SPICEDOUBLE_CELL(200000)
    cover=ckcov(ScanPlatformCK,ckids[0],False,"INTERVAL",0.0,"SCLK",cover)
    print(wncard(cover))
    ticks=[]
    ets=[]
    sclkstrs=[]
    for i in range(wncard(cover)):
        tick0,tick1=wnfetd(cover,i)
        assert tick0==tick1,"Not a discrete pointing kernel like expected"
        et0=sct2e(sclkid,tick0)
        sclkstr=scdecd(sclkid, tick0)
        #assert tick0==sce2t(sclkid,et0),"Didn't round-trip"
        print(f"{i:7d},{tick0:20.7f},{sclkstr},{et0:20.7f},{etcal(et0)}")
        ticks.append(tick0)
        sclkstrs.append(scdecd(sclkid,tick0))
        ets.append(et0)
    return ticks,sclkstrs,ets


def mtx_to_euler(R, verbose=False):
    """
    Calculate the Euler angles in the Proper ZYZ Extrinsic form, compatible
    with co-elevation, co-azimuth and twist for the Voyager spacecraft.
    See https://omoikane.kwansystems.org/wiki/index.php/Matrix_to_Euler_Angle

    :param R: rotation matrix. If this is not a proper rotation matrix, you
              will silently get a wrong answer.
    :return: a tuple (twist,co-elevation,co-azimuth) of Euler angles

    The proper ZYZ extrinsic form means that from the home position (~90deg pitch up
    from squared-up Von Karman position) you:

    * rotate around reference Z axis by twist
    * rotate around reference Y axis by co-elevation
    * rotate around reference Z axis again by co-azimuth
    """
    ce = R[2, 2]
    se = np.sqrt(1 - ce ** 2)
    ca = R[0, 2] / se
    sa = R[1, 2] / se
    ct = -R[2, 0] / se
    st = R[2, 1] / se

    # Check Pythagorean identity
    if verbose:
        print(ca ** 2 + sa ** 2, " Pythag(a) should be 1")
        print(ct ** 2 + st ** 2, " Pythag(t) should be 1")

    # Check other components
    if verbose:
        print(R[0, 0], " r00 should be ", ca * ce * ct - sa * st)
        print(R[0, 1], " r01 should be ", -ca * ce * st - ct * sa)
        print(R[1, 0], " r10 should be ", ca * st + ce * ct * sa)
        print(R[1, 1], " r11 should be ", ca * ct - ce * sa * st)
    twist = np.arctan2(st, ct)
    coelevation = np.arctan2(se, ce)
    coazimuth = np.arctan2(sa, ca)

    # Check that the result is a rotation matrix with no reflection
    R = rot_z(c=ca, s=sa) @ rot_y(c=ce, s=se) @ rot_z(c=ct, s=st)
    if verbose:
        print(np.linalg.det(R), " det should be 1")
    if np.linalg.det(R) < 0:
        raise ValueError("We hit one")
    return (twist, coelevation, coazimuth)


def get_eulers(sc:int,
               ets:list[float],
               ticks:list[float],
               sclkstrs:list[str],
               et0:float=None,
               et1:float=None):
    vgr=-sc%10
    if et0 is None:
        et0=ets[0]
    if et1 is None:
        et1=ets[-1]
    twists = []  # Twists
    coelevations = []  # Coelevations
    coazimuths = []  # Coazimuths
    eet = []  # Teph's
    j2000 = datetime(2000, 1, 1, 12, 0, 0)
    dts = []  # Timestamps
    for (i, (et, tick, sclk)) in enumerate(zip(ets, ticks, sclkstrs)):
        if et < et0 or et > et1:
            continue
        print(i, etcal(et))
        try:
            M = pxform(f"VG{vgr}_SCAN_PLATFORM", "VG{vgr}_AZ_EL", et)
            print("No spice error")
        except:
            print("Spice error (probably no super_ck)")
            MT,tickout=ckgp(sc*1000-100,tick,0,f"VG{vgr}_AZ_EL")
            M=MT.T
            #ls_spice(verbose=True)
            #continue
        twist, coelevation, coazimuth = mtx_to_euler(M)
        twists.append(twist)
        coelevations.append(coelevation)
        coazimuths.append(coazimuth)
        eet.append(et)
        dts.append(et_to_dt(et))
    twists = np.degrees(twists)
    coelevations = np.degrees(coelevations)
    coazimuths = np.degrees(coazimuths)

    plt.plot(dts, twists, 'r+', label='twist')
    plt.plot(dts, coelevations, 'g+', label='coelevation')
    plt.plot(dts, coazimuths, 'b+', label='coazimuth')
    plt.legend()
    plt.pause(0.001)
    return eet,twists,coelevations,coazimuths,dts



def make_slews(ets:list[float],coelevations:list[float],coazimuths:list[float],
               allowed_slew_rates:tuple[float]=(0.08, 0.33, 1.0),settle_time:float=5.0):
    """
    Given a bunch of discrete elevations and azimuths,
    generate continuous slews that hit those targets

    :param ets: ephemeris time of each discrete pointing
    :param coelevations: coelevation at each time
    :param coazimuths: coazimuth at each time
    :return: list of tuples, one for each segment of motion at constant speed
      * ET of segment start
      * ET of segment end
      * coelevation at start of segment (deg)
      * coazimuth at start of segment (deg)
      * elevation speed (deg/s)
      * azimuth speed (deg/s)

    Note:
      Input data might generate extremely short segments. Segments shorter than
      one Voyager clock tick are dropped, with the preceding segment extended to
      cover the dropped segment.

      As a result, the time coverage will be continuous, but the pointing may
      have discontinuities. These will generally be small.
    """

    slew_rates = []
    dt_seg=[]
    et_seg=[]
    coel0_seg=[]
    coaz0_seg=[]
    dcoel_seg=[]
    dcoaz_seg=[]
    for i_nac,(et0,coel0,coaz0,et1,coel1,coaz1) in enumerate(zip(ets[:-1],coelevations[:-1],coazimuths[:-1],
                                                                 ets[1: ],coelevations[1: ],coazimuths[1: ])):
        #distances to slew in each direction
        d_el = np.abs(coel1 - coel0)
        d_az = np.abs(coaz1 - coaz0)
        d = d_el + d_az
        #time available to do the slew
        d_t = np.abs(et1 - et0)
        time_available = d_t - settle_time
        #required speed of slew
        req_spd = d / time_available
        found = False
        for i_slew, slew_rate in enumerate(allowed_slew_rates):
            if slew_rate >= req_spd:
                found = True
                break
        if not found:
            #if the slew is too fast, slew as fast as needed (but no faster)
            #in order to hit the marks
            print(f"Warning -- Slew {i_nac:5d} at {et_to_dt(et0)} is too fast: {d_el:7.3f} el and {d_az:7.3f} az in {d_t:7.3f}sec")
            slew_rate = req_spd
        slew_rates.append(slew_rate)
        #time to do each slew
        y_el = d_el / slew_rate
        y_az = d_az / slew_rate
        # Set up first segment - non-motion from previous time to current time-(x+y_el+y_az)
        t00 = et0 - settle_time
        t01 = et1 - (settle_time + y_el + y_az)
        if t00 > t01:
            raise ValueError("Time travel?")
        if(t01-t00<0.06):
            print("less than 1tick")
        et_seg.append(t00)
        dt_seg.append(et_to_dt(t00))
        coel0_seg.append(coel0)
        coaz0_seg.append(coaz0)
        dcoel_seg.append(0)
        dcoaz_seg.append(0)
        # Set up second segment - motion around elevation axis
        t10 = t01
        t11 = et1 - (settle_time + y_az)
        if t10 > t11:
            raise ValueError("Time travel?")
        if(t11-t10<0.06):
            print("less than 1tick")
        et_seg.append(t10)
        dt_seg.append(et_to_dt(t10))
        coel0_seg.append(coel0)
        coaz0_seg.append(coaz0)
        dcoel_seg.append(np.copysign(slew_rate,(coel1-coel0)))
        dcoaz_seg.append(0)
        # Set up third segment - motion around azimuth axis
        t20 = t11
        t21 = et1 - (settle_time)
        if t20 > t21:
            raise ValueError("Time travel?")
        if(t21-t20<0.06):
            print("less than 1tick")
        et_seg.append(t20)
        dt_seg.append(et_to_dt(t20))
        coel0_seg.append(coel1)
        coaz0_seg.append(coaz0)
        dcoel_seg.append(0)
        dcoaz_seg.append(np.copysign(slew_rate, (coaz1 - coaz0)))
    plt.figure(2)
    plt.plot(dt_seg,coel0_seg,label='coel')
    plt.plot(dt_seg,coaz0_seg,label='coaz')
    plt.plot(dt_seg[::3],linterp(0,-180,1.2,180,np.array(slew_rates)),'*',label='slew_rate')
    for rate in allowed_slew_rates:
        plt.plot([dt_seg[0],dt_seg[-1]],np.array((1,1))*linterp(0, -180, 1.2, 180, rate), '-', label='{rate=:.02f}')
    plt.legend()
    plt.pause(0.001)
    return et_seg,dt_seg,coel0_seg,coaz0_seg,dcoel_seg,dcoaz_seg


def merge_slews(et_seg, dt_seg, coel0_seg, coaz0_seg, dcoel_seg, dcoaz_seg):
    et_res=et_seg[0:2]
    dt_res=dt_seg[0:2]
    coel0_res=coel0_seg[0:2]
    coaz0_res=coaz0_seg[0:2]
    dcoel_res=dcoel_seg[0:2]
    dcoaz_res=dcoaz_seg[0:2]
    for i_slew,(et,dt,coel0,coaz0,dcoel,dcoaz) in enumerate(
            zip(et_seg[2:], dt_seg[2:], coel0_seg[2:], coaz0_seg[2:], dcoel_seg[2:], dcoaz_seg[2:]),start=2):
        et_end=et
        et_mid=et_res[-1]
        et_beg=et_res[-2]
        coel_end = coel0
        coel_mid = coel0_res[-1]
        coel_beg = coel0_res[-2]
        coaz_end = coaz0
        coaz_mid = coaz0_res[-1]
        coaz_beg = coaz0_res[-2]
        dt_end=et_end-et_mid
        if dt_end<0.06:
            # calculate the slope necessary to hit the *end* of the current segment
            # from the *beginning* of the previous segment
            dt_twoseg=et_end-et_beg
            dcoel=coel_end-coel_beg
            dcoaz=coaz_end-coaz_beg
            # change the previous segment so that it runs the whole beginning to end
            et_res[-1]=et_end
            dt_res[-1]=dt
            coel0_res[-1]=coel_end
            coaz0_res[-1]=coaz_end
            dcoel_res[-2]=dcoel/dt_twoseg
            dcoaz_res[-2]=dcoaz/dt_twoseg
        else:
            # add this segment
            et_res.append(et_end)
            dt_res.append(dt)
            coel0_res.append(coel_end)
            coaz0_res.append(coaz_end)
            dcoel_res.append(dcoel)
            dcoaz_res.append(dcoaz)
    return et_res,dt_res,coel0_res,coaz0_res,dcoel_res,dcoaz_res


def msopck(et_seg: list[float], dt_seg:list[datetime],
           coel0_seg:list[float], coaz0_seg:list[float],
           dcoel_seg:list[float], dcoaz_seg:list[float],
           sc=-32, planet="Neptune"):
    msopck_txt = rf"""
    Voyager {-sc%10} {planet} encounter scan platform kernel

    This kernel represents the orientation of the scan platform
    relative to the az/el frame. It is based upon two kernels:

    * One is derived from the ISS SEDR, whose raw data gives the
      position of the scan platform at each NAC image, but no
      data in between. This has been processed using optical 
      navigation, and therefore is relative to the star catalog
      and not the orientation of the spacecraft. This gives pointing
      only at discrete points in time.
    * The other is the orientation of the spacecraft relative
      to the stars. This has the measured position at a period
      of as low as 48s over almost the entire mission, including
      the {planet} near encounter.

    The present kernel is constructed by figuring the orientation
    of the scan platform relative to the az/el frame, and then
    resolving that orientation into azimuth, elevation, and twist
    (around the instrument boresight) angles. Since Voyager doesn't
    have a twist actuator, the twist angle should be zero. Experience
    shows that it is near zero almost all of the time.

    Once these angles are obtained, we do the following:

    * Force the twist angle to zero. This will not change the 
    direction of the boresight, but will change the orientation
    of the image around the boresight.
    * Figure out the slew. The spacecraft always used one axis 
    of motion at a time, in any of three slew rates. It always 
    hits its mark in time to take the image. We will simulate this by:
      - Figure the distance to move in azimuth and elevation
      - Assume that the scan platform hits its mark X seconds before
        the exposure (time of data point in NAC kernel). We will try
        x=5sec for now.
      - Before the scan platform reaches its mark, it slews in one axis
        at a time at the slowest rate allowed by the distance and time
        between two NAC kernel data points. We add up the distance to
        slew in each direction, and divide by each slew rate to get
        the time required to do the slew at that rate. The longest
        slew time which fits in between the NAC kernel points is used --
        this would put minimum wear on the mechanisms in real life,
        and generates slower slews which are more likely to be seen in
        an animation. We call the slew time Y, with Y_el the elevation
        portion and Y_az the azimuth portion.
      - Generate a zero-motion position from the previous NAC point to 
        X+Y_el+Y_az from the current NAC time
      - Generate an elevation slew starting at the previous NAC point
        with a constant rate around the Y axis, starting at X+Y_el+Y_az
        before the current NAC time
      - Generate an azimuth slew starting at the azimuth of the previous
        NAC point and the elevation of the current point, with a constant
        rate around the Z axis, starting at X+Y_az before the current 
        NAC time. This will result in hitting the mark X seconds before
        the NAC time.
      - Generate a zero-motion position from X seconds before the current
        NAC time to the NAC time.

    Since all of this is done with Euler azimuth/elevation angles, the result
    will be in the az/el frame of the spacecraft.

    Times are in ET. I would have preferred to specify in SCLK, but
    I am not sure how to interpolate. At all of the ETs of the photos,
    the ET matches exactly the SCLK in the source kernel given the 
    SCLK kernel listed below, raw from 

    \begindata
    INPUT_DATA_TYPE = 'EULER ANGLES'
    INPUT_TIME_TYPE = 'ET'
    ANGULAR_RATE_PRESENT= 'YES'
    CK_TYPE = 2
    INSTRUMENT_ID={sc*1000-100}
    REFERENCE_FRAME_NAME='VG{-sc%10}_AZ_EL'
    PRODUCER_ID='C. Jeppesen, Kwan Systems'
    FRAMES_FILE_NAME='data/spice/fk/vg{-sc%10}_v02.tf'
    SCLK_FILE_NAME='data/spice/sclk/vg{-sc%10}00046.tsc'
    LSK_FILE_NAME='data/spice/lsk/naif0012.tls'
    EULER_ANGLE_UNITS='DEGREES'
    EULER_ROTATIONS_ORDER=('Z','Y','Z')
    EULER_ROTATIONS_TYPE='SPACE'
    ANGULAR_RATE_FRAME='REFERENCE'
    """
    with open(f'data/scratch/v{-sc%10}{planet[0].lower()}_msopck_header.txt', 'w') as ouf:
        print(msopck_txt, file=ouf)
    # The inputs include speed around elevation and azimuth axis,
    # of which one is always zero. The file wants rotation axis and angle,
    # ie a vector with direction parallel to the rotation axis and length
    # equal to the rotation speed in deg/s. In order to do this, we need
    # some trigonometry because the elevation axis rotates with azimuth.
    with open(f'data/scratch/v{-sc%10}{planet[0].lower()}_msopck_data.txt', 'w') as ouf:
        for et0, et1, coel0, coaz0, dcoel, dcoaz in zip(et_seg[:-1], et_seg[1:], coel0_seg[:-1], coaz0_seg[:-1],
                                                        dcoel_seg[:-1], dcoaz_seg[:-1]):
            x=-np.sin(np.radians(coaz0))*dcoel
            y= np.cos(np.radians(coaz0))*dcoel
            print(f"{et0:25.15e} {et1:25.15e} 0 {coel0:25.15e} {coaz0:25.15e} {x:25.15e} {y:25.15e} {dcoaz:25.15e}",file=ouf)
    try:
        os.remove(f"data/spice/ck/v{-sc%10}{planet[0].lower()}_slew.bc")
    except FileNotFoundError:
        pass  # no error, file is already not present
    subprocess.call(f"~/bin/msopck data/scratch/v{-sc%10}{planet[0].lower()}_msopck_header.txt data/scratch/v{-sc%10}{planet[0].lower()}_msopck_data.txt data/spice/ck/v{-sc%10}{planet[0].lower()}_slew.bc",
                shell=True)


def main():
    encounters={
        "V1J":("vg1_jup_version1_type1_iss_sedr.bc",-31,"Jupiter"),
        "V2N":("vg2_nep_version1_type1_iss_sedr.bc",-32,"Neptune"),
        "V2U":("vg2_ura_version1_type1_iss_sedr.bc",-32,"Uranus"),
    }
    encounter="V1J"
    ScanPlatformCK = "data/spice/ck/"+encounters[encounter][0]
    furnish_kernels(1,ScanPlatformCK)
    ticks,sclkstrs,ets=get_pointing_instances(ScanPlatformCK,sclkid=encounters[encounter][1])
    ets,twists,coels,coazs,dts=get_eulers(sc=encounters[encounter][1],
                                          ets=ets,
                                          ticks=ticks,
                                          sclkstrs=sclkstrs,
                                          et0=str2et("1979-03-01 00:00:00 TDB"),
                                          et1=str2et("1979-03-07 00:00:00 TDB"))
    et_seg,dt_seg,coel0_seg,coaz0_seg,dcoel_seg,dcoaz_seg=make_slews(ets,coels,coazs)
    et_seg,dt_seg,coel0_seg,coaz0_seg,dcoel_seg,dcoaz_seg=merge_slews(et_seg,dt_seg,coel0_seg,coaz0_seg,dcoel_seg,dcoaz_seg)
    msopck(et_seg,dt_seg,coel0_seg,coaz0_seg,dcoel_seg,dcoaz_seg,planet=encounters[encounter][2],sc=encounters[encounter][1])
    plt.show()


if __name__=="__main__":
    main()
