/* Test the Great Red Spot location.

The Project Pluto site references a calculator on Sky and Telescope web site at
https://skyandtelescope.org/observing/interactive-sky-watching-tools/transit-times-of-jupiters-great-red-spot/
This is perfect. It states that the spot had a definite location (71degW) on 2025-03-01, with a drift rate of
xdeg/month. It can calculate the transit time for any date, so we might as well choose 2025-03-01 so we don't
have to worry about the drift rate changing. From that we get transit times of:

Enter a date (mm/dd/yyyy):
03/01/2025
Universal Times of Red Spot transits
centered on date:
03/01/2025 @ 05:52 UT

03/01/2025 @ 15:47 UT

03/02/2025 @ 01:43 UT
Corresponding local dates & times of Red Spot transits:
02/28/2025 @ 11:52 pm

03/01/2025 @ 09:47 am

03/01/2025 @ 07:43 pm
Note: Local times are based on the time zone offset of from UT given by your Web browser. Hour(s):
-6

Note that this is definitely as-observed from Earth, so light time must definitely be included.
*/

#furnsh "data/spice/spk/de440s.bsp"
#furnsh "data/spice/spk/jup230l.bsp"
#furnsh "data/spice/lsk/naif0012.tls"
#furnsh "data/spice/pck/pck00011.tpc"
#furnsh "data/spice/pck/jupiter_system2.tpc"
#include "pov/KwanMath.inc"
#include "pov/SpiceQuat.inc"
#declare ET=str2et("2025-03-01 05:52:00 UTC")+clock*86400;
PrintNumber("ET: ", ET)
// In all of these, LT+S means where does the target appear to be when observed at given time,
// so position is where the planet "was" LT in the past. Aim a telescope at this point to see it.
#declare JupiterPos=spkezr("599",ET,"J2000","LT+S","399"); // Pos of Jupiter
#declare SunPos=spkezr("10",ET,"J2000","LT+S","399");
#declare c=299792.458; //km/s
#declare LT=vlength(JupiterPos)/c;
PrintVector("JupiterPos: ",JupiterPos)
PrintNumber("LT:         ",LT)

#declare PlanetA=71516;
#declare PlanetB=PlanetA;
#declare PlanetC=66871;
#local Q=pxform("IAU_JUPITER","J2000",ET-LT); //pxform in this version makes a quaternion compatible with QuatTrans

sphere {
  0,1
  scale <PlanetA,PlanetB,PlanetC>
  pigment {
    image_map {
      png "data/textures/JupiterMap.png"
      map_type spheroid
      flatness 1-PlanetC/PlanetA
    }
    scale <-1,1,1>
    scale PlanetA
    rotate x*90
    rotate z*256 //This puts the center of the Great Red Spot on the System II prime meridian
    rotate -z*71 //This puts the center at 71degW, correct for 2025-03-01
  }
  QuatTrans(Q,JupiterPos) // transform with given quaternion, then translate so origin moves to given position
}

camera {
  up y
  right -4/3*x
  sky z
  location <0,0,0>
  look_at JupiterPos
  angle 0.02
}

light_source {
  SunPos
  color rgb <1,1,1>
}