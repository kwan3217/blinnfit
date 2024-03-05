#declare CelestialSphereRad=8e5;
//#declare StarRatio=1440*45/CamAngle;
#declare LimitMag=6.0;
#declare LimitStars=4000;
//#declare ConstellSelect="Sgr";
#include "StarsRight.inc"
#include "SpiceQuat.inc"

#declare RefFrame="ECLIPB1950";
#declare ET=0;

//From VoyagerNeptuneB, frame 2000
#declare CamLat=-8.2484919199613;
#declare CamLon=290.743486161355;
#declare CamAngle=47.9465502489589;
#declare CamTwist=-0.00687845413119703;
#declare RightDenom=2.91142742556547;
#declare MasterSky=z;


camera {
  #declare LookAt=LLR2XYZ(radians(CamLat),radians(CamLon),1);
  PrintVector("LookAt: ",LookAt)
  up y
  right -x*4/RightDenom
  #declare SSky=vnormalize(vcross(LookAt,MasterSky));
  #declare CSky=vnormalize(vcross(SSky,LookAt));
  sky cos(radians(CamTwist))*CSky+sin(radians(CamTwist))*SSky
  angle CamAngle
  location <0,0,0>
  look_at LookAt
}

object {
  Stars
  QuatTrans(pxform("J2000",RefFrame,ET),<0,0,0>)
}

#declare

object {
  Voyager
}

