#declare NACAngle=0.424;
#declare CelestialSphereRad=8e5;
#declare StarRatio=1440*45/NACAngle;
#declare LimitMag=6.0;
#declare LimitStars=10000;
#include "StarsRight.inc"
#include "SpiceQuat.inc"
#declare CameraFrame="VG2_ISSNA";

#declare RefFrame="ECLIPB1950";
#declare ET=0;
#furnsh "data/spice/fk/vg2_v02.tf"
#furnsh "data/spice/lsk/naif0012.tls"
#furnsh "data/spice/sclk/vg200045.tsc"
#furnsh "data/spice/spk/vgr2_nep097.bsp"
#declare Type1=0;
#if(Type1)
  #furnsh "data/spice/ck/vg2_nep_version1_type1_iss_sedr.bc"
  #declare NWnd=ckcov("data/spice/ck/vg2_nep_version1_type1_iss_sedr.bc",-32100,"INTERVAL",0,"TDB",30000);
  PrintNumber("NWnd: ",NWnd)
  /*
  #local IWnd=0;
  #while(IWnd<NWnd)
    #local Timout=timout(ckgetcov(IWnd,0),"YYYY-MM-DDTHR:MN:SC.######Z::UTC");
    #debug concat(str(IWnd,-5,0)," ",Timout,"\n")
    #local IWnd=IWnd+1;
  #end
  */
  #declare ET=ckgetcov(19282,0);
#else
  #furnsh "data/spice/ck/v2n_slew.bc"
  #furnsh "data/spice/ck/vgr2_super.bc"
  #declare ET=-326758504.439144;
#end
PrintNumber("ET: ",ET)
#debug concat(etcal(ET),"\n")

#declare ToCameraFrame=pxform("J2000",CameraFrame,ET);
PrintQuat("ToCameraFrame: ",ToCameraFrame)

#declare vel=<0,0,0>;
#declare pos=spkezr("801",ET,"J2000","NONE","-32",vel);
PrintVector("Pos: ",pos)
PrintVector("Vel: ",vel)

object {
  Stars
  QuatTrans(ToCameraFrame,<0,0,0>)
}

sphere {
  pos,1352.6
  QuatTrans(ToCameraFrame,<0,0,0>)
  pigment {color rgb <1,1,1>}
  finish {ambient 1}
}

camera {
  up y
  right -x
  sky y
  location <0,0,0>
  look_at z
  angle NACAngle
}

