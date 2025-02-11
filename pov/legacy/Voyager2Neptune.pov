#version 3.7
global_settings{assumed_gamma 1.0}
#include "KwanMath.inc"
#declare RefFrame="ECLIPB1950";
#furnsh "data/spice/fk/vg2_v02.tf"
#furnsh "data/spice/lsk/naif0012.tls"
#furnsh "data/spice/sclk/vg200045.tsc"
#furnsh "data/spice/spk/vgr2_nep097.bsp"
#furnsh "data/spice/ck/v2n_slew.bc"
#furnsh "data/spice/ck/vgr2_super.bc"
#declare ET0=str2et("1989-08-24 16:25:00 UTC");
#declare ET=ET0+86400*clock;
PrintNumber("frame_number: ",frame_number)
PrintNumber("clock: ",clock)
PrintNumber("ET: ",ET)
#debug concat(timout(ET,"YYYY-MM-DD HR:MN:SC.###"),"\n")

#declare NeptuneVel=<0,0,0>;
#declare NeptunePos=spkezr("899",ET,RefFrame,"LT+S","-32",NeptuneVel);
PrintVector("NeptunePos: ",NeptunePos)
PrintVector("NeptuneVel: ",NeptuneVel)
#declare SunVel=<0,0,0>;
#declare SunPos=spkezr("899",ET,RefFrame,"LT+S","-32",SunVel);
PrintVector("SunPos: ",SunPos)
PrintVector("SunVel: ",SunVel)
#declare TritonVel=<0,0,0>;
#declare TritonPos=spkezr("801",ET,RefFrame,"LT+S","-32",TritonVel);
PrintVector("TritonPos: ",TritonPos)
PrintVector("TritonVel: ",TritonVel)

#declare Target=vnormalize(TritonPos+z*290);
PrintVector("Target:     ",Target)
#declare StageRight=vnormalize(vcross(Target,z));
PrintVector("StageRight: ",StageRight)
#declare StageUp=vnormalize(vcross(StageRight,Target));
PrintVector("StageUp: ",StageUp)

#include "VoyagerSimple.inc"
#include "SpiceQuat.inc"

union {
  union {
    object {ScanPlatform rotate x*90 rotate z*90}
    cylinder {0,x* 36,1 pigment{color rgb x} finish {ambient 0.5}}
    cylinder {0,y* 36,1 pigment{color rgb y} finish {ambient 0.5}}
    cylinder {0,z*144,1 pigment{color rgb z} finish {ambient 0.5}}
    QuatTrans(pxform("VG2_ISSNA","VG2_SC_BUS",ET),ScanPivot)
  }
  union {
    cylinder {0,x* 36,1 pigment{color <1,1,1>} finish {ambient 0.5}}
    sphere   {  x* 36,3 pigment{color rgb x  } finish {ambient 0.5}}
    cylinder {0,y* 36,1 pigment{color <1,1,1>} finish {ambient 0.5}}
    sphere   {  y* 36,3 pigment{color rgb y  } finish {ambient 0.5}}
    cylinder {0,z* 36,1 pigment{color <1,1,1>} finish {ambient 0.5}}
    sphere   {  z* 36,3 pigment{color rgb z  } finish {ambient 0.5}}
    QuatTrans(pxform("VG2_AZ_EL","VG2_SC_BUS",ET),ScanPivot)
  }
  object {VoyagerBody}
  scale 0.0254
  QuatTrans(pxform("VG2_SC_BUS",RefFrame,ET),Target*20+StageRight*5+StageUp*2)
}


light_source {
  -vnormalize(NeptunePos)*1000
  color rgb <1,1,1>
}

sphere {
  NeptunePos,25000
  pigment {color rgb <0,0,1>}
}

sphere {
  TritonPos,1250
  pigment {color rgb <1,0,1>}
}


camera {
  up y
  right -x*4/3
  sky z
  location <0,0,0>
  look_at Target
  angle 48.5
}

PrintNumber("frame_number: ",frame_number)
PrintNumber("clock: ",clock)
PrintNumber("ET: ",ET)
#debug concat(timout(ET,"YYYY-MM-DD HR:MN:SC.###"),"\n")
