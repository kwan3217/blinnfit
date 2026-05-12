#furnsh "data/spice/vgr2.tm"

#include "VoyagerSimple.inc"

cylinder {0,x*100,0.1 pigment {color rgb x}}
cylinder {0,y*100,0.1 pigment {color rgb y}}
cylinder {0,z*100,0.1 pigment {color rgb z}}
//cylinder {-z*2,-z*3,144.93/2/39.37 pigment {color rgb <1,1,1>}}


#declare ET=str2et("1986-01-22 20:45:00")+clock*600;
#declare RefFrame="J2000";
union {
  Voyager(ET)
  PrintQuat("VG2_SC_BUS: ",pxform("VG2_SC_BUS",RefFrame,ET))
 // QuatTrans(pxform("VG2_SC_BUS",RefFrame,ET),<0,0,0>)
}

light_source {
  <-20,20,-20>*1000
  color rgb <1,1,1>
}

camera {
  up y
  right -x*4/3
  sky -z
  angle 45
  location <-300,00,0>*1.5/39.37
  look_at <0,0,0>
}
