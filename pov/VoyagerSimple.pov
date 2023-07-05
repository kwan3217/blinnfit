#include "VoyagerSimple.inc"

object {
  VoyagerBody
}

object {
  ScanPlatform
  //transform
  //rotate -z*90
  translate ScanPivot
}

cylinder {0,x*100,0.1 pigment {color rgb x}}
cylinder {0,y*100,0.1 pigment {color rgb y}}
cylinder {0,z*100,0.1 pigment {color rgb z}}

light_source {
  <20,20,-20>*1000
  color rgb <1,1,1>
}

camera {
  up y
  right -x*4/3
  sky y
  angle 15
  location <300,500,-100>*1.5
  look_at ScanPivot
}
