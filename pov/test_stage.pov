#include "Freehand.inc"

Stage(<0,0,0>,<10,0,0>,-z,4/3,2*degrees(atan(2/3)),0.25,0)

camera {
  up -y
  right CamRight
  location Location
  look_at LookAt
  sky CamSky
  angle CamAngle
}

sphere {
  0,1
  pigment {color rgb <1,0,0>}
  finish {ambient 1 diffuse 0}
  translate <10,0,0>
}

union {
  text {
    ttf "data/fonts/UbuntuMono-R.ttf" concat("Frame ",str(frame_number,4,0)) 0,0
    scale <1,-1,1>*0.025
    translate  y*0.425
    translate -0.5*CamRight
    pigment {color rgb <1,1,1>}
    finish {ambient 1 diffuse 0}
  }
  cylinder {
    -0.5*CamRight-0.5*y,0.5*CamRight-0.5*y,0.001
    pigment {color rgb <1,0,1>}
    finish {ambient 1 diffuse 0}
  }
  cylinder {
    -0.5*CamRight+0.5*y,0.5*CamRight+0.5*y,0.001
    pigment {color rgb <0,1,0>}
    finish {ambient 1 diffuse 0}
  }
  cylinder {
    -0.5*CamRight-0.5*y,-0.5*CamRight+0.5*y,0.001
    pigment {color rgb <0,1,1>}
    finish {ambient 1 diffuse 0}
  }
  cylinder {
    +0.5*CamRight-0.5*y,+0.5*CamRight+0.5*y,0.001
    pigment {color rgb <1,0,0>}
    finish {ambient 1 diffuse 0}
  }
  translate z
  transform {HudMatrix}
  no_shadow
}
