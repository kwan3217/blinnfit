#include "pov/KwanMath.inc"
#include "pov/SpiceQuat.inc"
#declare NACFovSize=0.424;

// Furnish the kernels
#furnsh "data/spice/vgr2.tm"

// Time of first FOV for each sequence for each body. This defines the
// time that the plane through the center is defined. If exactly zero,
// don't draw the plane
#declare FOVEt=array[800];
#declare FOVEt[704]=scs2e(-32,"3/26836:22:768");
#declare FOVEt[701]=0;
#declare FOVEt[702]=0;
#declare FOVEt[703]=0;
#declare FOVEt[705]=0;
#declare FOVEt[799]=0;

#declare BodyNames=array[800];
#declare BodyNames[799]="Uranus";
#declare BodyNames[701]="Ariel";
#declare BodyNames[702]="Umbriel";
#declare BodyNames[703]="Titania";
#declare BodyNames[704]="Oberon";
#declare BodyNames[705]="Miranda";

// Times of FOV images. First number is body to project on, second is start time, third is stop time
#declare VOBEST_ETm=-439784344.977951; // Time that zoom-in finishes
#declare VOBEST_ET1=-439784344.977951; // Time that VOBEST frames disappear
#declare NAExpEt=array[4][6] {
//body  start                    end            R G B
//VOBEST
  {704,scs2e(-32,"3/26836:22:768"),VOBEST_ET1,  1,1,1}, //C2683623, clear filter
  {704,scs2e(-32,"3/26836:24:736"),VOBEST_ET1,  1,1,1}, //C2683625, clear filter
  {704,scs2e(-32,"3/26836:26:736"),VOBEST_ET1,0.5,0,1}, //C2683627, violet filter
  {704,scs2e(-32,"3/26836:28:736"),VOBEST_ET1,  0,1,0}, //C2683629, green filter
}

#include "pov/fov.inc"

// Position of spacecraft in target body frame at moment of first image
// in VxBEST sequence. This is the time at which the plane catching
// the off-body FOV is created.

#declare SCBFPos0=spkezr("-32",VOBEST_ET1-1,"IAU_OBERON","NONE","704");
FOVBody(-32,"Narrow",704,VOBEST_ET1-1,NACFovSize,0.03*NACFovSize,NAExpEt)

camera {
  //location (y*5-z*2)*PlanetRp
  location SCBFPos0
  look_at <0,0,0>
  sky z
  right -x*image_width/image_height
  angle 2.3
}

light_source {
  spkezr("10",FOVEt[704],"IAU_OBERON","NONE","704")
  color rgb <1,1,1>
}

