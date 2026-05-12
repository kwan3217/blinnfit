#include "KwanMath.inc"

// Furnish the kernels
#furnsh "data/spice/vgr2.tm"
#furnsh "data/spice/lsk/naif0012.tls"

// Use time of closest approach as one possible time reference
#declare CalCA="1986-01-24 17:58:51 UTC";
#declare scsCA="3/26847:50:451";
#declare ETCA_Cal=str2et(CalCA);
PrintNumber("ETCA_Cal: ",ETCA_Cal)
#declare ETCA_scs=scs2e(-32,scsCA);
PrintNumber("ETCA_scs: ",ETCA_scs)
PrintNumber("Delta: ",ETCA_Cal-ETCA_scs)

// Angle keyframe table. First we just make the table where we mark
// whether each corner is a corner or is eased.
#declare AngleTable=array[2][6] {
//Frame Angle Ease Body cx cy scr scx scy
/*
  {   0, 97.8,   0,799,}, //Show all orbits
  { 100, 97.8,   0}, //Start zoom in to standard zoom
  { 220, 59.0, -20}, //Finish zoom to standard
  { 401, 59.0,   0}, //Start zoom into SigSag
  { 520,  6.8, -90}, //Finish zoom in
  { 910,  6.8,  90}, //Start zoom out
  {1033, 59.0,   0}, //Finish zoom out
  {1292, 59.0,   0}, //Start zoom into 704, VOBEST
  {1380,  2.3, -10}, //Finish zoom in
  {1555,  2.3,  10}, //Start zoom out
  {1622, 59.0,   0}, //Finish zoom out
  {2240, 59.0,   0}, //Start zoom into 702, VUBEST
  {2325,  2.3,   0}, //Finish zoom in
  {2480,  2.3,  20}, //Start zoom out
  {2570, 59.0,   0}, //Finish zoom out
  {2988, 59.0,   0}, //Start zoom into 703, VTBEST
  {3075,  2.3, -10}, //Finish zoom in
  {3230,  2.3,   0}, //Start zoom out
  {3318, 59.0,   0}, //Finish zoom out
  {3740, 59.0,   0}, //Start zoom into 701, VABEST
  {3830,  2.3, -15}, //Finish zoom in
  {4040,  2.3,  10}, //Start zoom out
  {4100, 59.0,   0}, //Finish zoom out
  {4190, 59.0,   0}, //Start zoom into 705, VMBEST
  {4240,  2.3, -15}, //Finish zoom in. Note that distance to Miranda changes significantly over sequence
  {4578,  2.3,  15}, //Start zoom out
  {4682, 59.0,   0}, //Finish zoom out */
  {4790, 59.0,   0,799,686,426} //
  {6906, 59.0,   0,799,}  //Constant to the end of the video
}


// Times of FOV images. First number is body to project on, second is start time, third is stop time
#declare VOBEST_ETm=FET[1382]; // Time that zoom-in finishes
#declare VOBEST_ET1=-439784344.977951; // Time that VOBEST frames disappear
//At least VOBEST is off on Blinn's video, by about an hour. By the time the video
//zooms in on Oberon, the observations have already been taken. Rather than
//re-doing the cinematography, we will shift when each exposure is shown. The
//geometry of each exposure will be from when it was actually taken, but
//it will be shown when it was shown in Blinn's work.
#declare NAExpEt=array[4][7] {
//body  actual                      start                  end         R G B
//VOBEST
  {704,scs2e(-32,"3/26836:22:768"),ETCA_Cal-30093.689+  0,VOBEST_ET1,  1,1,1}, //C2683623, clear filter
  {704,scs2e(-32,"3/26836:24:736"),ETCA_Cal-30093.689+ 60,VOBEST_ET1,  1,1,1}, //C2683625, clear filter
  {704,scs2e(-32,"3/26836:26:736"),ETCA_Cal-30093.689+120,VOBEST_ET1,0.5,0,1}, //C2683627, violet filter
  {704,scs2e(-32,"3/26836:28:736"),ETCA_Cal-30093.689+180,VOBEST_ET1,  0,1,0}, //C2683629, green filter
}
//Print the start and end times for each row
#declare I=0;
#while(I<dimension_size(NAExpEt,1))
  PrintNumber("Exposure ",I)
  PrintNumber("  ET: ",NAExpEt[I][1])
  #debug concat("  TDB: ",timout(NAExpEt[I][1],"YYYY-MM-DD HR:MN:SC.### ::TDB"),"\n")
  #debug concat("  UTC: ",timout(NAExpEt[I][1],"YYYY-MM-DD HR:MN:SC.### ::UTC"),"\n")
  #declare I=I+1;
#end
PrintNumber("  VOBEST ET1: ",VOBEST_ET1)
#debug concat("  TDB: ",timout(VOBEST_ET1,"YYYY-MM-DD HR:MN:SC.### ::TDB"),"\n")
#debug concat("  UTC: ",timout(VOBEST_ET1,"YYYY-MM-DD HR:MN:SC.### ::UTC"),"\n")
#declare ET=FET[frame_number];
PrintNumber("  ET: ",ET)
#debug concat("  TDB: ",timout(ET,"YYYY-MM-DD HR:MN:SC.### ::TDB"),"\n")
#debug concat("  UTC: ",timout(ET,"YYYY-MM-DD HR:MN:SC.### ::UTC"),"\n")

#include "VoyagerScene.inc"
