#include "data/scratch/pointing_VoyagerUranusHD.inc"
#include "data/scratch/patch_pointing_VoyagerUranusHD_VOBEST.inc"
#include "data/scratch/patch_pointing_VoyagerUranusHD_VMBEST.inc"
#include "KwanMath.inc"
#declare Planet=799;
#declare NACFovSize=0.424;

// Furnish the kernels
#furnsh "data/spice/vgr2.tm"

#declare CalCA="1986-01-24 17:58:51 UTC";
#declare scsCA="3/26847:50:451";
#declare ETCA_Cal=str2et(CalCA);
PrintNumber("ETCA_Cal: ",ETCA_Cal)
#declare ETCA_scs=scs2e(-32,scsCA);
PrintNumber("ETCA_scs: ",ETCA_scs)
PrintNumber("Delta: ",ETCA_Cal-ETCA_scs)


// Times of FOV images. First number is body to project on, second is start time, third is stop time
#declare VOBEST_ETm=FET[1382]; // Time that zoom-in finishes
#declare VOBEST_ET1=-439784344.977951; // Time that VOBEST frames disappear
#declare VMBEST_ET0=FET[4310];
#declare VMBEST_dET=FET[4340]-VMBEST_ET0;
#declare VMBEST_ET1=FET[4670];
//At least VOBEST is off on Blinn's video, by about an hour. By the time the video
//zooms in on Oberon, the observations have already been taken. Rather than
//re-doing the cinematography, we will shift when each exposure is shown. The
//geometry of each exposure will be from when it was actually taken, but
//it will be shown when it was shown in Blinn's work.
#declare NAExpEt=array[12][7] {
//body  actual                      start                  end         R G B
//VOBEST
  {704,scs2e(-32,"3/26836:22:768"),ETCA_Cal-30093.689+  0,VOBEST_ET1,  1,1,1}, //C2683623, clear filter
  {704,scs2e(-32,"3/26836:24:736"),ETCA_Cal-30093.689+ 60,VOBEST_ET1,  1,1,1}, //C2683625, clear filter
  {704,scs2e(-32,"3/26836:26:736"),ETCA_Cal-30093.689+120,VOBEST_ET1,0.5,0,1}, //C2683627, violet filter
  {704,scs2e(-32,"3/26836:28:736"),ETCA_Cal-30093.689+180,VOBEST_ET1,  0,1,0}, //C2683629, green filter
//VMBEST
  {705,scs2e(-32,"3/26846:07:768"),VMBEST_ET0+VMBEST_dET*0,VMBEST_ET1,  1,1,1}, //C2684608, clear filter
  {705,scs2e(-32,"3/26846:10:768"),VMBEST_ET0+VMBEST_dET*1,VMBEST_ET1,  1,1,1}, //C2684611, clear filter
  {705,scs2e(-32,"3/26846:13:768"),VMBEST_ET0+VMBEST_dET*2,VMBEST_ET1,  1,1,1}, //C2684614, clear filter
  {705,scs2e(-32,"3/26846:16:768"),VMBEST_ET0+VMBEST_dET*3,VMBEST_ET1,  1,1,1}, //C2684617, clear filter
  {705,scs2e(-32,"3/26846:19:768"),VMBEST_ET0+VMBEST_dET*4,VMBEST_ET1,  1,1,1}, //C2684620, clear filter
  {705,scs2e(-32,"3/26846:22:768"),VMBEST_ET0+VMBEST_dET*5,VMBEST_ET1,  1,1,1}, //C2684623, clear filter
  {705,scs2e(-32,"3/26846:25:768"),VMBEST_ET0+VMBEST_dET*6,VMBEST_ET1,  1,1,1}, //C2684626, clear filter
  {705,scs2e(-32,"3/26846:28:768"),VMBEST_ET0+VMBEST_dET*7,VMBEST_ET1,  1,1,1}, //C2684629, clear filter
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

#ifndef (Fpixc_sc[frame_number])
#declare Fpixc_sc[frame_number]=Linterp(4250,Fpixc_sc[4250],4590,Fpixc_sc[4590],frame_number);
#end

#include "VoyagerScene.inc"
