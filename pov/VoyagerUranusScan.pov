// Always point at the spacecraft so we can see where the scan platform is pointing

#include "data/scratch/pointing_VoyagerUranusHD.inc"
#include "data/scratch/patch_pointing_VoyagerUranusHD_VOBEST.inc"
#include "KwanMath.inc"
#declare Planet=799;
#declare NACFovSize=0.424;

#declare FCopy=1048;
#declare FLat[frame_number]=FLat[FCopy];
#declare FLon[frame_number]=FLon[FCopy];
#declare FAngle[frame_number]=20;
#declare FTwist[frame_number]=FTwist[FCopy];
#declare FRight[frame_number]=FRight[FCopy];
#declare Fpixc_sc[frame_number]=<640,360,200>;


// Furnish the kernels
#furnsh "data/spice/vgr2.tm"


// Time of closest approach. For Voyager 2 at Uranus, we have both a UTC time and a SCLK time
#declare CalCA="1986-01-24 17:58:51 UTC";
#declare scsCA="3/26847:50:451";
#declare ETCA_Cal=str2et(CalCA);
PrintNumber("ETCA_Cal: ",ETCA_Cal)
#declare ETCA_scs=scs2e(-32,scsCA);
PrintNumber("ETCA_scs: ",ETCA_scs)
PrintNumber("Delta: ",ETCA_Cal-ETCA_scs)


#declare VOBEST_ETm=FET[1382]; // Time that zoom-in finishes
#declare VOBEST_ET1=-439784344.977951; // Time that VOBEST frames disappear
#declare NAExpEt=array[4][7] {
//body  actual                      start                  end         R G B
//VOBEST
  {704,scs2e(-32,"3/26836:22:768"),ETCA_Cal-30093.689+  0,VOBEST_ET1,  1,1,1}, //C2683623, clear filter
  {704,scs2e(-32,"3/26836:24:736"),ETCA_Cal-30093.689+ 60,VOBEST_ET1,  1,1,1}, //C2683625, clear filter
  {704,scs2e(-32,"3/26836:26:736"),ETCA_Cal-30093.689+120,VOBEST_ET1,0.5,0,1}, //C2683627, violet filter
  {704,scs2e(-32,"3/26836:28:736"),ETCA_Cal-30093.689+180,VOBEST_ET1,  0,1,0}, //C2683629, green filter
}


#include "VoyagerScene.inc"
