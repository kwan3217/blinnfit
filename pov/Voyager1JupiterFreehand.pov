#include "KwanMath.inc"
#include "Freehand.inc"

// Furnish the kernels
#furnsh "data/spice/vgr1.tm"
#furnsh "data/spice/lsk/naif0012.tls"
#furnsh "data/spice/pck/pck00011.tpc"
#furnsh "data/spice/pck/jupiter_system2.tpc"

// Use time of closest approach as one possible time reference
#declare CalCA="1979-03-05 12:05:26 TDB";  // from Voyager Hyperbolic Elements
#declare ET_CA=str2et(CalCA);
PrintNumber("ETCA_Cal: ",ET_CA)

#declare RubberClock=array[12][2] {
{1238,	-657299160.0}
{1346,	-657292812.046332}
{1387,	-657290822.175032}
{1388,	-657290343.397684}
{1389,	-657290284.620335}
{1390,	-657290225.842986}
{1391,	-657290167.065637}
{1701,	-657272606.672883}
{2701,	-657220199.519541}
{3756,	-657159019.772765}
{4515,	-657115048.349255}
{4797,	-657098669.14052}
}

#declare JupiterPos=array[29][3] {
{ 330,	218.086650516856,	229.128544196091},
{ 430,	207.991758097495,	223.355790534115},
{ 530,	199.897673147538,	214.631523301627},
{ 630,	197.925496001986,	207.942872090036},
{ 670,	199.580115082222,	206.025989560958},
{ 700,	194.313484215603,	201.72075979344 },
{ 717,	185.640994657917,	198.365917743208},
{ 800,	178.05830242027 ,	191.462338009167},
{ 894,	195.999814421613,	194.361074695563},
{ 900,	196.543810731952,	194.39902701721 },
{ 979,	213.307043783177,	197.025811555916},
{1000,	217.421302205295,	198.885076811364},
{1100,	237.885670711657,	203.281969384618},
{1200,	258.499711844441,	206.641237640298},
{1238,	268.172951591134,	208.74939357576 },
{1300,	279.278662184253,	209.835532426236},
{1346,	289.821192681776,	213.05068735074 },
{1387,	297.306934827983,	214.165522399942},
{1400,	299.053536481189,	215.09862107381 },
{1500,	321.712009715949,	220.01519533549 },
{1600,	340.127919393118,	222.037790289735},
{1700,	353.636735817922,	227.268348575425},
{1701,	359.11554427203 ,	226.675356969349},
{1800,	371.40406732944 ,	234.33982446882 },
{1900,	374.125967649316,	233.957801093284},
{2000,	338.56836422701 ,	252.934236465958},
{3756,	366.960729754548,	208.896668393555},
{4515,	356.986465711134,	213.815450059114},
{4797,	333.791771516755,	226.936599745126}
}

#declare VGR=1;
#declare ET=Arrayterp(RubberClock,0,1,frame_number);
PrintNumber(concat("ET (",etcal(ET),"): "),ET)

PlanetMod(599,UnivFrame,-30-VGR,ET,256-60)
PlanetMod(501,UnivFrame,-30-VGR,ET,0)
PlanetMod(502,UnivFrame,-30-VGR,ET,0)
PlanetMod(503,UnivFrame,-30-VGR,ET,0)
PlanetMod(504,UnivFrame,-30-VGR,ET,0)
#if(image_width>1000)
PlanetMod(505,UnivFrame,-30-VGR,ET,0)
#else
sphere {
  0,1
  pigment {color rgb <1,0,1>}
  scale 1000
  translate spkezr("505",ET,UnivFrame,"NONE",str(-30-VGR,0,0))
}
#end
PlanetMod( 10,UnivFrame,-30-VGR,ET,0)

#declare xp=Arrayterp(JupiterPos,0,1,frame_number);
#declare yp=Arrayterp(JupiterPos,0,2,frame_number);
#declare xd=Linterp(0,-0.5,640,0.5,xp);
#declare yd=Linterp(0,-0.5,480,0.5,yp);

#declare UniverseFrame="ECLIPB1950";
#declare JupiterPos=spkezr("599",ET,UniverseFrame,"LT+S",str(-30-VGR,0,0));
#declare IoPos=spkezr("501",ET,UniverseFrame,"LT+S",str(-30-VGR,0,0));
#declare GanymedePos=spkezr("503",ET,UniverseFrame,"LT+S",str(-30-VGR,0,0));
#declare CallistoPos=spkezr("504",ET,UniverseFrame,"LT+S",str(-30-VGR,0,0));
#switch(frame_number)
  #range(0,2185)
    Stage(<0,0,0>,JupiterPos,z,4/3,18,xd,yd)
    #declare DrawVoy=1;
    #break
  #range(2186,2279)
    // Do Stage first, so that we have a comparison
    //#debug "Stage\n"
    //Stage(<0,0,0>,JupiterPos,z,4/3,18,xd,yd)
    // Then do doublestage
    //#debug "DoubleStage\n"
    DoubleStage(<0,0,0>,
                JupiterPos,xd,yd,
                IoPos,xd,yd,
                Linterp(2186,0,2279,1,frame_number),
                z,4/3,18)
    #declare DrawVoy=1;
    #break
  #range(2280,2872)
    Stage(<0,0,0>,IoPos,z,4/3,18,xd,yd)
    #declare DrawVoy=0;
    #break
  #range(2873,2966)
    DoubleStage(<0,0,0>,
                IoPos,xd,yd,
                GanymedePos,xd,yd,
                Linterp(2873,0,2966,1,frame_number),
                z,4/3,18)
    #declare DrawVoy=1;
    #break
  #range(2967,3372)
    Stage(<0,0,0>,GanymedePos,z,4/3,18,xd,yd)
    #declare DrawVoy=0;
    #break
  #range(3373,3466)
    DoubleStage(<0,0,0>,
                GanymedePos,xd,yd,
                JupiterPos,xd,yd,
                Linterp(3373,0,3466,1,frame_number),
                z,4/3,18)
    #declare DrawVoy=1;
    #break
  #range(3466,3780)
    Stage(<0,0,0>,JupiterPos,z,4/3,18,xd,yd)
    #declare DrawVoy=0;
    #break
  #range(3781,3873)
    DoubleStage(<0,0,0>,
                JupiterPos,xd,yd,
                CallistoPos,xd,yd,
                Linterp(3781,0,3873,1,frame_number),
                z,4/3,18)
    #declare DrawVoy=1;
    #break
  #range(3784,4279)
    Stage(<0,0,0>,CallistoPos,z,4/3,18,xd,yd)
    #declare DrawVoy=0;
    #break
  #range(4280,4371)
    DoubleStage(<0,0,0>,
                CallistoPos,xd,yd,
                JupiterPos,xd,yd,
                Linterp(4280,0,4371,1,frame_number),
                z,4/3,18)
    #declare DrawVoy=1;
    #break
  #range(4371,9999)
    Stage(<0,0,0>,JupiterPos,z,4/3,18,xd,yd)
    #declare DrawVoy=0;
    #break
#end

union {
  union {
    text {
      ttf "data/fonts/UbuntuMono-R.ttf" timout(ET,"YYYY-MM-DD HR:MN:SC.###::TDB TDB") 0,0
      scale <1,-1,1>*0.025
      translate  y*0.5
    }
    text {
      ttf "data/fonts/UbuntuMono-R.ttf" timout(ET,"YYYY-MM-DD HR:MN:SC.###::UTC UTC") 0,0
      scale <1,-1,1>*0.025
      translate  y*0.475
    }
    text {
      ttf "data/fonts/UbuntuMono-R.ttf" concat("Frame ",str(frame_number,4,0)) 0,0
      scale <1,-1,1>*0.025
      translate  y*0.425
    }
    #if((ET-ET_CA)<0)
      #declare DT=ET_CA-ET;
      #declare Utime="J-";
    #else
      #declare DT=ET-ET_CA;
      #declare Utime="J+";
    #end
    #declare HR=floor(DT/3600);
    #declare MN=floor(mod(DT,3600)/60);
    #declare SC=mod(DT,60);
    #declare Utime=concat(Utime,str(HR,-2,0),":",str(MN,-2,0),":",str(SC,-6,3)," (",str(DT,10,3),")");
    text {
      ttf "data/fonts/UbuntuMono-R.ttf" Utime 0,0
      scale <1,-1,1>*0.025
      translate  y*0.45
    }
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

camera {
  up -y
  right CamRight
  location Location
  look_at LookAt
  sky CamSky
  angle CamAngle
}

#declare SunPos=spkezr("10",ET,UniverseFrame,"LT+S",str(-30-VGR,0,0));

light_source {
  SunPos
  color rgb <1,1,1>
}

#declare StarRatio=1440*45/CamAngle;
#declare CelestialSphereRad=1e7;
#include "StarsRight.inc"
object {
  Stars
  QuatTrans(pxform("J2000",UnivFrame,ET),<0,0,0>)
}

#if(DrawVoy)
  #include "VoyagerSimple.inc"
  #furnsh "data/spice/ck/vg1_jup_version1_type1_iss_sedr.bc"
  #declare VoyRight=StageRight*CamRight;
  #declare VoyDown=StageDown;
  #declare VoyOut=StageOut*DirectionScale;
  #declare VoyagerPos=vnormalize(VoyOut+VoyRight*0.25+VoyDown*0.25)*50;
  union {
    Voyager(VGR,ET)
    #declare ToCameraFrame=pxform(concat("VG",str(VGR,0,0),"_SC_BUS"),UnivFrame,ET);
    QuatTrans(ToCameraFrame,VoyagerPos)
  }
#end

PrintNumber("frame_number: ",frame_number)

