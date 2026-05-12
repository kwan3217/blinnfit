// Clean sheet Voyager
#include "KwanMath.inc"
#include "SpiceQuat.inc"
#include "VoyagerSimple.inc"

#declare BodyNames=array[800];
#declare BodyNames[599]="Jupiter";
#declare BodyNames[501]="Io";
#declare BodyNames[502]="Europa";
#declare BodyNames[503]="Ganymede";
#declare BodyNames[504]="Callisto";
#declare BodyNames[505]="Amalthea";
#declare BodyNames[799]="Uranus";
#declare BodyNames[701]="Ariel";
#declare BodyNames[702]="Umbriel";
#declare BodyNames[703]="Titania";
#declare BodyNames[704]="Oberon";
#declare BodyNames[705]="Miranda";

#declare UnivFrame="ECLIPB1950";

// Name is now looked up from BodyNames so not passed in
#macro PlanetMod(SpiceID,UnivFrame,Center,ET)
  // :param SpiceId: Spice ID of planet in question
  // :param UnivFrame: Frame used for the universe
  // :param Center: Object which is at the origin of the universe frame
  // :param ET: Ephemeris time for body
  #local Name=BodyNames[SpiceID];
  #local found=0;
  #local PlanetA=gdpool(concat("BODY",str(SpiceID,3,0),"_RADII"),0,found);
  #local PlanetB=gdpool(concat("BODY",str(SpiceID,3,0),"_RADII"),1,found);
  #local PlanetC=gdpool(concat("BODY",str(SpiceID,3,0),"_RADII"),2,found);
  #declare SpiceIDStr=str(SpiceID,3,0);
  #declare BodyFrameStr=concat("IAU_",strupr(Name));
  #local R=spkezr(SpiceIDStr,ET,UnivFrame,"NONE",str(Center,3,0));
  PrintVector("R: ",R)
  #local Q=pxform(BodyFrameStr,UnivFrame,ET);
  PrintQuat("Q: ",Q)
  union {
    sphere {
      0,1
      scale <PlanetA,PlanetB,PlanetC>
      texture {
        pigment {color rgb <1,0,1>} // Intentionally ugly color to show where we have no map
        finish {ambient 0 diffuse 1}
      }
      #if(file_exists(concat("data/textures/",Name,"Map.png")))
      texture {
        pigment {
          image_map {
            // Maps are transparent where there is no data
            png concat("data/textures/",Name,"Map.png")
            map_type spheroid
            flatness 1-PlanetC/PlanetA
          }
          scale <-1,1,1>
          rotate x*90
          scale PlanetA
        }
      }
      #end
    }
    QuatTrans(Q,R)
  }
#end

PrintNumber("frame_number: ",frame_number)
// We no longer look at FET -- define ET yourself in the calling file
//#declare ET=FET[frame_number];
#debug concat(etcal(ET),"\n")

#include "fov.inc"

#if(Planet=599)
  PlanetMod(599,UnivFrame,-30-VGR,ET)
  PlanetMod(501,UnivFrame,-30-VGR,ET)
  PlanetMod(502,UnivFrame,-30-VGR,ET)
  PlanetMod(503,UnivFrame,-30-VGR,ET)
  PlanetMod(504,UnivFrame,-30-VGR,ET)
  PlanetMod(505,UnivFrame,-30-VGR,ET)
#end

#if(Planet=799)
  PlanetMod(799,"Uranus" ,UnivFrame,-32,ET)
  PlanetMod(701,"Ariel"  ,UnivFrame,-32,ET)
  PlanetMod(702,"Umbriel",UnivFrame,-32,ET)
  PlanetMod(703,"Titania",UnivFrame,-32,ET)
  FOVBody(-32,"Narrow",704,ET,NACFovSize,0.03*NACFovSize,NAExpEt,UnivFrame)
  FOVBody(-32,"Narrow",705,ET,NACFovSize,0.03*NACFovSize,NAExpEt,UnivFrame)
#end

#declare ToSun=spkezr("10",ET,UnivFrame,"None",str(-30-VGR,0,0));
PrintVector("ToSun: ",ToSun)
#declare DirToSun=vnormalize(ToSun);
PrintVector("DirToSun: ",DirToSun)

#declare Location=<0,0,0>;
//FLat and FLon no longer an array, look it up in calling code
#ifdef(FLat)
#declare LookAt=LLR2XYZ(radians(FLat[frame_number]),radians(FLon[frame_number]),1);
#end
#ifdef(FTwist)
#declare SSky=vnormalize(vcross(LookAt,z));
#declare CSky=vnormalize(vcross(SSky,LookAt));
#declare Sky=cos(radians(FTwist[frame_number]))*CSky+sin(radians(FTwist[frame_number]))*SSky;
#end
#ifdef(FRight)
#declare CamRight=FRight[frame_number];
#declare CamAngle=FAngle[frame_number];
#end

// Used for calculating stage directions, not directly in camera.
#declare Direction=0.5*CamRight/tan(radians(CamAngle)/2);

#declare StageOut=-vnormalize(Location-LookAt);
#declare StageRight=vnormalize(vcross(StageOut,Sky));
#declare StageDown=vnormalize(vcross(StageOut,StageRight));
PrintVector(concat("Fpixc_sc[",str(frame_number,4,0),"]: "),Fpixc_sc[frame_number])
#declare VoyagerOut=Direction;
#declare VoyagerRight=Linterp(0,-0.5,NominalImageWidth,0.5,Fpixc_sc[frame_number].x)*CamRight;
#declare VoyagerDown=Linterp(0,-0.5,NominalImageHeight,0.5,Fpixc_sc[frame_number].y);
PrintVector("VoyagerPosCam: ",<VoyagerRight,VoyagerDown,VoyagerOut>)
#declare VoyagerPos=(VoyagerOut  *StageOut
                    +VoyagerRight*StageRight
                    +VoyagerDown *StageDown)
                    *1300/Fpixc_sc[frame_number].z;
PrintVector("VoyagerPos: ",VoyagerPos)

camera {
  up y
  right -x*CamRight
  sky Sky
  location Location
  look_at LookAt
  angle CamAngle
}

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
    #if((ET-str2et("1986-01-24 17:58:51 UTC"))<0)
      #declare DT=str2et("1986-01-24 17:58:51 UTC")-ET;
      #declare Utime="U-";
    #else
      #declare DT=ET-str2et("1986-01-24 17:58:51 UTC");
      #declare Utime="U+";
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
    text {
      ttf "data/fonts/UbuntuMono-R.ttf" concat("Frame ",str(frame_number,4,0)) 0,0
      scale <1,-1,1>*0.025
      translate  y*0.425
    }
    translate z
    translate -x*0.5*CamRight
  }
  matrix <StageRight.x,StageRight.y,StageRight.z,
          StageDown.x,StageDown.y,StageDown.z
          StageOut.x*Direction,StageOut.y*Direction,StageOut.z*Direction,
          Location.x,Location.y,Location.z>
  pigment {color rgb <1,1,1>}
  finish {ambient 1 diffuse 0}
}

union {
  Voyager(ET)
  #declare ToCameraFrame=pxform("VG2_SC_BUS",UnivFrame,ET);
  QuatTrans(ToCameraFrame,VoyagerPos)
}


light_source {
  ToSun
  color rgb <1,1,1>
}

#declare StarRatio=1440*45/CamAngle;
#include "StarsRight.inc"
object {
  Stars
  QuatTrans(pxform("J2000",UnivFrame,ET),<0,0,0>)
}

PrintNumber("frame_number: ",frame_number)
#debug concat(etcal(ET)," ",timout(ET,"YYYY-MM-DD HR:MN:SC.### ::UTC"),"\n")

