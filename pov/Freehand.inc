// Clean sheet Voyager
#include "KwanMath.inc"
#include "SpiceQuat.inc"
#include "slerp.inc"

#declare BodyNames=array[800];
#declare BodyNames[ 10]="Sun";
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

// Default Universe frame. Use this in your own code, and override with care if you
// want to change it. Recommended to override it first, before any macro calls.
#declare UnivFrame="ECLIPB1950";

// Name is now looked up from BodyNames so not passed in
#macro PlanetMod(SpiceID,UnivFrame,Center,ET,ExtraRot)
  // :param SpiceId: Spice ID of planet in question
  // :param UnivFrame: Frame used for the universe
  // :param Center: Object which is at the origin of the universe frame
  // :param ET: Ephemeris time for body
  #local Name=BodyNames[SpiceID];
  #debug concat("Body: ",Name,"\n")
  #local found=0;
  #local PlanetA=gdpool(concat("BODY",str(SpiceID,0,0),"_RADII"),0,found);
  PrintNumber("PlanetA: ",PlanetA)
  #local PlanetB=gdpool(concat("BODY",str(SpiceID,0,0),"_RADII"),1,found);
  PrintNumber("PlanetB: ",PlanetB)
  #local PlanetC=gdpool(concat("BODY",str(SpiceID,0,0),"_RADII"),2,found);
  PrintNumber("PlanetC: ",PlanetC)
  #declare SpiceIDStr=str(SpiceID,3,0);
  #declare BodyFrameStr=concat("IAU_",strupr(Name));
  #local R=spkezr(SpiceIDStr,ET,UnivFrame,"NONE",str(Center,3,0));
  #local Q=pxform(BodyFrameStr,UnivFrame,ET);
  union {
    sphere {
      0,1
      scale <PlanetA,PlanetB,PlanetC>
      #if(SpiceID=10)
        #debug "Sun special model\n"
        no_shadow
        texture {
          pigment {color rgb <1,1,0>}
          finish {ambient 1 diffuse 0}
        }
        // Due to astronomical scaling, scale down the model and distance by a factor of 100
        #local R=R/100;
        scale 0.01
      #else
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
      #end
      rotate z*ExtraRot
    }
    PrintVector("R: ",R)
    PrintQuat("Q: ",Q)
    QuatTrans(Q,R)
  }
#end

#macro LocalStage(CameraRu,ActorRu, Sky, RightScale, Angle, x_d,y_d)
/*
Calculate *a* stage.
:param CameraRu: Offset from POV-Ray origin to camera, often zero. TODO - implement nonzero
:param ActorRu: Spice ID of background actor that is at given position in camera. This is often a planet
                     or moon in
:param Sky: Sky vector in POV-Ray(=SpiceFrame) frame
:param RightScale: Length of right vector. Right vector will be x*Right, and up vector will be -y. Setting Right
              equal to rendering aspect ratio will result in square pixels
:param Angle: horizontal full FOV angle, deg
:return: An array of vectors:
  [0] - StageRight - unit vector camera +x in universe frame
  [1] - StageDown  - unit vector camera +y in universe frame
  [2] - StageOut   - unit vector camera +z in universe frame
  [3] - Location   - location of camera in universe. Intended to be used as `location Location` in camera statement. UNTESTED for Location!=<0,0,0>
  [4] - LookAt     - LookAt vector for camera in universe. Intended to be used as `look_at LookAt` in camera statement.
  [5] - CamSky     - sky vector for camera in universe. Intended to be used as `sky CamSky` in camera statement.
  [6] - CamRight   - Right vector for camera in universe, indended to be used as `up -y right CamRight` in camera statement.
  [7] - CamAngle (in result[6].x) - Horizontal field of view in degrees, intended to be used as `angle CamAngle.x` in camera statement.
  [8] - DirectionScale (in result[7].x) - Length of direction vector, used for HUD scaling
*/
  #local Result=array[9];
  // Calculate direction from angle and aspect ratio
  #local DirectionScale=0.5*RightScale/tan(radians(Angle)/2);
  #local Result[8]=<DirectionScale,0,0>;
  #local CamAngle=Angle;
  #local Result[7]=<CamAngle,0,0>;
  #local CamRight=x*RightScale;
  #local Result[6]=CamRight;
  PrintNumber("DirectionScale: ",DirectionScale)
  // Create actor vector in stage frame
  #local z_d=1;
  #local r_d=<x_d,y_d,z_d>;
  PrintVector("r_d: ",r_d)
  // Transform actor vector into camera orthonomral frame
  #local r_n=r_d*<RightScale,1,DirectionScale>;
  PrintVector("r_n: ",r_n)
  // Implement Point-Toward. We have R=[p_r, s_r, u_r] in reference (universe) space, and
  // BT=[p_b^T;s_b^T;u_b^T] in body (camera orthonormal) space.
  #local Pr=vnormalize(ActorRu-CameraRu);
  PrintVector("Pr: ",Pr)
  #local Pb=vnormalize(r_n);
  PrintVector("Pb: ",Pb)
  #local Tr=vnormalize(SkyVec);
  PrintVector("Tr: ",Tr)
  #local Tb=-y;
  PrintVector("Tb: ",Tb)
  #local Sr=vnormalize(vcross(Pr,Tr));
  PrintVector("Sr: ",Sr)
  #local Sb=vnormalize(vcross(Pb,Tb));
  PrintVector("Sb: ",Sb)
  #local Ur=vcross(Pr,Sr);
  PrintVector("Ur: ",Ur)
  #local Ub=vcross(Pb,Sb);
  PrintVector("Ub: ",Ub)
  #local R=array[3][3] {{Pr.x,Sr.x,Ur.x},
                        {Pr.y,Sr.y,Ur.y},
                        {Pr.z,Sr.z,Ur.z}};
  PrintMatrix("R: ",R)
  #local BT=array[3][3] {{Pb.x,Pb.y,Pb.z},
                         {Sb.x,Sb.y,Sb.z},
                         {Ub.x,Ub.y,Ub.z}};
  PrintMatrix("BT: ",BT)
  #local M_rb=MxM(R,BT)
  PrintMatrix("M_rb: ",M_rb)
  // Stage vectors in camera universe space. These are appropriate for placing foreground actors. They are all unit
  // length so you will have to take into account z distance, direction, and aspect to place objects.
  #local StageRight=MxV(M_rb,x);
  #local Result[0]=StageRight;
  #local StageDown=MxV(M_rb,y);
  #local Result[1]=StageDown;
  #local StageOut=MxV(M_rb,z);
  #local Result[2]=StageOut;
  // Assign camera variables
  #local Location=CameraRu;
  #local Result[3]=Location;
  #local LookAt=Location+StageOut;
  #local Result[4]=LookAt;
  #local CamSky=vnormalize(SkyVec);
  #local Result[5]=CamSky;
  Result
#end

#macro Stage(CameraRu,ActorRu,SkyVec,RightScale,Angle,x_d,y_d)
/*
Stage the background actor.

This function sets:
  * (/) StageRight (unit vector in POV-Ray frame)
  * (/) StageDown (unit vector in POV-Ray frame)
  * (/) StageOut (unit vector in POV-Ray frame)
  * (/) Location (vector, position in POV-Ray frame)
  * (/) LookAt (vector, position in POV-Ray frame)
  * (/) CamSky (unit vector in POV-Ray frame)
  * (/) CamRight (vector in POV-Ray frame)
  * (/) CamAngle (scalar, deg full horizontal FOV)
  * (/) DirectionScale (scalar, relative direction length)
  * (/) HudMatrix (transformation)
These are used to create the camera in the Camera() macro, and to place the HUD
in the HUD() macro, so this should be called first (or at least some routine that
sets these vars)
:param CameraRu: Offset from POV-Ray origin to camera, often zero. TODO - implement nonzero
:param ActorRu: Spice ID of background actor that is at given position in camera. This is often a planet
                     or moon in
:param Sky: Sky vector in POV-Ray(=SpiceFrame) frame
:param Right: Length of right vector. Right vector will be x*Right, and up vector will be -y. Setting Right
              equal to rendering aspect ratio will result in square pixels
:param Angle: horizontal full FOV angle, deg
:param x_d: X display coordinate of actor, -0.5=left edge of screen, +0.5=right edge
:param y_d: Y display coordinate of actor, -0.5=top edge of screen, +0.5=bottom edge
*/
  #local ThisLocalStage=LocalStage(CameraRu,ActorRu, SkyVec, RightScale, Angle, x_d,y_d);
  /*
  [0] - StageRight
  [1] - StageDown
  [2] - StageOut
  [3] - Location
  [4] - LookAt
  [5] - CamSky
  [6] - CamRight
  [7] - CamAngle (in result[7].x)
  [8] - DirectionScale (in result[8].x)
  [9] - HUDRight
*/
  // Calculate direction from angle and aspect ratio
  #declare DirectionScale=ThisLocalStage[8].x;
  #declare CamAngle=ThisLocalStage[7].x;
  #declare CamRight=ThisLocalStage[6];
  #declare StageRight=ThisLocalStage[0];
  #declare StageDown=ThisLocalStage[1];
  #declare StageOut=ThisLocalStage[2];
  PrintVector(" StageRight: ",StageRight)
  PrintVector(" StageDown:  ",StageDown )
  PrintVector(" StageOut:   ",StageOut  )
  PrintMatrix("M: ",M3b1b(StageRight,StageDown,StageOut))
  // Assign camera variables
  #declare Location=ThisLocalStage[3];
  #declare LookAt=ThisLocalStage[4];
  #declare CamSky=ThisLocalStage[5];
  // Transformation appropriate for HUD. This properly projects things in the plane z=DirectionScale
  // to the screen. The top edge is at y=-0.5, bottom at y=+0.5, left at x=-CamRight*0.5, right at x=+CamRight*0.5
  // Note that this is uniform in the XY plane, so it won't distort text etc.
  // Note that POV-Ray convention looks like a transpose of my convention so
  //  M=[ M00 M01 M02 Tx]
  //    [ M10 M11 M12 Ty]
  //    [ M20 M21 M22 Tz]
  //    [   0   0   0  1]
  // is represented as:
  // transform{matrix <M00,M10,M20,
  //                   M01,M11,M21,
  //                   M02,M12,M22,
  //                   Tx ,Ty ,Tz  >}
  #declare HudMatrix=transform{
    MatrixVerbose("      Stage HudMatrix: ",
                  StageRight.x,StageRight.y,StageRight.z,
                  StageDown.x,StageDown.y,StageDown.z
                  StageOut.x*DirectionScale,StageOut.y*DirectionScale,StageOut.z*DirectionScale,
                  Location.x,Location.y,Location.z)
  };
#end

#macro DoubleStage(CameraRu, //Position of camera, currently required to be origin of universe space
       ActorRu0,x0_d,y0_d, //position of primary background actor in primary stage
       ActorRu1,x1_d,y1_d, //position of secondary background actor in secondary stage
       T, //interpolation parameter. When t==0, camera is aligned with primary stage. When t==1, camera is aligned to secondary stage
       SkyVec,Right,Angle //sky vector, aspect ratio scalar, and horizontal field of view angle scalar in deg
)
  // Primary stage
  #local Stage0=LocalStage(CameraRu,ActorRu0, SkyVec, Right, Angle, x0_d,y0_d);
  #local Stage0Right=Stage0[0];
  #local Stage0Down =Stage0[1];
  #local Stage0Out  =Stage0[2];
  #local M0=M3b1b(Stage0Right,Stage0Down,Stage0Out);
  PrintVector("Stage0Right: ",Stage0Right)
  PrintVector("Stage0Down:  ",Stage0Down )
  PrintVector("Stage0Out:   ",Stage0Out  )
  PrintMatrix("M0: ",M0)
  // Secondary stage
  #local Stage1=LocalStage(CameraRu,ActorRu1, SkyVec, Right, Angle, x1_d,y1_d);
  #local Stage1Right=Stage1[0];
  #local Stage1Down =Stage1[1];
  #local Stage1Out  =Stage1[2];
  #local M1=M3b1b(Stage1Right,Stage1Down,Stage1Out);
  PrintMatrix("M1: ",M1)
  // inTerpolated stage (T for parameter $t$)
  #local MT=Slerp(M0,M1,T);
  PrintNumber("T:  ",T)
  PrintMatrix("MT: ",MT)
  #local StageTRight=<MT[0][0],MT[1][0],MT[2][0]>;
  #local StageTDown =<MT[0][1],MT[1][1],MT[2][1]>;
  #local StageTOut  =<MT[0][2],MT[1][2],MT[2][2]>;
  PrintVector("Stage0Right: ",Stage0Right)
  PrintVector("StageTRight: ",StageTRight)
  PrintVector("Stage0Down:  ",Stage0Down)
  PrintVector("StageTDown:  ",StageTDown)
  PrintVector("Stage0Out:   ",Stage0Out)
  PrintVector("StageTOut:   ",StageTOut)
  // Camera vectors, follows interpolation
  #declare Location=CameraRu;
  #declare LookAt=Location+StageTOut;
  #declare CamSky=vnormalize(-StageTDown);
  // Foreground Stage vectors, follows stage 0
  #declare StageRight=Stage0Right;
  #declare StageDown =Stage0Down;
  #declare StageOut  =Stage0Out;
  // Camera parameters, constant between frames so use stage 0. Note that angle etc might change over time,
  // but will not change between stage 0 and stage 1. Any angle etc change will occur to both stages
  // in lockstep, so it's only specified once.
  #declare DirectionScale=Stage0[8].x;
  #declare CamAngle=Stage0[7].x;
  #declare CamRight=Stage0[6];
  // HUD Matrix, follow interpolation
  #declare HudMatrix=transform{
    MatrixVerbose("DoubleStage HudMatrix: ",
            StageTRight.x,StageTRight.y,StageTRight.z,
            StageTDown.x,StageTDown.y,StageTDown.z
            StageTOut.x*DirectionScale,StageTOut.y*DirectionScale,StageTOut.z*DirectionScale,
            Location.x,Location.y,Location.z)
  };
#end


