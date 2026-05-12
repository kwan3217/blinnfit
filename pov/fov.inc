#macro Pyramid(Angle)
  // Infinite square pyramid stretching from origin along +z axis
  intersection {
    plane { x,0 rotate  y*Angle/2}
    plane {-x,0 rotate -y*Angle/2}
    plane { y,0 rotate -x*Angle/2}
    plane {-y,0 rotate  x*Angle/2}
  }
#end

#macro HollowPyramid(Angle,DAngle,ExtraRot)
  difference {
    Pyramid(Angle+DAngle/2)
    Pyramid(Angle-DAngle/2)
    rotate ExtraRot
  }
#end

#macro FOVtex(InTex,OutTex)
  texture {
    object {
      FOVs
      texture {OutTex}
      texture {InTex}
    }
  }
#end

#macro FOVat(SCSpiceID,Cam,BodySpiceID,ET,FOVSize,LineWidth)
  /***
  :param SCSpiceID:     Spice ID of spacecraft (int)
  :param Cam:           Narrow 'N' or Wide 'W' angle camera
  :param BodySpiceID:   Spice ID of target body (int)
  :param ET:            Time of FOV
  :param FOVSize:       Angular field of view size
  :param LineWidth:     Angular field of view line width
  ***/
  #local VgrId=mod(-SCSpiceID,10);
  #local CamFrame=concat("VG",str(VgrId,1,0),"_ISS",strupr(substr(Cam,1,1)),"A");
  #local BodyFrame=concat("IAU_",strupr(BodyNames[BodySpiceID]));
  #local Q=pxform(CamFrame,BodyFrame,ET);
  #local SCBFPos=spkezr(str(SCSpiceID,3,0),ET,BodyFrame,"NONE",str(BodySpiceID,3,0));
  #debug "Field of view cast"
  PrintNumber("  ET: ",ET)
  #debug concat("  TDB: ",timout(ET,"YYYY-MM-DD HR:MN:SC.### ::TDB"),"\n")
  #debug concat("  UTC: ",timout(ET,"YYYY-MM-DD HR:MN:SC.### ::UTC"),"\n")
  PrintQuat  ("  Q:       ",Q)
  PrintVector("  SCBFPos: ",SCBFPos)
  union {
    HollowPyramid(FOVSize,LineWidth,<0,-0.1,0>)
    QuatTrans(Q,SCBFPos)
  }
#end

#declare BodyCol=0;
#declare ActualCol=1;
#declare StartCol=2;
#declare EndCol=3;
#declare RCol=4;
#declare GCol=5;
#declare BCol=6;

#macro FOVBody(SCSpiceID,Cam,BodySpiceID,ET,FOVSize,LineWidth,FOVTable,UnivFrame)
  /***
  Draw a body with several FOVs projected onto it

  :param SCSpiceID:     Spice ID of spacecraft (int)
  :param Cam:           Narrow 'N' or Wide 'W' angle camera
  :param BodySpiceID:   Spice ID of target body (int)
  :param ET:            Time of FOV
  :param FOVSize:       Angular field of view size
  :param LineWidth:     Angular field of view line width
  :param FOVTable:      2D array, one row for each exposure. Exposure consists of elements:
    * Body to project onto
    * Start ET when exposure is visible
    * Stop ET when exposure is no longer visible
    * Red, green, and blue color components for this FOV (not currently used)
  :param UnivFrame:     Spice name of Universe frame
  ***/
  // Body stuff that we are using regardless
  #local found=0;
  #local BodyA=gdpool(concat("BODY",str(BodySpiceID,3,0),"_RADII"),0,found);
  #local BodyB=gdpool(concat("BODY",str(BodySpiceID,3,0),"_RADII"),1,found);
  #local BodyRe=(BodyA+BodyB)/2;
  #local BodyRp=gdpool(concat("BODY",str(BodySpiceID,3,0),"_RADII"),2,found);
  #local SCSpiceIDStr=str(SCSpiceID,3,0);
  #local BodyIDStr=str(BodySpiceID,3,0);
  #local BodyName=BodyNames[BodySpiceID];
  #local VgrId=mod(-SCSpiceID,10);
  #local CamFrame=concat("VG",str(VgrId,1,0),"_ISS",strupr(substr(Cam,1,1)),"A");
  #local BodyFrame=concat("IAU_",strupr(BodyNames[BodySpiceID]));

  // Figure time range of FOVs
  #local ET0= 9e99; //Effectively +infinity
  #local ET1=-9e99; //Effectively -infinity
  #local I=0;
  #local N=dimension_size(FOVTable,1);
  PrintNumber("N: ",N)
  #local Found0=0;
  #local Found1=0;
  #while(I<N)
    #if(BodySpiceID=FOVTable[I][BodyCol])
      #if(FOVTable[I][StartCol]<ET0)
        #local ET0=FOVTable[I][StartCol];
        #local Found0=1;
      #end
      #if(FOVTable[I][EndCol]>ET1)
        #local ET1=FOVTable[I][EndCol];
        #local Found1=1;
      #end
    #end
    #local I=I+1;
  #end
  PrintNumber(concat("ET0 (",timout(ET0,"YYYY-MM-DD HR:MN:SC.### ::TDB"),"): "),ET0)
  PrintNumber(concat("ET1 (",timout(ET1,"YYYY-MM-DD HR:MN:SC.### ::TDB"),"): "),ET1)
  #local SCBFPos0=spkezr(SCSpiceIDStr,ET0,BodyFrame,"NONE",BodyIDStr);
  // Background plane
  #local BackgroundPlane=plane {vnormalize(SCBFPos0),0}
  // Body shape
  #local BodyShape=sphere {
    0,1
    scale <BodyRe,BodyRp,BodyRp>
  }

  //Define layered planet texture
  #local BodyTex=texture {
    pigment {color rgb <1,0,1>} // Intentionally ugly color to show where we have no map
    finish {ambient 0 diffuse 1}
  } texture {
    pigment {
      image_map {
        // Maps are transparent where there is no data
        png concat("data/textures/",BodyName,"Map.png")
        map_type spheroid
        flatness 1-BodyRp/BodyRe
      }
      scale <-1,1,1>
      rotate x*90
      scale BodyRe
    }
  }

  union {
    // Draw everything in the body-centered body-fixed frame then rotate and translate it to the universe frame
    // Determine if we are drawing __any__ FOVs
    #if(ET>=ET0 & ET<=ET1)
      // Yes, we are
      #local FOVs=union {
        #local I=0;
        #while(I<N)
          #local DrawI=true; // Draw unless we decide not to
          PrintNumber(concat("FOVTable[",str(I,0,0),"][Start] (",timout(FOVTable[I][StartCol],"YYYY-MM-DD HR:MN:SC.### ::TDB"),"): "),FOVTable[I][StartCol])
          PrintNumber(concat("  ET (",timout(ET,"YYYY-MM-DD HR:MN:SC.### ::TDB"),"): "),ET)
          PrintNumber(concat("FOVTable[",str(I,0,0),"][End  ] (",timout(FOVTable[I][  EndCol],"YYYY-MM-DD HR:MN:SC.### ::TDB"),"): "),FOVTable[I][  EndCol])
          #if(BodySpiceID!=FOVTable[I][BodyCol])
            #debug concat("Wrong body -- we are on ",str(BodySpiceID,0,0),", but this row is for ",str(FOVTable[I][BodyCol],0,0),"\n")
            #local DrawID=false;
          #end
          #if(FOVTable[I][StartCol]>ET)
            #debug concat("Not time yet\n")
            #local DrawI=false;
          #end
          #if(FOVTable[I][  EndCol]<ET)
            #debug concat("Past time\n")
            #local DrawI=false;
          #end
          #if(DrawI)
            // Draw the FOV
            FOVat(SCSpiceID,Cam,BodySpiceID,FOVTable[I][ActualCol],FOVSize,LineWidth)
          #end
          #local I=I+1;
        #end
      }
      object {
        BackgroundPlane
        texture {
          object{FOVs
          // Blank space (not on the FOV lines)
          texture{pigment{color rgbt <1,1,1,1>}}
          // White (on the lines)
          texture{pigment{color rgb <1,1,1>}finish{ambient 1 diffuse 0}}
          }
        }
      }
      object {
        BodyShape
        texture {
          object {BackgroundPlane
            // Front side of plane (not "inside")
            texture {
              object {FOVs
                // Body surface (not on the FOV lines)
                texture{BodyTex}
                // White (on the lines)
                texture{pigment{color rgb <1,1,1>}finish{ambient 1 diffuse 0}}
              }
            }
            // Back side of plane
            texture {BodyTex}
          }
        }
      }
    #else
      // Just draw the body
      object {
        BodyShape
        texture{BodyTex}
      }
    #end
    // Now that we've drawn everything, rotate it to the universe frame
    // and translate it to it's spacecraft-relative position
    #local R=spkezr(BodyIDStr,ET,UnivFrame,"NONE",SCSpiceIDStr);
    PrintVector("FOVBody R: ",R)
    #debug concat("BodyFrame: ",BodyFrame,"\n")
    #debug concat("UnivFrame: ",UnivFrame,"\n")
    #local Q=pxform(BodyFrame,UnivFrame,ET);
    PrintQuat("FOVBody Q: ",Q)
    QuatTrans(Q,R)
  }
#end