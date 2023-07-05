#declare NSides=10;
#declare FA=360/NSides;
#declare BusR=38;
#declare BusH=17;

#declare Decagon=prism {
   linear_spline
   0,1,10
   <1,0>
   #local I=1;
   #while(I<NSides)
     ,<cos(radians(I*FA)),sin(radians(I*FA))>
     #local I=I+1;
   #end
   rotate x*90
   rotate z*18
}

#declare Louver=union {
  #declare LouverD=1;
  #declare LouverInsetL=3;
  #declare LouverInsetR=LouverInsetL;
  #declare LouverInsetT=1;
  #declare LouverInsetB=LouverInsetT;
  #declare BladeN=8;
  #declare LouverR= sin(radians(FA/2))*BusR-LouverInsetR;
  #declare LouverL=-sin(radians(FA/2))*BusR+LouverInsetL;
  #declare LouverW=LouverR-LouverL;
  #declare LouverFrame=0.5;
  #declare LouverX=cos(radians(FA/2))*BusR;
  #declare BladePitch=(LouverW-(2*LouverFrame))/BladeN;
  #declare BladeOpen=30;
  box {
    <LouverX        ,LouverL,-     LouverInsetT>,
    <LouverX+LouverD,LouverR,-BusH+LouverInsetB>
    pigment {color rgb <1,1,1>}
  }
  #declare I_Blade=0;
  #while(I_Blade<BladeN)
    box {
      <0  ,0         ,-LouverInsetT-LouverFrame>,
      <0.1,BladePitch,-BusH+LouverInsetB+LouverFrame>
      rotate -z*BladeOpen
      translate x*(LouverX+LouverD)
      translate y*(LouverL+LouverFrame+BladePitch* I_Blade)
      pigment {color rgb <1,1,1>}
    }
    #declare I_Blade=I_Blade+1;
  #end
}


#macro Tubes(V,E,R)
  #local I=0;
  #while(I<dimension_size(E,1))
    cylinder {
      V[E[I].x],V[E[I].y],R
    }
    #local I=I+1;
  #end
#end

#declare LaunchMount=union {
  #local LaunchMountV=array[8] {
    < 18, 25,0>,< 18,-25,0>,
    <-18,-25,0>,<-18, 25,0>,
    <  0, 20,44.5>,< 20,  0,44.5>,
    <  0,-20,44.5>,<-20,  0,44.5>,
  }
  #local LaunchMountE=array[8] {
    <0,5>,<1,5>,
    <1,6>,<2,6>,
    <2,7>,<3,7>,
    <3,4>,<0,4>
  }
  Tubes(LaunchMountV,LaunchMountE,1)
  pigment {color rgb <0.5,0.5,0.5>}
}

#declare RTGMount=union {
  #local RTGMountV=array[7] {
    < 21,-28,0>, //0, bottom outboard +x
    <-21,-28,0>, //1, bottom outboard -x
    < 21,-28,-BusH>, //2, top outboard +x
    <-21,-28,-BusH>, //3, top outboard -x
    <0,-BusR,0>,   //4,bottom inboard
    < 11,-60,-24>,  //5,pivot +x
    <-11,-60,-24>,  //6,pivot +x
  }
  #local RTGMountE=array[6] {
    <0,5>,
    <2,5>,
    <4,5>,
    <1,6>,
    <3,6>,
    <4,6>
  }
  Tubes(RTGMountV,RTGMountE,0.5)
  pigment {color rgb <1,1,1>}
}

#declare RTGTruss=union {
  #local RTGTrussV=array[8] {
    < 11, -60,-24>,  //0,pivot +x
    <-11, -60,-24>,  //1,pivot -x
    <  0, -84,-13>,  //2 field joint x0
    <  8, -84,-18>,  //3 field joint +x
    < -8, -84,-18>,  //4 field joint -x
    <  0,-104,-11>,  //5 RTG x0
    <  3,-104,  0>,  //6 RTG +x
    < -3,-104,  0>,  //7 RTG -x
  }
  #local RTGTrussE=array[14] {
    <0,2>,
    <0,3>,
    <1,2>,
    <1,4>,
    <1,3>,
    <2,3>,
    <3,4>,
    <4,2>,
    <2,6>,
    <2,7>,
    <3,5>,
    <3,6>,
    <4,5>,
    <4,7>
  }
  Tubes(RTGTrussV,RTGTrussE,0.5)
  pigment {color rgb <1,1,1>}
}

#declare RTG=union {
  cylinder {
    <0,0,0>,<0,-20,0>,6
  }
  #local I=0;
  #while(I<3)
    box {
      <-0.25,0,-8>,<0.25,-20,8>
      rotate y*120*I
    }
  #local I=I+1;
  #end
  pigment {color rgb <0.5,0.5,0.5>}
}

#declare RTGs=union {
  object {RTGTruss}
  #local I=0;
  #while(I<3)
    object {RTG translate -y*I*24 rotate x*5 translate <0,-103,-5>}
    #local I=I+1;
  #end
}

#declare GoldRecord=cylinder {
  <0,0,0>,<-3,0,0>,6
  translate -x*cos(radians(FA/2))*BusR
  translate +y*1
  translate -z*9.5
  pigment {color rgb <1,1,0>}
}


#declare MainBus=difference {
  object {Decagon
    scale <BusR,BusR,BusH>
    translate -z*BusH
  }
  object {Decagon
    scale <BusR-8,BusR-8,BusH+2>
    translate -z*(BusH+1)
  }
  translate y*1.12 //bus offset
  pigment {color rgb <0.5,0.5,0.5>}
}

#declare MagBoom=union {
  #local MagBoomLength=13*39.97;
  #local MagBoomBays=60;
  #local MagBayLength=MagBoomLength/MagBoomBays;
  #local MagBoomR=5;
  #local WhipR=0.2;
  cylinder {
    <0,0,0>,<0,-16,0>,MagBoomR
    pigment {color rgb <1,1,1>}
  }
  #local I=0;
  #while(I<3)
    cylinder {
      <0,-16,MagBoomR-WhipR>,<0,-16-MagBoomLength,MagBoomR-WhipR>,WhipR
      rotate y*(I*120+40)
      pigment {color rgb <1.0,0.75,0>}
    }
    #local I=I+1;
  #end
  #local I=0;
  #while(I<MagBoomBays)
    union {
      Tubes(array[6] {
        <sin(radians(0*120+40))*(MagBoomR-WhipR),0,cos(radians(0*120+40))*(MagBoomR-WhipR)>,
        <sin(radians(1*120+40))*(MagBoomR-WhipR),0,cos(radians(1*120+40))*(MagBoomR-WhipR)>,
        <sin(radians(2*120+40))*(MagBoomR-WhipR),0,cos(radians(2*120+40))*(MagBoomR-WhipR)>,
        <sin(radians(0*120+40))*(MagBoomR-WhipR),-MagBayLength,cos(radians(0*120+40))*(MagBoomR-WhipR)>,
        <sin(radians(1*120+40))*(MagBoomR-WhipR),-MagBayLength,cos(radians(1*120+40))*(MagBoomR-WhipR)>,
        <sin(radians(2*120+40))*(MagBoomR-WhipR),-MagBayLength,cos(radians(2*120+40))*(MagBoomR-WhipR)>,
      },array[6] {
        <0,4>,
        <1,5>,
        <2,3>,
        <3,4>,
        <4,5>,
        <5,3>
      },WhipR/2)
      pigment {color rgb <1,1,1>}
      translate -y*(16+I*MagBayLength)
    }
    #local I=I+1;
  #end
  rotate x*40
  translate <0,-52,-28>
}

#declare HGADish=union {
  //Measurements are technically of the outside surface
  #local HGA_R=144.93/2; //numerator is diameter, exact value from blueprint
  #local HGA_D=27; //depth measured from blueprint
  #local HGA_c=HGA_D/(HGA_R*HGA_R);
  #local CollarR=35;
  #local CollarD=HGA_c*CollarR*CollarR;
  quadric {
    <HGA_c,HGA_c,0>,<0,0,0>,<0,0,-1>,0
    clipped_by {
      cylinder {
        <0,0,0>,<0,0,HGA_D>,HGA_R
      }
    }
    scale <1,1,-1>
    translate -z*36
  }
  cylinder {
    <0,0,-33>,<0,0,-36-CollarD>,CollarR
    open
  }
  pigment {color rgb <1,1,1>}
}

#declare HGASupport=union {
  Tubes(array[9] {
    <-21,-28,-BusH>,// 0, connecting point between bay 7 and 8
    < 21,-28,-BusH>,// 1, connecting point between bay 9 and 10
    <-28, 10,-BusH>,// 2, connecting point between bay 5 and 6
    <-17,25,-BusH>, // 3, between bay 4 and 5
    < 28, 10,-BusH>,// 4, connecting point between bay 1 and 2
    < 17,25,-BusH>, // 5, between bay 2 and 3
    < 0,-34,-33>, //6, dish connection between bay 8 and 9
    <-29,17,-33>, //7, over bay 5,
    < 29,17,-33>, //7, over bay 2,
  },array[6] {
    <0,6>,
    <1,6>,
    <2,7>,
    <3,7>,
    <4,8>,
    <5,8>
  },0.5)
  pigment {color rgb <0.5,0.5,0.5>}
}

#declare HGAFeed=union {
  Tubes(array[6] {
    <13*cos(radians(0*120+90)),13*sin(radians(0*120+90)),-87>, //0
    <13*cos(radians(1*120+90)),13*sin(radians(1*120+90)),-87>, //1
    <13*cos(radians(2*120+90)),13*sin(radians(2*120+90)),-87>, //2
    <36*cos(radians(0*120+30)),36*sin(radians(0*120+30)),-36-HGA_c*36*36>, //3
    <36*cos(radians(1*120+30)),36*sin(radians(1*120+30)),-36-HGA_c*36*36>, //4
    <36*cos(radians(2*120+30)),36*sin(radians(2*120+30)),-36-HGA_c*36*36>, //5
  },array[6] {
    <3,0>,
    <3,2>,
    <4,0>,
    <4,1>,
    <5,1>,
    <5,2>
  },0.5)
  cone {
    <0,0,-87.5>,15,<0,0,-80>,0
  }
  cone {
    <0,0,-87.5>,4,<0,0,-97>,1.5
    pigment {color rgb <1,1,0>}
  }
  pigment {color rgb <1,1,1>}
}

#declare PWA=union {
  box {<-3,-3,-3>,<3,3,3> pigment {color rgb <1,1,1>}}
  cylinder {0,-y*10*39.37,0.2 pigment {color rgb <1,0.5,0>}}
  cylinder {0,-x*10*39.37,0.2 pigment {color rgb <1,0.5,0>}}
  rotate z*45
  rotate -x*50
  translate <0,-40,-20>
}

#declare ScienceTruss=union {
  Tubes(array[5] {
    <-17,25,-BusH>, // 0, between bay 4 and 5
    < 17,25,-BusH>, // 1, between bay 2 and 3
    <  0,37,-BusH>, // 2, between bay 3 and 4
    < 8,45,-43>, //3
    <-8,45,-43>, //4
  },array[5] {
    <0,3>,
    <2,3>,
    <1,4>,
    <2,4>,
    <1,3>
  },0.5)
  pigment {color rgb <0.5,0.5,0.5>}
}

#declare ScienceBoom=union {
  Tubes(array[10] {
    < 8, 0, 0>, //0
    <-8, 0, 0>, //1
    < 4, 5, 4>, //2
    < 4, 5,-4>, //3
    <-4, 5,-4>, //4
    <-4, 5, 4>, //5
    < 4,96, 4>, //6
    < 4,96,-4>, //7
    <-4,96,-4>, //8
    <-4,96, 4>, //9
  },array[12] {
    <0,2>,
    <0,3>,
    <1,4>,
    <1,5>,
    <2,3>,<3,4>,<4,5>,<5,2>,
    <2,6>,<3,7>,<4,8>,<5,9>
  },0.5)
  #local SciBoomNBays=7;
  #local SciBoomBayLen=12;
  #local IBay=0;
  #while(IBay<SciBoomNBays)
    union {
    Tubes(array[8] {
      < 4, 0, 4>, //0
      < 4, 0,-4>, //1
      <-4, 0,-4>, //2
      <-4, 0, 4>, //3
      < 4,SciBoomBayLen, 4>, //4
      < 4,SciBoomBayLen,-4>, //5
      <-4,SciBoomBayLen,-4>, //6
      <-4,SciBoomBayLen, 4>, //7
    },array[8] {
      <0,1>,<1,2>,<2,3>,<3,0>,
      <0,7>,<1,4>,<2,5>,<3,6>
    },0.5)
    translate y*(IBay*SciBoomBayLen+5)
    }
    #local IBay=IBay+1;
  #end
  cylinder {
    <0,96,4.5>,<0,96,-4.5>,4.5
  }
  box {
    <-4.5,89,-4.5>,<4.5,96,4.5>
  }
  cylinder {
    <0,96,-8>,<0,96,26>,2
  }
  pigment {color rgb <0.5,0.5,0.5>}
  rotate -x*7 //deployed position
  //rotate x*80  //stowed position
  translate <0,45,-43>
}

#declare FuelTank=sphere {
  <0,0,-12.5>,14
  pigment {color rgb <1,1,1>}
}

#declare VoyagerBody=union {
  object {MainBus}
  object {Louver rotate -z*( 1-1)*FA}
  object {Louver rotate -z*( 2-1)*FA}
  object {Louver rotate -z*( 5-1)*FA}
  object {Louver rotate -z*( 7-1)*FA}
  object {GoldRecord}
  object {LaunchMount}
  object {RTGMount}
  object {RTGs}
  object {MagBoom}
//  object {HGADish}
  object {HGASupport}
  object {HGAFeed}
  object {PWA}
  object {ScienceTruss}
  object {FuelTank}
  object {ScienceBoom}
  no_shadow
}

#declare ISSNA=union {
  cylinder {
    <0,18,0>,<0,-9,0>,4.5
  }
  cylinder {
    <0,-9,0>,<0,-18,0>,3
  }
  translate -x*13
  pigment {color rgb <0.5,0.5,0.5>}
}

#declare ISSWA=union {
  cylinder {
    <0,15,0>,<0,5,0>,1.5
  }
  cylinder {
    <0,5,0>,<0,-4,0>,3
  }
  translate <-13,0,-9>
  pigment {color rgb <0.5,0.5,0.5>}
}

#declare PPS=union {
  cylinder {
    <0,-2,0>,<0,-10,0>,4
  }
  cylinder {
    <0,-10,0>,<0,-15,0>,3
  }
  cylinder {
    <0,18,0>,<0,-10,0>,4
    open
    clipped_by { plane {y,0 rotate x*70 translate y*5 rotate -y*45}}
  }
  pigment {color rgb <0.5,0.5,0.5>}
  translate <-13,0,9>
}

#declare IRIS=union {
  difference {
    union {
      cylinder {
        <0,11,0>,<0,-8,0>,11
      }
      cone {
        <0,-8,0>,11,<0,-12,0>,3.5
      }
    }
    cylinder {
      <0,12,0>,<0,0.1,0>,10.5
    }
    cylinder {
      <0,1,0>,<0,0,0>,10.4
      pigment {color rgb <1,1,0>}
    }
  }
  cylinder {
    <0,0,0>,<0,14,0>,1.5
  }
  cone {
    <0,14,0>,3.5,<0,17,0>,2.25
  }
  pigment {color rgb <0.5,0.5,0.5>}
  translate <7,0,14>
}

#declare UVS=box {
  <17,10,4>,<22,-10,10>
  pigment {color rgb <0.5,0.5,0.5>}
}

#declare ScanPlatform=union {
  object {ISSNA}
  object {ISSWA}
  object {PPS}
  object {IRIS}
  object {UVS}
  box {
    <-3,8,-16>,<-8,-8,10>
  }
  box {
    <4,6,-7>,<16,-2,3>
  }
  cylinder {
    <-5,0,0>,<22,0,0>,2
  }
  pigment {color rgb <0.5,0.5,0.5>}
  no_shadow
}

#declare ScanPivot=<0,143,-29.5>;
