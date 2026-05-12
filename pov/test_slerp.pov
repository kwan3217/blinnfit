#include "KwanMath.inc"
#include "slerp.inc"

#include "slerp_ref.inc"

#macro Get(I)
  array[3][3] {{TestPoints[I][0][0],TestPoints[I][0][1],TestPoints[I][0][2]},
               {TestPoints[I][1][0],TestPoints[I][1][1],TestPoints[I][1][2]},
               {TestPoints[I][2][0],TestPoints[I][2][1],TestPoints[I][2][2]}}
#end

#declare I=0;
#declare N=dimension_size(TestPoints,1);
#declare M0=Get(0);
PrintMatrix("M0: ",M0)
#declare M1=Get(N-1);
PrintMatrix("M1: ",M1)
#while(I<N)
  #declare T=I/(N-1);
  PrintNumber("T: ",T)
  #declare Mt=Slerp(M0,M1,T);
  PrintMatrix("Mt: ",Mt)
  PrintMatrix("ref:",Get(I))
  #declare mGetI=MxS(Get(I),-1)
  PrintMatrix("dif:",MpM(Mt,mGetI))
  #declare I=I+1;
#end