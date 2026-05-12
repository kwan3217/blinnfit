#include "matrix.inc"

#macro Slerp(M0,M1,T)
/*  Spherical linear interpolation of SO(3) rotations. Most naive translation of the given
    algorithm, with no caching, currying, etc.

    https://jupyter.kwansystems.org/lab/tree/MekwanicalUniverse/Linear%20Algebra/Rotation/Slerp%202.ipynb

    :param M0: Initial orientation, in the form of a matrix transforming camera to universe
               at T=0, itself in the form of a 3x3 array
    :param M1: Final orientation, in same form, for T=1
    :param T: interpolation parameter
*/
    #local P10=MTxM(M1,M0);
    PrintMatrix("P10: ",P10)
    #local Trace=P10[0][0]+P10[1][1]+P10[2][2];
    PrintNumber("Trace: ",Trace)
    #local Theta1=acos((Trace-1)/2); // Note, theta matches units with return value of arccos, so radians
    PrintNumber("Theta1/deg: ",degrees(Theta1))
    #local A=<P10[2][1]-P10[1][2],
              P10[0][2]-P10[2][0],
              P10[1][0]-P10[0][1]>;
    #local A=A/vlength(A);
    PrintVector("A: ",A)
    #local Ax=array[3][3] {{   0,-A.z, A.y},
                           { A.z,   0,-A.x},
                           {-A.y, A.x,   0}};
    PrintMatrix("Ax: ",Ax)
    #local Ax2=MxM(Ax,Ax);
    PrintMatrix("Ax2: ",Ax2)
    #local Theta=T*Theta1;
    PrintNumber("Theta/deg: ",degrees(Theta))
    #local Pt0=Identity(3);
    PrintMatrix("Pt0: ",Pt0)
    Mfma(Pt0,Ax ,  sin(Theta))
    Mfma(Pt0,Ax2,1-cos(Theta))
    #local Mt=MxMT(M0,Pt0)
    Mt
#end



