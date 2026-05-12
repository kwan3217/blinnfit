// Matrix math stuff

#macro MxM(A,B)
  /* Matrix Multiplication
  :param A: left matrix, mxn
  :param B: right matrix, nxp
  :return: Matrix product A@B
  */
  #local M=dimension_size(A,1);
  #local N=dimension_size(A,2);
  #local P=dimension_size(B,2);
  #local C=array[M][P];
  #local I=0;
  #while(I<M)
    #local K=0;
    #while(K<P)
      #local J=0;
      #local C[I][K]=0;
      #while(J<N)
        #local C[I][K]=C[I][K]+A[I][J]*B[J][K];
        #local J=J+1;
      #end
      #local K=K+1;
    #end
    #local I=I+1;
  #end
  C
#end

#macro MTxM(A,B)
  /* Matrix transpose Multiplication
  :param A: left matrix, mxn
  :param B: right matrix, nxp
  :return: Matrix product A.T@B
  */
  #local M=dimension_size(A,2);
  #local N=dimension_size(A,1);
  #local P=dimension_size(B,2);
  #local C=array[M][P];
  #local I=0;
  #while(I<M)
    #local K=0;
    #while(K<P)
      #local J=0;
      #local C[I][K]=0;
      #while(J<N)
        #local C[I][K]=C[I][K]+A[J][I]*B[J][K]; // Note that order of A subscripts are swapped, to give A.T
        #local J=J+1;
      #end
      #local K=K+1;
    #end
    #local I=I+1;
  #end
  C
#end

#macro MxMT(A,B)
  /* Matrix transpose Multiplication
  :param A: left matrix, mxn
  :param B: right matrix, nxp
  :return: Matrix product A@B.T
  */
  #local M=dimension_size(A,1);
  #local N=dimension_size(A,2);
  #local P=dimension_size(B,1);
  #local C=array[M][P];
  #local I=0;
  #while(I<M)
    #local K=0;
    #while(K<P)
      #local J=0;
      #local C[I][K]=0;
      #while(J<N)
        #local C[I][K]=C[I][K]+A[I][J]*B[K][J]; // Note that order of B subscripts are swapped, to give B.T
        #local J=J+1;
      #end
      #local K=K+1;
    #end
    #local I=I+1;
  #end
  C
#end


#macro MxV(A,V)
  /* Matrix Vector Multiplication M@V (where V is treated as a 3-element column vector)
  :param A: Matrix, must be 3x3
  :param V: vector, only x, y, and z are used
  */
  #local C=<A[0][0]*V.x+A[0][1]*V.y+A[0][2]*V.z,
            A[1][0]*V.x+A[1][1]*V.y+A[1][2]*V.z,
            A[2][0]*V.x+A[2][1]*V.y+A[2][2]*V.z>;
  C
#end

#macro MTxV(A,V)
  /* Matrix Vector Multiplication M.T@V
  :param A: Matrix, must be 3x3
  :param V: vector, only x, y, and z are used
  */
  #local C=<A[0][0]*V.x+A[1][0]*V.y+A[2][0]*V.z,
            A[0][1]*V.x+A[1][1]*V.y+A[2][1]*V.z,
            A[0][2]*V.x+A[1][2]*V.y+A[2][2]*V.z>;
  C
#end

#macro strmul(S,N)
  #local I=0;
  #local Result="";
  #while(I<N)
    #local Result=concat(Result,S)
    #local I=I+1;
  #end
  Result
#end

#macro PrintMatrix(Label,A)
  #local M=dimension_size(A,1);
  #local N=dimension_size(A,2);
  #local I=0;
  #while(I<M)
    #local J=0;
    #if(I=0)
      #debug Label
    #else
      #debug strmul(" ",strlen(Label))
    #end
    #debug "["
    #while(J<N)
      #debug concat(str(A[I][J],10,6),",")
      #local J=J+1;
    #end
    #debug "]\n"
    #local I=I+1;
  #end
#end

#macro Identity(N)
  #local Result=array[N][N];
  #local I=0;
  #while(I<N)
    #local J=0;
    #while(J<N)
      #local Result[I][J]=0;
      #local J=J+1;
    #end
    #local Result[I][I]=1;
    #local I=I+1;
  #end
  Result
#end

#macro MpM(A,B)
  /* Matrix addition
  :param A: left matrix, mxn
  :param B: right matrix, nxp
  :return: Matrix sum A+B
  */
  #local M=dimension_size(A,1);
  #local N=dimension_size(A,2);
  #local C=array[M][N];
  #local I=0;
  #while(I<M)
    #local J=0;
    #while(J<N)
      #local C[I][J]=A[I][J]+B[I][J];
      #local J=J+1;
    #end
    #local I=I+1;
  #end
  C
#end

#macro MxS(A,S)
  /* Matrix Scaling
  :param A: left matrix, mxn
  :param S: scalar
  :return: Scaled Matrix S*A
  */
  #local M=dimension_size(A,1);
  #local N=dimension_size(A,2);
  #local C=array[M][N];
  #local I=0;
  #while(I<M)
    #local J=0;
    #while(J<N)
      #local C[I][J]=A[I][J]*S;
      #local J=J+1;
    #end
    #local I=I+1;
  #end
  C
#end

#macro Mpp(A,B)
  /* Matrix addition accumulation
  :param A: left matrix, mxn
  :param B: right matrix, nxp
  :Side effect: A<-A+B
  */
  #local M=dimension_size(A,1);
  #local N=dimension_size(A,2);
  #local I=0;
  #while(I<M)
    #local J=0;
    #while(J<N)
      #declare A[I][J]=A[I][J]+B[I][J];
      #local J=J+1;
    #end
    #local I=I+1;
  #end
#end

#macro Mfma(A,B,S)
  /* Matrix addition accumulation
  :param A: left matrix, mxn
  :param B: right matrix, nxp
  :param S: scalar for B
  :Side effect: A<-A+B*s
  */
  #local M=dimension_size(A,1);
  #local N=dimension_size(A,2);
  #local I=0;
  #while(I<M)
    #local J=0;
    #while(J<N)
      #declare A[I][J]=A[I][J]+B[I][J]*S;
      #local J=J+1;
    #end
    #local I=I+1;
  #end
#end

#macro M3b1b(e0,e1,e2)
  /* Create a transformation matrix using the 3blue1brown method:
     For each column of the matrix, put down where the corresponding
     basis vector of the *from* space lands in *to* space. This will
     transform a vector in *from* space to *to* space.
     :param e0: location where basis vector 0 (x) lands in *to* space
     :param e1:                             1 (y)
     :param e2:                             2 (z)
     */
  array[3][3] {{e0.x,e1.x,e2.x},
               {e0.y,e1.y,e2.y},
               {e0.z,e1.z,e2.z}}
#end

#macro MatrixVerbose(Label,
                     Ix,Iy,Iz,
                     Jx,Jy,Jz,
                     Kx,Ky,Kz,
                     Tx,Ty,Tz)
  // Print a matrix transformation as we make it,
  // because its elements are not available
  // once they are in the transformation
  // :param ix,iy,iz: components of i-hat of *body* (from) frame in *reference* (to) frame
  // :param jx,jy,jz: components of j-hat of *body* (from) frame in *reference* (to) frame
  // :param kx,ky,kz: components of k-hat of *body* (from) frame in *reference* (to) frame
  // :param tx,ty,tz: location of *body* (from) origin in *reference* (to) frame
  #debug concat(Label                    ,"<",str(Ix,9,6),",",str(Iy,9,6),",",str(Iz,9,6),",\n")
  #debug concat(strmul(" ",strlen(Label))," ",str(Jx,9,6),",",str(Jy,9,6),",",str(Jz,9,6),",\n")
  #debug concat(strmul(" ",strlen(Label))," ",str(Kx,9,6),",",str(Ky,9,6),",",str(Kz,9,6),",\n")
  #debug concat(strmul(" ",strlen(Label))," ",str(Tx,9,6),",",str(Ty,9,6),",",str(Tz,9,6),">\n")
  matrix <Ix,Iy,Iz,
          Jx,Jy,Jz,
          Kx,Ky,Kz,
          Tx,Ty,Tz>
#end