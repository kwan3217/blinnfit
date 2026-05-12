"""

This is a clean sheet derivation of slerp with minimal help from Grok. I thought of almost all these ideas myself and used Grok as a duck to bounce my ideas off of.

$$
\def\q#1{{\mathbf{#1}}}
\def\M#1{{[\mathbf{#1}]}}
\def\MM#1#2{{[\mathbf{#1}{#2}]}}
\def\E{\operatorname{E}}
\def\cov{\operatorname{cov}}
\def\T{^\mathsf{T}}
\def\C{^*}
$$

# Problem Statement
A camera is at the origin of a universe frame $U$. At interpolation parameter (to first order we can think of this as scaled time) $t=0$ the camera is oriented in frame $C_0$, by which I mean the camera principal vectors like boresight, detector right, detector down are parallel to the unit-length basis vectors of $C_0$. At $t=1$, the camera is oriented in frame $C_1$.

Both of these frames can be quantified by the rotation matrix required to transform a vector in frame $C_n$ to $U$. This is specified to be a pure rotation in SO(3), since we imagine swinging a physical camera from one pointing to another. This matrix $\MM{M}{_{U,c_n}}$ follows our usual convention of subscripting the matrix such that its *to* frame is first, and its *from* frame is second. This allows us to write:

$$\vec{v}_U=\MM{M}{_{U,c_n}}\vec{v}_{c_n}$$

We note that the $U$ subscripts and $c_n$ subscripts are each adjacent. This notation will help us out shortly.

# Pan Matrix
There is *some* pan matrix $\M{P}$ which relates $\MM{M}{_{U,c_0}}$ and $\MM{M}{_{U,c_1}}$. We will construct an equation including this matrix and assign its subscripts by adjacency, and then interpret its physical meaning. We have for any arbitrary vector:

$$\begin{eqnarray}
\vec{v}_U&=&\MM{M}{_{U,c_0}}\vec{v}_{c_0} \\
 &=&\MM{M}{_{U,c_1}}\MM{P}{_{yx}}\vec{v}_{c_0}
 \end{eqnarray}$$

Here we are constrained by the subscripts to choose to use $\MM{M}{_{U,c_1}}$ *on the left* because there is no other way to get a $U$ to match the $U$ subscript on the vector on the left of the equals. So now we match up adjacent subscripts on $\M{P}$:

$$\begin{eqnarray}
c_1&=&y \\
c_0&=&x \\
\MM{P}{_{yx}}&=&\MM{P}{_{c_1,c_0}}\\
\vec{v}_U&=&\MM{M}{_{U,c_1}}\MM{P}{_{c_1,c_0}}\vec{v}_{c_0}
\end{eqnarray}$$

Since this holds for all vectors $\vec{v}$, we identify $\MM{M}{_{U,c_1}}\MM{P}{_{c_1,c_0}}$ with $\MM{M}{_{U,c_0}}$ and from there we can solve for $\M{P}$:

$$\begin{eqnarray}
\MM{M}{_{U,c_1}}\MM{P}{_{c_1,c_0}}&=&\MM{M}{_{U,c_0}} \\
\MM{M}{_{U,c_1}}^{-1}\MM{M}{_{U,c_1}}\MM{P}{_{c_1,c_0}}&=&\MM{M}{_{U,c_1}}^{-1}\MM{M}{_{U,c_0}} \\
\MM{P}{_{c_1,c_0}}&=&\MM{M}{_{U,c_1}}^{-1}\MM{M}{_{U,c_0}} \\
 &=&\MM{M}{_{U,c_1}}\T\MM{M}{_{U,c_0}} & \mbox{ since in SO(3), }\M{M}^{-1}=\M{M}\T\\
 &=&\MM{M}{_{c_1,U}}\MM{M}{_{U,c_0}} & \mbox{ since in SO(3), }\MM{M}{_{ba}}\T=\MM{M}{_{ab}}\\
\end{eqnarray}$$

This last line makes perfect sense. A transform from $C_0$ to any arbitrary frame, then from that arbitrary frame to $C_1$ is the same as a transform straight from $C_0$ to $C_1$. The pan matrix by itself therefore acts thusly:

$$\vec{v}_{c_1}=\MM{P}{_{c_1,c_0}}\vec{v}_{c_0}$$

Now like all SO(3) matrices, $\M{P}$ can be represented as an axis and angle. Vectors on the axis of any rotation are not changed by that rotation. That's what axis means, and coincidentally that's also what eigenvector means. Any vector parallel or antiparallel to the axis is an eigenvector, and the eigenvalue is 1. Similarly, the angle of rotation is directly related to the trace (sum of diagonal elements) of $\M{P}$.

We *could* put down the derivation of these two facts, or we can trust the results printed by Grok and originally found by Euler and Rodrigues:

$$\begin{eqnarray}
\operatorname{trace}(\MM{P}{_{c_1,c_0}})&=&1+2\cos(\theta)\\
\vec{a}&=&\operatorname{vnormalize}\left(\begin{bmatrix}
P_{21}-P_{12} \\
P_{02}-P_{20} \\
P_{10}-P_{01}
\end{bmatrix}\right)
\end{eqnarray}$$

We have these corner cases, which won't be a practical issue in the camera slewing problem:

* If $\theta=0$, then the vector inside the normalize function will be zero length, and therefore have no well-defined direction to normalize. In this case, *any* vector is an axis and it doesn't matter since the angle is 0.
* If $\theta=\pi=180^\circ$, then apparently the antisymmetric method is unstable and we need another method.

As I stated, Euler and Rodrigues figure this out centuries ago, and we have [this formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Matrix_notation) with Rodrigues' name on it:

$$\M{P}(\vec{a},\theta)=\M{1}+\sin\theta\MM{a}{_\times}+(1-\cos\theta)\MM{a}{_\times}^2$$

where $\MM{a}{_\times}$ is the skew-symmetric matrix:

$$\MM{a}{_\times}=\begin{bmatrix}
0 & -a_z & a_y\\
a_z & 0 & -a_x \\
-a_y & a_x & 0
\end{bmatrix}$$

and $\MM{a}{_\times}^2=\MM{a}{_\times}\MM{a}{_\times}$.

Note that since the matrices depend on $\vec{a}$ and not $\theta$, they don't need to be recomputed for a different angle.

WARNING: The following is an independent idea that happens to NOT WORK. Do not use it.

>Another way to do this that might make sense if there are a lot of angles is to compute the matrices when rotated $0^\circ$ and $90^\circ$ and get the weighted sum where the weight is the cosine and sine of theta:
>
>$$\begin{eqnarray}
\M{P}(\theta=0)&=&\M{1}\\
\M{P}(\theta=90^\circ)&=&\M{K}\\
 &=&\M{1}+\sin\theta\MM{a}{_\times}+(1-\cos\theta)\MM{a}{_\times}^2\\
 &=&\M{1}+1\MM{a}{_\times}+(1-0)\MM{a}{_\times}^2\\
\M{P}(\theta)&=&\cos\theta\M{1}+\sin\theta\M{K}
\end{eqnarray}$$
>
>In any case, since the parameter $t$ varies from 0 to 1, the linear interpolation of $\theta$ is super-easy -- it's just $\theta=t\theta_1$.

Either way, once we have $\M{P}$, we can calculate the intermediate orientation of the body $\MM{M}{_U,c_t}$:
$$\begin{eqnarray}
\MM{M}{_{U,c_t}}\MM{P}{_{c_t,c_0}}&=&\MM{M}{_{U,c_0}} \\
\MM{M}{_{U,c_t}}\MM{P}{_{c_t,c_0}}\MM{P}{_{c_t,c_0}}\T&=&\MM{M}{_{U,c_0}}\MM{P}{_{c_t,c_0}}\T \\
\MM{M}{_{U,c_t}}&=&\MM{M}{_{U,c_0}}\MM{P}{(\theta)_{c_t,c_0}}\T \\
 &=&\MM{M}{_{U,c_0}}\MM{P}{(\theta)_{c_0,c_t}} \\
\end{eqnarray}$$
This gives a transformation matrix to the universe frame from the interpolated camera frame, which is isomorphic to and can be used as the camera orientation.

# Algorithm
All of this leads to the following algorithm:
* Input:
   * $\MM{M}{_{U,c_0}}$ - Camera orientation at $t=0$ in the form of a matrix that transforms vectors in camera space to vectors in universe space
   * $\MM{M}{_{U,c_1}}$ - Camera orientation at $t=1$
   * $t$ - interploation parameter. Normally covers $0\le t \le 1$ but it is well-behaved and works as expected for any $-\infty \lt t \lt \infty$.
* Output:
   * $\MM{M}{_{U,c_t}}$ - Camera orientation at $t$
* Algorithm:
   * Pre-calculate the following, valid for given matrices $\M{M}$ and *any* $t$:$$\begin{eqnarray}
\MM{P}{_{c_1,c_0}}&=&\MM{M}{_{U,c_1}}\T\MM{M}{_{U,c_0}}\\
\theta_1&=&\operatorname{acos}\left(\frac{\operatorname{trace}(\MM{P}{_{c_1,c_0}})-1}{2}\right) \\
\vec{a}&=&\operatorname{vnormalize}\left(\begin{bmatrix}
P_{21}-P_{12} \\
P_{02}-P_{20} \\
P_{10}-P_{01}
\end{bmatrix}\right) \\
\MM{a}{_\times}&=&\begin{bmatrix}
0 & -a_z & a_y\\
a_z & 0 & -a_x \\
-a_y & a_x & 0
\end{bmatrix}\\
\MM{a}{_\times}^2&=&\MM{a}{_\times}\MM{a}{_\times}
\end{eqnarray}$$
   * For given $t$, calculate the following:$$\begin{eqnarray}
\theta&=&t\theta_1 \\
\MM{P}{_{c_t,c_0}}&=&\M{1}+\sin\theta\MM{a}{_\times}+(1-\cos\theta)\MM{a}{_\times}^2 \\
\MM{M}{_{U,c_t}}&=&\MM{M}{_{U,c_0}}\MM{P}{_{c_t,c_0}}\T \\
\end{eqnarray}$$

"""

import numpy as np
from typing import Callable, Optional
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def slerp_precalc(M0,M1,verbose:bool=False):
    P10=M1.T@M0
    if verbose:
        print(f"{P10=}")
    trace=P10[0,0]+P10[1,1]+P10[2,2]
    if verbose:
        print(f"{trace=}")
    theta1=np.arccos((trace-1)/2) # Note, theta matches units with return value of arccos, so radians
    if verbose:
        print(f"{theta1=} rad")
    a=np.array([[P10[2,1]-P10[1,2]],
                [P10[0,2]-P10[2,0]],
                [P10[1,0]-P10[0,1]]])
    a/=np.linalg.norm(a)
    if verbose:
        print(f"{a=}")
    ax=np.array([[      0,-a[2,0], a[1,0]],
                 [ a[2,0],      0,-a[0,0]],
                 [-a[1,0], a[0,0],      0]])
    if verbose:
        print(f"{ax=}")
    ax2=ax@ax
    if verbose:
        print(f"{ax2=}")
    return ax,ax2,theta1


def slerp_core(*,M0,ax,ax2,theta1,t,verbose:bool=False):
    theta=t*theta1
    if verbose:
        print(theta.shape)
    Pt0=np.eye(3)
    if verbose:
        print(Pt0.shape)
    Pt0=Pt0+np.sin(theta)*ax
    if verbose:
        print(Pt0.shape)
    Pt0=Pt0+(1-np.cos(theta))*ax2
    if verbose:
        print(Pt0.shape)
    # We are doing this instead of Pt0.T because of the difference with bundles:
    # np.zeros((10,20,30,40)).T reverses all dimensions so shape (10,20,30,40) turns into
    # (40,30,20,10) which is wrong for bundles of matrices.
    # np.matrix_transpose(np.zeros(10,20,30,40)) on the other hand gives
    # shape (10,20,40,30), just reversing the last two, following the de facto
    # convention for matrix bundles.
    Pt0T=np.matrix_transpose(Pt0)
    Mt=M0@Pt0T
    if verbose:
        print(Mt.shape)
    return Mt


def slerp_naive(M0:np.ndarray,M1:np.ndarray,t:float|np.ndarray,verbose:bool=False)->np.ndarray:
    """
    Spherical linear interpolation of SO(3) rotations. Most naive translation of the given
    algorithm, with no caching, currying, etc.

    :param M0: Initial orientation, in the form of a matrix transforming camera to universe at t=0
    :param M1: Final orientation, in same form, for t=1
    :param t: interpolation parameter. Can be an array but must broadcast against a 3x3 array, so (...,1,1)
    """
    ax,ax2,theta1=slerp_precalc(M0,M1,verbose=verbose)
    return slerp_core(M0=M0,ax=ax,ax2=ax2,theta1=theta1,t=t,verbose=verbose)


MInterp=Callable[[float|np.ndarray,bool],np.ndarray]


def slerp(M0:np.ndarray,M1:np.ndarray,verbose:bool=False)->MInterp:
    """
    Currying function, used when we have one set of matrices but many interpolation points. Use it like this:
    s=slerp(M0,M1)
    Mta=s(ta)
    Mtb=s(tb)
    Mtc=slerp(M0c,M1c)(tc) # if you only want one

    :param M0: Initial orientation, in the form of a matrix transforming camera to universe at t=0
    :param M1: Final orientation, in same form, for t=1
    :return: function which takes the interpolation parameter and returns the interpolated orientation
    """
    ax,ax2,theta1=slerp_precalc(M0,M1,verbose=verbose)
    def inner(t:float|np.ndarray,verbose:bool=False)->np.ndarray:
        """
        :param t: interpolation parameter. To be multi-dimensional I think it has to be shape (...,1,1)
        """
        return slerp_core(M0=M0,ax=ax,ax2=ax2,theta1=theta1,t=t,verbose=verbose)
    return inner


def plot_slerp(s:MInterp, title:str='Slerp Visualization: Basis Vector Interpolation'):

    # Set up 3D plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot unit sphere (wireframe)
    u = np.linspace(0, 2 * np.pi, 20)
    v = np.linspace(0, np.pi, 20)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_wireframe(x, y, z, color='gray', alpha=0.2)

    # Basis vectors for plotting
    basis = np.eye(3)  # Standard basis vectors [e1, e2, e3]
    colors = ['r', 'g', 'b']  # Colors for x, y, z axes
    labels = ['x-axis', 'y-axis', 'z-axis']

    # Plot basis vectors for M0 (t=0)
    M = s(0)
    for i in range(3):
        v = M @ basis[:, i]  # Rotate basis vector
        ax.quiver(0, 0, 0, v[0], v[1], v[2], color=colors[i], linewidth=2, label=f'M0 {labels[i]} (t=0)')

    # Plot basis vectors for M1 (t=1)
    M = s(1)
    for i in range(3):
        v = M @ basis[:, i]
        ax.quiver(0, 0, 0, v[0], v[1], v[2], color=colors[i], linewidth=2, linestyle='--', label=f'M1 {labels[i]} (t=1)')

    # Plot interpolation paths for each basis vector
    t_values = np.linspace(0, 1, 50)
    for i in range(3):
        path = []
        for t in t_values:
            M = s(t)
            v = M @ basis[:, i]
            path.append(v)
        path = np.array(path)
        ax.plot(path[:, 0], path[:, 1], path[:, 2], color=colors[i], linestyle=':', label=f'{labels[i]} path')

    # Set plot properties
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    ax.set_zlim([-1.5, 1.5])
    ax.legend()
    ax.set_title(title)

    # Ensure equal aspect ratio
    ax.set_box_aspect([1, 1, 1])

    plt.show()


# Function to generate random SO(3) matrix
def random_so3_matrix() -> np.ndarray:
    """
    Generate a random 3x3 rotation matrix in SO(3) using uniform random quaternions.
    :return: A 3x3 rotation matrix (orthogonal, det=1).
    """
    theta1, theta2 = np.random.uniform(0, 2 * np.pi, size=2)
    z = np.random.uniform(0, 1)
    s1, s2 = np.sqrt(1 - z), np.sqrt(z)
    w = s1 * np.cos(theta1)
    x = s1 * np.sin(theta1)
    y = s2 * np.cos(theta2)
    z = s2 * np.sin(theta2)
    q = np.array([w, x, y, z])
    q /= np.linalg.norm(q)
    w, x, y, z = q
    R = np.array([
        [1 - 2 * (y ** 2 + z ** 2), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x ** 2 + z ** 2), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x ** 2 + y ** 2)]
    ])
    return R


def generate_ref_inc(s:MInterp,
                     t_values:Optional[np.ndarray]=None,
                    filename:str="slerp_ref.inc",
                    POVVarName:str="TestPoints",
                    verbose:bool=False):
    if t_values is None:
        t_values= np.linspace(0, 1, 51)
    output=f"#declare {POVVarName}=array[{len(t_values)}][3][3] {{ //first index is test point, second is row, third is column\n"
    Ms = s(t_values.reshape(-1,1,1))
    for i,t in enumerate(t_values):
        M=Ms[i,:,:]
        output+=f"/*t={t:7.4f}*/ {{{{{M[0][0]:10.6f},{M[0][1]:10.6f},{M[0][2]:10.6f}}},\n"
        output+=f"               {{{M[1][0]:10.6f},{M[1][1]:10.6f},{M[1][2]:10.6f}}},\n"
        output+=f"               {{{M[2][0]:10.6f},{M[2][1]:10.6f},{M[2][2]:10.6f}}}}},\n"
    output+="}\n"
    if verbose:
        print(output)
    with open(filename,"wt") as ouf:
        print(output,file=ouf)


