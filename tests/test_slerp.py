"""
Describe purpose of this script here

Created: 8/6/25
"""
import numpy as np

from slerp import slerp_naive, slerp, plot_slerp, random_so3_matrix, generate_ref_inc


def test_slerp_broadcast():
    # Test broadcasting - M0 is aligned to U, M1 is a 90deg right-handed rotation around the z axis.
    # This points C1 +x at U +y, C1 +y at U -x, and C1 +z at U +z. Using the 3b1b matrix building,
    # we take each basis vector of C1, record where it lands in U, as the column of the matrix converting C1 to U:
    M0 = np.eye(3)
    M1 = np.array([[0.0, -1.0, 0.0],
                   [1.0, 0.0, 0.0],
                   [0.0, 0.0, 1.0]])
    # Trivial cases: for t=0, we should get M0, and for t=1, we should get M1
    print(slerp_naive(M0, M1, 0.0))
    assert np.allclose(slerp_naive(M0, M1, 0.0), M0)
    print(slerp_naive(M0, M1, 1.0))
    assert np.allclose(slerp_naive(M0, M1, 1.0), M1)
    # Interesting case: for t=0.5, should have 0.707 in the upper left 2x2. Specifically something like:
    # x_ct is between +x_u and +y_u, so [.707,.707,0]^T
    # y_ct is between +y_u and -x_u, so [-.707,.707,0]^T
    # z_ct is still z_u
    print(slerp_naive(M0, M1, 0.5))
    s = slerp(M0, M1)
    print(s(0.5))

    # Test broadcasting
    t = np.arange(0, 1, 0.01).reshape(-1, 1, 1)  # Shape (100, 1, 1)
    Mt = slerp_naive(M0, M1, t)
    Mt = s(t)  # Shape (100, 3, 3)
    print("Testing slerp broadcasting:")
    print("Shape of Mt:", Mt.shape)
    print("s(0) = M0:\n", Mt[0, :, :])
    print("s(0.5):\n", Mt[50, :, :])
    print("s(1) = M1:\n", Mt[-1, :, :])

    plot_slerp(s)


def test_identity_to_random():
    # Set random seed for reproducibility
    np.random.seed(3217)

    M0 = np.identity(3)
    M1 = random_so3_matrix()
    s = slerp(M0, M1)

    plot_slerp(s, title="Identity to random SO(3)")


def test_generate_ref_inc():
    """
    Despite being called a test, this *generates* slerp_ref.inc
    :return:
    """
    # Set random seed for reproducibility
    np.random.seed(3217)

    M0 = np.identity(3)
    M1 = random_so3_matrix()
    s = slerp(M0, M1)

    generate_ref_inc(s)


def test_frame2186():
    # Test case Blinn - M0 and M1 are from the Jupiter Voyager 1 flyby
    M0 = np.array([[-0.607490, -0.381916, 0.696488],
                   [-0.585954, -0.376535, -0.717551],
                   [0.536296, -0.844015, 0.004957]])
    M1 = np.array([[-0.456767, -0.392419, -0.798355],
                   [0.611991, 0.512712, -0.602157],
                   [0.645624, -0.763632, 0.005967]])

    s = slerp(M0, M1)
    plot_slerp(s, title="Blinn matrix interpolations")
    generate_ref_inc(s, filename="jupiter_io_pan.inc")


def test_frame2186_t0():
    # Test case Blinn - M0 and M1 are from the Jupiter Voyager 1 flyby
    M0 = np.array([[-0.607490, -0.381916, 0.696488],
                   [-0.585954, -0.376535, -0.717551],
                   [0.536296, -0.844015, 0.004957]])
    M1 = np.array([[-0.456767, -0.392419, -0.798355],
                   [0.611991, 0.512712, -0.602157],
                   [0.645624, -0.763632, 0.005967]])

    s = slerp(M0, M1, verbose=True)
    Mt=s(t=0.0,verbose=True)
    plot_slerp(s, title="Blinn matrix interpolations")
    generate_ref_inc(s, filename="jupiter_io_pan.inc")


