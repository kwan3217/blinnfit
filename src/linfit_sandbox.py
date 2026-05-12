"""
Describe purpose of this script here

Created: 8/6/25
"""

import numpy as np
from matplotlib import pyplot as plt


def simple_lin_regression(xs:np.ndarray,ys:np.ndarray)->tuple[float,float]:
    sy=np.sum(ys)
    sxx=np.sum(xs*xs)
    sx=np.sum(xs)
    sxy=np.sum(xs*ys)
    sx2=sx**2
    n=len(xs)
    #common denominator
    d=(n*sxx-sx2)
    # intercept coefficient alpha
    alpha=(sy*sxx-sx*sxy)/d
    # slope coefficient beta
    beta=(n*sxy-sx*sy)/d
    # change to familiar y=mx+b
    m=beta
    b=alpha
    return m,b




def main():
    rubber_clock=np.array([[1238,	-657299160.0     ],
                           [1346,	-657292812.046332],
                           [1387,	-657290822.175032],
                           [1388,	-657290343.397684],
                           [1389,	-657290284.620335],
                           [1390,	-657290225.842986],
                           [1391,	-657290167.065637],
                           [1701,	-657272606.672883],
                           [2701,	-657220199.519541],
                           [3756,	-657159019.772765],
                           [4515,	-657115048.349255],
                           [4797,	-657098669.14052 ]])
    rubber_clock_frame=rubber_clock[:,0]
    rubber_clock_et=rubber_clock[:,1]
    m,b=simple_lin_regression(xs=rubber_clock_frame,ys=rubber_clock_et)
    lin_et=m*rubber_clock_frame+b
    plt.figure("rubber clock")
    plt.subplot(2,1,1)
    plt.plot(rubber_clock_frame,(rubber_clock_et-rubber_clock_et[0])/3600,'*')
    plt.plot(rubber_clock_frame,(lin_et-rubber_clock_et[0])/3600,'-')
    plt.subplot(2,1,2)
    plt.plot(rubber_clock_frame,(rubber_clock_et-lin_et)/3600,'*')
    plt.plot(rubber_clock_frame,0*rubber_clock_frame,'-')

    JupiterPos=np.array([[ 330,	218.086650516856,	229.128544196091],
                         [ 430,	207.991758097495,	223.355790534115],
                         [ 530,	199.897673147538,	214.631523301627],
                         [ 630,	197.925496001986,	207.942872090036],
                         [ 670,	199.580115082222,	206.025989560958],
                         [ 700,	194.313484215603,	201.72075979344 ],
                         [ 717,	185.640994657917,	198.365917743208],
                         [ 800,	178.05830242027 ,	191.462338009167],
                         [ 894,	195.999814421613,	194.361074695563],
                         [ 900,	196.543810731952,	194.39902701721 ],
                         [ 979,	213.307043783177,	197.025811555916],
                         [1000,	217.421302205295,	198.885076811364],
                         [1100,	237.885670711657,	203.281969384618],
                         [1200,	258.499711844441,	206.641237640298],
                         [1238,	268.172951591134,	208.74939357576 ],
                         [1300,	279.278662184253,	209.835532426236],
                         [1346,	289.821192681776,	213.05068735074 ],
                         [1387,	297.306934827983,	214.165522399942],
                         [1400,	299.053536481189,	215.09862107381 ],
                         [1500,	321.712009715949,	220.01519533549 ],
                         [1600,	340.127919393118,	222.037790289735],
                         [1700,	353.636735817922,	227.268348575425],
                         [1701,	359.11554427203 ,	226.675356969349],
                         [1800,	371.40406732944 ,	234.33982446882 ],
                         [1900,	374.125967649316,	233.957801093284],
                         [2000,	338.56836422701 ,	252.934236465958],
                         [3756,	366.960729754548,	208.896668393555],
                         [4515,	356.986465711134,	213.815450059114],
                         [4797,	333.791771516755,	226.936599745126]
                        ])
    jupiter_pos_frames=JupiterPos[:,0]
    jupiter_pos_xs=JupiterPos[:,1]
    jupiter_pos_ys=JupiterPos[:,2]
    plt.figure("Jupiter pos")
    plt.subplot(2,2,1)
    plt.plot(jupiter_pos_xs,jupiter_pos_ys)
    for x,y,s in zip(jupiter_pos_xs,jupiter_pos_ys,[str(frame) for frame in jupiter_pos_frames]):
        plt.text(x,y,s)
    plt.subplot(2,2,2)
    plt.plot(jupiter_pos_frames,jupiter_pos_ys,'*')
    plt.subplot(2,2,3)
    plt.plot(jupiter_pos_xs,jupiter_pos_frames,'*')
    plt.show()


if __name__ == "__main__":
    main()
