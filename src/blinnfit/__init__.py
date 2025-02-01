"""
Describe purpose of this script here

Created: 1/31/25
"""
from matplotlib import pyplot as plt


def main():
    from blinnfit.gui import CameraMount
    CameraMount(casename="VoyagerUranusHD")
    plt.show()


if __name__ == "__main__":
    main()


