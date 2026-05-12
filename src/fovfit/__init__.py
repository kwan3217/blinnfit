"""
Describe purpose of this script here

Created: 1/31/25
"""
import argparse

from matplotlib import pyplot as plt


def main():
    from fovfit.gui import FovFit
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Do field-of-view fits for the given case")

    # Add arguments
    parser.add_argument("case_name", default="SuperTrajectoryB",type=str, help="The name of the case")
    parser.add_argument("-f", "--frame", type=int, default=1600, help="Initial frame number")
    parser.add_argument("-b", "--seqid", type=str, default='VOBEST', help="Spice body id")
    parser.add_argument("-i", "--imgid", type=int, default=0, help="Spice body id")
    parser.add_argument("-N", "--narrow_angle", action='store_true', help="Use the Narrow Angle Camera field of view")

    # Parse arguments
    args = parser.parse_args()

    # Use the arguments
    case_name = args.case_name
    start_frame = args.frame

    # Example of how to use these in your script
    print(f"Running simulation for case: {case_name}")
    print(f"Starting from frame: {start_frame}")
    #CameraMount(casename="VoyagerUranusHD")
    FovFit(casename=case_name,initframe=start_frame,seqid=args.seqid,imgid=args.imgid)
    plt.show()


if __name__ == "__main__":
    main()


