"""
Experiment with edge detection

Created: 2/5/25
"""
import cv2
import numpy as np
from matplotlib import pyplot as plt, image as mpimg


def main():
    plt.figure()
    plt.subplot(2,1,1)
    framenum=2000
    infn = f"data/frames/VoyagerUranusHD/frame{framenum:04d}.png"
    img = mpimg.imread(infn)
    img=(img[:,:,2]*255).astype(np.uint8) # Take just the blue channel
    plt.imshow(img)
    # Apply Gaussian blur to reduce noise
    #img = cv2.GaussianBlur(img, (5, 5), 0)

    # Canny edge detection
    edges = cv2.Canny(img, 50, 150)  # Adjust thresholds
    plt.subplot(2,1,2)
    plt.imshow(edges)
    plt.show()


if __name__ == "__main__":
    main()
