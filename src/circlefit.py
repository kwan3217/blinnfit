import cv2
import numpy as np
import matplotlib.pyplot as plt

# Load the image
img = cv2.imread('data/frames/SuperTrajectoryB/frame0600.png')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
img2=gray*0
img2[305:340,595:630]=gray[305:340,595:630]
# Apply Hough Transform
circles = cv2.HoughCircles(img2, cv2.HOUGH_GRADIENT, 1, 20, param1=50, param2=30, minRadius=1, maxRadius=50)

plt.imshow(img2)
# Ensure circles were found
if circles is not None:
    circles = np.uint16(np.around(circles))
    q = np.arange(0, np.pi * 2, 0.01)
    c = np.cos(q)
    s = np.sin(q)
    for i in circles[0,:]:
        plt.plot(c * i[2] + i[0], s * i[2] + i[1])
        plt.plot(i[0],i[1],'+')
else:
    print("No circles found")
plt.show()