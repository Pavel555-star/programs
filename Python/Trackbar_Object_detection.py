import cv2
import numpy as np

def nothing(x):
    pass
    
cv2.namedWindow('Mask')
cv2.createTrackbar('LH', 'Mask', 0, 255, nothing)
cv2.createTrackbar('LS', 'Mask', 0, 255, nothing)
cv2.createTrackbar('LV', 'Mask', 0, 255, nothing)
cv2.createTrackbar('HH', 'Mask', 255, 255, nothing)
cv2.createTrackbar('HS', 'Mask', 255, 255, nothing)
cv2.createTrackbar('HV', 'Mask', 255, 255, nothing)


while True:
    img = cv2.imread('input.png', -1);
    
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)

    l_h = cv2.getTrackbarPos('LH', 'Mask')
    l_s = cv2.getTrackbarPos('LS', 'Mask')
    l_v = cv2.getTrackbarPos('LV', 'Mask')
    h_h = cv2.getTrackbarPos('HH', 'Mask')
    h_s = cv2.getTrackbarPos('HS', 'Mask')
    h_v = cv2.getTrackbarPos('HV', 'Mask')
    
    l_b = np.array([l_h, l_s, l_v])
    h_b = np.array([h_h, h_s, h_v])
    
    mask = cv2.inRange(hsv, l_b, h_b)
    
    res = cv2.bitwise_and(img, img, mask = mask)
    
    cv2.imshow('image', img)
    cv2.imshow('mask', mask)
    cv2.imshow('res', res)
    
    k = cv2.waitKey(1) & 0xFF
    if k == 27:
        break

cv2.imwrite('output.png', res)
cv2.destroyAllWindows()

