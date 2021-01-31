import cv2
import sys

infile = "input.png"
outfile = "output.png"
n = len(sys.argv)
if  n >= 2:
    infile = sys.argv[1]
if  n >= 3:
    outfile = sys.argv[2]
    
img = cv2.imread("input.png")
img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY) # Přete obrázek input.png

threshold = cv2.adaptiveThreshold(img, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 15, 1)
threshold = cv2.bitwise_not(threshold)
cv2.imwrite("output.png", threshold) # Uloží výstup do obrázku output.png

