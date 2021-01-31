import cv2
import sys

infile = "input.png"
outfile = "output.png"
n = len(sys.argv)
if  n >= 2:
    infile = sys.argv[1]
if  n >= 3:
    outfile = sys.argv[2]

img = cv2.imread(infile) # Přete obrázek input.png
canny = cv2.Canny(img, 60, 100)

canny = cv2.bitwise_not(canny)
cv2.imwrite(outfile, canny) # Uloží výstup do obrázku output.png

