import cv2

img = cv2.imread('input.png')
imgGrey = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
_, thresh = cv2.threshold(imgGrey, 240, 255, cv2.THRESH_BINARY)
contours, _ = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

f = open("output.csv", "w")
i = 0

for contour in contours:
    i = i + 1
    approx = cv2.approxPolyDP(contour, 0.01* cv2.arcLength(contour, True), True)
    cv2.drawContours(img, [approx], 0, (0, 0 ,0), 5)
    x = approx.ravel()[0]
    y = approx.ravel()[1]
    if len(approx) == 3:
        cv2.putText(img, 'Triangle ' + str(i), (x, y), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 
        (0, 0, 0))
        f.write( 'Triangle ' + str(i) + ' [x, y]\n')
        
    elif len(approx) == 4:
        (x, y, w, h) = cv2.boundingRect(approx)
        aspectRatio = float(w)/h
        if aspectRatio >= 0.95 and aspectRatio <= 1.05:
            cv2.putText(img, 'Square ' + str(i), (x, y), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 
            (0, 0, 0))
            f.write( 'Square ' + str(i) + ' [x; y]\n')
        else:
            cv2.putText(img, 'Rectangle ' + str(i), (x, y), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 
            (0, 0, 0))
            f.write( 'Rectangle ' + str(i) + ' [x; y]\n')
            
    elif len(approx) == 5:
        cv2.putText(img, 'Pentagon ' + str(i), (x, y), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 
        (0, 0, 0))
        f.write('Pentagon ' + str(i) + ' [x; y]\n')
        
    elif len(approx) == 6:
        cv2.putText(img, 'Hexagon '  + str(i), (x, y), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 
        (0, 0, 0))
        f.write('Hexagon '  + str(i) + ' [x; y]\n')
        
    elif len(approx) == 7:
        cv2.putText(img, 'Heptagon ' + str(i), (x, y), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 
        (0, 0 ,0))
        f.write('Heptagon ' + str(i) + ' [x; y]\n')
        
    elif len(approx) == 8:
        cv2.putText(img, 'Octagon ' + str(i), (x, y), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 
        (0, 0, 0))
        f.write('Octagon ' + str(i) + ' [x; y]\n')
        
    else:
        (x, y, w, h) = cv2.boundingRect(approx)
        aspectRatio = float(w)/h
        if aspectRatio >= 0.95 and aspectRatio <= 1.05:
            cv2.putText(img, 'Circle ' + str(i), (x, y), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 
            (0, 0, 0))
            f.write('Circle [x; y]' + str(i) + '\n')
        else:
            cv2.putText(img, 'Elipse ' + str(i), (x, y), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 
            (0, 0, 0))
            f.write('Elipse [x; y]' + str(i) + '\n')
            
    r = contour.ravel()
    for position in range(len(r)):
        if position % 2 == 0:
            f.write(str(r[position]) + '; '  + str(r[position + 1]) + '\n')
            
cv2.imwrite('output.png', img)
f.close()
cv2.imshow('shapes', img)
cv2.waitKey(0)
cv2.destroyAllWindows()
