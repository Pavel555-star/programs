from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import csv
import sys

fig = plt.figure()
ax = fig.gca(projection='3d')

for i in range (1, len(sys.argv)):
    with open(sys.argv[i], 'r') as csvfile: # Načte a zpracuje vstupní soubory uvedené v parametrech příkazu
        plots= csv.reader(csvfile, delimiter=';')
        line_count = 0
        x=[]
        y_RS=[]
        y_BC=[]
        z=[]
        for row in plots:
            if line_count == 0:
                line_count += 1
                x_str = row[0]
                y_RS_str = row[1]
                y_BC_str = row[2]
            else:
                x.append(float(row[0]))
                y_RS.append(float(row[1]))
                y_BC.append(float(row[2]))
                z.append(float(i))
                line_count += 1
        if i == 1:
            ax.plot(x, z, y_RS, 'b-',  label=y_RS_str)
            ax.plot(x, z, y_BC, 'r-', label=y_BC_str)
        else:
            ax.plot(x, z, y_RS, 'b-')
            ax.plot(x, z, y_BC, 'r-')
                
ax.legend(loc="upper left")
ax.set_title('Fractal analysis')
ax.set_xlabel(x_str)
ax.set_ylabel('number of file')

plt.show()
