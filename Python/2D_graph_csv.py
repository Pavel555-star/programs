import matplotlib.pyplot as plt
import csv

x=[]
y_RS=[]
y_BC=[]

with open('output.csv', 'r') as csvfile: # Načte a zpracuje soubor output.csv
    plots= csv.reader(csvfile, delimiter=';')
    line_count = 0
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
            line_count += 1

plt.plot(x, y_RS,'b-',  label=y_RS_str)
plt.plot(x, y_BC, 'r-', label=y_BC_str)

plt.legend(loc="upper left")
plt.title('Fractal analysis')
plt.xlabel(x_str)
plt.ylabel(' ')
plt.show()
