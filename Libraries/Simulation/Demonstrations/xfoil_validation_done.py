import matplotlib.pyplot as plt
import xfoil_interface 
from importing import *

class Rfoil:
    def __init__(self, name, chord, xcor, ycor, polars = [], Re_list = []):
        self.name = name
        self.chord = chord
        self.xcor = xcor[:]
        self.ycor = ycor[:]
        self.polars = polars[:]
        self.Re_list = Re_list[:]
            
def pack_airfoils():
    with open(r"C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil_validation\soarTech8_data\ALL.DAP", 'r') as inp:
        x_dist, y_dist, foil_list = [], [], []
        first = True
        for line in inp:
            if line == '\n':
                temp = Rfoil(name, chord, x_dist, y_dist)
                foil_list.append(temp)
                first = True
                x_dist.clear()
                y_dist.clear()
                continue
            contents = line.split()
            if first:
                name = contents[0]
                chord = float(contents[-1])
                first = False
                continue
            x_dist.append(float(contents[0]))
            y_dist.append(float(contents[1]))

    with open(r"C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil_validation\soarTech8_data\ALL.PD", 'r') as inp:
        while True:
            polars = []
            Re_list = []
            test = inp.readline()
            if test == "":
                break
            elif test == "\n":
                continue
            name = test.split()[1]
            inp.readline()
            curves_num = int(inp.readline())
            for _ in range(curves_num):
                cl_dist, cd_dist, alpha_dist = [], [], []
                Reynolds = float(inp.readline())
                Re_list.append(Reynolds)
                points_num = int(inp.readline())
                for _ in range(points_num):
                    contents = inp.readline().split()
                    alpha_dist.append(float(contents[2]))
                    cd_dist.append(float(contents[1]))
                    cl_dist.append(float(contents[0]))
                polars.append([alpha_dist, cl_dist, cd_dist])
            inp.readline()
            inp.readline()
            for foil in foil_list:
                if foil.name == name:
                    foil.Re_list = Re_list[:]
                    foil.polars = polars[:]
                    
    return foil_list

def write_foil_default(foil, coord_file = 'airfoils\\airfoil.txt'):
    x_dist = foil.xcor
    y_dist = foil.ycor
    length = len(x_dist)
    try:
        with open(coord_file, 'w') as inp:
            for i in range(length):
                inp.write(f'{x_dist[i]} {y_dist[i]}\n')
        return True
    except:
        return False

foils = pack_airfoils()

for foil in foils:
    write_foil_default(foil)
    a, cl, cd = xfoil_interface.get_curve_com_default(foil.Re_list[0], -2, 10, 0.5, timeout = 30)
    plt.plot(cl, cd, label = f"XFoil Default")
    plt.plot(foil.polars[0][1][:], foil.polars[0][2][:], label = f'{foil.Re_list[0]} experimental')
    plt.legend()
    plt.title(foil.name)
    plt.show()