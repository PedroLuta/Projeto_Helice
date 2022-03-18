import numpy as nup

def read_foil(airfoil_file, header_lines):
    x, y = nup.loadtxt(airfoil_file, skiprows = header_lines, unpack = True)
    return x, y

def count_header_lines(file):
    with open(file, 'r') as inp:
        i = 0
        while True:
            try:
                content = inp.readline()
                cont = content.split()
                _ = float(cont[0])
                break
            except:
                i += 1
    return i

airfoil = "airfoils\\sunnysky.dat"

x, y = read_foil(airfoil_file = airfoil, header_lines = count_header_lines(airfoil))
x_delta = nup.max(x) - nup.min(x)
y_delta = nup.max(y) - nup.min(y)

V = nup.array([[0, 0]])

i = 1
invert = False
while i < len(x) - 1:
    tanvec_x = (x[i + 1] - x[i - 1])/2
    tanvec_y = ((y[i + 1] - y[i - 1])/2)*(x_delta/y_delta)
    if invert:
        tanvec_x = -tanvec_x
        tanvec_y = -tanvec_y
    normvec_x = -tanvec_y
    normvec_y = tanvec_x
    placeholder = ((normvec_x**2) + (normvec_y**2))**0.5
    normvec_x /= placeholder
    normvec_y /= placeholder
    V = nup.append(V, [[normvec_x, normvec_y]], axis = 0)
    i += 1

x_offset = nup.array([])
y_offset = nup.array([])

i = 1
offset = 0.01
while i < len(x) - 1:
    x_offset = nup.append(x_offset, x[i] + V[i][0]*offset)
    y_offset = nup.append(y_offset, y[i] + V[i][1]*offset/(x_delta/y_delta))
    i += 1 