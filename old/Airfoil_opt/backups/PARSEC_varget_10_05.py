import numpy as nup
import matplotlib.pyplot as plt

def read_foil(airfoil_file, header_lines):
    with open(airfoil_file, 'r') as inp:
        x, y = nup.loadtxt(airfoil_file, skiprows = header_lines, unpack = True)
    return x, y

def split_upper_lower(x, y):
    i = 0
    diff_remember = 999
    while i < len(x):
        if x[i] > diff_remember:
            if len(x)%2 == 0:
                x_upper = x[:i-1]
                x_lower = x[i-1:]
                y_upper = y[:i-1]
                y_lower = y[i-1:]
                x_upper = nup.flip(x_upper)
                y_upper = nup.flip(y_upper)
                #x_upper = nup.array([0])
                #y_upper = nup.array([0])
                #x_lower = nup.array([0])
                #y_lower = nup.array([0])
                #ii = 0
                #while ii < len(x_upperb):
                #    #print(ii)
                #    x_upper = nup.append(x_upper, x_upperb[ii]) #x_upper.append(x_upperb[ii])
                #    y_upper = nup.append(y_upper, y_upperb[ii]) #y_upper.append(y_upperb[ii])
                #    x_lower = nup.append(x_lower, x_lowerb[ii]) #x_lower.append(x_lowerb[ii])
                #    y_lower = nup.append(y_lower, y_lowerb[ii]) #y_lower.append(y_lowerb[ii])
                #    ii += 1
                break
            else:
                x_upper = x[:i]
                x_lower = x[i-1:]
                y_upper = y[:i]
                y_lower = y[i-1:]
                x_upper = nup.flip(x_upper)
                y_upper = nup.flip(y_upper)
                break
        diff_remember = x[i]
        i += 1 
    return x_lower, y_lower, x_upper, y_upper

x, y = read_foil('airfoil.txt', header_lines = 0)
#x, y = read_foil('airfoil_sunnysky.dat', header_lines = 2)
#x, y = read_foil('airfoil_clarky.txt', header_lines = 0)
x_lo, y_lo, x_u, y_u = split_upper_lower(x, y)
#print(x_lo)
#print(y_lo)
#print(x_u)
#print(y_u)
#print(len(x_lo))
#print(len(y_lo))
#print(len(x_u))
#print(len(y_u))
i = 0
camber_line = nup.array([])
while i < len(y_u):
    camber_line = nup.append(camber_line, (y_u[i] + y_lo[i])/2) #camber_line.append((y_u[i] + y_lo[i])/2)
    i += 1
#print(camber_line)
thick_u = y_u - camber_line
thick_lo = camber_line - y_lo

#Zup1, Xup1 = y_u[nup.argmax(thick_u)], x_u[nup.argmax(thick_u)] # max(y_u - camber_line),  #y_u.index(max(y_u))
#Zlo1, Xlo1 = y_lo[nup.argmax(thick_lo)], x_lo[nup.argmax(thick_lo)] #y_lo.index(max(y_lo))
#print(Zup1)
#print(Xup1)
#print(Zlo1)
#print(Xlo1)

Zup2, Xup2 = y_u[nup.argmax(y_u)], x_u[nup.argmax(y_u)] # max(y_u - camber_line),  #y_u.index(max(y_u))
Zlo2, Xlo2 = y_lo[nup.argmin(y_lo)], x_lo[nup.argmin(y_lo)] #y_lo.index(max(y_lo))
#print(Zup2)
#print(Xup2)
#print(Zlo2)
#print(Xlo2)

fig, ax = plt.subplots()
ax.plot(x_u, camber_line)
ax.plot(x_u, y_u)
ax.plot(x_lo, y_lo)
ax.xaxis.grid(True, which='major')
plt.show()