from importing import *
from PARSEC_functions import *

def apply_thickness_y(y, thick):
    for i in range(len(y)):
        y[i] *= thick
    return y

def apply_chord_x(x, chord):
    for i in range(len(x)):
        x[i] *= chord
    return x

def apply_offset_x(x, offset):
    for i in range(len(x)):
        x[i] += offset
    return x

def apply_offset_y(y, offset):
    for i in range(len(y)):
        y[i] += offset
    return y

def apply_rotation_deg(x, y, angle):
    return apply_rotation_rad(x, y, math.radians(angle))

def apply_rotation_rad(x, y, angle):
    rot_mat = [[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]]
    x_rot = nup.array([])
    y_rot = nup.array([])
    for i in range(len(x)):
        case = nup.array([[x[i]], [y[i]]])
        rotated = nup.matmul(rot_mat, case)
        x_rot = nup.append(x_rot, rotated[0])
        y_rot = nup.append(y_rot, rotated[1])
    return x_rot, y_rot

max_thickness_location = 0.14
chord = 2
extra_thickness = 1 #Make the airfoil more thick than it is

airfoil_file = r'airfoils\\sunnysky.dat'
x, y =  read_foil(airfoil_file, header_lines = count_header_lines(airfoil_file))
plt.plot(x, y)
y = apply_thickness_y(y, extra_thickness)
x = apply_chord_x(x, chord)
x = apply_offset_x(x, -max_thickness_location*chord) #OFFSET ALWAYS AFTER EVERYTHING
plt.plot(x, y)
x, y = apply_rotation_deg(x, y, 30)
plt.plot(x, y)
plt.show()
