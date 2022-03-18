import numpy as nup
import matplotlib.pyplot as plt

def read_foil(airfoil_file, header_lines):
    x, y = nup.loadtxt(airfoil_file, skiprows = header_lines, unpack = True)
    return x, y

def split_upper_lower(x, y):
    i = 0
    diff_remember = 999
    while i < len(x):
        if x[i] > diff_remember:
            if len(x)%2 == 0:
                x_lower = x[i-1:]
                y_lower = y[i-1:]
                x_upper = nup.flip(x[:i-1])
                y_upper = nup.flip(y[:i-1])
                break
            else:
                x_lower = x[i-1:]
                y_lower = y[i-1:]
                x_upper = nup.flip(x[:i])
                y_upper = nup.flip(y[:i])
                break
        diff_remember = x[i]
        i += 1 
    #while True:
    #    #print(x_lower)
    #    #print(x_upper)
    #    print(int(round(len(x_lower)/2)))
    #    print(int(round(len(x_upper)/2)))
    #    if len(x_upper) > len(x_lower):
    #        x_upper = nup.delete(x_upper, int(round(len(x_upper)/2)), 0)
    #        y_upper = nup.delete(y_upper, int(round(len(y_upper)/2)), 0)
    #        #i += 1
    #    elif len(x_upper) < len(x_lower):
    #        x_lower = nup.delete(x_lower, int(round(len(x_lower)/2)), 0)
    #        y_lower = nup.delete(y_lower, int(round(len(y_lower)/2)), 0)
    #        #i += 1
    #    else:
    #        break
    return x_lower, y_lower, x_upper, y_upper

#def radians_between_vectors(a, b):
#    unit_vector_1 = a/nup.linalg.norm(a)
#    unit_vector_2 = b/nup.linalg.norm(b)
#    dot_product = nup.dot(unit_vector_1, unit_vector_2)
#    angle_rad = nup.arccos(dot_product)
#    return angle_rad
#
#def degrees_between_vectors(a, b):
#    return nup.degrees(radians_between_vectors(a, b))

def get_PARSEC_x_mat(x_curve):
    i = 0
    x_mat = nup.ones((x_curve.size, 6))
    while i < len(x_curve):
        ii = 0
        while ii < 6:
            x_mat[i][ii] = x_curve[i]**(ii + 0.5)
            ii += 1
        i += 1
    return x_mat

def get_PARSEC_curve_coefficients(x_lo, x_u, y_lo, y_u):
    x_mat_lo = get_PARSEC_x_mat(x_lo)
    x_mat_u = get_PARSEC_x_mat(x_u)
    y_mat_lo = nup.transpose(y_lo)
    y_mat_u = nup.transpose(y_u)
    x_mat_lot = nup.transpose(x_mat_lo)
    x_mat_ut = nup.transpose(x_mat_u)
    x_times_xt_lo = nup.matmul(x_mat_lot, x_mat_lo)
    x_times_xt_u = nup.matmul(x_mat_ut, x_mat_u)
    y_times_xt_lo = nup.matmul(x_mat_lot, y_mat_lo)
    y_times_xt_u = nup.matmul(x_mat_ut, y_mat_u)
    x_inverse_lo = nup.linalg.inv(x_times_xt_lo)
    x_inverse_u = nup.linalg.inv(x_times_xt_u)
    answer_lo = nup.matmul(y_times_xt_lo, x_inverse_lo)
    answer_u = nup.matmul(y_times_xt_u, x_inverse_u)
    return answer_lo, answer_u

def get_camber_line(y_lo, y_u, x_u):
    i = 0
    camber_line = nup.array([])
    #while True:
    #    if len(y_u) > len(y_lo):
    #        y_u = nup.delete(y_u, int(round(len(y_u)/2)), 0)
    #        x_u = nup.delete(x_u, int(round(len(x_u)/2)), 0)
    #    elif len(y_u) < len(y_lo):
    #        y_lo = nup.delete(y_lo, int(round(len(y_lo)/2)), 0)
    #    else:
    #        break
    while i < len(y_u):
        camber_line = nup.append(camber_line, (y_u[i] + y_lo[i])/2) 
        i += 1
    return camber_line, x_u

def write_PARSEC_zcoor(x_curve, a1, a2, a3, a4, a5, a6):
    i = 0
    y_poli = nup.array([])
    while i < len(x_curve):
        y_poli = nup.append(y_poli, a1*(x_curve[i]**0.5) + a2*(x_curve[i]**1.5) + a3*(x_curve[i]**2.5) + a4*(x_curve[i]**3.5) + a5*(x_curve[i]**4.5) + a6*(x_curve[i]**5.5))
        i += 1
    return y_poli

def get_PARSEC_parameters(answer_lo, answer_u):
    pass




#x, y = read_foil('airfoil.txt', header_lines = 0)
x, y = read_foil('airfoil_s1223.txt', header_lines = 1)
#x, y = read_foil('airfoil_clarky.txt', header_lines = 0)
#x, y = read_foil('airfoil_sunnysky.dat', header_lines = 2)
x_lo, y_lo, x_u, y_u = split_upper_lower(x, y)
camber_line, x_camb = get_camber_line(y_lo, y_u, x_u)


answer_lo, answer_u = get_PARSEC_curve_coefficients(x_lo, x_u, y_lo, y_u)
R_le_lo = (answer_lo[0]**2)/2
R_le_u = (answer_u[0]**2)/2
zte = (nup.sum(answer_lo) + nup.sum(answer_u))/2 
dzte = nup.sum(answer_u) - nup.sum(answer_lo)

y_poli_lo = write_PARSEC_zcoor(x_lo, answer_lo[0], answer_lo[1], answer_lo[2], answer_lo[3], answer_lo[4], answer_lo[5]) 
y_poli_u = write_PARSEC_zcoor(x_u, answer_u[0], answer_u[1], answer_u[2], answer_u[3], answer_u[4], answer_u[5]) 





fig, ax = plt.subplots()
ax.plot(x_camb, camber_line)
ax.plot(x_u, y_u)
ax.plot(x_lo, y_lo)
ax.plot(x_lo, y_poli_lo)
ax.plot(x_u, y_poli_u)
ax.xaxis.grid(True, which='major')
plt.show()
