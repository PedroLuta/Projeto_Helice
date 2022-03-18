import numpy as nup
import matplotlib.pyplot as plt



def PARSEC_main_func(ans, x):
    return ans[0]*(x**0.5) + ans[1]*(x**1.5) + ans[2]*(x**2.5) + ans[3]*(x**3.5) + ans[4]*(x**4.5) + ans[5]*(x**5.5)

def PARSEC_derivative1(ans, x):
    return 0.5*ans[0]*(x**(-0.5)) + 1.5*ans[1]*(x**0.5) + 2.5*ans[2]*(x**1.5) + 3.5*ans[3]*(x**2.5) + 4.5*ans[4]*(x**3.5) + 5.5*ans[5]*(x**4.5)

def PARSEC_derivative2(ans, x):
    return -0.25*ans[0]*(x**(-1.5)) + 0.75*ans[1]*(x**(-0.5)) + 3.75*ans[2]*(x**0.5) + 8.75*ans[3]*(x**1.5) + 15.75*ans[4]*(x**2.5) + 24.75*ans[5]*(x**3.5) 

#def PARSEC_derivative3(ans, x):
#    return 0.375*ans[0]*(x**(-2.5)) - 0.375*ans[1]*(x**(-1.5)) + 1.875*ans[2]*(x**(-0.5)) + 13.125*ans[3]*(x**0.5) + 39.375*ans[4]*(x**1.5) + 86.625*ans[5]*(x**2.5)

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
    return x_lower, y_lower, x_upper, y_upper

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
    while i < len(y_u):
        camber_line = nup.append(camber_line, (y_u[i] + y_lo[i])/2) 
        i += 1
    return camber_line, x_u

def write_PARSEC_zcoor(x_curve, ans):
    i = 0
    y_poli = nup.array([])
    while i < len(x_curve):
        y_poli = nup.append(y_poli, PARSEC_main_func(ans, x_curve[i]))
        i += 1
    return y_poli

#def newt_rhap_PARSEC_main(ans, tolerance = 0.0001, init = 0.01):
#    xi = init
#    while True:
#        delta = PARSEC_main_func(ans, xi)/PARSEC_derivative1(ans, xi)
#        xn = xi - delta
#        if abs(xn - xi) < tolerance:
#            break
#        xi = xn
#    
#    return xn

def newt_rhap_PARSEC_first(ans, tolerance = 0.0001, init = 0.01):
    xi = init
    #i = 0
    while True:
        delta = PARSEC_derivative1(ans, xi)/PARSEC_derivative2(ans, xi)
        xn = xi - delta
        if abs(xn - xi) < tolerance:
            break
        xi = xn
        #i += 1
    #print(i)
    return xn

#def newt_rhap_PARSEC_second(ans, tolerance = 0.0001, init = 0.01):
#    xi = init
#    while True:
#        delta = PARSEC_derivative2(ans, xi)/PARSEC_derivative3(ans, xi)
#        xn = xi - delta
#        if abs(xn - xi) < tolerance:
#            break
#        xi = xn
#    
#    return xn

def get_PARSEC_parameters(answer_lo, answer_u):
    rlo = (answer_lo[0]**2)/2
    rup = (answer_u[0]**2)/2
    zte = (nup.sum(answer_lo) + nup.sum(answer_u))/2 
    dzte = nup.sum(answer_u) - nup.sum(answer_lo)
    xlo = newt_rhap_PARSEC_first(answer_lo, tolerance = 0.00000001, init = 0.0001)
    xup = newt_rhap_PARSEC_first(answer_u, tolerance = 0.00000001, init = 0.0001)
    zlo = PARSEC_main_func(answer_lo, xlo)
    zup = PARSEC_main_func(answer_u, xup)
    zxxlo = PARSEC_derivative2(answer_lo, xlo)
    zxxup = PARSEC_derivative2(answer_u, xup)
    ate = (nup.degrees(nup.arctan(PARSEC_derivative1(answer_u, 1))) + nup.degrees(nup.arctan(PARSEC_derivative1(answer_lo, 1))))/2
    bte = nup.degrees(nup.arctan(PARSEC_derivative1(answer_lo, 1))) - nup.degrees(nup.arctan(PARSEC_derivative1(answer_u, 1)))

    return [ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo]




#x, y = read_foil('airfoil.txt', header_lines = 0)
#x, y = read_foil('airfoil_s1223.txt', header_lines = 1)
#x, y = read_foil('airfoil_clarky.txt', header_lines = 0)
x, y = read_foil('airfoil_sunnysky.dat', header_lines = 2)
x_lo, y_lo, x_u, y_u = split_upper_lower(x, y)
#camber_line, x_camb = get_camber_line(y_lo, y_u, x_u)


answer_lo, answer_u = get_PARSEC_curve_coefficients(x_lo, x_u, y_lo, y_u)
y_poli_lo = write_PARSEC_zcoor(x_lo, answer_lo) 
y_poli_u = write_PARSEC_zcoor(x_u, answer_u) 
PARSEC_params = get_PARSEC_parameters(answer_lo, answer_u)
print(PARSEC_params)





fig, ax = plt.subplots()
#ax.plot(x_camb, camber_line)
ax.plot(x_u, y_u, 'r', label = "original")
ax.plot(x_lo, y_lo, 'r') 
ax.plot(x_lo, y_poli_lo, 'b--', label = "PARSEC")
ax.plot(x_u, y_poli_u, 'b--') 
ax.legend()
ax.xaxis.grid(True, which='major')
plt.show()
