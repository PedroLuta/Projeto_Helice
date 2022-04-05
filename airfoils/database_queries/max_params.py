import numpy as nup
import matplotlib.pyplot as plt
import os



def PARSEC_main_func(ans, x):
    #print("main")
    if x < 0:
        return -1337
    return ans[0]*(x**0.5) + ans[1]*(x**1.5) + ans[2]*(x**2.5) + ans[3]*(x**3.5) + ans[4]*(x**4.5) + ans[5]*(x**5.5)

def PARSEC_derivative1(ans, x):
    #print("der1")
    if x < 0:
        return -1337
    return 0.5*ans[0]*(x**(-0.5)) + 1.5*ans[1]*(x**0.5) + 2.5*ans[2]*(x**1.5) + 3.5*ans[3]*(x**2.5) + 4.5*ans[4]*(x**3.5) + 5.5*ans[5]*(x**4.5)

def PARSEC_derivative2(ans, x):
    #print("der1")
    if x < 0:
        return -1337
    return -0.25*ans[0]*(x**(-1.5)) + 0.75*ans[1]*(x**(-0.5)) + 3.75*ans[2]*(x**0.5) + 8.75*ans[3]*(x**1.5) + 15.75*ans[4]*(x**2.5) + 24.75*ans[5]*(x**3.5) 

def read_foil(airfoil_file, header_lines):
    #print("read_foil")
    try:
        x, y = nup.loadtxt(airfoil_file, skiprows = header_lines, unpack = True)
    except:
        x, y = nup.array([0, 0]), nup.array([0, 0])
    return x, y

def split_upper_lower(x, y):
    #print("split")
    i = 0
    diff_remember = 999
    while i < len(x):
        if x[i] > diff_remember:
            index0 = i - 1
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
    delta = y[index0]
    y_upper = y_upper - delta
    y_lower = y_lower - delta
    delta = x[index0]
    x_upper = x_upper - delta
    x_lower = x_lower - delta
    return x_lower, y_lower, x_upper, y_upper

def get_PARSEC_x_mat(x_curve):
    #print("x_mat")
    i = 0
    x_mat = nup.ones((x_curve.size, 6))
    while i < len(x_curve):
        if x_curve[i] < 0:
            x_curve[i] = 0
        ii = 0
        while ii < 6:
            x_mat[i][ii] = x_curve[i]**(ii + 0.5)
            ii += 1
        i += 1
    return x_mat

def get_PARSEC_curve_coefficients(x_lo, x_u, y_lo, y_u):
    #print("coeffs")
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

def newt_rhap_PARSEC_first(ans, tolerance = 0.0001, init = 0.01, cutbreak = False):
    #print("newt")
    xi = init
    retry = False
    while True:
        if xi < 0 or xi > 1:
            retry = True
            break
        delta = PARSEC_derivative1(ans, xi)/PARSEC_derivative2(ans, xi)
        xn = xi - delta
        if abs(xn - xi) < tolerance:
            break
        xi = xn
    if retry:
        if cutbreak:
            return -1337
        xn = newt_rhap_PARSEC_first(ans, tolerance = 0.0001, init = 0.3, cutbreak = True)
    return xn

def get_PARSEC_parameters(answer_lo, answer_u):
    #print("parameters")
    rlo = (answer_lo[0]**2)/2
    rup = (answer_u[0]**2)/2
    zte = (nup.sum(answer_lo) + nup.sum(answer_u))/2 
    dzte = nup.sum(answer_u) - nup.sum(answer_lo)
    xlo = newt_rhap_PARSEC_first(answer_lo, tolerance = 0.00000001, init = 0.0001)
    xup = newt_rhap_PARSEC_first(answer_u, tolerance = 0.00000001, init = 0.0001)
    if (xlo == -1337) or (xup == -1337):
        return [0, 1, 0, 0, 0, 0.5, 0.1, 0, 0, 0.5, -0.1, 2]
    zlo = PARSEC_main_func(answer_lo, xlo)
    zup = PARSEC_main_func(answer_u, xup)
    zxxlo = PARSEC_derivative2(answer_lo, xlo)
    zxxup = PARSEC_derivative2(answer_u, xup)
    ate = (nup.degrees(nup.arctan(PARSEC_derivative1(answer_u, 1))) + nup.degrees(nup.arctan(PARSEC_derivative1(answer_lo, 1))))/2
    bte = nup.degrees(nup.arctan(PARSEC_derivative1(answer_lo, 1))) - nup.degrees(nup.arctan(PARSEC_derivative1(answer_u, 1)))

    return [ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo]

def count_header_lines(file):
    #print("headers")
    with open(file, 'r') as inp:
        i = 0
        while True and (i < 20):
            try:
                content = inp.readline()
                cont = content.split()
                _, _, _ = float(cont[-1]), float(cont[0]), float(cont[1])
                break
            except:
                i += 1
    if i == 20:
        return 1
    return i

def curve_cross(aup, alo):
    for x in nup.linspace(0, 1, 101):
        zu = PARSEC_main_func(aup, x)
        zlo = PARSEC_main_func(alo, x)
        if zlo > zu:
            return True
    return False

list_at = os.listdir(r"C:\Users\PEDRO\Desktop\IC DE HÉLICE\git_control\git_prop\airfoils\query_airfoiltools")
list_uiuc = os.listdir(r"C:\Users\PEDRO\Desktop\IC DE HÉLICE\git_control\git_prop\airfoils\query_UIUC")

#print(len(list_at))
#print(len(list_uiuc))


ate_max, ate_min, bte_max, zte_min, zte_max, dzte_max, xup_max, xup_min, xlo_max, xlo_min, zup_max, zlo_min, \
     rup_max, rlo_max, zxxup_max, zxxup_min, zxxlo_max, zxxlo_min, zup_min, zlo_max = 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, -10, 0, 0, 10, 10, -10

i = 0

for airfoil in list_at:
    try:
        #print(airfoil)
        airfoil_file = r'C:\Users\PEDRO\Desktop\IC DE HÉLICE\git_control\git_prop\airfoils\query_airfoiltools\\' + airfoil
        x, y = read_foil(airfoil_file, header_lines = count_header_lines(airfoil_file))
        if (nup.all(x == 0)) and (nup.all(y == 0)):
            #print(airfoil)
            continue
        #    continue
        x_lo, y_lo, x_u, y_u = split_upper_lower(x, y)
        answer_lo, answer_u = get_PARSEC_curve_coefficients(x_lo, x_u, y_lo, y_u)
        #if curve_cross(answer_u, answer_lo):
        #    continue
        PARSEC_params = get_PARSEC_parameters(answer_lo, answer_u)
    except:
        continue
    #if airfoil == "fx79w660a-il.dat" or airfoil == 'dbln526-il.dat' or airfoil == 'ah93w480b-il.dat':
    #    continue
    i += 1
    #print(i)
    ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo = PARSEC_params[0], PARSEC_params[1], PARSEC_params[2], PARSEC_params[3], PARSEC_params[4], \
        PARSEC_params[5], PARSEC_params[6], PARSEC_params[7], PARSEC_params[8], PARSEC_params[9], PARSEC_params[10], PARSEC_params[11]
    #print([ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo])
    if ate > ate_max:
        ate_max = ate
        print(airfoil + f" {ate} ate (aft)")
    if ate < ate_min:
        ate_min = ate
        print(airfoil + f" {ate} ate (aft)")
    if (bte > bte_max) and (dzte >= 0) and (bte <= 45):
        bte_max = bte
        print(airfoil + f" {bte} bte (aft)")
    if zte < zte_min:
        zte_min = zte
        print(airfoil + f" {zte} zte (aft)")
    if zte > zte_max:
        zte_max = zte
        print(airfoil + f" {zte} zte (aft)")
    if dzte > dzte_max:
        dzte_max = dzte
        print(airfoil + f" {dzte} dzte (aft)")
    if xup > xup_max:
        xup_max = xup
        print(airfoil + f" {xup} xup (aft)")
    if xup < xup_min:
        xup_min = xup
        print(airfoil + f" {xup} xup (aft)")
    if (xlo > xlo_max) and (xlo < 1):
        xlo_max = xlo
        print(airfoil + f" {xlo} xlo (aft)")
    if xlo < xlo_min:
        xlo_min = xlo
        print(airfoil + f" {xlo} xlo (aft)")
    if zup > zup_max:
        zup_max = zup
        print(airfoil + f" {zup} zup (aft)")
    if (zlo < zlo_min) and (zlo < 1):
        zlo_min = zlo
        print(airfoil + f" {zlo} zlo (aft)")
    if zup < zup_min:
        zup_min = zup
        print(airfoil + f" {zup} zup (aft)")
    if (zlo > zlo_max) and (zlo < 1):
        zlo_max = zlo
        print(airfoil + f" {zlo} zlo (aft)")
    if rup > rup_max:
        rup_max = rup
        print(airfoil + f" {rup} rup (aft)")
    if rlo > rlo_max:
        rlo_max = rlo
        print(airfoil + f" {rlo} rlo (aft)")
    if zxxup > zxxup_max:
        zxxup_max = zxxup
        print(airfoil + f" {zxxup} zxxup (aft)")
    if zxxup < zxxup_min:
        zxxup_min = zxxup
        print(airfoil + f" {zxxup} zxxup (aft)")
    if zxxlo > zxxlo_max:
        zxxlo_max = zxxlo
        print(airfoil + f" {zxxlo} zxxlo (aft)")
    if zxxlo < zxxlo_min:
        zxxlo_min = zxxlo
        print(airfoil + f" {zxxlo} zxxlo (aft)")

for airfoil in list_uiuc:
    try:
        #print(airfoil)
        airfoil_file = r'C:\Users\PEDRO\Desktop\IC DE HÉLICE\git_control\git_prop\airfoils\query_UIUC\\' + airfoil
        x, y = read_foil(airfoil_file, header_lines = count_header_lines(airfoil_file))
        if (nup.all(x == 0)) and (nup.all(y == 0)):
            #print(airfoil)
            continue
        #    continue
        x_lo, y_lo, x_u, y_u = split_upper_lower(x, y)
        answer_lo, answer_u = get_PARSEC_curve_coefficients(x_lo, x_u, y_lo, y_u)
        #if curve_cross(answer_u, answer_lo):
        #    continue
        PARSEC_params = get_PARSEC_parameters(answer_lo, answer_u)
    except:
        continue
    i += 1
    #print(i)
    ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo = PARSEC_params[0], PARSEC_params[1], PARSEC_params[2], PARSEC_params[3], PARSEC_params[4], \
        PARSEC_params[5], PARSEC_params[6], PARSEC_params[7], PARSEC_params[8], PARSEC_params[9], PARSEC_params[10], PARSEC_params[11]
    #print([ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo])
    if ate > ate_max:
        ate_max = ate
        print(airfoil + " ate (uiuc)")
    if ate < ate_min:
        ate_min = ate
        print(airfoil + " ate (uiuc)")
    if (bte > bte_max) and (dzte >= 0) and (bte <= 45):
        bte_max = bte
        print(airfoil + " bte (uiuc)")
    if zte < zte_min:
        zte_min = zte
        print(airfoil + " zte (uiuc)")
    if zte > zte_max:
        zte_max = zte
        print(airfoil + " zte (uiuc)")
    if dzte > dzte_max:
        dzte_max = dzte
        print(airfoil + " dzte (uiuc)")
    if xup > xup_max:
        xup_max = xup
        print(airfoil + " xup (uiuc)")
    if xup < xup_min:
        xup_min = xup
        print(airfoil + " xup (uiuc)")
    if (xlo > xlo_max) and (xlo < 1):
        xlo_max = xlo
        print(airfoil + " xlo (uiuc)")
    if xlo < xlo_min:
        xlo_min = xlo
        print(airfoil + " xlo (uiuc)")
    if zup > zup_max:
        zup_max = zup
        print(airfoil + " zup (uiuc)")
    if (zlo < zlo_min) and (zlo < 1):
        zlo_min = zlo
        print(airfoil + " zlo (uiuc)")
    if zup < zup_min:
        zup_min = zup
        print(airfoil + " zup (uiuc)")
    if (zlo > zlo_max) and (zlo < 1):
        zlo_max = zlo
        print(airfoil + " zlo (uiuc)")
    if rup > rup_max:
        rup_max = rup
        print(airfoil + " rup (uiuc)")
    if rlo > rlo_max:
        rlo_max = rlo
        print(airfoil + " rlo (uiuc)")
    if zxxup > zxxup_max:
        zxxup_max = zxxup
        print(airfoil + " zxxup (uiuc)")
    if zxxup < zxxup_min:
        zxxup_min = zxxup
        print(airfoil + " zxxup (uiuc)")
    if zxxlo > zxxlo_max:
        zxxlo_max = zxxlo
        print(airfoil + " zxxlo (uiuc)")
    if zxxlo < zxxlo_min:
        zxxlo_min = zxxlo
        print(airfoil + " zxxlo (uiuc)")

min_max = [ate_max, ate_min, bte_max, zte_min, zte_max, dzte_max, xup_max, \
    xup_min, xlo_max, xlo_min, zup_max, zup_min, zlo_max, zlo_min, \
        rup_max, rlo_max, zxxup_max, zxxup_min, zxxlo_max, zxxlo_min]
print('[ate_max, ate_min, bte_max, zte_min, zte_max, dzte_max, xup_max, xup_min, xlo_max, xlo_min, zup_max, zup_min, zlo_max, zlo_min, rup_max, rlo_max, zxxup_max, zxxup_min, zxxlo_max, zxxlo_min]')
print(min_max)
print(f'Airfoils calculated: {i}')


