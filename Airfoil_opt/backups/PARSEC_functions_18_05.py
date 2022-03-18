import random
import numpy as nup
import matplotlib.pyplot as plt
import os

from numpy.linalg.linalg import LinAlgError

#general functions
def PARSEC_main_func(ans, x):
    return ans[0]*(x**0.5) + ans[1]*(x**1.5) + ans[2]*(x**2.5) + ans[3]*(x**3.5) + ans[4]*(x**4.5) + ans[5]*(x**5.5)

def PARSEC_derivative1(ans, x):
    return 0.5*ans[0]*(x**(-0.5)) + 1.5*ans[1]*(x**0.5) + 2.5*ans[2]*(x**1.5) + 3.5*ans[3]*(x**2.5) + 4.5*ans[4]*(x**3.5) + 5.5*ans[5]*(x**4.5)

def PARSEC_derivative2(ans, x):
    return -0.25*ans[0]*(x**(-1.5)) + 0.75*ans[1]*(x**(-0.5)) + 3.75*ans[2]*(x**0.5) + 8.75*ans[3]*(x**1.5) + 15.75*ans[4]*(x**2.5) + 24.75*ans[5]*(x**3.5)

def newt_rhap_PARSEC_first(ans, tolerance = 0.0001, init = 0.01, cutbreak = False):
    xi = init
    retry = False
    while True:
        if xi < 0:
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
        xn = newt_rhap_PARSEC_first(ans, tolerance = tolerance, init = 0.3, cutbreak = True)
    return xn

def curve_cross(aup, alo):
    for x in nup.linspace(0, 1, 101):
        zu = PARSEC_main_func(aup, x)
        zlo = PARSEC_main_func(alo, x)
        if zlo > zu:
            return True
    return False

def plot_foil(x_u, x_lo, y_u, y_lo):
    fig, ax = plt.subplots()
    ax.plot(x_u, y_u, 'r')
    ax.plot(x_lo, y_lo, 'r')
    ax.xaxis.grid(True, which='major')
    plt.show()

def plot_file(file):
    x, y = read_foil(file, header_lines = count_header_lines(file))
    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.xaxis.grid(True, which='major')
    plt.show()

def plot_foilxy(x, y):
    fig, ax = plt.subplots()
    ax.plot(x, y, 'r')
    ax.xaxis.grid(True, which='major')
    plt.show()

def get_camber_line(y_lo, y_u, x_u):
    i = 0
    camber_line = nup.array([])
    while i < len(y_u):
        camber_line = nup.append(camber_line, (y_u[i] + y_lo[i])/2) 
        i += 1
    return camber_line, x_u









#generate foil functions

def PARSEC_randvargen():
    ate = random.uniform(-45, 45) #-30, 40
    bte = random.uniform(0, 180) #0, 140
    zte = random.uniform(-0.1, 0.1) #-0.07, 0.09
    dzte = random.uniform(0, 0.05) #0, 0.26
    xup = random.uniform(0, 0.5) #0.07, 0.59
    xlo = random.uniform(0, 0.5) #0, 1
    zup = random.uniform(0.0001, 0.35) #0.02, 0.35
    zlo = random.uniform(-0.35, -0.0001) #-0.32, 0.08
    rup = random.uniform(0, 0.28) 
    rlo = random.uniform(0, 0.18) 
    zxxup = random.uniform(-5, 0) 
    zxxlo = random.uniform(0, 100) #0, 24
    return nup.array([ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo])

#def PARSEC_randvargen():
#    ate = random.uniform(-30, 11) #-30, 40
#    bte = random.uniform(0, 45) #0, 140
#    zte = random.uniform(-0.02, 0.04) #-0.07, 0.09
#    dzte = random.uniform(0, 0.05) #0, 0.26
#    xup = random.uniform(0.11, 0.55) #0.07, 0.59
#    xlo = random.uniform(0.008, 0.6) #0, 1
#    zup = random.uniform(0.01, 0.35) #0.02, 0.35
#    zlo = random.uniform(-0.33, 0.08) #-0.32, 0.08
#    rup = random.uniform(0, 0.28) 
#    rlo = random.uniform(0, 0.18) 
#    zxxup = random.uniform(-6, 0) 
#    zxxlo = random.uniform(0, 74) #0, 24
#    return nup.array([ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo])

def get_curve_coeffs(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo):
    Cup = nup.zeros((6, 6))
    Cup[0] = [1, 1, 1, 1, 1, 1]
    Cup[1] = [xup**0.5, xup**1.5, xup**2.5, xup**3.5, xup**4.5, xup**5.5]
    Cup[2] = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    Cup[3] = [0.5*(xup**(-0.5)), 1.5*(xup**0.5), 2.5*(xup**1.5), 3.5*(xup**2.5), 4.5*(xup**3.5), 5.5*(xup**4.5)]
    Cup[4] = [-0.25*(xup**(-1.5)), 0.75*(xup**(-0.5)), 3.75*(xup**0.5), 8.75*(xup**1.5), 15.75*(xup**2.5), 24.75*(xup**3.5)]
    Cup[5] = [1, 0, 0, 0, 0, 0]

    bup = nup.array(([zte + (dzte/2)], [zup], [nup.tan(nup.radians(ate - (bte/2)))], [0], [zxxup], [(2*rup)**0.5]))

    Clo = nup.zeros((6, 6))
    Clo[0] = [1, 1, 1, 1, 1, 1]
    Clo[1] = [xlo**0.5, xlo**1.5, xlo**2.5, xlo**3.5, xlo**4.5, xlo**5.5]
    Clo[2] = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    Clo[3] = [0.5*(xlo**(-0.5)), 1.5*(xlo**0.5), 2.5*(xlo**1.5), 3.5*(xlo**2.5), 4.5*(xlo**3.5), 5.5*(xlo**4.5)]
    Clo[4] = [-0.25*(xlo**(-1.5)), 0.75*(xlo**(-0.5)), 3.75*(xlo**0.5), 8.75*(xlo**1.5), 15.75*(xlo**2.5), 24.75*(xlo**3.5)]
    Clo[5] = [1, 0, 0, 0, 0, 0]

    blo = nup.array(([zte - (dzte/2)], [zlo], [nup.tan(nup.radians(ate + (bte/2)))], [0], [zxxlo], [-(2*rlo)**0.5]))
    try:
        aup = nup.squeeze(nup.transpose(nup.matmul(nup.linalg.inv(Cup), bup)))
        alo = nup.squeeze(nup.transpose(nup.matmul(nup.linalg.inv(Clo), blo)))
    except LinAlgError:
        return [1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]

    return aup, alo

def airfoilgen(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo, coord_file, np = 100, check_cross = False, even_spacing = False):
    try:
        aup, alo = get_curve_coeffs(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo)
        if check_cross:
            if curve_cross(aup, alo):
                return False

        with open(coord_file, 'w') as inp:
            if even_spacing:
                for x in nup.linspace(1, 0, np + 1):
                    inp.write(f'{x} {PARSEC_main_func(aup, x)}\n')
                for x in nup.linspace(0, 1, np + 1):
                    inp.write(f'{x} {PARSEC_main_func(alo, x)}\n')
            else:
                for i in range(1, np + 1):
                    theta = (180.0/(np - 1))*(np - i)
                    x = 0.5 - (0.5*nup.cos(theta*nup.pi/180.0))
                    inp.write(f'{x} {PARSEC_main_func(aup, x)}\n')
                for i in range(np, 0, -1):
                    theta = (180.0/(np - 1))*(np - i)
                    x = 0.5 - (0.5*nup.cos(theta*nup.pi/180.0))
                    inp.write(f'{x} {PARSEC_main_func(alo, x)}\n')

        return True
    except:
        return False

def get_curve_coeffs_params(PARSEC_params):
    ate = PARSEC_params[0]
    bte = PARSEC_params[1]
    zte = PARSEC_params[2]
    dzte = PARSEC_params[3]
    rup = PARSEC_params[4]
    xup = PARSEC_params[5]
    zup = PARSEC_params[6]
    zxxup = PARSEC_params[7]
    rlo = PARSEC_params[8]
    xlo = PARSEC_params[9]
    zlo = PARSEC_params[10]
    zxxlo = PARSEC_params[11]
    return get_curve_coeffs(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo)

def airfoilgen_params(PARSEC_params, coord_file, np = 100):
    ate = PARSEC_params[0]
    bte = PARSEC_params[1]
    zte = PARSEC_params[2]
    dzte = PARSEC_params[3]
    rup = PARSEC_params[4]
    xup = PARSEC_params[5]
    zup = PARSEC_params[6]
    zxxup = PARSEC_params[7]
    rlo = PARSEC_params[8]
    xlo = PARSEC_params[9]
    zlo = PARSEC_params[10]
    zxxlo = PARSEC_params[11]
    return airfoilgen(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo, coord_file = coord_file, np = np)

#def check_wavy(aup, alo):
#    up_lim = 2
#    lo_lim = 3
#    upc1 = False
#    upc2 = False
#    loc1 = False
#    loc2 = False
#    loc3 = False
#    zu_remember = 0
#    zlo_remember = 0
#    for x in nup.linspace(0, 1, 101):
#        zu = PARSEC_main_func(aup, x)
#        zlo = PARSEC_main_func(alo, x)
#
#        if (zu < zu_remember) and (not upc1):
#            upc1 = True
#        elif (zu > zu_remember) and (not upc2):
#            upc2 = True
#        elif (zu < zu_remember) and upc2:
#            return True
#
#        if (zlo > zlo_remember) and (not loc1):
#            loc1 = True
#        elif (zlo < zlo_remember) and (not loc2):
#            loc2 = True
#        elif (zlo > zlo_remember) and (not loc3):
#            loc3 = True
#        elif (zlo < zlo_remember) and loc3:
#            return True
#
#        zlo_remember = zlo
#        zu_remember = zu
#    return False

def check_wavy(aup, alo, tolerance = 0.0001):
    up_lim = 2
    lo_lim = 3
    u_counter = 0
    lo_counter = 0
    for x in nup.linspace(0, 1, 501):
        zprimeu = PARSEC_derivative1(aup, x)
        zprimelo = PARSEC_derivative1(alo, x)
        if abs(zprimeu) < tolerance:
            u_counter += 1
        if abs(zprimelo) < tolerance:
            lo_counter += 1
    
    print(u_counter)
    print(lo_counter)

    if lo_counter > lo_lim:
        return True
    if u_counter > up_lim:
        return True

    return False

def check_valid(aup, alo, params):
    if nup.all(aup == 1) or nup.all(alo == 1):
        return False
    ate = params[0]
    bte = params[1]
    zte = params[2]
    dzte = params[3]
    rup = params[4]
    xup = params[5]
    zup = params[6]
    zxxup = params[7]
    rlo = params[8]
    xlo = params[9]
    zlo = params[10]
    zxxlo = params[11]

    if curve_cross(aup, alo):
        #print("Discarded cross")
        return False
    if PARSEC_derivative2(aup, xup) > 0:
        #print("Discarded zxxup > 0")
        return False
    if PARSEC_derivative2(alo, xlo) < 0:
        #print("Discarded zxxlo < 0")
        return False
    for i in nup.linspace(0, 1, 26):
        if PARSEC_main_func(aup, i) > PARSEC_main_func(aup, xup):
            #print("Discarded z > zup")
            return False
        if PARSEC_main_func(alo, i) < PARSEC_main_func(alo, xlo):
            #print("Discarded z < zlo")
            return False
    #if check_wavy(aup, alo):
    #    print("Discarded wavy")
    #    return False
    return True








#utilize existing airfoils

def read_foil(airfoil_file, header_lines):
    x, y = nup.loadtxt(airfoil_file, skiprows = header_lines, unpack = True)
    return x, y

def count_header_lines(file):
    with open(file, 'r') as inp:
        i = 0
        while True:
            try:
                content = inp.readline()
                cont = content.split() #.replace("\n", '')
                _ = float(cont[0])
                break
            except:
                i += 1
                print(i)
    return i

def split_upper_lower(x, y):
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
    i = 0
    x_mat = nup.ones((x_curve.size, 6))
    while i < len(x_curve):
        ii = 0
        while ii < 6:
            x_mat[i][ii] = x_curve[i]**(ii + 0.5)
            ii += 1
        i += 1
    return x_mat 

def write_PARSEC_zcoor(ans, x_curve = 0, np = 50):
    #i = 0
    y_poli = nup.array([])
    x_dist = nup.array([])
    #while i < len(x_curve):
    #    y_poli = nup.append(y_poli, PARSEC_main_func(ans, x_curve[i]))
    #    i += 1
    
    for i in range(np, 0, -1):
        theta = (180.0/(np - 1))*(np - i)
        x = 0.5 - (0.5*nup.cos(theta*nup.pi/180.0))
        x_dist = nup.append(x_dist, x)
        y_poli = nup.append(y_poli,PARSEC_main_func(ans, x))
    return y_poli, x_dist

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
    return answer_u, answer_lo

def get_PARSEC_parameters(answer_u, answer_lo):
    rlo = (answer_lo[0]**2)/2
    rup = (answer_u[0]**2)/2
    zte = (nup.sum(answer_lo) + nup.sum(answer_u))/2 
    dzte = nup.sum(answer_u) - nup.sum(answer_lo)
    #print(11)
    xlo = newt_rhap_PARSEC_first(answer_lo, tolerance = 0.00000001, init = 0.0001)
    #print(22)
    xup = newt_rhap_PARSEC_first(answer_u, tolerance = 0.00000001, init = 0.0001)
    #print(33)
    zlo = PARSEC_main_func(answer_lo, xlo)
    zup = PARSEC_main_func(answer_u, xup)
    zxxlo = PARSEC_derivative2(answer_lo, xlo)
    zxxup = PARSEC_derivative2(answer_u, xup)
    ate = (nup.degrees(nup.arctan(PARSEC_derivative1(answer_u, 1))) + nup.degrees(nup.arctan(PARSEC_derivative1(answer_lo, 1))))/2
    bte = nup.degrees(nup.arctan(PARSEC_derivative1(answer_lo, 1))) - nup.degrees(nup.arctan(PARSEC_derivative1(answer_u, 1)))

    return [ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo]












#VARGET WRITTEN BELOW-----------------------------------
#airfoil_file = r'C:\Users\PEDRO\Desktop\IC DE HÉLICE\git_control\git_prop\airfoils\airfoil.txt'
#airfoil_file = r'C:\Users\PEDRO\Desktop\IC DE HÉLICE\git_control\git_prop\airfoils\query_uiuc\ua79sff-il.dat'
#airfoil_file = r'C:\Users\PEDRO\Desktop\IC DE HÉLICE\git_control\git_prop\airfoils\query_airfoiltools\b707b-il.dat'
#x, y = read_foil(airfoil_file, header_lines = count_header_lines(airfoil_file))
##print(1)
#x_lo, y_lo, x_u, y_u = split_upper_lower(x, y)
##print(2)
#answer_lo, answer_u = get_PARSEC_curve_coefficients(x_lo, x_u, y_lo, y_u)
##print(3)
#y_poli_lo = write_PARSEC_zcoor(x_lo, answer_lo) 
##print(4)
#y_poli_u = write_PARSEC_zcoor(x_u, answer_u) 
#PARSEC_params = get_PARSEC_parameters(answer_lo, answer_u)
##print(5)
#print(PARSEC_params)
#fig, ax = plt.subplots()
#ax.plot(x_u, y_u, 'r', label = "original")
#ax.plot(x_lo, y_lo, 'r') 
#ax.plot(x_lo, y_poli_lo, 'b--', label = "PARSEC")
#ax.plot(x_u, y_poli_u, 'b--') 
#ax.legend()
#ax.xaxis.grid(True, which='major')
#plt.show()


#list_at = os.listdir(r"C:\Users\pedro\Desktop\Iniciação_Científica\git_prop\airfoils\query_airfoiltools")
#list_uiuc = os.listdir(r"C:\Users\pedro\Desktop\Iniciação_Científica\git_prop\airfoils\query_UIUC")
#
#for airfoil in list_at:
#    airfoil_file = r'C:\Users\pedro\Desktop\Iniciação_Científica\git_prop\airfoils\query_airfoiltools\\' + airfoil
#    print(airfoil)
#    x, y = read_foil(airfoil_file, header_lines = count_header_lines(airfoil_file))
#    x_lo, y_lo, x_u, y_u = split_upper_lower(x, y)
#    answer_lo, answer_u = get_PARSEC_curve_coefficients(x_lo, x_u, y_lo, y_u)
#    y_poli_lo = write_PARSEC_zcoor(x_lo, answer_lo) 
#    y_poli_u = write_PARSEC_zcoor(x_u, answer_u) 
#    PARSEC_params = get_PARSEC_parameters(answer_lo, answer_u)
#    #print(PARSEC_params)
#    fig, ax = plt.subplots()
#    ax.plot(x_u, y_u, 'r', label = "original")
#    ax.plot(x_lo, y_lo, 'r') 
#    ax.plot(x_lo, y_poli_lo, 'b--', label = "PARSEC")
#    ax.plot(x_u, y_poli_u, 'b--') 
#    ax.legend()
#    ax.xaxis.grid(True, which='major')
#    plt.show()
#for airfoil in list_uiuc:
#    airfoil_file = r'C:\Users\pedro\Desktop\Iniciação_Científica\git_prop\airfoils\query_UIUC\\' + airfoil
#    x, y = read_foil(airfoil_file, header_lines = count_header_lines(airfoil_file))
#    x_lo, y_lo, x_u, y_u = split_upper_lower(x, y)
#    answer_lo, answer_u = get_PARSEC_curve_coefficients(x_lo, x_u, y_lo, y_u)
#    y_poli_lo = write_PARSEC_zcoor(x_lo, answer_lo) 
#    y_poli_u = write_PARSEC_zcoor(x_u, answer_u) 
#    PARSEC_params = get_PARSEC_parameters(answer_lo, answer_u)
#    print(PARSEC_params)
#    fig, ax = plt.subplots()
#    ax.plot(x_u, y_u, 'r', label = "original")
#    ax.plot(x_lo, y_lo, 'r') 
#    ax.plot(x_lo, y_poli_lo, 'b--', label = "PARSEC")
#    ax.plot(x_u, y_poli_u, 'b--') 
#    ax.legend()
#    ax.xaxis.grid(True, which='major')
#    plt.show()

#RANDVARGEN WRITTEN BELOW---------------------------------
#coord_file = r'C:\Users\PEDRO\Desktop\IC DE HÉLICE\git_control\git_prop\airfoils\airfoil.txt'
#while True:
#    params = PARSEC_randvargen()
#    run_check = airfoilgen_params(params, coord_file = coord_file, np = 50) 
#    aup, alo = get_curve_coeffs_params(params)
#    #plot_file(coord_file)
#    if check_valid(aup, alo, params):
#        break
#print('[ate, bte, zte, dzte, \nrup, xup, zup, zxxup, \nrlo, xlo, zlo, zxxlo]')
#print(params[:4])
#print(params[4:8])
#print(params[8:])
#plot_file(coord_file)

#run_check = airfoilgen_params(params, coord_file = coord_file, np = 50) 
#plot_file(coord_file)
#aup, alo = get_curve_coeffs_params(params)
#valid_check = check_valid(aup, alo, params)


#coord_file = r'C:\Users\PEDRO\Desktop\IC DE HÉLICE\git_control\git_prop\airfoils\airfoil.txt'
#while True:
#    params = PARSEC_randvargen()
#    run_check = airfoilgen_params(params, coord_file = coord_file, np = 50) 
#    aup, alo = get_curve_coeffs_params(params)
#    x = check_wavy(aup, alo, tolerance = 0.01)
#    plot_file(coord_file)