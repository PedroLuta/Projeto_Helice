import random
import numpy as nup
import matplotlib.pyplot as plt

from numpy.linalg.linalg import LinAlgError

class RegFoil:
    def __init__(self, name, xcor, ycor, polars = [], Re_list = []):
        self.name = name
        self.xcor = xcor[:]
        self.ycor = ycor[:]
        self.polars = polars[:]
        self.Re_list = Re_list[:]

class PFoil:
    def __init__(self, random = False, selig_file = "", lednicer_file = "", plot_diff = False):
        if selig_file != "":
            self.load_from_selig(selig_file, plot_diff)
            self.Cl0 = 0
            self.Cd0 = 0
        elif lednicer_file != "":
            self.Cl0 = 0
            self.Cd0 = 0
        elif random:
            temp = rand_params_gen()
            self.set_parameters_vec(temp)
            #self.name = ""
            #self.ate = temp[0]
            #self.bte = temp[1]
            #self.zte = temp[2]
            #self.dzte = temp[3]
            #self.rup = temp[4]
            #self.xup = temp[5]
            #self.zup = temp[6]
            #self.zxxup = temp[7]
            #self.rlo = temp[8]
            #self.xlo = temp[9]
            #self.zlo = temp[10]
            #self.zxxlo = temp[11]
            #self.params = [self.ate, self.bte, self.zte,\
            #     self.dzte, self.rup, self.xup, self.zup, self.zxxup, self.rlo, self.xlo, \
            #         self.zlo, self.zxxlo]
            #self.params = temp[:]
            #self.aup, self.alo = coeffs_from_params_vec(self.params)
            self.Cl0 = 0
            self.Cd0 = 0
        else:
            self.name = ""
            self.set_parameters_vec([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], coeffs = False)
            self.aup = []
            self.alo = []
            self.Cl0 = 0
            self.Cd0 = 0
    
    def set_name(self, name):
        self.name = name

    def set_parameters_vec(self, parvec, coeffs = True):
        self.ate = parvec[0]
        self.bte = parvec[1]
        self.zte = parvec[2]
        self.dzte = parvec[3]
        self.rup = parvec[4]
        self.xup = parvec[5]
        self.zup = parvec[6]
        self.zxxup = parvec[7]
        self.rlo = parvec[8]
        self.xlo = parvec[9]
        self.zlo = parvec[10]
        self.zxxlo = parvec[11]
        self.params = parvec[:]
        if coeffs:
            self.aup, self.alo = coeffs_from_params_vec(parvec)

    def set_parameters(self, ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo):
        self.set_parameters_vec([ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo])

    def set_coeffs(self, aup, alo):
        self.aup = aup
        self.alo = alo

    def get_params_vec(self):
        return [self.ate, self.bte, self.zte, self.dzte, self.rup, self.xup, self.zup, self.zxxup, \
            self.rlo, self.xlo, self.zlo, self.zxxlo]
    
    def get_coeffs(self):
        return self.aup, self.alo

    def load_from_selig(self, file, plot):
        x, y = read_foil(file, header_lines = count_header_lines(file))
        x_lo, y_lo, x_u, y_u = split_upper_lower(x, y)
        self.aup, self.alo = coeffs_from_curves(x_lo, x_u, y_lo, y_u)
        parvec = parameters_from_coeffs(self.aup, self.alo)
        self.set_parameters_vec(parvec)
        #self.ate = parvec[0]
        #self.bte = parvec[1]
        #self.zte = parvec[2]
        #self.dzte = parvec[3]
        #self.rup = parvec[4]
        #self.xup = parvec[5]
        #self.zup = parvec[6]
        #self.zxxup = parvec[7]
        #self.rlo = parvec[8]
        #self.xlo = parvec[9]
        #self.zlo = parvec[10]
        #self.zxxlo = parvec[11]
        #self.params = parvec[:]
        if plot:
            y_poli_up, x_dist_up = write_zcoor(self.aup)
            y_poli_lo, x_dist_lo = write_zcoor(self.alo)
            _, ax = plt.subplots()
            ax.plot(x_u, y_u, 'r')
            ax.plot(x_lo, y_lo, 'r', label = "original")
            ax.plot(x_dist_up, y_poli_up, 'b--')
            ax.plot(x_dist_lo, y_poli_lo, 'b--', label = "PARSEC aproximation")
            ax.xaxis.grid(True, which='major')
            ax.legend()
            plt.show()
            



    def plot(self):
        y_poli_up, x_dist_up = write_zcoor(self.aup)
        y_poli_lo, x_dist_lo = write_zcoor(self.alo)
        fig, ax = plt.subplots()
        ax.plot(x_dist_up, y_poli_up, 'r')
        ax.plot(x_dist_lo, y_poli_lo, 'r')
        ax.xaxis.grid(True, which='major')
        plt.show()

    def plot_from_params(self):
        plot_from_params_vec([self.ate, self.bte, self.zte, self.dzte, self.rup, self.xup, self.zup, self.zxxup, \
            self.rlo, self.xlo, self.zlo, self.zxxlo])




















#general functions
def main_func(ans, x):
    return ans[0]*(x**0.5) + ans[1]*(x**1.5) + ans[2]*(x**2.5) + ans[3]*(x**3.5) + ans[4]*(x**4.5) + ans[5]*(x**5.5)

def derivative1(ans, x):
    return 0.5*ans[0]*(x**(-0.5)) + 1.5*ans[1]*(x**0.5) + 2.5*ans[2]*(x**1.5) + 3.5*ans[3]*(x**2.5) + 4.5*ans[4]*(x**3.5) + 5.5*ans[5]*(x**4.5)

def derivative2(ans, x):
    return -0.25*ans[0]*(x**(-1.5)) + 0.75*ans[1]*(x**(-0.5)) + 3.75*ans[2]*(x**0.5) + 8.75*ans[3]*(x**1.5) + 15.75*ans[4]*(x**2.5) + 24.75*ans[5]*(x**3.5)















#find_out methods

def newt_rhap_first(ans, tolerance = 0.0001, init = 0.01, cutbreak = False):
    xi = init
    retry = False
    while True:
        if xi < 0:
            retry = True
            break
        delta = derivative1(ans, xi)/derivative2(ans, xi)
        xn = xi - delta
        if abs(xn - xi) < tolerance:
            break
        xi = xn
    if retry:
        if cutbreak:
            return -1337
        xn = newt_rhap_first(ans, tolerance = tolerance, init = 0.3, cutbreak = True)
    return xn










#validity checkers

def check_wavy(aup, alo, up_lim = 2, lo_lim = 3):
    zu_remember1 = 0
    zlo_remember1 = 0
    zu_remember2 = 0
    zlo_remember2 = 0
    lo_count = 0
    up_count = 0
    first = True
    for x in nup.linspace(0, 1, 101, endpoint = False):
        zu = main_func(aup, x)
        zlo = main_func(alo, x)
        if x == 0:
            continue
        if first:
            zu_remember1 = zu_remember2
            zlo_remember1 = zlo_remember2
            zlo_remember2 = zlo
            zu_remember2 = zu
            first = False
            continue

        if (nup.sign(zlo_remember2 - zlo_remember1) != nup.sign(zlo - zlo_remember2)):
            lo_count += 1

        if (nup.sign(zu_remember2 - zu_remember1) != nup.sign(zu - zu_remember2)):
            up_count += 1

        zu_remember1 = zu_remember2
        zlo_remember1 = zlo_remember2
        zlo_remember2 = zlo
        zu_remember2 = zu

    if up_count > up_lim:
        return True
    if lo_count > lo_lim:
        return True
    return False

def check_wavy_der(aup, alo, up_lim = 2, lo_lim = 3):
    lo_count = 0
    up_count = 0
    x_remember = 0
    first = True
    for x in nup.linspace(0, 1, 101, endpoint = False):
        if x == 0:
            continue
        if first:
            x_remember = x
            first = False
            continue
        if (nup.sign(derivative1(alo, x_remember)) != nup.sign(derivative1(alo, x))):
            lo_count += 1
        if (nup.sign(derivative1(aup, x_remember)) != nup.sign(derivative1(aup, x))):
            up_count += 1

        x_remember = x

    if up_count > up_lim:
        return True
    if lo_count > lo_lim:
        return True
    return False

def curve_cross(aup, alo):
    for x in nup.linspace(0, 1, 101):
        zu = main_func(aup, x)
        zlo = main_func(alo, x)
        if zlo > zu:
            return True
    return False

def check_choke(aup, alo):
    delta_remember = 0
    changed = False
    for x in nup.linspace(0, 1, 101):
        zu = main_func(aup, x)
        zlo = main_func(alo, x)
        delta = zu - zlo
        if (delta < delta_remember) and not changed:
            changed = True
        elif changed:
            if delta > delta_remember:
                return True
        delta_remember = delta
    return False

def check_max(aup, alo, xup, xlo):
    for i in nup.linspace(0, 1, 26):
        if main_func(aup, i) > main_func(aup, xup):
            return True
        if main_func(alo, i) < main_func(alo, xlo):
            return True
    return False

def check_valid(aup, alo, params, checkWavy = False):
    xup = params[5]
    xlo = params[9]

    if nup.all(aup == 1) or nup.all(alo == 1):
        return False
    if curve_cross(aup, alo):
        return False
    if check_choke(aup, alo):
        return False
    if check_max(aup, alo, xup, xlo):
        return False
    if derivative2(aup, xup) > 0:
        return False
    if derivative2(alo, xlo) < 0:
        return False
    if checkWavy:
        if check_wavy(aup, alo):
            return False
    return True







#auxiliary functions

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

def get_x_mat(x_curve):
    i = 0
    x_mat = nup.ones((x_curve.size, 6))
    while i < len(x_curve):
        ii = 0
        while ii < 6:
            x_mat[i][ii] = x_curve[i]**(ii + 0.5)
            ii += 1
        i += 1
    return x_mat 

def write_zcoor(ans, np = 50):
    y_poli = nup.array([])
    x_dist = nup.array([])
    
    for i in range(np, 0, -1):
        theta = (180.0/(np - 1))*(np - i)
        x = 0.5 - (0.5*nup.cos(theta*nup.pi/180.0))
        x_dist = nup.append(x_dist, x)
        y_poli = nup.append(y_poli, main_func(ans, x))
    return y_poli, x_dist








#generate foil functions

def rand_params_gen():
    ate = random.uniform(-30, 15) #-30, 40, -30, 11
    bte = random.uniform(0, 45) #0, 140, 0, 45
    zte = random.uniform(-0.05, 0.05) #-0.07, 0.09, -0.02, 0.04
    dzte = random.uniform(0, 0.05) #0, 0.26, 0, 0.05
    rup = random.uniform(0, 0.1) 
    xup = random.uniform(0, 0.5) #0.07, 0.59, 0.11, 0.55
    zup = random.uniform(0.0001, 0.35) #0.02, 0.35, 0.01, 0.35
    zxxup = random.uniform(-5, 0) #-6, 0
    rlo = random.uniform(0, 0.1)
    xlo = random.uniform(0, 0.3) #0, 1, 0.008, 0.6
    zlo = random.uniform(-0.35, -0.0001) #-0.32, 0.08, -0.33, 0.08
    zxxlo = random.uniform(0, 10) #0, 24, 0, 74
    return nup.array([ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo])

def limits_from_foil(foil_params, factor):
    limits = nup.zeros((12, 2))
    i = 0
    while i < len(foil_params):
        limits[i][0] = foil_params[i] - (foil_params[i]*factor)
        limits[i][1] = foil_params[i] + (foil_params[i]*factor)
        i += 1
    return limits

def flex_params_gen(limits):
    ate = random.uniform(limits[0][0], limits[0][1]) 
    bte = random.uniform(limits[1][0], limits[1][1]) 
    zte = random.uniform(limits[2][0], limits[2][1]) 
    dzte = random.uniform(limits[3][0], limits[3][1]) 
    rup = random.uniform(limits[4][0], limits[4][1]) 
    xup = random.uniform(limits[5][0], limits[5][1]) 
    zup = random.uniform(limits[6][0], limits[6][1]) 
    zxxup = random.uniform(limits[7][0], limits[7][1]) 
    rlo = random.uniform(limits[8][0], limits[8][1]) 
    xlo = random.uniform(limits[9][0], limits[9][1]) 
    zlo = random.uniform(limits[10][0], limits[10][1]) 
    zxxlo = random.uniform(limits[11][0], limits[11][1])
    return nup.array([ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo])

def coeffs_from_params(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo):
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

def write_from_params(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo, coord_file, np = 100, even_spacing = False):
    try:
        aup, alo = coeffs_from_params(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo)

        with open(coord_file, 'w') as inp:
            if even_spacing:
                for x in nup.linspace(1, 0, np + 1):
                    inp.write(f'{x} {main_func(aup, x)}\n')
                for x in nup.linspace(0, 1, np + 1):
                    inp.write(f'{x} {main_func(alo, x)}\n')
            else:
                for i in range(1, np + 1):
                    theta = (180.0/(np - 1))*(np - i)
                    x = 0.5 - (0.5*nup.cos(theta*nup.pi/180.0))
                    inp.write(f'{x} {main_func(aup, x)}\n')
                for i in range(np, 0, -1):
                    theta = (180.0/(np - 1))*(np - i)
                    x = 0.5 - (0.5*nup.cos(theta*nup.pi/180.0))
                    inp.write(f'{x} {main_func(alo, x)}\n')

        return True
    except:
        return False

def write_from_coeffs(aup, alo, coord_file, np = 100, even_spacing = False):
    try:
        with open(coord_file, 'w') as inp:
            if even_spacing:
                for x in nup.linspace(1, 0, np + 1):
                    inp.write(f'{x} {main_func(aup, x)}\n')
                for x in nup.linspace(0, 1, np + 1):
                    inp.write(f'{x} {main_func(alo, x)}\n')
            else:
                for i in range(1, np + 1):
                    theta = (180.0/(np - 1))*(np - i)
                    x = 0.5 - (0.5*nup.cos(theta*nup.pi/180.0))
                    inp.write(f'{x} {main_func(aup, x)}\n')
                for i in range(np, 0, -1):
                    theta = (180.0/(np - 1))*(np - i)
                    x = 0.5 - (0.5*nup.cos(theta*nup.pi/180.0))
                    inp.write(f'{x} {main_func(alo, x)}\n')

        return True
    except:
        return False

def coeffs_from_curves(x_lo, x_u, y_lo, y_u):
    x_mat_lo = get_x_mat(x_lo)
    x_mat_u = get_x_mat(x_u)
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

def parameters_from_coeffs(answer_u, answer_lo):
    rlo = (answer_lo[0]**2)/2
    rup = (answer_u[0]**2)/2
    zte = (nup.sum(answer_lo) + nup.sum(answer_u))/2 
    dzte = nup.sum(answer_u) - nup.sum(answer_lo)
    xlo = newt_rhap_first(answer_lo, tolerance = 0.00000001, init = 0.0001)
    xup = newt_rhap_first(answer_u, tolerance = 0.00000001, init = 0.0001)
    zlo = main_func(answer_lo, xlo)
    zup = main_func(answer_u, xup)
    zxxlo = derivative2(answer_lo, xlo)
    zxxup = derivative2(answer_u, xup)
    ate = (nup.degrees(nup.arctan(derivative1(answer_u, 1))) + nup.degrees(nup.arctan(derivative1(answer_lo, 1))))/2
    bte = nup.degrees(nup.arctan(derivative1(answer_lo, 1))) - nup.degrees(nup.arctan(derivative1(answer_u, 1)))

    return [ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo]









#PLOTTING
def plot_from_split(x_u, x_lo, y_u, y_lo):
    fig, ax = plt.subplots()
    ax.plot(x_u, y_u, 'r')
    ax.plot(x_lo, y_lo, 'r')
    ax.xaxis.grid(True, which='major')
    plt.show()

def plot_from_file(file):
    x, y = read_foil(file, header_lines = count_header_lines(file))
    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.xaxis.grid(True, which='major')
    plt.show()

def plot_from_xy(x, y):
    fig, ax = plt.subplots()
    ax.plot(x, y, 'r')
    ax.xaxis.grid(True, which='major')
    plt.show()

def plot_from_coeffs(aup, alo):
    y_poli_u, x_u = write_zcoor(aup)
    y_poli_lo, x_lo = write_zcoor(alo)
    plot_from_split(x_u, x_lo, y_poli_u, y_poli_lo)

def plot_from_params_vec(params_vec):
    aup, alo = coeffs_from_params_vec(params_vec)
    plot_from_coeffs(aup, alo)








#comodity functions
def coeffs_from_params_vec(PARSEC_params):
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
    return coeffs_from_params(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo)

def write_from_params_vec(PARSEC_params, coord_file, np = 100):
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
    return write_from_params(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo, coord_file = coord_file, np = np)



#myFoil = PFoil(selig_file = r'C:\Users\PEDRO\Desktop\Airfoils\s1223.txt', plot_diff = True)







