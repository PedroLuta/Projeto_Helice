import optimization_NSGA2 as optim
import xfoil_interface
import induction
from importing import *

class Chord_Dist:
    def __init__(self, xi, Croot, chrom = [], coeffs = []):
        self.xi = xi
        self.Croot = Croot
        self.params = chrom
        self.coeffs = coeffs
    def set_params(self, params):
        self.params = params
    def get_params(self):
        return self.params
    def set_coeffs(self, coeffs):
        self.coeffs = coeffs
    def get_coeffs(self):
        return self.coeffs
    def plot_show(self):
        dist = []
        x_dist = nup.linspace(self.xi, 1, 51)
        for x in x_dist:
            dist.append(chord_dist_main(self.coeffs, x))
        plt.plot(x_dist, dist)
        plt.show()
def chord_dist_main(ans, x):
    return ans[0]*x**5 + ans[1]*x**4 + ans[2]*x**3 + ans[3]*x**2 + ans[4]*x + ans[5] 
def chord_dist_first(ans, x):
    return 5*ans[0]*x**4 + 4*ans[1]*x**3 + 3*ans[2]*x**2 + 2*ans[3]*x + ans[4]
def chord_dist_second(ans, x):
    return 20*ans[0]*x**3 + 12*ans[1]*x**2 + 6*ans[2]*x + 2*ans[3]
def true_from_params(params):
    xm = Auxiliary.linear_interpolate(0, 1, xi, 1, params[0]) 
    Cmax = params[1]
    xp = Auxiliary.linear_interpolate(0, 1, xm, 1, params[2])
    Ctip = Auxiliary.linear_interpolate(0, 1, 0, Cmax, params[4])
    Cplat = Auxiliary.linear_interpolate(0, 1, Ctip, Cmax, params[3])
    return xm, Cmax, xp, Cplat, Ctip
def true_from_dist(dist):
    xi = dist.xi
    Croot = dist.Croot
    chrom = dist.get_params()
    xm = Auxiliary.linear_interpolate(0, 1, xi, 1, chrom[0]) 
    Cmax = chrom[1]
    xp = Auxiliary.linear_interpolate(0, 1, xm, 1, chrom[2])
    Ctip = Auxiliary.linear_interpolate(0, 1, 0, Cmax, chrom[4])
    Cplat = Auxiliary.linear_interpolate(0, 1, Ctip, Cmax, chrom[3])
    return xi, Croot, xm, Cmax, xp, Cplat, Ctip
def coeffs_from_params(xi, Croot, params):
    xm, Cmax, xp, Cplat, Ctip = true_from_params(params)
    A = nup.zeros((6, 6))
    A[0] = [xi**5, xi**4, xi**3, xi**2, xi, 1]
    A[1] = [xm**5, xm**4, xm**3, xm**2, xm, 1]
    A[2] = [xp**5, xp**4, xp**3, xp**2, xp, 1]
    A[3] = [1, 1, 1, 1, 1, 1]
    A[4] = [5*xm**4, 4*xm**3, 3*xm**2, 2*xm, 1, 0]
    A[5] = [20*xp**3, 12*xp**2, 6*xp, 2, 0, 0]
    b = nup.array([Croot, Cmax, Cplat, Ctip, 0, 0])
    try:
        coeffs = nup.squeeze(nup.transpose(nup.matmul(nup.linalg.inv(A), b)))
    except LinAlgError:
        return [1, 1, 1, 1, 1, 1]
    return coeffs
def check_valley(coeffs, xi):
    changed = False
    z_remember = 0
    for x in nup.linspace(xi, 1, 100):
        z_new = chord_dist_main(coeffs, x)
        if z_new < z_remember:
            changed = True
        if changed:
            if z_new > z_remember:
                return True
        z_remember = z_new
    return False
def check_max_chord(dist):
    coeffs = dist.coeffs
    params = dist.params
    xi, Croot, xm, Cmax, xp, Cplat, Ctip = true_from_dist(dist)
    for x in nup.linspace(xi, 1, 100):
        if chord_dist_main(coeffs, x) > Cmax:
            return True
    return False
def check_valid(dist):
    coeffs = dist.coeffs
    xi, Croot, xm, Cmax, xp, Cplat, Ctip = true_from_dist(dist)
    if chord_dist_second(coeffs, xm) >= 0:
        return False
    if chord_dist_first(coeffs, xi) < 0:
        return False
    if chord_dist_first(coeffs, xp) > 0:
        return False
    if chord_dist_first(coeffs, 1) >= 0:
        return False
    if check_valley(coeffs, xi):
        return False
    if check_max_chord(dist):
        return False
    return True
def calculate_T_Q(vi, rpm, Blades, R, dist, sections = 21, airfoil = 'airfoils\\sunnysky.dat', rho = 1.225, dvisc = 1.8/100000, alphas = [0, 10, 1]):
    xi = dist.xi
    coeffs = dist.get_coeffs()
    #print(coeffs)
    a1, a2, astep = alphas[0], alphas[1], alphas[2]
    kvisc = dvisc/rho
    Beta_vector = []
    r_vector = []
    dT_vector = []
    dQ_vector = []
    radps = (rpm*2*pi)/60
    for rr in nup.linspace(xi*R, R, num = sections):
        chord = R*chord_dist_main(coeffs, rr/R)
        #print(chord)
        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc) 
        alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_comMethod(Re, a1, a2, astep, afile = airfoil)
        #alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_fileMethod(Re, a1 = a1, a2 = a2, astep = astep, afile = airfoil)
        alpha, Cl, Cd = xfoil_interface.calculate_most_eff_alpha_simple(alpha_c, Cl_c, Cd_c)
        try:
            dT, dQ, phi, _, _ = induction.jitted_induction_Ftip(radps, rr, Cl, Cd, Blades, rho, R, chord, vi)
        except:
            dT, dQ, phi, _, _ = 0, 0, 0, 0, 0 
        Beta_vector.append(alpha + math.degrees(phi))
        r_vector.append(rr)
        dQ_vector.append(dQ)
        dT_vector.append(dT)
    #print(dT_vector)
    #print(dQ_vector)
    Thrust = Auxiliary.area_under_curve(r_vector, dT_vector)
    Torque = Auxiliary.area_under_curve(r_vector, dQ_vector)
    #print(Thrust, Torque)
    return Thrust, Torque

def create1stgen(n, R = 0.3, xi = 0, Croot = 0.05, limits = [[0, 1], [0, -1], [0, 1], [0, 1], [0, 1]]):
    if limits[1][1] == -1:
        limits[1][1] = R/2
    population = []
    while len(population) < n:
        chromossome = optim.generate_chrom(limits)
        temp = optim.Individual(chromossome)
        tester = Chord_Dist(xi, Croot, chromossome, coeffs = coeffs_from_params(xi, Croot, chromossome))
        if check_valid(tester):
            population.append(temp)
    return population
def evaluate_population(population, vi, rpm, Blades, R, sections = 21):
    for dist in population:
        temp = Chord_Dist(xi, Croot, chrom = dist.get_chrom(), coeffs = coeffs_from_params(xi, Croot, dist.get_chrom()))
        if check_valid(temp):
            T, Q = calculate_T_Q(vi, rpm, Blades, R, temp, sections = sections)
        else:
            T, Q = 0, 1000
        dist.set_ObjVal([-T, Q])
    #population = optim.assign_fronts(population)
    #population = optim.assign_crowding(population)

if __name__ == '__main__':    
    start = time.time()
    vi = 10
    rpm = 2300
    Blades = 4
    R = 0.3
    sections = 5
    limits = [[0, 1], [0, R/2], [0, 1], [0, 1], [0, 1]]

    xi = 0.1
    Croot = 0.03
    #[xm, Cmax, xp, Cplat, Ctip] 

    n_ind = 30
    n_gen = 50

    n_objvals = 2

    generation = 1
    while generation <= n_gen:
        print(f'Geração: {generation}')
        if generation == 1:
            current_pop = create1stgen(n_ind, xi = xi, Croot = Croot, R = R, limits = limits)
            evaluate_population(current_pop, vi, rpm, Blades, R, sections = sections)
            for dist in current_pop:
                objvals = dist.get_ObjVal()
                #print(objvals)
                Thrust, Torque = objvals[0], objvals[1]
                plt.scatter(Torque, Thrust, 4, "r") 
        new_gen = optim.create_nextGen(current_pop, limits, n_objvals)
        evaluate_population(new_gen, vi, rpm, Blades, R, sections = sections)
        current_pop = optim.reinsert(current_pop, new_gen)
        for dist in current_pop:
            objvals = dist.get_ObjVal()
            Thrust, Torque = objvals[0], objvals[1]
            if Torque == 1000:
                continue
            plt.scatter(Torque, Thrust, 4, "r")
        generation += 1
    end = time.time()
    print(end - start)
    #plt.show()
#
    #for dist in current_pop:
    #    temp = Chord_Dist(xi, Croot, chrom = dist.get_chrom(), coeffs = coeffs_from_params(xi, Croot, dist.get_chrom()))
    #    temp.plot_show()