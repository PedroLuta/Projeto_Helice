import optimization_NSGA2 as optim
import xfoil_interface
import induction
from importing import *

def chord_dist_main(ans, x):
    return ans[0]*x**5 + ans[1]*x**4 + ans[2]*x**3 + ans[3]*x**2 + ans[4]*x

def chord_dist_first(ans, x):
    return 5*ans[0]*x**4 + 4*ans[1]*x**3 + 3*ans[2]*x**2 + 2*ans[3]*x + ans[4]

def chord_dist_second(ans, x):
    return 20*ans[0]*x**3 + 12*ans[1]*x**2 + 6*ans[2]*x + 2*ans[3]



def coeffs_from_params(params):
    xmax, Zmax, Zxxmax, Ztip, Zxtip = params[0], params[1], params[2], params[3], params[4] 
    A = nup.zeros((5, 5))
    A[0] = [xmax**5, xmax**4, xmax**3, xmax**2, xmax]
    A[1] = [5*xmax**4, 4*xmax**3, 3*xmax**2, 2*xmax, 1]
    A[2] = [20*xmax**3, 12*xmax**2, 6*xmax, 2, 0]
    A[3] = [1, 1, 1, 1, 1]
    A[4] = [5, 4, 3, 2, 1]
    
    b = nup.array([Zmax, 0, Zxxmax, Ztip, Zxtip])
    try:
        coeffs = nup.squeeze(nup.transpose(nup.matmul(nup.linalg.inv(A), b)))
    except LinAlgError:
        return [1, 1, 1, 1, 1, 1]
    return coeffs

#def check_valley(coeffs):
#    changed = False
#    z_remember = 0
#    for x in nup.linspace(0, 1, 100):
#        z_new = chord_dist_main(coeffs, x)
#        if z_new < z_remember:
#            changed = True
#        if changed:
#            if z_new > z_remember:
#                return True
#        z_remember = z_new
#    return False

#def check_max_chord(coeffs, xm):
#    max_c = chord_dist_main(coeffs, xm)
#    for x in nup.linspace(xi, 1, 100):
#        if chord_dist_main(coeffs, x) > max_c:
#            return True
#    return False

def check_minval(coeffs, Croot, minval = 0.01):
    for x in nup.linspace(0, 1, 100, endpoint = False):
        if chord_dist_main(coeffs, x) + Croot <= minval:
            return True
    return False

def check_maxval(coeffs, max_val):
    for x in nup.linspace(0, 1, 100):
        if chord_dist_main(coeffs, x) > max_val:
            return True
    return False

def check_valid(params, Croot, maxval):
    coeffs = coeffs_from_params(params)
    #xmax, Zmax, Zxxmax, Ztip, Zxtip = params[0], params[1], params[2], params[3], params[4] 
    if chord_dist_first(coeffs, 0) < 0:
        return False
    if check_minval(coeffs, Croot):
        return False
    if check_maxval(coeffs, maxval):
        return False
    #if check_valley(coeffs, xi):
    #    return False
    #if check_max_chord(dist):
    #    return False
    return True




def calculate_T_Q(vi, rpm, Blades, R, xi, Croot, params, airfoil = 'airfoils\\sunnysky.dat', sections = 21, rho = 1.225, dvisc = 1.8/100000, alphas = [0, 10, 1]):
    coeffs = coeffs_from_params(params)
    a1, a2, astep = alphas[0], alphas[1], alphas[2]
    kvisc = dvisc/rho
    Beta_vector = []
    r_vector = []
    dT_vector = []
    dQ_vector = []
    radps = (rpm*2*pi)/60

    for rr in nup.linspace(xi*R, R, num = sections):
        get_chord_radius = Auxiliary.linear_interpolate(xi, 1, 0, 1, rr/R)
        chord = R*chord_dist_main(coeffs, get_chord_radius) + Croot
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
    Thrust = Auxiliary.area_under_curve(r_vector, dT_vector)
    Torque = Auxiliary.area_under_curve(r_vector, dQ_vector)
    return Thrust, Torque

def create1stgen(n, R = 0.3, Croot = 0.05, limits = [[0, 1], [-1, -1], [-10, -0.0001], [-1, -1], [0, 1]]):
    if limits[1][0] == -1:
        limits[1][0] = -Croot
    if limits[1][1] == -1:
        limits[1][1] = R/2
    if limits[3][0] == -1:
        limits[3][0] = -Croot
    if limits[3][1] == -1:
        limits[3][1] = R/2
    population = []
    while len(population) < n:
        chromossome = optim.generate_chrom(limits)
        if check_valid(chromossome, Croot, maxval = R/2):
            population.append(optim.Individual(chromossome))
    return population

def evaluate_population(population, xi, Croot, vi, rpm, Blades, R, sections = 21, maxval = -1):
    if maxval == -1:
        maxval = R/2
    for dist in population:
        if check_valid(dist.chrom, Croot, maxval = maxval):
            T, Q = calculate_T_Q(vi, rpm, Blades, R, xi, Croot, dist.chrom, sections = sections)
        else:
            T, Q = 0, 1000
        dist.set_ObjVal([-T, Q])
    #population = optim.assign_fronts(population)
    #population = optim.assign_crowding(population)

def plot_dist(individual):
    coeffs = coeffs_from_params(individual.chrom)
    x_dist = nup.linspace(0, 1, 100)
    y_dist = []
    for x in x_dist:
        y_dist.append(chord_dist_main(coeffs, x))
    plt.plot(x_dist, y_dist)
    plt.show()



vi = 10
rpm = 2300
Blades = 4
R = 0.3
sections = 5

xi = 0.1
Croot = 0.03

#xmax, Zmax, Zxxmax, Ztip, Zxtip
limits = [[0, 1], [0, R/2], [-10, -0.0001], [-Croot, R/2], [0, 1]]

n_ind = 10
n_gen = 5

n_objvals = 2

#current_pop = create1stgen(n_ind, Croot = Croot, R = R, limits = limits)
#for dist in current_pop:
#    plot_dist(dist)
generation = 1
while generation <= n_gen:
    print(f'Geração: {generation}')
    if generation == 1:
        current_pop = create1stgen(n_ind, Croot = Croot, R = R, limits = limits)
        evaluate_population(current_pop, xi, Croot, vi, rpm, Blades, R, sections = sections)
        for dist in current_pop:
            objvals = dist.get_ObjVal()
            #print(objvals)
            Thrust, Torque = objvals[0], objvals[1]
            plt.scatter(Torque, Thrust, 4, "r") 
    new_gen = optim.create_nextGen(current_pop, limits, n_objvals)
    evaluate_population(new_gen, xi, Croot, vi, rpm, Blades, R, sections = sections)
    current_pop = optim.reinsert(current_pop, new_gen)
    for dist in current_pop:
        objvals = dist.get_ObjVal()
        #print(dist.get_chrom())
        #print(objvals)
        Thrust, Torque = objvals[0], objvals[1]
        if Torque == 1000:
            continue
        plt.scatter(Torque, Thrust, 4, "r")
    generation += 1
plt.show()

#for dist in current_pop:
#    temp = Chord_Dist(xi, Croot, chrom = dist.get_chrom(), coeffs = coeffs_from_params(xi, Croot, dist.get_chrom()))
#    temp.plot_show()