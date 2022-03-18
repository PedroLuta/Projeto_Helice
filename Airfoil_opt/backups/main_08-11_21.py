#from induction_qprop import induction_qprop
from importing import *

#functions import
import PARSEC_functions
import optimization_NSGA2 as optim
import optimization_GA as optimGA
#import optimization_GA as optimGA
import Auxiliary
import simulate

def evaluate_foil(chrom):
    global V 
    global Pmax 
    global radps 
    global Blades 
    global p 
    global R
    foil = PARSEC_functions.PFoil(params_vec = chrom)
    foil.write_file()
    converged = False
    for _ in range(15):
        T_dist, Q_dist, Beta_dist, r_dist = simulate.vortex(V[0], radps, Blades, p, R, sections = sections)
        T = Auxiliary.area_under_curve(r_dist, T_dist)
        Q = Auxiliary.area_under_curve(r_dist, Q_dist)
        if abs(Q*radps - Pmax) > 10:
            radps = ((Pmax/Q) + radps)/2
        else:
            converged = True
            break
    T1 = (T/(Q*radps))*Pmax
    if converged:
        converged = False
        for _ in range(15):
            T_dist, Q_dist, r_vector, _, _, _, _ = simulate.fixed_pitch_qprop(V[1], radps, Blades, p, R, Beta_dist, sections = sections)
            T = Auxiliary.area_under_curve(r_vector, T_dist)
            Q = Auxiliary.area_under_curve(r_vector, Q_dist)
            if abs(Q*radps - Pmax) > 10:
                radps = ((Pmax/Q) + radps*2)/3
            else:
                converged = True
                break
        T2 = (T/(Q*radps))*Pmax
        if converged:
            print("done simulating")
            return [-T1, -T2]
    print("discarded")
    return [0, 0]

def evaluate_foil2(chrom):
    global V 
    global Pmax 
    global radps 
    global Blades 
    global p 
    global R
    foil = PARSEC_functions.PFoil(params_vec = chrom)
    foil.write_file()
    T_dist, Q_dist, Beta_dist, r_dist = simulate.vortex(V[0], radps_vec[0], Blades, p, R, sections = sections)
    T1 = Auxiliary.area_under_curve(r_dist, T_dist)
    Q1 = Auxiliary.area_under_curve(r_dist, Q_dist)
    T_dist, Q_dist, r_vector, _, _, _, _ = simulate.fixed_pitch_qprop(V[1], radps_vec[1], Blades, p, R, Beta_dist, sections = sections)
    T2 = Auxiliary.area_under_curve(r_vector, T_dist)
    Q2 = Auxiliary.area_under_curve(r_vector, Q_dist)
    print("Done")
    return [-T1/Q1, -T2/Q2]

def evaluate_foil_single(chrom):
    global V 
    global Pmax 
    global radps 
    global Blades 
    global p 
    global R
    foil = PARSEC_functions.PFoil(params_vec = chrom)
    foil.write_file()
    converged = False
    for _ in range(20):
        T_dist, Q_dist, Beta_dist, r_dist = simulate.vortex(V[0], radps, Blades, p, R, sections = sections)
        T = Auxiliary.area_under_curve(r_dist, T_dist)
        Q = Auxiliary.area_under_curve(r_dist, Q_dist)
        if abs(Q*radps - Pmax) > 10:
            radps = ((Pmax/Q) + radps)/2
        else:
            converged = True
            break
    T1 = (T/(Q*radps))*Pmax
    if converged:
        print("Done converged")
        return -T1
    print("fail")
    return 0

def is_valid(chrom):
    foil = PARSEC_functions.PFoil(params_vec = chrom)
    return PARSEC_functions.check_valid(foil.aup, foil.alo, foil.get_params_vec())

#absolute_limits = [[-30, 15], [0.0001, 45], [-0.05, 0.05], [0.0001, 0.05], [0.0001, 0.1], [0.0001, 0.8], [0.0001, 0.35], [-5, 0.0001], [0.0001, 0.1], [0.0001, 0.5], [-0.35, -0.0001], [0.0001, 10]]
absolute_limits = [[-30, 15], [0.0001, 45], [0, 0], [0.0001, 0.05], [0.0001, 0.1], [0.0001, 0.8], [0.0001, 0.35], [-5, 0.0001], [0.0001, 0.1], [0.0001, 0.5], [-0.35, -0.0001], [0.0001, 10]]

n_ind = 25
n_gen = 30
convergence = 15 #generations without progress for convergence
search = 0.2 #percentage of parameters
n_objvals = 1
mut_rate = 2/len(absolute_limits)
t_size = 2
#elitism_ratio = 0.01
V = [0, 6]
Pmax = 580
rpm = [2300, 2500]
rpm_init = 2300
radps_vec = [x*2*pi/60 for x in rpm]
radps = rpm_init*2*pi/60
Blades = 4
p = 0.1
R = 0.3
sections = 5

if __name__ == '__main__':
    foil_init = PARSEC_functions.PFoil(selig_file = "airfoils\\sunnysky.dat")
    init_chrom = foil_init.get_params_vec()
    Tinit = evaluate_foil(init_chrom)
    print(Tinit)
    search_limits = PARSEC_functions.limits_from_foil(init_chrom, factor = search)
    search_limits[2] = [0, 0]

    optimized_pop = optim.run_GA_convergence_test(n_ind, search_limits, evaluate_foil, valid_func = is_valid, mut_rate = mut_rate, t_size = t_size, gwcfc = convergence)

    with open("Solutions\\out.txt", 'w') as out:
       out.write("Solutions found by the NSGA2 for conditions below:\n")
       for v in V:
           out.write(f"{v}\n")
       out.write("\n\n")
       for ind in optimized_pop:
           out.write(f'{ind.get_chrom()}\n')
           for x in range(len(ind.get_ObjVal())):
               out.write(f'Thrust_{x}: {ind.get_ObjVal(x)}\n')
           out.write('\n')

