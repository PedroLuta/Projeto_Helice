#from induction_qprop import induction_qprop
from importing import *

#functions import
#import xfoil_interface
#import induction
import PARSEC_functions
#import optimization_old
import optimization_NSGA2 as optim
import optimization_GA as optimGA
import Auxiliary
import simulate

def evaluate_foil(foil, V, radps, Blades, p, R, sections, Pmax):
    foil.write_file2()
    tries = 0
    converged = False
    while tries < 15:
        T_dist, Q_dist, Beta_dist, r_dist = simulate.vortex(V, radps, Blades, p, R, sections = sections)
        #T_dist, Q_dist, Beta_dist, r_dist = simulate.momentum_Ftip(V, radps, Blades, p, R, sections = sections)
        #T_dist, Q_dist, Beta_dist, r_dist = simulate.momentum_simple(V, radps, Blades, p, R, sections = sections)
        T = Auxiliary.area_under_curve(r_dist, T_dist)
        Q = Auxiliary.area_under_curve(r_dist, Q_dist)
        if abs(Q*radps - Pmax) < 1: #(Q*radps > Pmax + 1) or (Pmax - Q*radps > 10):
            
            radps = ((Pmax/Q) + radps*2)/3
            #radps = radps_new
            #radps = radps + ((Pmax - Q*radps)/Q)
            #print(f'{radps*Q} {Q} {T}')
            #print(f"iterating, try number {tries}")
            tries += 1
        else:
            converged = True
            break
            #print("not iterating")
    #print("Done")
    if converged:
        return (T/(Q*radps))*Pmax#, (Q/(Q*radps))*Pmax
    return -0.5#, 1000

def evaluate_population(pop, V, radps, Blades, p, R, sections, Pmax):
    for individual in pop:
        foil = PARSEC_functions.PFoil()
        foil.set_parameters_vec(individual.get_chrom())
        if PARSEC_functions.check_valid(foil.aup, foil.alo, foil.get_params_vec()):
            T= evaluate_foil(foil, V, radps, Blades, p, R, sections, Pmax) #, Q 
        else:
            T = -0.5#, 1000 #, Q
        individual.set_ObjVal([-T])


absolute_limits = [[-30, 15], [0.0001, 45], [-0.05, 0.05], [0.0001, 0.05], [0.0001, 0.1], [0.0001, 0.8], [0.0001, 0.35], [-5, 0.0001], [0.0001, 0.1], [0.0001, 0.5], [-0.35, -0.0001], [0.0001, 10]]
#[ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo]

n_ind = 15
n_gen = 10
search = 0.2 #percentage of parameters
n_objvals = 1
mut_rate = -1
t_size = 2

V = 10
Pmax = 600
rpm_init = 2300
radps = 2*pi*rpm_init/60
Blades = 4
p = 0.1
R = 0.3
sections = 5

foil_init = PARSEC_functions.PFoil(selig_file = "airfoils\\sunnysky.dat")
init_chrom = foil_init.get_params_vec()
size_chrom = len(init_chrom)
search_limits = PARSEC_functions.limits_from_foil(init_chrom, factor = search)
#print(search_limits)
generation = 1
while generation <= n_gen:
    print(f"Geração {generation}")
    if generation == 1:
        individuals = 0 
        current_pop = []
        while individuals < n_ind:
            #chromossome = optim.generate_chrom(absolute_limits)
            #chromossome = optim.generate_chrom(search_limits)
            #chromossome = optim.generate_chrom(absolute_limits)
            chromossome = optimGA.generate_chrom(search_limits)
            #current_pop.append(optim.Individual(chrom = chromossome))
            #individuals += 1
            temp = PARSEC_functions.PFoil()
            temp.set_parameters_vec(chromossome)
            if PARSEC_functions.check_valid(temp.aup, temp.alo, temp.get_params_vec()):
                #print(chromossome)
                current_pop.append(optimGA.Individual(chrom = chromossome))
                individuals += 1
        evaluate_population(current_pop, V, radps, Blades, p, R, sections, Pmax)
        #for foil in current_pop:
        #    if foil.get_ObjVal(1) == 1000:
        #        continue
        #    plt.scatter(foil.get_ObjVal(1), foil.get_ObjVal(0), 4, "r")
    new_gen = optimGA.create_nextGen(current_pop, search_limits, mut_rate = mut_rate, t_size = t_size)
    print("Analisando Descendentes")
    evaluate_population(new_gen, V, radps, Blades, p, R, sections, Pmax) 
    #for foil in new_gen:
    #    #print(foil.get_chrom())
    #    if foil.get_ObjVal(1) == 1000:
    #        continue
    #    plt.scatter(foil.get_ObjVal(1), foil.get_ObjVal(0), 4,'r') 
    current_pop = optimGA.reinsert(current_pop, new_gen)
    generation += 1

T1 = evaluate_foil(foil_init, V, radps, Blades, p, R, sections, Pmax) #, Q1
#plt.scatter(Q1, -T1, 4, 'b')
#plt.show()

for foil in current_pop:
    print(foil.get_Front())
    print(foil.get_Crowding())
    print(foil.get_Dominated_counter())
    #print(foil.get_Dominated_counter())
    print(foil.get_ObjVal())
    temp = PARSEC_functions.PFoil()
    temp.set_parameters_vec(foil.get_chrom())
    x1, y1, x2, y2 = temp.dists()
    x3, y3, x4, y4 = foil_init.dists()
    plt.plot(x1, y1, 'b', label = "optimized")
    plt.plot(x2, y2, 'b')
    plt.plot(x3, y3, 'r', label = "initial")
    plt.plot(x4, y4, 'r')
    plt.ylim(-0.5, 0.5)
    plt.xlim(0, 1)
    plt.legend()
    plt.show()
