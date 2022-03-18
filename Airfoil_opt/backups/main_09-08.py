#from induction_qprop import induction_qprop
from importing import *

#functions import
import PARSEC_functions
import optimization_NSGA2 as optim
import optimization_GA as optimGA
import Auxiliary
import simulate


absolute_limits = [[-30, 15], [0.0001, 45], [-0.05, 0.05], [0.0001, 0.05], [0.0001, 0.1], [0.0001, 0.8], [0.0001, 0.35], [-5, 0.0001], [0.0001, 0.1], [0.0001, 0.5], [-0.35, -0.0001], [0.0001, 10]]
#[ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo]

n_ind = 10
n_gen = 10
search = 0.2 #percentage of parameters
n_objvals = 1
mut_rate = 2/len(absolute_limits)
t_size = 2

V = 10
Pmax = 600
rpm_init = 2300
Blades = 4
p = 0.1
R = 0.3
sections = 5

class state:
    def __init__(self, V, rpm, Pmax):
        self.V = V
        self.rpm = rpm
        self.radps = rpm*2*pi/60
        self.Pmax = Pmax

class prop:
    def __init__(self, Blades, p, R, sections):
        self.Blades = Blades
        self.p = p
        self.R = R
        self.sections = sections

def evaluate_foil(chrom):#, state, prop):
    foil = PARSEC_functions.PFoil(params_vec = chrom)
    foil.write_file()
    radps = state1.radps
    converged = False
    for _ in range(50):
        T_dist, Q_dist, Beta_dist, r_dist = simulate.vortex(state1.V, state1.radps, prop1.Blades, prop1.p, prop1.R, sections = prop1.sections)
        T = Auxiliary.area_under_curve(r_dist, T_dist)
        Q = Auxiliary.area_under_curve(r_dist, Q_dist)
        if abs(Q*radps - state1.Pmax) > 1:
            radps = ((state1.Pmax/Q) + radps)/2
        else:
            converged = True
            break
    if converged:
        return (T/(Q*radps))*state1.Pmax
    return 0

def is_valid(chrom):
    foil = PARSEC_functions.PFoil(params_vec = chrom)
    return PARSEC_functions.check_valid(foil.aup, foil.alo, foil.get_params_vec())

#def evaluate_population(pop):#, state, prop):
#    for individual in pop:
#        foil = PARSEC_functions.PFoil(params_vec = individual.get_chrom())
#        if PARSEC_functions.check_valid(foil.aup, foil.alo, foil.get_params_vec()):
#            T = evaluate_foil(foil, state, prop)
#            if T == infinite:
#                individual.set_valid(False)
#            else:
#                individual.set_ObjVal(-T)
#        else:
#            individual.set_valid(False)
        

state1 = state(V, rpm_init, Pmax)
prop1 = prop(Blades, p, R, sections)

foil_init = PARSEC_functions.PFoil(selig_file = "airfoils\\sunnysky.dat")
T1 = evaluate_foil(foil_init.get_params_vec())#, state1, prop1)
init_chrom = foil_init.get_params_vec()
search_limits = PARSEC_functions.limits_from_foil(init_chrom, factor = search)

optimized_pop = optimGA.run_GA(n_gen, n_ind, search_limits, evaluate_foil, valid_func = is_valid, mut_rate = mut_rate, t_size = t_size)

#individuals = 0 
#current_pop = []
#while individuals < n_ind:
#    chromossome = optimGA.generate_chrom(search_limits)
#    temp = PARSEC_functions.PFoil(params_vec = chromossome)
#    #temp.set_parameters_vec(chromossome)
#    if PARSEC_functions.check_valid(temp.aup, temp.alo, temp.get_params_vec()):
#        current_pop.append(optimGA.Individual(chrom = chromossome))
#        individuals += 1
#print("Geração 1")
#evaluate_population(current_pop, state1, prop1)

#generation = 1
#for generation in range(2, n_gen + 1): 
#    print(f"Geração {generation}")
#    new_gen = optimGA.create_offspring(current_pop, search_limits, mut_rate = mut_rate, t_size = t_size)
#    print("Analisando Descendentes")
#    evaluate_population(new_gen, state1, prop1) 
#    current_pop = optimGA.reinsert_immortal_individuals(current_pop, new_gen)

#x3, y3, x4, y4 = foil_init.dists()

#for foil in current_pop:
#    temp = PARSEC_functions.PFoil(params_vec = foil.get_chrom())
#    #temp.set_parameters_vec(foil.get_chrom())
#    x1, y1, x2, y2 = temp.dists()
#    plt.plot(x1, y1, 'b', label = f"optimized Thrust = {foil.get_ObjVal()}")
#    plt.plot(x2, y2, 'b')
#    plt.plot(x3, y3, 'r', label = f"initial Thrust = {T1}")
#    plt.plot(x4, y4, 'r')
#    plt.ylim(-0.5, 0.5)
#    plt.xlim(0, 1)
#    plt.legend()
#    plt.show()
