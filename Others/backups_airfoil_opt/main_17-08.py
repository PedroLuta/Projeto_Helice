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

n_ind = 5
n_gen = 10
search = 0.2 #percentage of parameters
n_objvals = 1
mut_rate = 2/len(absolute_limits)
t_size = 2
elitism_ratio = 0.01

V = 10
Pmax = 600
rpm_init = 2300
radps = rpm_init*2*pi/60
Blades = 4
p = 0.1
R = 0.3
sections = 5

#class state:
#    def __init__(self, V, rpm, Pmax):
#        self.V = V
#        self.rpm = rpm
#        self.radps = rpm*2*pi/60
#        self.Pmax = Pmax
#
#class prop:
#    def __init__(self, Blades, p, R, sections):
#        self.Blades = Blades
#        self.p = p
#        self.R = R
#        self.sections = sections

def evaluate_foil(chrom):#, state, prop):
    foil = PARSEC_functions.PFoil(params_vec = chrom)
    foil.write_file()
    converged = False
    for _ in range(30):
        T_dist, Q_dist, Beta_dist, r_dist = simulate.vortex(V, radps, Blades, p, R, sections = sections)
        T = Auxiliary.area_under_curve(r_dist, T_dist)
        Q = Auxiliary.area_under_curve(r_dist, Q_dist)
        if abs(Q*radps - Pmax) > 1:
            radps = ((Pmax/Q) + radps)/2
        else:
            converged = True
            break
    if converged:
        return (T/(Q*radps))*Pmax
    return 0

def is_valid(chrom):
    foil = PARSEC_functions.PFoil(params_vec = chrom)
    return PARSEC_functions.check_valid(foil.aup, foil.alo, foil.get_params_vec())

#state1 = state(V, rpm_init, Pmax)
#prop1 = prop(Blades, p, R, sections)

foil_init = PARSEC_functions.PFoil(selig_file = "airfoils\\sunnysky.dat")
T1 = evaluate_foil(foil_init.get_params_vec())
init_chrom = foil_init.get_params_vec()
search_limits = PARSEC_functions.limits_from_foil(init_chrom, factor = search)

optimized_pop = optimGA.run_GA(n_gen, n_ind, search_limits, evaluate_foil, valid_func = is_valid, mut_rate = mut_rate, t_size = t_size, elitism_ratio = elitism_ratio)

for ind in optimized_pop:
    foil = PARSEC_functions.PFoil(params_vec = ind.get_chrom())
    dists = foil.dists()
    plt.plot(dists[0], dists[1], 'b')
    plt.plot(dists[2], dists[3], 'b', label = "Optimized")
    dists2 = foil_init.dists()
    plt.plot(dists2[0], dists2[1], 'r')
    plt.plot(dists2[2], dists2[3], 'r', label = "original")
    plt.ylim(-0.5, 0.5)
    plt.xlim(0, 1)
    plt.legend()
    plt.show()
