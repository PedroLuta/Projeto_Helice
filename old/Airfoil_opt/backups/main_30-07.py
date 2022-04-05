#from induction_qprop import induction_qprop
from importing import *

#functions import
#import xfoil_interface
#import induction
import PARSEC_functions
import optimization_old
import optimization_NSGA2
import Auxiliary
import simulate

def evaluate_foil(foil, V, rpm, Blades, p, R, sections):
    foil.write_file()
    T_dist, Q_dist, Beta_dist, r_dist = simulate.vortex(V, rpm, Blades, p, R, sections = sections)
    #T_dist, Q_dist, Beta_dist, r_dist = simulate.momentum_Ftip(V, rpm, Blades, p, R, sections = sections)
    #T_dist, Q_dist, Beta_dist, r_dist = simulate.momentum_simple(V, rpm, Blades, p, R, sections = sections)
    T = Auxiliary.area_under_curve(r_dist, T_dist)
    Q = Auxiliary.area_under_curve(r_dist, Q_dist)
    foil.set_TQ(T, Q)

def evaluate_population(pop, V, rpm, Blades, p, R, sections):
    for foil in pop:
        if PARSEC_functions.check_valid(foil.aup, foil.alo, foil.get_params_vec()):
            evaluate_foil(foil, V, rpm, Blades, p, R, sections)
        else:
            foil.set_TQ(-0.5, 1000)
    copy = pop[:]
    copy = optimization_old.assign_fronts(copy)
    copy = optimization_old.assign_crowding(copy)
    return copy

def plot_population_TQ(population):
    for foil in population:
        plt.plot(foil.Torque, -foil.Thrust, "ro")
    plt.show()

def plot_comparative(population, init_foil):
    for foil in population:
        x1, y1, x2, y2 = foil.dists()
        x3, y3, x4, y4 = init_foil.dists()
        plt.plot(x1, y1, 'b', label = "optimized")
        plt.plot(x2, y2, 'b')
        plt.plot(x3, y3, 'r', label = "initial")
        plt.plot(x4, y4, 'r')
        plt.legend()
        plt.show()


absolute_limits = [[-30, 15], [0.0001, 45], [-0.05, 0.05], [0.0001, 0.05], [0.0001, 0.1], [0.0001, 0.8], [0.0001, 0.35], [-5, 0.0001], [0.0001, 0.1], [0.0001, 0.5], [-0.35, -0.0001], [0.0001, 10]]
#[ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo]

n_ind = 5
n_gen = 5
search = 0.2

V = 10
rpm = 2300
Blades = 4
p = 0.15
R = 0.3
sections = 5

foil_init = PARSEC_functions.PFoil(selig_file = "airfoils\\sunnysky.dat")
evaluate_foil(foil_init, V, rpm, Blades, p, R, sections = 5)
#print(foil_init.Thrust)
#print(foil_init.Torque)
init_chrom = foil_init.get_params_vec()
size_chrom = len(init_chrom)

generation = 1
#best = []
while generation <= n_gen:
    print(f"Geração {generation}")
    if generation == 1:
        current_pop, limits = optimization_old.create_1stGen_from_airfoil(n_ind, init_chrom, factor = search)
        current_pop = evaluate_population(current_pop, V, rpm, Blades, p, R, sections)
        for foil in current_pop:
            plt.scatter(foil.Torque, -foil.Thrust, 4, "r") 
    #new_gen = create_nextGen(current_pop, limits, mut_rate = 2/size_chrom)
    new_gen = optimization_old.create_nextGen_Fronts(current_pop, limits, mut_rate = 2/size_chrom)
    new_gen = evaluate_population(new_gen, V, rpm, Blades, p, R, sections) 
    #current_pop = reinsert_Thrust(current_pop, new_gen)
    current_pop = optimization_old.reinsert_Fronts(current_pop, new_gen)
    for foil in current_pop:
        if foil.Torque == 1000:
            continue
        plt.scatter(foil.Torque, -foil.Thrust, 4,'r') 
    #best.append(current_pop[0])
    generation += 1
plt.scatter(foil_init.Torque, -foil_init.Thrust, 4, 'b')
plt.show()
plot_comparative(current_pop, foil_init)