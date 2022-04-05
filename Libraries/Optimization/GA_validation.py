from importing import *

#functions import
import PARSEC_functions
import optimization_NSGA2 as optim
import optimization_GA as optimGA

from pymoo.factory import get_problem
from pymoo.util.plotting import plot


absolute_limits = [[-32.768, 32.768], [-32.768, 32.768]]

n_ind = 20
n_gen = 30
mut_rate = 0.1
t_size = 2

def evaluate(chrom):
    a = 20
    b = 1/5
    c = 2*pi
    sum_squared = 0
    sum_cos = 0
    for x in chrom:
        sum_squared += x**2
        sum_cos += math.cos(c*x)
    f = -a*euler**(-b*(sum_squared/len(chrom))**0.5) - euler**(sum_cos/len(chrom)) + a + euler
    return f

def is_valid(chrom):
    return True

optimized_pop = optimGA.run_GA(n_gen, n_ind, absolute_limits, evaluate, mut_rate = mut_rate, t_size = t_size, valid_func = is_valid) 


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
for individual in optimized_pop:
    chrom = individual.get_chrom()
    func = evaluate(chrom)
    #print(f'chrom: {individual.get_chrom()} Front: {individual.get_Front()} Crowding: {individual.get_Crowding()} Dominated_counter: {individual.get_Dominated_counter()}')
    ax.scatter(chrom[0], chrom[1], func)
plt.show()
#problem = get_problem("bnh")
#plt.title("População final")
#plot(problem.pareto_front(), no_fill = True)
#plt.show()