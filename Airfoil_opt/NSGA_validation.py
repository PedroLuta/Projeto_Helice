from importing import *

#functions import
import PARSEC_functions
import optimization_NSGA2 as optim
import optimization_GA as optimGA

from pymoo.factory import get_problem
from pymoo.util.plotting import plot




def evaluate(chrom):
    x1, x2 = chrom[0], chrom[1]
    f1 = 4*(x1**2) + 4*(x2**2)
    f2 = (x1 - 5)**2 + (x2 - 5)**2
    return [f1, f2]

def is_valid(chrom):
    x1, x2 = chrom[0], chrom[1]
    C1 = (x1 - 5)**2 + x2**2 #≤ 25
    C2 = (x1 - 8)**2 + (x2 + 3)**2 #≥ 7.7
    if (C1 > 25) or C2 < 7.7:
        return False
    return True

from pymoo.factory import get_problem
from pymoo.util.plotting import plot

absolute_limits = [[0, 5], [0, 3]]

n_ind = 300
archive_size = 600
mut_rate = 0.0
t_size = 2
convergence = 4

optimized_pop = optim.run_GA_convergence_new(n_ind, absolute_limits, evaluate, mut_rate = mut_rate, t_size = t_size, valid_func = is_valid, gwcfc = convergence, archive_size = archive_size)

for individual in optimized_pop:
    funcs = evaluate(individual.get_chrom())
    plt.plot(funcs[0], funcs[1], 'bo')
plt.plot(funcs[0], funcs[1], 'bo', label = "Individual")
plt.plot(funcs[0], funcs[1], 'r', label = "Pareto-front")
problem = get_problem("bnh")
plt.title("NSGA-II Convergence test")
plt.xlabel("f1")
plt.ylabel("f2")
plot(problem.pareto_front(), no_fill = True, labels = 'Pareto-front')
#plt.legend()
#plt.show()