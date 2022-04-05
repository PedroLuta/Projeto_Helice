from importing import *
import PARSEC_functions
import main

file = "Solutions\\10.txt"
n_header_lines = 5
n_objval = 2

pop = []
with open(file, 'r') as inp:
    Tinit = inp.readline().replace('[', '').replace(']', '').split(',')
    Tinit1 = float(Tinit[0])
    Tinit2 = float(Tinit[1])
    for _ in range(n_header_lines):
        inp.readline()
    while True:
        line = inp.readline()
        if line == '':
            break
        if line == '\n':
            continue
        line = line.replace('[', '').replace(']', '').split(',')
        chrom_vec = []
        for value in line:
            chrom_vec.append(float(value))
        objval_vec = []
        for _ in range(n_objval):
            content = inp.readline().split()
            value = float(content[-1])
            objval_vec.append(value)
        pop.append([chrom_vec, objval_vec])




foil_init = PARSEC_functions.PFoil(selig_file = "airfoils\\sunnysky.dat")
init_chrom = foil_init.get_params_vec()
#Tinit = main.evaluate_foil_single(init_chrom)

for ind in pop:
   foil = PARSEC_functions.PFoil(params_vec = ind[0])
   dists = foil.dists()
   plt.plot(dists[0], dists[1], 'b')
   plt.plot(dists[2], dists[3], 'b', label = f"Optimized")
   dists2 = foil_init.dists()
   plt.plot(dists2[0], dists2[1], 'r')
   plt.plot(dists2[2], dists2[3], 'r', label = f"original")
   plt.ylim(-0.5, 0.5)
   plt.xlim(0, 1)
   plt.legend()
   plt.show()

# for ind in pop:
#     plt.plot(ind[1][0], ind[1][1], 'bo', ms = 3)
# if Tinit1 != 0:
#     plt.plot(Tinit1, Tinit2, 'ro', ms = 5, label = "Original")
# plt.plot(pop[0][1][0], pop[0][1][1], 'bo', ms = 3, label = "Optimized")
# plt.xlabel("-Thrust at condition 1")
# plt.ylabel("-Thrust at condition 2")
# plt.legend()
# plt.show()