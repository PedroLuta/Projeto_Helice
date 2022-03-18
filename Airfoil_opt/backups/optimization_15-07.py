from PARSEC_functions import *
import math

#Genetic Algorythm
def create_1stGen_from_airfoil(num, foil_params, factor = 0.2, clean_foils = True):
    limits = limits_from_foil(foil_params, factor)
    population = []
    while len(population) < num:
        params_vec = params_gen(limits)
        temp = PFoil()
        temp.set_parameters_vec(params_vec)
        if clean_foils:
            if check_valid(temp.aup, temp.alo, temp.params):
                population.append(temp)
        else:
            population.append(temp)
    return population, limits

def create_1stGen_random(num, clean_foils = True):
    population = []
    while len(population) < num:
        temp = PFoil(random = True)
        if clean_foils:
            if check_valid(temp.aup, temp.alo, temp.params):
                population.append(temp)
        else:
            population.append(temp)
    return population

#def assign_random_probability(population):
#    prob_list = []
#    for _ in population:
#        prob_list.append(random.random())
#    i = 0
#    for foil in population:
#        foil.Probability = prob_list[i]/sum(prob_list)
#        i += 1
#
#def assign_uniform_probability(population):
#    num_foils = len(population)
#    for foil in population:
#        foil.Probability = 1/num_foils
#
#def assign_sorted_probability(population):
#    sum_Thrust = 0
#    for foil in population:
#        if math.isnan(foil.Thrust):
#            pass
#        else:
#            sum_Thrust += foil.Thrust
#    for foil in population:
#        if math.isnan(foil.Thrust):
#            foil.Probability = 0
#        else:
#            foil.Probability = foil.Thrust/sum_Thrust

def crossover_all_genes(mother, father):
    son = PFoil()
    mother_genes = mother.get_params_vec()
    father_genes = father.get_params_vec()
    son_genes = []
    i = 0
    while i < len(mother_genes):
        chosen = random.random()
        if chosen > 0.5:
            son_genes.append(mother_genes[i])
        else:
            son_genes.append(father_genes[i])
        i += 1
    son.set_parameters_vec(son_genes)
    return son

def choose_random_foil_intrinsic_probability(population):
    random.seed()
    chosen = random.random()
    accumulator = 0
    for foil in population:
        accumulator += foil.Probability
        if chosen < accumulator:
            chosen_foil = foil
            break
    return chosen_foil

def choose_foil_tournament(population, tournament_size = 2):
    index_list = []
    for _ in range(tournament_size):
        rand_index = random.randint(0, len(population) - 1)
        while rand_index in index_list:
            rand_index = random.randint(0, len(population) - 1)
        index_list.append(rand_index)
    best = PFoil()
    for index in index_list:
        if population[index].Thrust > best.Thrust:
            best = population[index]
    return best

def create_nextGen(population, limits, mut_rate = -1, mut_range = -1):
    if mut_rate == -1:
        mut_rate == 1/len(population[0].get_params_vec())
    if mut_range == -1:
        mut_range = 0.2
    next_gen = []
    while len(next_gen) < len(population):
        #mother = choose_random_foil_intrinsic_probability(population)
        mother = choose_foil_tournament(population, tournament_size = 2)
        #father = choose_random_foil_intrinsic_probability(population)
        father = choose_foil_tournament(population, tournament_size = 2)
        son = crossover_all_genes(mother, father)
        son = mutate(son, mut_rate, mut_range, limits)
        #if check_valid(son.aup, son.alo, son.get_params_vec()):
        next_gen.append(son)
    return next_gen

def mutate(foil, mut_rate, mut_range, limits):
    i = 0
    chromossome = foil.get_params_vec()
    while i < len(chromossome):
        mut = random.random()
        if mut < mut_rate:
            chromossome[i] = random.uniform(chromossome[i] - chromossome[i]*mut_range, chromossome[i] + chromossome[i]*mut_range)
            if chromossome[i] < limits[i][0]:
                chromossome[i] = limits[i][0]
            if chromossome[i] > limits[i][1]:
                chromossome[i] = limits[i][1]
        i += 1
    foil.set_parameters_vec(chromossome)
    return foil

def reinsert(pop_ini, pop_new):
    pop_final = []
    pop_ini.extend(pop_new)
    while len(pop_final) < len(pop_new):
        best = PFoil()
        i = 0
        hold = 0
        while i < len(pop_ini):
            if pop_ini[i].Thrust > best.Thrust:
                best = pop_ini[i]
                hold = i
            i += 1
        pop_final.append(best)
        pop_ini.pop(hold)
    return pop_final

#pop = [PFoil(Thrust = 1), PFoil(Thrust = 2), PFoil(Thrust = 3), PFoil(Thrust = 4), PFoil(Thrust = 5)]
#while True:
#    best = choose_foil_tournament(pop, tournament_size = 3)
#    print(best.Thrust)

