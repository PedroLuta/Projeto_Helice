import PARSEC_functions
from importing import *

#First Generation
#DONE GENERALIZING
def create_1stGen_from_airfoil(num, foil_params, factor = 0.2, clean_foils = True):
    limits = PARSEC_functions.limits_from_foil(foil_params, factor)
    population = []
    while len(population) < num:
        params_vec = PARSEC_functions.params_gen(limits)
        temp = PARSEC_functions.PFoil()
        temp.set_parameters_vec(params_vec)
        if clean_foils:
            if PARSEC_functions.check_valid(temp.aup, temp.alo, temp.params):
                population.append(temp)
        else:
            population.append(temp)
    return population, limits



#Further Generations
#DONE GENERALIZING
def create_nextGen(population, limits, mut_rate = -1):
    if mut_rate == -1:
        mut_rate == 1/len(population[0].get_params_vec())
    next_gen = []
    while len(next_gen) < len(population):
        mother = choose_foil_tournament_Thrust(population, tournament_size = 2)
        father = choose_foil_tournament_Thrust(population, tournament_size = 2)
        son = crossover_all_genes(mother, father)
        son = mutate(son, mut_rate, limits)
        next_gen.append(son)
    return next_gen

def create_nextGen_Fronts(population, limits, mut_rate = -1):
    population = assign_fronts(population)
    population = assign_crowding(population)
    if mut_rate == -1:
        mut_rate == 1/len(population[0].get_params_vec())
    next_gen = []
    while len(next_gen) < len(population):
        mother = choose_foil_tournament_Fronts(population, tournament_size = 2)
        father = choose_foil_tournament_Fronts(population, tournament_size = 2)
        son = crossover_all_genes(mother, father)
        son = mutate(son, mut_rate, limits)
        next_gen.append(son)
    return next_gen



#Reinsert
#DONE GENERALIZING
def reinsert_Thrust(pop_ini, pop_new):
    pop_final = []
    pop_ini.extend(pop_new)
    while len(pop_final) < len(pop_new):
        best = PARSEC_functions.PFoil()
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

def reinsert_Fronts(pop_ini, pop_new):
    pop_final = []
    pop_ini.extend(pop_new)
    while len(pop_final) < len(pop_new):
        best = PARSEC_functions.PFoil()
        best.set_Front(999)
        best.set_Crowding(0)
        i = 0
        hold = 0
        while i < len(pop_ini):
            if pop_ini[i].get_Front() < best.get_Front():
                best = pop_ini[i]
                hold = i
            elif pop_ini[i].get_Front() == best.get_Front():
                if pop_ini[i].get_Crowding() > best.get_Crowding():
                    best = pop_ini[i]
                    hold = i
            i += 1
        pop_final.append(best)
        pop_ini.pop(hold)
    return pop_final



#Crossover functions
#DONE GENERALIZING
def crossover_single_point(mother, father):
    son1 = PARSEC_functions.PFoil()
    son2 = PARSEC_functions.PFoil()
    mother_genes = mother.get_params_vec()
    father_genes = father.get_params_vec()
    crossover_point = random.randint(0, len(mother_genes) - 1)
    son1_genes = mother_genes[:crossover_point]
    son1_genes.extend(father_genes[crossover_point:])
    son2_genes = father_genes[:crossover_point]
    son2_genes.extend(mother_genes[crossover_point:])
    son1.set_parameters_vec(son1_genes)
    son2.set_parameters_vec(son2_genes)
    return son1, son2

def crossover_two_points(mother, father):
    son1 = PARSEC_functions.PFoil()
    son2 = PARSEC_functions.PFoil()
    mother_genes = mother.get_params_vec()
    father_genes = father.get_params_vec()
    crossover_point1 = random.randint(0, len(mother_genes) - 1)
    crossover_point2 = random.randint(0, len(mother_genes) - 1)
    if (crossover_point2 == crossover_point1) or (abs(crossover_point1 - crossover_point2) == len(mother_genes) - 1):
        return mother, father
    if crossover_point2 < crossover_point1:
        temp = crossover_point1
        crossover_point1 = crossover_point2
        crossover_point2 = temp
    son1_genes = mother_genes[:crossover_point1]
    son1_genes.extend(father_genes[crossover_point1:crossover_point2])
    son1_genes.extend(mother_genes[crossover_point2:])
    son2_genes = father_genes[:crossover_point1]
    son2_genes.extend(mother_genes[crossover_point1:crossover_point2])
    son2_genes.extend(father_genes[crossover_point2:])
    son1.set_parameters_vec(son1_genes)
    son2.set_parameters_vec(son2_genes)
    return son1, son2

def crossover_all_genes(mother, father):
    son = PARSEC_functions.PFoil()
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



#Mutation
#DONE GENERALIZING
def mutate(foil, mut_rate, limits):
    i = 0
    chromossome = foil.get_params_vec()
    while i < len(chromossome):
        mut = random.random()
        if mut < mut_rate:
            chromossome[i] = PARSEC_functions.params_gen_single(limits[i])
        i += 1
    foil.set_parameters_vec(chromossome)
    return foil

def mutate_generalized(chrom, mut_rate, limits):
    i = 0
    while i < len(chrom):
        mut = random.random()
        if mut < mut_rate:
            chrom[i] = generate_random_parameter(limits[i])
        i += 1
    return chrom



#Chromossome Selection
#DONE GENERALIZING
def choose_foil_tournament_Thrust(population, tournament_size = 2):
    index_list = []
    for _ in range(tournament_size):
        rand_index = random.randint(0, len(population) - 1)
        while rand_index in index_list:
            rand_index = random.randint(0, len(population) - 1)
        index_list.append(rand_index)
    best = PARSEC_functions.PFoil()
    for index in index_list:
        if population[index].Thrust > best.Thrust:
            best = population[index]
    return best

def choose_foil_tournament_Fronts(population, tournament_size = 2):
    index_list = []
    for _ in range(tournament_size):
        rand_index = random.randint(0, len(population) - 1)
        while rand_index in index_list:
            rand_index = random.randint(0, len(population) - 1)
        index_list.append(rand_index)
    best = PARSEC_functions.PFoil()
    best.set_Front(999)
    best.set_Crowding(0)
    for index in index_list:
        if population[index].get_Front() < best.get_Front():
            best = population[index]
        elif population[index].get_Front() == best.get_Front():
            if population[index].get_Crowding() > best.get_Crowding():
                best = population[index]
    return best

def assign_domination(population):
    i = 0
    if len(population) == 1:
        population[i].set_Dominated_counter(0)
        return
    elif len(population) == 0:
        return
    while i < len(population):
        dominated_counter = 0
        j = 0
        while j < len(population):
            if j == i:
                j += 1
                continue
            if ((population[j].Thrust >= population[i].Thrust) and (population[j].Torque <= population[i].Torque)) \
                and ((population[j].Thrust > population[i].Thrust) or (population[j].Torque < population[i].Torque)):
                dominated_counter += 1
            j += 1
        population[i].set_Dominated_counter(dominated_counter)
        i += 1

def assign_fronts(population):
    copy = population[:]
    assign_domination(copy)
    new_pop = []
    front = 1
    while len(copy) > 0:
        i = 0
        while i < len(copy):
            if copy[i].get_Dominated_counter() == 0:
                copy[i].set_Front(front)
                new_pop.append(copy[i])
                copy.pop(i)
            else:
                i += 1
        front += 1
        assign_domination(copy)
    return new_pop

def assign_crowding(population):
    copy = sort_Thrust(population)
    copy[0].set_Crowding(9999)
    copy[-1].set_Crowding(9999)
    maxT = copy[0].Thrust
    minT = copy[-1].Thrust
    i = 1
    while i < len(copy) - 1:
        dT = abs(PARSEC_functions.linear_interpolate(minT, maxT, 0, 1, copy[i + 1].Thrust) - PARSEC_functions.linear_interpolate(minT, maxT, 0, 1, copy[i - 1].Thrust))
        copy[i].set_Crowding(dT)
        i += 1
    copy = sort_Torque(copy)
    copy[0].set_Crowding(9999)
    copy[-1].set_Crowding(9999)
    minQ = copy[0].Torque
    maxQ = copy[-1].Torque
    i = 1
    while i < len(copy) - 1:
        dQ = abs(PARSEC_functions.linear_interpolate(minQ, maxQ, 0, 1, copy[i + 1].Torque) - PARSEC_functions.linear_interpolate(minQ, maxQ, 0, 1, copy[i - 1].Torque))
        copy[i].set_Crowding(copy[i].get_Crowding() + dQ)
        i += 1
    return copy


#Auxiliary
#DONE GENERALIZING (ALWAYS MINIMIZE)
def sort_Thrust(population):
    new_pop = []
    copy = population[:]
    while len(copy) > 0:
        maxT = -1
        i = 0
        hold = 0
        while i < len(copy):
            if copy[i].Thrust > maxT:
                maxT = copy[i].Thrust
                hold = i
            i += 1
        new_pop.append(copy[hold])
        copy.pop(hold)
    return new_pop

def sort_Torque(population):
    new_pop = []
    copy = population[:]
    while len(copy) > 0:
        minQ = 1001
        i = 0
        hold = 0
        while i < len(copy):
            if copy[i].Torque < minQ:
                minQ = copy[i].Torque
                hold = i
            i += 1
        new_pop.append(copy[hold])
        copy.pop(hold)
    return new_pop

def calculate_domination(indA, indB, n_objvals):
    i = 0
    while i < n_objvals:
        if i == 0:
            boolA = indA.get_ObjVal(i) <= indB.get_ObjVal(i)
            boolB = indA.get_ObjVal(i) < indB.get_ObjVal(i)
            i += 1
        boolA = boolA and (indA.get_ObjVal(i) <= indB.get_ObjVal(i))
        boolB = boolB or (indA.get_ObjVal(i) < indB.get_ObjVal(i))
        i += 1
    return (boolA and boolB) #if True, A dominates B, if False, A does not dominate B

def generate_random_parameter(limits):
    parameter = random.uniform(limits[0], limits[1])
    return parameter
