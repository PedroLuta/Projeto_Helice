from importing import *
import Auxiliary

class Individual:
    def __init__(self, chrom = [], ObjVal = 0, valid = True):
        self.chrom = chrom
        self.ObjVal = ObjVal
        self.valid = valid

    def set_chrom(self, chrom):
        self.chrom = chrom
    def get_chrom(self):
        return self.chrom

    def set_valid(self, valid):
        self.valid = valid 
    def get_valid(self):
        return self.valid

    def set_ObjVal(self, value):
        self.ObjVal = value
    def get_ObjVal(self):
        return self.ObjVal


def run_GA_convergence_test(n_ind, gene_limits, evaluate_func, valid_func = Auxiliary.return_true, mut_rate = 0.05, t_size = 2, elitism_ratio = 0.05, gwcfc = 10):
    print("Geração 1")
    current_pop = create_first_gen(n_ind, evaluate_func, valid_func, gene_limits)
    current_pop = sort_ObjVal(current_pop)
    best_ind = current_pop[0]
    stillness_count = 0
    generation = 2
    while stillness_count < gwcfc:
        print(f"Geração {generation}")
        offspring = create_offspring(current_pop, gene_limits, evaluate_func, valid_func, mut_rate = mut_rate, t_size = t_size)
        current_pop = reinsert_immortal_individuals(current_pop, offspring)
        current_pop = sort_ObjVal(current_pop)
        if best_ind.get_ObjVal() > current_pop[0].get_ObjVal():
            best_ind = current_pop[0]
            stillness_count = 0
        else: 
            stillness_count += 1
        #print(best_ind.get_chrom())
        print(f'O melhor indivíduo reina há {stillness_count} gerações')
        generation += 1
    return current_pop


def run_GA_validate(n_gen, n_ind, gene_limits, evaluate_func, valid_func = Auxiliary.return_true, mut_rate = 0.05, t_size = 2, elitism_ratio = 0.05):
    print("Geração 1")
    current_pop = create_first_gen(n_ind, evaluate_func, valid_func, gene_limits)
    current_pop = sort_ObjVal(current_pop)
    #plt.ion()
    plt.plot(1, current_pop[0].get_ObjVal(), 'ro')
    #plt.draw()
    #fig = plt.figure()
    #ax = fig.add_subplot(projection='3d')   
    for generation in range(2, n_gen + 1): 
        print(f"Geração {generation}")
        #for individual in current_pop:
        #    chrom = individual.get_chrom()
        #    func = evaluate_func(chrom)
        #    ax.scatter(chrom[0], chrom[1], func)
        offspring = create_offspring(current_pop, gene_limits, evaluate_func, valid_func, mut_rate = mut_rate, t_size = t_size)
        #current_pop = reinsert_basic_elitism(current_pop, offspring, elitism_ratio = elitism_ratio)
        #plt.draw()
        current_pop = reinsert_immortal_individuals(current_pop, offspring)
        current_pop = sort_ObjVal(current_pop)
        plt.plot(generation, current_pop[0].get_ObjVal(), 'ro')
    plt.show()
    return current_pop


#main functions
def run_GA(n_gen, n_ind, gene_limits, evaluate_func, valid_func = Auxiliary.return_true, mut_rate = 0.05, t_size = 2, elitism_ratio = 0.05):
    print("Geração 1")
    current_pop = create_first_gen(n_ind, evaluate_func, valid_func, gene_limits)
    for generation in range(2, n_gen + 1): 
        print(f"Geração {generation}")
        offspring = create_offspring(current_pop, gene_limits, evaluate_func, valid_func, mut_rate = mut_rate, t_size = t_size)
        current_pop = reinsert_immortal_individuals(current_pop, offspring)
    return current_pop

def create_first_gen(n_ind, evaluate_func, valid_func, limits):
    individuals = 0 
    chrom_pop = []
    first_gen = []
    while individuals < n_ind:
        chromossome = generate_chrom(limits)
        if valid_func(chromossome) and chromossome not in chrom_pop:
            chrom_pop.append(chromossome)
            individuals += 1
    for chrom in chrom_pop:
        first_gen.append(Individual(chrom = chrom))
    evaluate_population(first_gen, evaluate_func, valid_func)
    return first_gen

def create_offspring(population, limits, evaluate_func, valid_func, mut_rate = 0.05, t_size = 2):
    chrom_pop = []
    init_chrom_pop = []
    length = len(population)
    for individual in population:
        init_chrom_pop.append(individual.get_chrom())
    population = strip_rejected(population)
    #population = strip_equal(population)
    while len(chrom_pop) < length:
        mother = select_tournament(population, tournament_size = t_size)
        father = select_tournament(population, tournament_size = t_size)
        son_chrom = crossoverFull(mother.get_chrom(), father.get_chrom())
        son_chrom = mutate(son_chrom, mut_rate, limits)
        if valid_func(son_chrom) and (son_chrom not in init_chrom_pop):
            chrom_pop.append(son_chrom)
    offspring = []
    for chrom in chrom_pop:
        offspring.append(Individual(chrom = chrom))
    evaluate_population(offspring, evaluate_func, valid_func)
    return offspring

def evaluate_population(population, evaluate_func, valid_func):
    for individual in population:
        valid = valid_func(individual.get_chrom())
        if not valid:
            individual.set_valid(False)
            continue
        obj_vals = evaluate_func(individual.get_chrom())
        individual.set_ObjVal(obj_vals)




#Reinsertion
def reinsert_immortal_individuals(pop_ini, pop_new):
    length = len(pop_ini)
    pop_final = []
    pop_new = strip_rejected(pop_new)
    pop_ini.extend(pop_new)
    while len(pop_final) < length:
        best = pop_ini[0]
        i = 1
        hold = 0
        while i < len(pop_ini):
            if best.get_ObjVal() > pop_ini[i].get_ObjVal():
                hold = i
                best = pop_ini[i]
            i += 1
        pop_final.append(best)
        pop_ini.pop(hold)
    return pop_final

def reinsert_basic_elitism(pop_ini, pop_new, elitism_ratio = 0.05):
    length = len(pop_ini)
    n_pass = math.ceil(length*elitism_ratio)
    pop_final = []
    pop_ini = sort_ObjVal(pop_ini[:])
    pop_new = sort_ObjVal(pop_new[:])
    for i in range(n_pass):
        pop_final.append(pop_ini[i])
    for i in range(length - n_pass):
        pop_final.append(pop_new[i])
    return pop_final




#Crossover
def crossover1Point(mother, father):
    crossover_point = random.randint(0, len(mother) - 1)
    son1_genes = mother[:crossover_point]
    son1_genes.extend(father[crossover_point:])
    son2_genes = father[:crossover_point]
    son2_genes.extend(mother[crossover_point:])
    return son1_genes, son2_genes

def crossover2Point(mother, father):
    crossover_point1 = random.randint(0, len(mother) - 1)
    crossover_point2 = random.randint(0, len(mother) - 1)
    if (crossover_point2 == crossover_point1) or (abs(crossover_point1 - crossover_point2) == len(mother) - 1):
        return mother, father
    if crossover_point2 < crossover_point1:
        temp = crossover_point1
        crossover_point1 = crossover_point2
        crossover_point2 = temp
    son1_genes = mother[:crossover_point1]
    son1_genes.extend(father[crossover_point1:crossover_point2])
    son1_genes.extend(mother[crossover_point2:])
    son2_genes = father[:crossover_point1]
    son2_genes.extend(mother[crossover_point1:crossover_point2])
    son2_genes.extend(father[crossover_point2:])
    return son1_genes, son2_genes

def crossoverFull(mother, father):
    son_genes = []
    length = len(mother)
    for i in range(length):
        chosen = random.random()
        if chosen > 0.5:
            son_genes.append(mother[i])
        else:
            son_genes.append(father[i])
    return son_genes




#Mutation
def mutate(chrom, mut_rate, limits):
    length = len(chrom)
    for i in range(length):
        mut = random.random()
        if mut < mut_rate:
            chrom[i] = Auxiliary.generate_random_parameter(limits[i])
    return chrom





#Selection
def select_tournament(population, tournament_size = 2):
    index_list = []
    for _ in range(tournament_size):
        rand_index = random.randint(0, len(population) - 1)
        while rand_index in index_list:
            rand_index = random.randint(0, len(population) - 1)
        index_list.append(rand_index)
    best = population[index_list[0]]
    for i in range(1, tournament_size):
        new = population[index_list[i]]
        best = battle(best, new)
    return best





#Auxiliary
def battle(indA, indB):
    if indA.get_ObjVal() < indB.get_ObjVal():
        return indA
    return indB

def strip_rejected(population):
    i = 0
    while i < len(population): 
        if not population[i].get_valid():
            population.pop(i)
            i -= 1
        i += 1
    return population

def strip_equal(population):
    res = []
    for i in population:
        if i not in res:
            res.append(i)
    return res

def sort_ObjVal(population):
    new_pop = []
    copy = population[:]
    while len(copy) > 0:
        min_obj_val = copy[0].get_ObjVal()
        i = 1
        hold = 0
        while i < len(copy):
            obj_val = copy[i].get_ObjVal()
            if obj_val < min_obj_val:
                min_obj_val = obj_val
                hold = i
            i += 1
        new_pop.append(copy[hold])
        copy.pop(hold)
    return new_pop

def generate_chrom(limits):
    chrom = []
    for limit in limits:
        chrom.append(Auxiliary.generate_random_parameter(limit))
    return chrom