from importing import *

class Individual:
    def __init__(self, chrom = [], front = -1, crowd = 0, domcount = infinite, ObjVal = [], valid = True):
        self.chrom = chrom
        self.front = front
        self.crowd = crowd
        self.domcount = domcount
        self.ObjVal = ObjVal
        self.valid = valid

    def set_valid(self, boolean):
        self.valid = boolean
    def get_valid(self):
        return self.valid

    def set_chrom(self, chrom):
        self.chrom = chrom
    def get_chrom(self):
        return self.chrom

    def set_Front(self, front):
        self.front = front
    def get_Front(self):
        return self.front

    def set_Crowding(self, crowd):
        self.crowd = crowd
    def get_Crowding(self):
        return self.crowd

    def set_Dominated_counter(self, counter):
        self.domcount = counter
    def get_Dominated_counter(self):
        return self.domcount

    def set_ObjVal(self, value, index = -1):
        if index == -1:
            self.ObjVal = value
        else:
            self.ObjVal[index] = value
    def get_ObjVal(self, index = -1):
        if index == -1:
            return self.ObjVal
        return self.ObjVal[index]

    def set_ID(self, ID):
        self.ID = ID
    def get_ID(self):
        return self.ID

class NSGA2:
    def __init__(self, n_ind = 0, n_gen = 0, mut_rate = 0, t_size = 0, a_size = 0, decimal_places = 4):
        self.n_ind = n_ind
        self.n_gen = n_gen
        self.mut_rate = mut_rate
        self.t_size = t_size
        self.a_size = a_size
        self.decimal_places = decimal_places

    def set_population(self, population, limits):
        y = {}
        for key in population:
            y[key] = population[key].copy()
        self.population = y
        self.pop_length = len(y)
        self.limits = limits[:]

    def set_functions(self, evaluate, validate):
        self.evaluate = evaluate
        self.validate = validate

    def run_v1(self):
        print("Geração 1")
        self.current_pop = self.create_first_gen()
        for generation in range(2, self.n_gen + 1): 
            print(f"Geração {generation}")
            offspring = create_offspring(current_pop, gene_limits, evaluate_func, valid_func, mut_rate = mut_rate, t_size = t_size)
            current_pop = reinsert(current_pop, offspring)
        return current_pop

    def create_first_gen(self):
        individuals = 0 
        chrom_pop = []
        first_gen = []
        while individuals < self.n_ind:
            chromossome = self.generate_chrom()
            if self.validate(chromossome) and chromossome not in chrom_pop:
                chrom_pop.append(chromossome)
                individuals += 1
        for chrom in chrom_pop:
            first_gen.append(Individual(chrom = chrom))
        self.evaluate_population(first_gen)
        first_gen = assign_fronts(first_gen)
        first_gen = assign_crowding(first_gen)
        return first_gen

    def evaluate_population(self, population):
        for individual in population:
            valid = self.validate(individual.get_chrom())
            if not valid:
                individual.set_valid(False)
                continue
            obj_vals = self.evaluate(individual.get_chrom())
            individual.set_ObjVal(obj_vals)

    def generate_chrom(self):
        chrom = []
        for limit in self.limits:
            param = Auxiliary.generate_random_parameter(limit, self.decimal_places)
            chrom.append(param)
        return chrom




















#Convergence through Pareto accumulator
def run_GA_convergence_v2(n_ind, gene_limits, evaluate_func, valid_func = Auxiliary.return_true, mut_rate = 0.0, t_size = 2, gwcfc = 10, archive_size = -1, decimal_places = 4):
    if archive_size == -1:
        archive_size == n_ind*2
    print("Geração 1")
    archive = []
    current_pop = create_first_gen(n_ind, evaluate_func, valid_func, gene_limits, decimal_places = decimal_places)
    counter = 0
    generation = 2
    while counter < gwcfc: 
        print(f"Geração {generation}")
        offspring = create_offspring(current_pop, gene_limits, evaluate_func, valid_func, mut_rate = mut_rate, t_size = t_size, decimal_places = decimal_places)
        current_pop = reinsert(current_pop, offspring)
        archive, counter = update_archive(archive, current_pop, archive_size, counter)
        print(counter)
        generation += 1
    return archive

#Convergence through area under front
def run_GA_convergence_v1(n_ind, gene_limits, evaluate_func, valid_func = Auxiliary.return_true, mut_rate = 0.0, t_size = 2, gwcfc = 10, tolerance_area = 0.01):
    print("Geração 1")
    current_pop = create_first_gen(n_ind, evaluate_func, valid_func, gene_limits)
    area = area_under_Front(current_pop)
    stillness_count = 0
    generation = 2
    while stillness_count < gwcfc:
        print(f"Geração {generation}")
        offspring = create_offspring(current_pop, gene_limits, evaluate_func, valid_func, mut_rate = mut_rate, t_size = t_size)
        current_pop = reinsert(current_pop, offspring)
        new_area = area_under_Front(current_pop)
        if abs(area - new_area) < tolerance_area:
            stillness_count += 1
        else: 
            area = new_area
            stillness_count = 0
        #for ind in current_pop:
        #    plt.plot(ind.get_ObjVal(0), ind.get_ObjVal(1), 'ro')
        #plt.show()
        print(area)
        print(f'Não há alterações nos melhores indivíduos há {stillness_count} gerações')
        generation += 1
    return current_pop



#main functions
# def run_GA(n_gen, n_ind, gene_limits, evaluate_func, valid_func = Auxiliary.return_true, mut_rate = 0.0, t_size = 2):
#     print("Geração 1")
#     current_pop = create_first_gen(n_ind, evaluate_func, valid_func, gene_limits)
#     for generation in range(2, n_gen + 1): 
#         print(f"Geração {generation}")
#         offspring = create_offspring(current_pop, gene_limits, evaluate_func, valid_func, mut_rate = mut_rate, t_size = t_size)
#         current_pop = reinsert(current_pop, offspring)
#     return current_pop

# def create_first_gen(n_ind, evaluate_func, valid_func, limits, decimal_places = 4):
#     individuals = 0 
#     chrom_pop = []
#     first_gen = []
#     while individuals < n_ind:
#         chromossome = generate_chrom(limits, decimal_places)
#         if valid_func(chromossome) and chromossome not in chrom_pop:
#             chrom_pop.append(chromossome)
#             individuals += 1
#     for chrom in chrom_pop:
#         first_gen.append(Individual(chrom = chrom))
#     evaluate_population(first_gen, evaluate_func, valid_func)
#     first_gen = assign_fronts(first_gen)
#     first_gen = assign_crowding(first_gen)
#     return first_gen

def create_first_gen_new(n_ind, evaluate_func, valid_func, limits, decimal_places = 4):
    individuals = 0 
    first_gen = []
    limits = nup.ndarray(limits)
    sampling = LHS(xlimits = limits)
    while individuals < n_ind:
        placeholder = n_ind
        chrom_pop = sampling(n_ind - individuals)
        chromossome = generate_chrom(limits, decimal_places)
        if valid_func(chromossome) and chromossome not in chrom_pop:
            chrom_pop.append(chromossome)
            individuals += 1
    for chrom in chrom_pop:
        first_gen.append(Individual(chrom = chrom))
    evaluate_population(first_gen, evaluate_func, valid_func)
    first_gen = assign_fronts(first_gen)
    first_gen = assign_crowding(first_gen)
    return first_gen
    
def create_offspring(population, limits, evaluate_func, valid_func, mut_rate = 0.0, t_size = 2, decimal_places = 4):
    init_chrom_pop = []
    chrom_pop = []
    length = len(population)
    for individual in population:
        init_chrom_pop.append(individual.get_chrom())
    population = strip_rejected(population)
    #population = strip_equal(population)
    while len(chrom_pop) < length:
        mother = select_tournament(population, tournament_size = t_size)
        father = select_tournament(population, tournament_size = t_size)
        son_chrom = crossoverFull(mother.get_chrom(), father.get_chrom())
        son_chrom = mutate(son_chrom, mut_rate, limits, decimal_places = decimal_places)
        if valid_func(son_chrom) and (son_chrom not in init_chrom_pop): #and (son_chrom not in chrom_pop):
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

def area_under_Front(population):
    pareto_front = return_pareto_front(population)
    area = calculate_hyper_volume(pareto_front)
    return area
    




#Reinsertion
def reinsert(pop_ini, pop_new):
    length = len(pop_ini)
    pop_final = []
    pop_ini.extend(pop_new)
    pop_ini = strip_rejected(pop_ini)
    pop_ini = strip_equal(pop_ini)
    pop_ini = assign_fronts(pop_ini)
    pop_ini = assign_crowding(pop_ini)
    while len(pop_final) < length and len(pop_ini) > 0:
        best = pop_ini[0]
        i = 1
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

def update_archive(archive, pop_new, archive_size, counter):
    # if len(archive) + len(pop_new) <= archive_size:
    #     archive.extend(pop_new)
    #     print("passed")
    #     return archive, 0
    init_archive = archive[:]
    # for ind in archive:
    #     ind.set_ID(1)
    # for ind in pop_new:
    #     ind.set_ID(2)
    pop_final = []
    archive.extend(pop_new)
    archive = strip_rejected(archive)
    archive = strip_equal(archive)
    archive = assign_fronts(archive)
    archive = assign_crowding(archive)
    archive = strip_not_pareto(archive)
    while len(pop_final) < archive_size and len(archive) > 0:
        best = archive[0]
        i = 1
        hold = 0
        while i < len(archive):
            if archive[i].get_Front() < best.get_Front():
                best = archive[i]
                hold = i
            elif archive[i].get_Front() == best.get_Front():
                if archive[i].get_Crowding() > best.get_Crowding():
                    best = archive[i]
                    hold = i
            i += 1
        pop_final.append(best)
        archive.pop(hold)
    if compare_pops(init_archive, pop_final):
        counter += 1
    else:
        counter = 0
    # for x in pop_final:
    #     xx = x.get_chrom()
    #     print(xx[0])
    #     print(xx[1])
    #     print('')
    # time.sleep(0.5)
    return pop_final, counter



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
def mutate(chrom, mut_rate, limits, decimal_places = 4):
    copy = chrom[:]
    length = len(chrom)
    for i in range(length):
        mut = random.random()
        if mut < mut_rate:
            param = Auxiliary.generate_random_parameter(limits[i], decimal_places)
            copy[i] = param
    return copy



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
def write_pop(population, file_name):
    pass

def calculate_hyper_volume(front):
    front = sort_ObjVal(front, 0)
    lenght = len(front)
    if lenght == 0:
        return 0
    n_objvals = len(front[0].get_ObjVal())
    #for i in range(n_objvals):
    #    minval = sort_ObjVal(front, i)[0].get_ObjVal(i)
    #    for ind in front:
    #        ind.set_ObjVal(ind.get_ObjVal(i) - minval, index = i)
    main_matrix = []
    for i in range(n_objvals):
        obj_vector = []
        for ind in front:
            obj_vector.append(ind.get_ObjVal(index = i))
        main_matrix.append(obj_vector)
    area = 0
    for i in range(lenght):
        accumulator = 1
        if i == 0:
            for j in range(n_objvals):
                accumulator *= main_matrix[j][i]
        else:
            for j in range(n_objvals - 1):
                accumulator *= (main_matrix[j][i] - main_matrix[j][i - 1])
            accumulator *= main_matrix[-1][i]
        area += accumulator
    return area

def return_pareto_front(population):
    population = assign_fronts(population)
    pareto_front = []
    min_front = infinite
    for ind in population:
        if ind.get_Front() < min_front:
            min_front = ind.get_Front()
    for ind in population:
        if ind.get_Front() == min_front:
            pareto_front.append(ind)
    return pareto_front

def assign_domination(population, n_objvals):
    if len(population) == 1:
        population[0].set_Dominated_counter(0)
        return
    elif len(population) == 0:
        return
    length = len(population)
    for i in range(length):
        dominated_counter = 0
        for j in range(length):
            if j == i:
                continue
            if calculate_domination(population[j], population[i], n_objvals):
                dominated_counter += 1
        population[i].set_Dominated_counter(dominated_counter)

def assign_fronts(population):
    n_objvals = len(population[0].get_ObjVal())
    copy = population[:]
    new_pop = []
    front = 1
    while len(copy) > 0:
        assign_domination(copy, n_objvals)
        i = 0
        while i < len(copy):
            if copy[i].get_Dominated_counter() == 0:
                copy[i].set_Front(front)
                new_pop.append(copy[i])
                copy.pop(i)
            else:
                i += 1
        front += 1
    return new_pop

def assign_crowding(population):
    n_objvals = len(population[0].get_ObjVal())
    length = len(population)
    reset_crowding(population)
    for i in range(n_objvals):
        copy = sort_ObjVal(population, i)
        copy[0].set_Crowding(infinite)
        copy[-1].set_Crowding(infinite)
        minval = copy[0].get_ObjVal(i)
        maxval = copy[-1].get_ObjVal(i)
        for j in range(1, length - 1):
            crowding = abs((copy[j + 1].get_ObjVal(i) - copy[j - 1].get_ObjVal(i))/(maxval - minval)) #abs(Auxiliary.linear_interpolate(minval, maxval, 0, 1, copy[j + 1].get_ObjVal(i)) - Auxiliary.linear_interpolate(minval, maxval, 0, 1, copy[j - 1].get_ObjVal(i)))
            copy[j].set_Crowding(copy[j].get_Crowding() + crowding)
    return copy

def reset_crowding(population):
    for ind in population:
        ind.set_Crowding(0)

def battle(indA, indB):
    if indA.get_Front() < indB.get_Front():
        return indA
    elif indA.get_Front() == indB.get_Front():
        if indA.get_Crowding() > indB.get_Crowding():
            return indA
    return indB

def sort_ObjVal(population, objval_index):
    new_pop = []
    copy = population[:]
    while len(copy) > 0:
        min_obj_val = infinite
        i = 0
        hold = 0
        while i < len(copy):
            obj_val = copy[i].get_ObjVal(objval_index)
            if obj_val < min_obj_val:
                min_obj_val = obj_val
                hold = i
            i += 1
        new_pop.append(copy[hold])
        copy.pop(hold)
    return new_pop

def calculate_domination(indA, indB, n_objvals):
    boolA = indA.get_ObjVal(0) <= indB.get_ObjVal(0)
    boolB = indA.get_ObjVal(0) < indB.get_ObjVal(0)
    for i in range(1, n_objvals):
        boolA = boolA and (indA.get_ObjVal(i) <= indB.get_ObjVal(i))
        boolB = boolB or (indA.get_ObjVal(i) < indB.get_ObjVal(i))
    return (boolA and boolB) #if True, A dominates B, if False, A does not dominate B

def strip_rejected(population):
    i = 0
    while i < len(population): 
        if not population[i].get_valid():
            population.pop(i)
            i -= 1
        i += 1
    return population

def strip_not_pareto(population):
    i = 0
    while i < len(population):
        if population[i].get_Front() != 1:
            population.pop(i)
        else:
            i += 1
    return population

def strip_equal(population):
    res = []
    chrom_vec = []
    chrom_vec.append(population[0].get_chrom())
    res.append(population[0])
    for i in range(1, len(population)):
        new_chrom = population[i].get_chrom()
        if new_chrom not in chrom_vec:
            chrom_vec.append(new_chrom)
            res.append(population[i])
    return res

def compare_pops(pop1, pop2): #returns False if they are different, True if they are equal
    pop1 = sort_ObjVal(pop1, 0)
    pop2 = sort_ObjVal(pop2, 0)
    chrom_pop1 = []
    chrom_pop2 = []
    for ind in pop1:
        chrom_pop1.append(ind.get_chrom())
    for ind in pop2:
        chrom_pop2.append(ind.get_chrom())
    for i in range(len(chrom_pop1)):
        if chrom_pop1[i] != chrom_pop2[i]:
            return False
    return True

# def generate_chrom(limits, decimal_places):
#     chrom = []
#     for limit in limits:
#         param = Auxiliary.generate_random_parameter(limit, decimal_places)
#         chrom.append(param)
#     return chrom

