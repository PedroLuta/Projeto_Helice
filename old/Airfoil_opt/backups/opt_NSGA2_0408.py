from importing import *
import Auxiliary

class Individual:
    def __init__(self, chrom = [], front = -1, crowd = 0, domcount = 9999, ObjVal = []):
        self.chrom = chrom
        self.front = front
        self.crowd = crowd
        self.domcount = domcount
        self.ObjVal = ObjVal

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

#First Generation and Chromossome Creation
def create_1stGenChroms(num, limits):
    population = []
    while len(population) < num:
        chrom = generate_chrom(limits)
        population.append(chrom)    
    return population

def generate_chrom(limits):
    chrom = []
    for limit in limits:
        param = Auxiliary.generate_random_parameter(limit)
        chrom.append(param)
    return chrom



#Subsequent Generations
def create_nextGenChroms(population, limits, n_objvals, mut_rate = -1):
    population = assign_fronts(population, n_objvals)
    population = assign_crowding(population, n_objvals)
    if mut_rate == -1:
        mut_rate == 1/len(population[0].get_chrom())
    next_gen = []
    while len(next_gen) < len(population):
        mother = select_tournament(population, tournament_size = 2)
        father = select_tournament(population, tournament_size = 2)
        son_chrom = crossoverFull(mother.get_chrom(), father.get_chrom())
        son_chrom = mutate(son_chrom, mut_rate, limits)
        next_gen.append(son_chrom)
    return next_gen

def create_nextGen(population, limits, n_objvals, mut_rate = -1):
    chrom_pop = create_nextGenChroms(population, limits, n_objvals, mut_rate = mut_rate)
    next_gen = []
    for chrom in chrom_pop:
        temp = Individual(chrom = chrom)
        next_gen.append(temp)
    return next_gen



#Reinsertion
def reinsert(pop_ini, pop_new, n_objvals):
    length = len(pop_ini)
    pop_final = []
    pop_ini.extend(pop_new)
    pop_ini = assign_fronts(pop_ini, n_objvals)
    pop_ini = assign_crowding(pop_ini, n_objvals)
    while len(pop_final) < length:
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
    i = 0
    while i < len(mother):
        chosen = random.random()
        if chosen > 0.5:
            son_genes.append(mother[i])
        else:
            son_genes.append(father[i])
        i += 1
    return son_genes



#Mutation
def mutate(chrom, mut_rate, limits):
    i = 0
    while i < len(chrom):
        mut = random.random()
        if mut < mut_rate:
            chrom[i] = Auxiliary.generate_random_parameter(limits[i])
        i += 1
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
    i = 1
    while i < len(index_list):
        new = population[index_list[i]]
        best = battle(best, new)
        i += 1
    return best



#Auxiliary
def assign_domination(population, n_objvals):
    if len(population) == 1:
        population[0].set_Dominated_counter(0)
        return
    elif len(population) == 0:
        return
    i = 0
    while i < len(population):
        dominated_counter = 0
        j = 0
        while j < len(population):
            if j == i:
                j += 1
                continue
            if calculate_domination(population[j], population[i], n_objvals):
                dominated_counter += 1
            j += 1
        population[i].set_Dominated_counter(dominated_counter)
        i += 1

def assign_fronts(population, n_objvals):
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

def assign_crowding(population, n_objvals):
    i = 0
    while i < n_objvals: #for i in range(n_objvals):
        copy = sort_ObjVal(population, i)
        copy[0].set_Crowding(infinite)
        copy[-1].set_Crowding(infinite)
        minval = copy[0].get_ObjVal(i)
        maxval = copy[-1].get_ObjVal(i)
        j = 1
        while j < len(copy) - 1:
            crowding = abs(Auxiliary.linear_interpolate(minval, maxval, 0, 1, copy[j + 1].get_ObjVal(i)) - Auxiliary.linear_interpolate(minval, maxval, 0, 1, copy[j - 1].get_ObjVal(i)))
            copy[i].set_Crowding(copy[i].get_Crowding() + crowding)
            j += 1
        i += 1
        return copy

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
