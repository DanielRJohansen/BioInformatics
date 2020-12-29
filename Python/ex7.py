import random
import matplotlib.pyplot as plt
import numpy as np
size = 100

class specimen():
    def __init__(self, g_type):
        self.genotype = g_type
        self.alleles = alleleFromGenotype(self.genotype)


def alleleFromGenotype(gt):
    if gt == 'P':
        return 'AA'
    elif gt == 'Q':
        return 'Aa'
    if gt == 'R':
        return 'aa'
    else:
        return False

def gtFromAllele(allele):
    if allele == 'AA':
        return 'P'
    elif allele == 'Aa':
        return 'Q'
    else:
        return 'R'

def gtToIndex(gt):
    if gt == 'P':
        return 0
    elif gt == 'Q':
        return 1
    if gt == 'R':
        return 2

def alToIndex(al):
    if al == 'A':
        return 0
    else:
        return 1

class ScoreBoard:
    def __init__(self):
        self.genotype_generatios = []
        self.allele_generations = []

    def addGeneration(self, population):
        gen = [0,0,0]
        al = [0,0]
        for i in range(len(population)):
            gen[gtToIndex(population[i].genotype)] += 1
            al[alToIndex(population[i].alleles[0])] += 1
            al[alToIndex(population[i].alleles[1])] += 1

        self.genotype_generatios.append(gen)
        self.allele_generations.append(al)

    def displayGraph(self):
        num_gens = len(self.genotype_generatios)
        g_generations = np.asarray(self.genotype_generatios)
        a_generations = np.asarray(self.allele_generations)
        genotypes = ['P', 'Q', 'R']
        alleles = ['A', 'a']
        plt.subplot(121)
        for i in range(3):
            plt.plot(range(num_gens), g_generations[:,i], label=genotypes[i])
        plt.legend(bbox_to_anchor=(.85, 0.95), loc='upper left', borderaxespad=0.)
        plt.title("Genotypes")

        plt.subplot(122)
        for i in range(2):
            plt.plot(range(num_gens), a_generations[:, i], label=alleles[i])
        plt.legend(bbox_to_anchor=(.85, 0.95), loc='upper left', borderaxespad=0.)
        plt.title("Alleles")

        plt.show()

class Environment:
    def __init__(self):
        self.population = self.initPop(0.3, 0.05, 0.65)
        self.population_head = size-1
        self.gt_hist = []
        self.scoreboard = ScoreBoard()
        self.scoreboard.addGeneration(self.population)


    def initPop(self, f1, f2, f3):
        population = []
        for i in range(size):
            if i <f1*size:
                gt = 'P'
            elif i <(f1+f2)*size:
                gt = 'Q'
            else:
                gt = 'R'

            population.append(specimen(gt))
        random.shuffle(population)
        return population

    def pick2(self, ):
        mates = random.sample(range(size), 2)
        return (self.population[mates[0]], self.population[mates[1]])

    def match2(self, mates):
        genotypes = []
        for i in range(2):
            for j in range(2):
                genotypes.append(mates[0].alleles[i] + mates[1].alleles[j])
        lucky_index = random.randrange(0,3)
        lucky_allele = genotypes[lucky_index]
        offspring = specimen(gtFromAllele(lucky_allele))
        return offspring

    def replaceOldest(self, offspring):
        self.population[self.population_head] = offspring
        self.population_head += 1
        if self.population_head == size:
            self.population_head = 0

    def newGen(self):
        mates = self.pick2()
        offspring = self.match2(mates)
        self.replaceOldest(offspring)
        self.scoreboard.addGeneration(self.population)



def ex7():
    env = Environment()
    for i in range(200):
        env.newGen()
    env.scoreboard.displayGraph()

