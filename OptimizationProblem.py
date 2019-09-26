# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 15:30:59 2018

@author: olechrs
"""

from random import uniform
import inspyred


class OptimizationProblem(object):

    def __init__(self):
        """ bla bla """
        self.popsize = 0
        self.maxgeneneration = 0
        self.goal = []
        self.lower_bound = []
        self.upper_bound = []

    def __del__(self):
        pass

    def generator(self, random, args):
        '''
        Initialize and create a random initial population with a uniform
        distribution.
        '''
        chromosome = []
        for lo, hi in zip(self.lower_bound, self.upper_bound):
            chromosome.append(uniform(lo, hi))
        return chromosome

    def evaluator_inductance(self, candidates, args):
        ''' bla bla '''
        fitness = []
        for c in candidates:
            c[:2] = [int(x) for x in c[:2]]
            fit = self.run_FEMM_inductance_prim(c, **args)
            values = []
            for f, g in zip(fit, self.goal):
                if g is None:
                    values.append(f)
                else:
                    values.append(abs(g-f))
            fitness.append(inspyred.ec.emo.Pareto(values))
        return fitness
