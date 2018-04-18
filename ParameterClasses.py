from enum import Enum
import numpy as np
import scipy.stats as stat
import math as math
import InputData as Data
import scr.MarkovClasses as MarkovCls
import scr.RandomVariantGenerators as Random
import scr.ProbDistParEst as Est

class HealthStats(Enum):
    """ health states of patients with HIV """
    WELL = 0
    STROKE = 1
    POST_STROKE = 2
    STROKE_DEATH = 3

class Therapies(Enum):
    """ mono vs. combination therapy """
    NONE = 0
    TREAT = 1

class _Parameters:

    def __init__(self, therapy):
        self._therapy = therapy  # selected therapy
        self._delta_t = Data.DELTA_T # simulation time step
        self._initialHealthState = HealthStats.WELL  # initial health state
        self._prob_matrix = [] # transition probability matrix of the selected therapy

    def get_initial_health_state(self):
        return self._initialHealthState

    def get_delta_t(self):
        return self._delta_t
er
    def get_transition_prob(self, state):
        return self._prob_matrix[state.value]

class ParametersFixed(_Parameters):
    def __init__(self, therapy):

        # initialize the base class
        _Parameters.__init__(self, therapy)

        # transition probabilities
        self._prob_matrix = Data.PROB_MATRIX

        # update the transition probability matrix if treatment
        if self._therapy == Therapies.TREAT:
            self._prob_matrix = Data.TREAT_PROB_MATRIX
