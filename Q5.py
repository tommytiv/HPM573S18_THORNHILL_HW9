import ParameterClasses as P
import MarkovModelClasses as MarkovCls
import SupportMarkovModel as SupportMarkov

# create a cohort
cohort = MarkovCls.Cohort(
    id=0,
    therapy=P.Therapies.TREAT)

# simulate the cohort
simOutputs = cohort.simulate()

# print the outcomes of this simulated cohort
SupportMarkov.print_outcomes(simOutputs, 'AntiCoagulant:')
