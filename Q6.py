
import ParameterClasses as P
import MarkovModelClasses as MarkovCls
import SupportMarkovModel as SupportMarkov

# simulating mono therapy
# create a cohort
cohort_none = MarkovCls.Cohort(
    id=0,
    therapy=P.Therapies.NONE)
# simulate the cohort
simOutputs_none = cohort_none.simulate()

# simulating combination therapy
# create a cohort
cohort_treat = MarkovCls.Cohort(
    id=0,
    therapy=P.Therapies.TREAT)
# simulate the cohort
simOutputs_treat = cohort_treat.simulate()

# print the estimates for the mean survival time and mean time to AIDS
SupportMarkov.print_outcomes(simOutputs_none, "No Therapy:")
SupportMarkov.print_outcomes(simOutputs_treat, "Anticoagulant Therapy:")

# print comparative outcomes
SupportMarkov.print_comparative_outcomes(simOutputs_none, simOutputs_treat)
