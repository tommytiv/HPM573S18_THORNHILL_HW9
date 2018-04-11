
# simulation settings
POP_SIZE = 2000     # cohort population size
SIM_LENGTH = 50    # length of simulation (years)
ALPHA = 0.05        # significance level for calculating confidence intervals
DISCOUNT = 0.03     # annual discount rate

DELTA_T = 1       # years

PSA_ON = False      # if probabilistic sensitivity analysis is on

PROB_MATRIX = [
    [0.75, 0.15, 0.00, 0.10],   # Well
    [0.00, 0.00, 1.00, 0.00],   # Stroke
    [0.00, 0.25, 0.55, 0.20],   # Post-Stroke
    [0.00, 0.00, 0.00, 1.00],   # Stroke-Death
    ]

# annual cost of each health state
ANNUAL_STATE_COST = [
    0,          # Well
    5000.00,    # Stroke
    0           # Post-Stroke
    ]

# treatment relative risk
TREATMENT_RR = 0.509
TREATMENT_RR_CI = 0.365, 0.71  # lower 95% CI, upper 95% CI

# annual probability of background mortality (number per year per 1,000 population)
ANNUAL_PROB_BACKGROUND_MORT = 8.15 / 1000

