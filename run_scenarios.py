'''
Run HPVsim scenarios
'''


# %% General settings

import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

# Standard imports
import numpy as np
import sciris as sc
import hpvsim as hpv

# Imports from this repository
import run_sims as rs


# What to run
debug = 0
n_seeds = [10, 1][debug]  # How many seeds to run per cluster


# %% Functions
def make_st(product='via', screen_coverage=0.15, treat_coverage=0.7, start_year=2018):
    """ Make screening & treatment intervention """

    age_range = [30, 50]
    len_age_range = (age_range[1]-age_range[0])/2
    model_annual_screen_prob = 1 - (1 - screen_coverage)**(1/len_age_range)

    screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | \
                                  (sim.t > (sim.people.date_screened + 5 / sim['dt']))
    screening = hpv.routine_screening(
        prob=model_annual_screen_prob,
        eligibility=screen_eligible,
        start_year=start_year,
        product=product,
        age_range=age_range,
        label='screening'
    )

    # Assign treatment
    screen_positive = lambda sim: sim.get_intervention('screening').outcomes['positive']
    assign_treatment = hpv.routine_triage(
        start_year=start_year,
        prob=1.0,
        annual_prob=False,
        product='tx_assigner',
        eligibility=screen_positive,
        label='tx assigner'
    )

    ablation_eligible = lambda sim: sim.get_intervention('tx assigner').outcomes['ablation']
    ablation = hpv.treat_num(
        prob=treat_coverage,
        annual_prob=False,
        product='ablation',
        eligibility=ablation_eligible,
        label='ablation'
    )

    excision_eligible = lambda sim: list(set(sim.get_intervention('tx assigner').outcomes['excision'].tolist() +
                                             sim.get_intervention('ablation').outcomes['unsuccessful'].tolist()))
    excision = hpv.treat_num(
        prob=treat_coverage,
        annual_prob=False,
        product='excision',
        eligibility=excision_eligible,
        label='excision'
    )

    radiation_eligible = lambda sim: sim.get_intervention('tx assigner').outcomes['radiation']
    radiation = hpv.treat_num(
        prob=treat_coverage/4,  # assume an additional dropoff in CaTx coverage
        annual_prob=False,
        product=hpv.radiation(),
        eligibility=radiation_eligible,
        label='radiation'
    )

    st_intvs = [screening, assign_treatment, ablation, excision, radiation]

    return st_intvs


def make_st_scenarios(products=None):
    if products is None: products = ['hpv', 'via']
    st_scenarios = dict()
    st_scenarios['No S&T'] = []
    for product in products:
        st_label = f'{product} screening'
        st_intvs = make_st(product=product)
        st_scenarios[st_label] = st_intvs
    return st_scenarios


def make_vx_scenarios(coverage_arr, product='bivalent', start_year=2018):

    routine_age = (13, 14)
    catchup_age = (14, 18)
    eligibility = lambda sim: (sim.people.doses == 0)

    vx_scenarios = dict()

    # Baseline
    vx_scenarios['No vx'] = []

    for cov_val in coverage_arr:
        label = f'{cov_val} vx coverage'

        routine_vx = hpv.routine_vx(
            prob=cov_val,
            start_year=start_year,
            product=product,
            age_range=routine_age,
            eligibility=eligibility,
            label='Routine vx'
        )

        catchup_vx = hpv.campaign_vx(
            prob=cov_val,
            years=start_year,
            product=product,
            age_range=catchup_age,
            eligibility=eligibility,
            label='Catchup vx'
        )
        vx_scenarios[label] = [routine_vx, catchup_vx]

    return vx_scenarios


def make_sims(calib_pars=None, vx_scenarios=None, st_scenarios=None):
    """ Set up all scenarios to run in parallel """

    all_msims = sc.autolist()
    for vx_name, vx_intv in vx_scenarios.items():
        for st_name, st_intv in st_scenarios.items():
            sims = sc.autolist()
            for seed in range(n_seeds):
                interventions = vx_intv + st_intv
                sim = rs.make_sim(calib_pars=calib_pars, debug=debug, interventions=interventions, end=2100, seed=seed)
                sim.label = f'{vx_name}, {st_name}'
                sims += sim
            all_msims += hpv.MultiSim(sims)

    msim = hpv.MultiSim.merge(all_msims, base=False)

    return msim


def run_scens(location=None, calib_pars=None, vx_scenarios=None, st_scenarios=None, verbose=0.2):
    """ Run the simulations """
    """ Run the simulations """
    msim = make_sims(location=location, calib_pars=calib_pars, vx_scenarios=vx_scenarios, st_scenarios=st_scenarios)
    # for sim in msim.sims:
    #     sim.run(verbose=verbose)
    msim.run(verbose=verbose)
    return msim


# %% Run as a script
if __name__ == '__main__':

    T = sc.timer()
    do_run = True
    do_save = False
    do_process = True

    # Define the specific scenarios, which are created above. This is where you can change the scenarios
    vx_scenarios = make_vx_scenarios(coverage_arr=[0.9], product='bivalent', start_year=2018)
    st_scenarios = make_st_scenarios(products=['hpv'])
    scen_labels = [vxl+', '+stl for vxl,stl in zip(vx_scenarios.keys(), st_scenarios.keys())]

    # Run scenarios (usually on VMs, runs n_seeds in parallel over M scenarios)
    if do_run:
        calib_pars = sc.loadobj('results/ethiopia_pars.obj')
        msim = run_scens(calib_pars=calib_pars, vx_scenarios=vx_scenarios, st_scenarios=st_scenarios)

        if do_save: msim.save(f'results/ethiopia_scenarios.msim')

        if do_process:

            metrics = ['year', 'asr_cancer_incidence', 'cancers', 'cancer_deaths']

            # Process results
            mlist = msim.split(chunks=len(scen_labels))
            msim_dict = sc.objdict()
            for si, scen_label in enumerate(scen_labels):
                reduced_sim = mlist[si].reduce(output=True)
                mres = sc.objdict({metric: reduced_sim.results[metric] for metric in metrics})
                msim_dict[scen_label] = mres

            sc.saveobj(f'results/ethiopia_scenarios.obj', msim_dict)

    print('Done.')
