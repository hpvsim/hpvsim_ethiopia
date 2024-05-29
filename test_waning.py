'''
Run simulations with/without waning immunity
'''

#%% Imports and settings
import sciris as sc
import hpvsim as hpv
import pylab as pl

do_plot = 1
do_run = 1

#%% Define the tests

def test_waning():
    sc.heading('Sim with and without waning')

    # Settings
    # Create and run the simulation
    pars = {
        'rand_seed': 1,
        'n_agents': 20e3,
        'start': 1970,
        'genotypes': [16, 18, 'hi5', 'ohr'],
        'burnin': 30,
        'end': 2060,
        'ms_agent_ratio': 100,
        'location': 'india',
        'verbose': 0.1
    }

    pars0 = sc.mergedicts(pars)
    pars1 = sc.mergedicts(pars, dict(use_waning=True, imm_decay={'form':'exp_decay', 'init_val':1, 'half_life':5}))

    routine_vx = hpv.routine_vx(product='bivalent', age_range=[9, 10], prob=0.9, start_year=2020)

    # Run the simulations and pull out the results
    s0 = hpv.Sim(pars0, interventions=routine_vx, label='No waning').run()
    s1 = hpv.Sim(pars1, interventions=routine_vx, label='Exponential waning').run()

    # Plot
    bi = 0  # indec for burnin
    fig, ax = pl.subplots(1, 1, figsize=(5, 4))
    res0 = s0.results
    res1 = s1.results
    what = 'hpv_incidence'
    pl.plot(res0['year'][bi:], res0[what][bi:], label='No waning')
    pl.plot(res0['year'][bi:], res1[what][bi:], color='r', label='Exponential waning')
    pl.legend()
    pl.title(what)
    fig.tight_layout()
    pl.show()

    return s0, s1


#%% Run as a script
if __name__ == '__main__':

    T = sc.tic()

    s0, s1 = test_waning()

    sc.toc(T)
    print('Done.')