# encoding: utf-8
import optimizer as opt
from plot_results import plot_fit
from extrapolate_dynamic_variables import extrapolate_var

optm = opt.Optimizer(False)
states = ['AC', 'RO', 'AM', 'RR', 'PA', 'AP', 'MA', 'PI', 'CE', 'RN', 'PB', 'PE', 'AL', 'SE']
array_num_points = [150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150]
sim_input = zip(states, array_num_points)
intervention_levels = [10, 9, 8]
num_intervation_levels = 3

n_iteration = 100
n_iteration_ref = 100
pop_gen = 30
pop_gen_ref = 30

for input_param in sim_input:
    ret = optm.run(input_param[0], num_intervation_levels, intervention_levels, input_param[1], n_iteration, pop_gen)
    ret_ref = optm.run_tune(input_param[0], ret, num_intervation_levels, intervention_levels, input_param[1],
                            n_iteration_ref, pop_gen_ref)
    ret_ref_ifr = optm.run_tune_ifr(input_param[0], ret, ret_ref.x[0], num_intervation_levels, intervention_levels,
                                    input_param[1], n_iteration_ref, pop_gen_ref)
    plot_fit(input_param[0], input_param[1])
    extrapolate_var(input_param[0], input_param[1])
    print('-'*60)
