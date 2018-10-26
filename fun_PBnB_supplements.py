from c_SubRegion import *

def fun_uni_sampler(c_subr, f_alpha, i_n_samp, i_n_rep, i_dim, fun_blackbox):
    import numpy as np
    import pandas as pd
    import copy

    def turn_to_power(list, power=1):
        if type(list) is list:
            return [number ** power for number in list]
        else:
            return [number ** power for number in [list]]
    c_subr.pd_sample_record = c_subr.pd_sample_record.sort_values(by="mean", ascending=True)  # sort before start
    c_subr.pd_sample_record = c_subr.pd_sample_record.reset_index(drop=True)  # reindex before start
    if len(c_subr.pd_sample_record) >= i_n_samp:  # if has enough number of sampling points
        for i in (i for i in range(0, len(c_subr.pd_sample_record)) if
                  c_subr.pd_sample_record.loc[i, '# rep'] < i_n_rep):  # check enough # reps or not, i is index
            # for j in range(int(record.loc[i,'# rep']),i_n_rep): # counting number of replications
            l_output = fun_blackbox(np.array(c_subr.pd_sample_record.loc[i, 'x1':'x' + str(i_dim)]),
                                    i_n_rep - int(c_subr.pd_sample_record.loc[i, '# rep']))
            c_subr.pd_sample_record.loc[i, 'mean'] = copy.copy(float(
                int(c_subr.pd_sample_record.loc[i, '# rep']) * c_subr.pd_sample_record.loc[i, 'mean'] + sum(
                    l_output)) / i_n_rep)
            c_subr.pd_sample_record.loc[i, 'SST'] = copy.copy(
                c_subr.pd_sample_record.loc[i, 'SST'] + sum(turn_to_power(l_output, 2)))
            c_subr.pd_sample_record.loc[i, 'var'] = copy.copy(float(
                c_subr.pd_sample_record.loc[i, 'SST'] - i_n_rep * pow(c_subr.pd_sample_record.loc[i, 'mean'], 2)) / (
                                                              i_n_rep - 1))
            c_subr.pd_sample_record.loc[i, '# rep'] = i_n_rep
    elif len(c_subr.pd_sample_record) < i_n_samp:  # if has not enough sampling points and replication
        if len(c_subr.pd_sample_record) >= 1:  # if already has sample points, first deal with them
            for i in (i for i in range(0, len(c_subr.pd_sample_record)) if c_subr.pd_sample_record.loc[
                i, '# rep'] < i_n_rep):  # check enough # reps or not for existing old sampling points
                l_output = fun_blackbox(np.array(c_subr.pd_sample_record.loc[i, 'x1':'x' + str(i_dim)]),
                                        i_n_rep - int(c_subr.pd_sample_record.loc[i, '# rep']))
                c_subr.pd_sample_record.loc[i, 'mean'] = copy.copy(float(
                    int(c_subr.pd_sample_record.loc[i, '# rep']) * c_subr.pd_sample_record.loc[i, 'mean'] + sum(
                        l_output)) / i_n_rep)
                c_subr.pd_sample_record.loc[i, 'SST'] = copy.copy(
                    c_subr.pd_sample_record.loc[i, 'SST'] + sum(turn_to_power(l_output, 2)))
                c_subr.pd_sample_record.loc[i, 'var'] = copy.copy(float(
                    c_subr.pd_sample_record.loc[i, 'SST'] - i_n_rep * pow(c_subr.pd_sample_record.loc[i, 'mean'],
                                                                          2)) / (i_n_rep - 1))
                c_subr.pd_sample_record.loc[i, '# rep'] = i_n_rep

        i_ini_length = len(c_subr.pd_sample_record)  # number of sampling point sofar in this subregion

        for i in range(i_ini_length, i_n_samp):  # create new rows for new sampling points
            c_subr.pd_sample_record.loc[i] = [1 for n in range(len(c_subr.pd_sample_record.columns))]
        for i in range(0, i_dim):  # create new sampling new point and add new sample point to dataframe
            a_new_sample = np.random.uniform(low=c_subr.l_coordinate_lower[i], high=c_subr.l_coordinate_upper[i],
                                             size=i_n_samp - i_ini_length)  # generate only for one dim
            index = [x for x in range(i_ini_length, i_n_samp)]
            c_subr.pd_sample_record.loc[i_ini_length:i_n_samp - 1, 'x' + str(i + 1)] = pd.Series(a_new_sample.tolist(),
                                                                                                 index)
            # c_subr.pd_sample_record.loc[index,'x'+str(i+1)]=pd.Series(a_new_sample.tolist(),index)
        for i in range(i_ini_length, i_n_samp):  # evalute new sampling points
            # select the input as arrat
            a_input = np.array(c_subr.pd_sample_record.loc[i, 'x1':'x' + str(i_dim)])
            l_output = fun_blackbox(a_input, i_n_rep)  # should be i_n_rep-.... and following is the same
            c_subr.pd_sample_record.loc[i, 'mean'] = np.mean(l_output)
            c_subr.pd_sample_record.loc[i, 'SST'] = sum(turn_to_power(l_output, 2))
            c_subr.pd_sample_record.loc[i, 'var'] = copy.copy(float(
                c_subr.pd_sample_record.loc[i, 'SST'] - i_n_rep * pow(c_subr.pd_sample_record.loc[i, 'mean'], 2)) / (
                                                              i_n_rep - 1))
            c_subr.pd_sample_record.loc[i, '# rep'] = i_n_rep

    # the following update i_min_sample, i_max_sample, f_min_diff_samplemean, and f_max_var
    c_subr.pd_sample_record = c_subr.pd_sample_record.sort_values(by="mean", ascending=True)
    c_subr.pd_sample_record = c_subr.pd_sample_record.reset_index(drop=True)  # reindex the sorted df
    if len(c_subr.pd_sample_record) > 0:
        c_subr.i_min_sample = c_subr.pd_sample_record.loc[0, 'mean']
        c_subr.i_max_sample = c_subr.pd_sample_record.loc[len(c_subr.pd_sample_record) - 1, 'mean']
    c_subr.f_min_diff_samplemean = min(c_subr.pd_sample_record['mean'].shift(-1) - c_subr.pd_sample_record['mean'])
    c_subr.f_max_var = max(c_subr.pd_sample_record.loc[:, 'var'])
    return c_subr


def fun_order_region(l_subr):
    import pandas as pd

    l_allsample_C = []
    for i in (i for i in l_subr if i.s_label == 'C'):
        l_allsample_C.append(i.pd_sample_record)
    pd_order_z = pd.concat(l_allsample_C)
    pd_order_z = pd_order_z.sort_values(by="mean", ascending=True)
    pd_order_z = pd_order_z.reset_index(drop=True)  # reindex the sorted df
    return pd_order_z


def fun_replication_update(l_subr, i_n_rep, f_alpha):
    import scipy.stats
    import math
    if list(i.f_min_diff_samplemean for i in l_subr if
            i.s_label == 'C' and i.b_activate == True) + [] == []:  # to prevent empty sequence
        f_d_star = 0.005
    elif min(i.f_min_diff_samplemean for i in l_subr if i.s_label == 'C' and i.b_activate == True) < 0.005:
        f_d_star = 0.005
    else:
        f_d_star = min(i.f_min_diff_samplemean for i in l_subr if i.s_label == 'C' and i.b_activate == True)
    f_var_star = max(i.f_max_var for i in l_subr if i.s_label == 'C' and i.b_activate == True)
    z = scipy.stats.norm.ppf(1 - f_alpha / 2)
    # to prevent the float NaN
    if math.isnan(z) == True or math.isnan(f_d_star) == True or math.isnan(f_var_star) == True:
        i_n_rep = i_n_rep
    else:
        i_n_rep = max(i_n_rep, 4 * int(math.ceil(pow(z, 2) * f_var_star / pow(f_d_star, 2))))
    return i_n_rep


def fun_order_region(l_subr):
    import pandas as pd

    l_allsample_C = []
    for i in (i for i in l_subr if i.s_label == 'C'):
        l_allsample_C.append(i.pd_sample_record)
    pd_order_z = pd.concat(l_allsample_C)
    pd_order_z = pd_order_z.sort_values(by="mean", ascending=True)
    pd_order_z = pd_order_z.reset_index(drop=True)  # reindex the sorted df
    return pd_order_z

def fun_pruning_indicator (c_subr, f_CI_u):
    if c_subr.i_min_sample > f_CI_u:
        c_subr.b_maintaining_indicator = False
        c_subr.b_pruning_indicator = True
    return c_subr

def fun_maintaining_indicator(c_subr, f_CI_l):
    if c_subr.i_max_sample < f_CI_l:
        c_subr.b_maintaining_indicator = True
        c_subr.b_pruning_indicator = False
    return c_subr

def fun_elite_indicator(c_subr, f_CI_l):
    if c_subr.i_max_sample < f_CI_l:
        c_subr.b_elite = True
        c_subr.b_worst = False
    return c_subr

def fun_worst_indicator (c_subr, f_CI_u):
    if c_subr.i_min_sample > f_CI_u:
        c_subr.b_elite=False
        c_subr.b_worst=True
    return c_subr

def fun_quantile_update (l_subr, f_delta):
    f_vol_C = sum(c.f_volumn for c in l_subr if c.s_label == 'C' and c.b_activate == True)
    f_vol_Pruning = sum(c.f_volumn for c in l_subr if c.b_pruning_indicator == True and c.b_activate == True)
    f_vol_Maintaining = sum(c.f_volumn for c in l_subr if c.b_maintaining_indicator == True and c.b_activate == True)
    f_delta = float(f_delta*f_vol_C-f_vol_Maintaining)/(f_vol_C-f_vol_Pruning-f_vol_Maintaining)
    return f_delta


def fun_CI_builder(l_subr, pd_order_z, f_delta_k, f_alpha_k, f_epsilon):
    from scipy.stats import binom
    import math
    f_vol_S = l_subr[0].f_volumn
    f_vol_C = sum(c.f_volumn for c in l_subr if c.s_label == 'C' and c.b_activate == True)
    f_vol_P = sum(c.f_volumn for c in l_subr if c.s_label == 'P' and c.b_activate == True)
    f_vol_M = sum(c.f_volumn for c in l_subr if c.s_label == 'M' and c.b_activate == True)
    f_delta_kl = f_delta_k - float(f_vol_P * f_epsilon) / (f_vol_S * f_vol_C)
    f_delta_ku = f_delta_k + float(f_vol_M * f_epsilon) / (f_vol_S * f_vol_C)
    f_max_r = binom.ppf(f_alpha_k / 2, len(pd_order_z), f_delta_kl)
    f_min_s = binom.ppf(1 - f_alpha_k / 2, len(pd_order_z), f_delta_ku)
    if math.isnan(f_max_r) == True:
        f_max_r = 0
    CI_l = pd_order_z.loc[f_max_r, 'mean']
    CI_u = pd_order_z.loc[f_min_s, 'mean']

    return [CI_u, CI_l]

def fun_pruning_labeler (c_subr):
    c_subr.s_label = 'P'
    return c_subr

def fun_maintaining_labeler (c_subr):
    c_subr.s_label = 'M'
    return c_subr


def plot2D(l_subr, l_ini_coordinate_lower, l_ini_coordinate_upper, i_k):
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(l_ini_coordinate_upper[0], l_ini_coordinate_upper[1])
    ax.plot(l_ini_coordinate_lower[0], l_ini_coordinate_lower[1])
    for i in (i for i in l_subr if i.b_activate == True):
        if i.s_label == 'M':
            alpha_value = 1
        elif i.s_label == 'P':
            alpha_value = 0.1
        else:  # i.s_label == 'C':
            alpha_value = 0.6
        ax.add_patch(
            patches.Rectangle(
                (i.l_coordinate_lower[0], i.l_coordinate_lower[1]),  # (x,y)
                i.l_coordinate_upper[0] - i.l_coordinate_lower[0],  # width
                i.l_coordinate_upper[1] - i.l_coordinate_lower[1],  # height
                alpha=alpha_value,
                edgecolor="black"
            )
        )
    fig.savefig('PBnB_plot_iteration ' + str(i_k) + '.png', dpi=200, bbox_inches='tight')

def fun_reg_branching(c_subr, i_dim, i_B):
    import numpy as np
    import copy
    # the following decides which diminesion should be divide
    my_list = (np.array(c_subr.l_coordinate_upper)-np.array(c_subr.l_coordinate_lower)).tolist()
    f_max_value = max(my_list)
    i_max_index = my_list.index(f_max_value)
    s_branching_dim = 'x'+str(i_max_index+1) # index
    l_subr = []
    # the following creates B subregions in the list of subregions
    for i in range(0,i_B):
        l_coordinate_lower = copy.deepcopy(c_subr.l_coordinate_lower)
        l_coordinate_upper = copy.deepcopy(c_subr.l_coordinate_upper)
        l_coordinate_lower[i_max_index] = float((c_subr.l_coordinate_upper[i_max_index]- c_subr.l_coordinate_lower[i_max_index])*i)/i_B+c_subr.l_coordinate_lower[i_max_index]
        l_coordinate_upper[i_max_index] = float((c_subr.l_coordinate_upper[i_max_index]- c_subr.l_coordinate_lower[i_max_index])*(i+1))/i_B+c_subr.l_coordinate_lower[i_max_index]
        l_new_branching_subr = c_SubRegion(i_dim, l_coordinate_lower, l_coordinate_upper)
        l_subr.append(l_new_branching_subr)
    # the following reallocate sthe sampling points
    for i in l_subr:
        i.pd_sample_record = c_subr.pd_sample_record[(c_subr.pd_sample_record[s_branching_dim] > i.l_coordinate_lower[i_max_index]) & (c_subr.pd_sample_record[s_branching_dim] < i.l_coordinate_upper[i_max_index]) ]
    for i in l_subr: # reindex the sampling points into 0 1 2...
        i.pd_sample_record = i.pd_sample_record.reset_index(drop=True)
        # update attributed based on datas
        if len(i.pd_sample_record) > 0:
            i.i_min_sample = min(i.pd_sample_record.loc[:,'mean'])
            i.i_max_sample = max(i.pd_sample_record.loc[:,'mean'])
            i.f_min_diff_samplemean = min(i.pd_sample_record['mean'].shift(-1) - i.pd_sample_record['mean'])
        if len(i.pd_sample_record) > 1:
            i.f_max_var = max(i.pd_sample_record.loc[:,'var'])
    return l_subr

def fun_order_subregion(c_subr):
    c_subr.pd_sample_record=c_subr.pd_sample_record.sort_values (by="mean", ascending=True)
    c_subr.pd_sample_record=c_subr.pd_sample_record.reset_index(drop=True) # reindex the sorted df
    if len(c_subr.pd_sample_record)>0:
        c_subr.i_min_sample = c_subr.pd_sample_record.loc[0,'mean']
        c_subr.i_max_sample = c_subr.pd_sample_record.loc[len(c_subr.pd_sample_record)-1,'mean']
    c_subr.f_min_diff_samplemean = min(c_subr.pd_sample_record['mean'].shift(-1) - c_subr.pd_sample_record['mean'])
    c_subr.f_max_var = max(c_subr.pd_sample_record.loc[:,'var'])
    return c_subr