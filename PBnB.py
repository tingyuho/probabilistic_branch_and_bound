# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 22:00:03 2017

@author: TingYu Ho

* To execute the file, type python PBnB.py
* All the parameter setting can be modified in m_intial_paramters_setting
* fun_PBnB_supplements contains all the functions used in PBnB
"""

import m_intial_paramters_setting as par
from fun_PBnB_supplements import *
from testfunction_rosenbrock import *
from testfunction_shifted_sin import *
import numpy as np
import sys

def PBnB (fun_blackbox): # PBnB(rosenbrock)

# step 0 intialization        

    i_dim = par.i_dim
    i_B = par.i_B
    f_delta_k = par.f_delta # <-- quantile setting
    f_alpha_k = float(par.f_alpha)/i_B
    f_epsilon_k = float(par.f_epsilon)/i_B
                
    i_R_k = par.i_replication
    i_k = 1
    i_k_c = par.i_k_b
    i_c_k = par.i_c
    i_stopping_maxK = par.i_stopping_maxK
    
    l_subr = [] # <-- list of subregion objects    
    l_subr.append(c_SubRegion(i_dim, par.l_coordinate_lower, par.l_coordinate_upper)) # <-- SIGMA_1={S}
    
# step 1: Sample total c_k points in current subregions with updated replication number

    while i_k <= i_stopping_maxK:
        print ('i_k:'+str(i_k))
        plot2D (l_subr, l_subr[0].l_coordinate_lower, l_subr[0].l_coordinate_upper,i_k)
        i_N_k = int(i_c_k/sum(1 for j in l_subr if j.s_label == 'C' and j.b_activate == True))
        for i in (i for i in l_subr if i.s_label == 'C' and i.b_activate == True):
            i = fun_uni_sampler(i, f_alpha_k, i_N_k, i_R_k, i_dim, fun_blackbox) # <-- c_subr, f_alpha, i_n_samp, i_n_rep, i_dim, fun_blackbox

# step 2: Order samples in each subregion and over all subregions byt estimated function values:
      
        for i in (i for i in l_subr if i.s_label == 'C' and i.b_activate == True and len(i.pd_sample_record)>0): 
            i = fun_order_subregion(i)
        
        #resampling
        i_R_k = fun_replication_update(l_subr, i_R_k, f_alpha_k)
        print ('i_R_k: '+str(i_R_k))
        for i in (i for i in l_subr if i.s_label == 'C' and i.b_activate == True):
            i = fun_uni_sampler(i, f_alpha_k, i_N_k, i_R_k, i_dim, fun_blackbox)
        
        #reorder
        for i in (i for i in l_subr if i.s_label == 'C' and i.b_activate == True and len(i.pd_sample_record)>0): 
            i = fun_order_subregion(i)

# step 3: Build CI of quantile

        pd_order_z = fun_order_region(l_subr)
        [f_CI_u,f_CI_l] = fun_CI_builder(l_subr, pd_order_z, f_delta_k, f_alpha_k, par.f_epsilon)
        print ('[f_CI_u,f_CI_l]: '+str([f_CI_u,f_CI_l]))
        
# step 4: Find the elite and worst subregions, and further sample with updated replications

        while (i_k_c < par.i_k_b or i_k == 1) and i_k<=i_stopping_maxK: # (3)
            print ('i_k:'+str(i_k))
            print ('i_k_c:'+str(i_k_c))
            plot2D(l_subr, l_subr[0].l_coordinate_lower, l_subr[0].l_coordinate_upper, i_k)
            for i in (i for i in l_subr if i.s_label == 'C' and i.b_activate == True):
                fun_elite_indicator(i,f_CI_l)
                fun_worst_indicator(i,f_CI_u)
                
            i_N_elite_worst = int(np.ceil(np.log(f_alpha_k)/np.log(1.-(par.f_epsilon/l_subr[0].f_volumn)))) # <-- number of sampling points for elite and worst subregions
            print ('i_N_elite_worst: '+str(i_N_elite_worst))
            for i in (i for i in l_subr if i.s_label == 'C' and i.b_activate == True and (i.b_worst == True or i.b_elite == True)): # <-- with label C and whose elite or worst==1
                i = fun_uni_sampler(i, f_alpha_k, i_N_elite_worst, i_R_k, i_dim, fun_blackbox)
            # reorder
            for i in (i for i in l_subr if i.s_label == 'C' and len(i.pd_sample_record)>0): 
                i = fun_order_subregion(i)
            
            # perform R_k^n-R_k
            i_R_elite_worst = fun_replication_update (l_subr, i_R_k, f_alpha_k) # <-- number of replication for all sampling points in elite and worst regions
            print ('i_R_elite_worst: '+str(i_R_elite_worst))
            for i in (i for i in l_subr if i.s_label == 'C' and i.b_activate == True and (i.b_worst == True or i.b_elite == True)): 
                i = fun_uni_sampler(i, f_alpha_k, i_N_elite_worst, i_R_elite_worst, i_dim, fun_blackbox)
                
# step 5: Maintain, prune, and branch

            for i in (i for i in l_subr if i.s_label == 'C' and i.b_activate == True and i.b_worst == True): # <-- whose  worst == 1    
                fun_pruning_indicator (i,f_CI_u)
            for i in(i for i in l_subr if i.s_label == 'C' and i.b_activate == True and i.b_elite == True): # <-- whose elite == 1    
                fun_maintaining_indicator (i,f_CI_l)
            # need to update the maintained set and the pruned set and reset all the indicator and label         
            f_delta_k = fun_quantile_update (l_subr, f_delta_k)
            
            for i in (i for i in l_subr if i.b_pruning_indicator == True and i.b_activate == True): # <-- whose  worst == 1    
                fun_pruning_labeler (i)
            for i in(i for i in l_subr if i.b_maintaining_indicator == True  and i.b_activate == True): # <-- whose elite == 1    
                fun_maintaining_labeler (i)
          
            # determine the next step
            b_stopping_branchable = all(i.b_branchable == False for i in (i for i in l_subr if i.b_activate == True)) or all(i.s_label != 'C' for i in (i for i in l_subr if i.b_activate == True)) # any subregion is branchable
            if b_stopping_branchable == True: # (1)
                for i in (i for i in l_subr if i.b_activate == True):
                    print ('l_coordinate_lower: '+str(i.l_coordinate_lower))
                    print ('l_coordinate_upper: '+str(i.l_coordinate_upper))
                    print ('activate: '+str(i.b_activate))
                    print ('label: '+str(i.s_label))
                    print ('worst: '+str(i.b_worst))
                    print ('elite: '+str(i.b_elite))
                print ('no branchable subregions')
                sys.exit() 
            else:    # (2) and (3)
                # beginning of branching-------------------------------
                l_temp_new_branching_subr = []
                for i in (i for i in l_subr if i.s_label == 'C' and i.b_activate == True and i.b_branchable == True):
                     i.b_activate = False # deacctivate the original subregion
                     l_temp_new_branching_subr += fun_reg_branching(i, i_dim, i_B)
                # begin: print the data sofar==================
                l_subr += l_temp_new_branching_subr# append the branching subregions to the list
                

                # end of branching-------------------------------
                if all(i.b_maintaining_indicator == False for i in (i for i in l_subr if i.b_activate == True)) and all(i.b_pruning_indicator == False for i in (i for i in l_subr if i.b_activate == True)):
                    i_k_c += 1
                else:
                    i_k_c = 0
                # reset pruning and maintaining updators
                for i in l_subr: 
                    i.b_maintaining_indicator = False
                    i.b_pruning_indicator = False
                # update others paramenters
                f_alpha_k = float(par.f_alpha)/pow(i_B,i_k)
                f_epsilon_k =  float(par.f_epsilon)/pow(i_B,i_k)
                i_c_k = par.i_c*i_k
                i_k +=1
                
                if i_k_c >= par.i_k_b:  # (2)
                    i_k_c = 1
                    break
         # the following save l_subr
                    #with open("test.dat", "wb") as f:
                    # pickle.dump(l_subr, f)
        # with open("test.dat", "rb") as f:
        #print pickle.load(f)
    print ('[f_CI_u,f_CI_l]: '+str([f_CI_u,f_CI_l]))
    for i in (i for i in l_subr if i.b_activate == True and i.s_label == 'P'):
        print ('l_coordinate_lower: '+str(i.l_coordinate_lower))
        print ('l_coordinate_upper: '+str(i.l_coordinate_upper))
        print ('activate: '+str(i.b_activate))
        print ('label: '+str(i.s_label))
        print ('worst: '+str(i.b_worst))
        print ('elite: '+str(i.b_elite))
        print ('i_min_sample: '+str(i.i_min_sample))
        print ('i_max_sample: '+str(i.i_max_sample))
        print ('f_min_diff_samplemean: '+str(i.f_min_diff_samplemean))
        print ('f_max_var: '+str(i.f_max_var))
        #print ('pd_sample_record: ')
        #print (i.pd_sample_record)
        print ('')
    for i in (i for i in l_subr if i.b_activate == True and i.s_label == 'M'):
        print ('l_coordinate_lower: '+str(i.l_coordinate_lower))
        print ('l_coordinate_upper: '+str(i.l_coordinate_upper))
        print ('activate: '+str(i.b_activate))
        print ('label: '+str(i.s_label))
        print ('worst: '+str(i.b_worst))
        print ('elite: '+str(i.b_elite))
        print ('i_min_sample: '+str(i.i_min_sample))
        print ('i_max_sample: '+str(i.i_max_sample))
        print ('f_min_diff_samplemean: '+str(i.f_min_diff_samplemean))
        print ('f_max_var: '+str(i.f_max_var))
    print ('[f_CI_u,f_CI_l]: '+str([f_CI_u,f_CI_l]))
    print ('reach the maximum number of iteration')
    plot2D (l_subr, l_subr[0].l_coordinate_lower, l_subr[0].l_coordinate_upper, 'final')
    sys.exit() 
    
class c_SubRegion(object):
    def __init__(self,dim,coordinate_lower,coordinate_upper):
        import pandas as pd #this is how I usually import pandas
        import numpy as np
        self.s_label = 'C' # C: undetermined, P:prune, M:maintain
        self.b_activate = True # become false if branching into two
        self.b_branchable = True
        self.b_elite = False
        self.b_worst = False
        self.b_maintaining_indicator = False # prepared to maintain in this time
        self.b_pruning_indicator = False # prepared to prune in this time
        self.l_sample = [] #list of sampling points in this subregion
        self.l_coordinate_lower = coordinate_lower # lowerbound of this region
        self.l_coordinate_upper = coordinate_upper # upperbound of this region
        self.f_volumn = np.prod(np.array(self.l_coordinate_upper)-np.array(self.l_coordinate_lower))    # volumn of the subregion
        self.l_coordinate_sample = []
        for i in range(0,dim): #sample coordinate
            self.l_coordinate_sample.append("x"+str(i+1))
        self.i_min_sample = 0 #minimum value of sampling points in this subregions
        self.i_max_sample = 0 #maximum value of sampling points in this subregions
        self.f_min_diff_samplemean = 0. #minimum value of difference of sorted sampling points in this subregions
        self.f_max_var = 0 #maximum value of variance of sampling points in this subregions
        self.l_sample_dataset = []
        self.pd_sample_record = pd.DataFrame([], columns=[])
        for i in range(0,dim): # create the dataset record the sampling data corordinate in this subregion
            self.pd_sample_record[self.l_coordinate_sample[i]] = pd.Series(float, index=self.pd_sample_record.index)
        self.pd_sample_record['# rep'] = pd.Series(0, index=self.pd_sample_record.index) # number of replication
        self.pd_sample_record['# rep'].astype(int)
        self.pd_sample_record['mean'] = pd.Series(0, index=self.pd_sample_record.index)
        self.pd_sample_record['var'] = pd.Series(0, index=self.pd_sample_record.index)
        self.pd_sample_record['SST'] = pd.Series(0, index=self.pd_sample_record.index)
        #self.pd_sample_record['iteration'] = pd.Series(0, index=self.pd_sample_record.index)

if __name__ == "__main__":
    fun_blackbox = testfunction_rosenbrock
    PBnB(fun_blackbox)