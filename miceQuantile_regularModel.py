""" DDM fit of mice response data using Chi-square quantile optimization procedure """

# For more details about the fitting procedure and Python package, see:

# [1] Ratcliff, R., and Tuerlinckx, F. (2002). Estimating parameters of the 
# diffusion model: Approaches to dealing with contaminant reaction times and 
# parameter variability. Psychonomic Bulletin & Review 9, 438–481 
# [2] Wiecki, T., Sofer, I., and Frank, M. (2013). HDDM: Hierarchical Bayesian 
# estimation of the Drift-Diffusion Model in Python. Frontiers in 
# Neuroinformatics 7, 14.

# Below, we also use ML fitting (not shown in paper).



# Load relevant packages:

import csv
import kabuki # 0.6.1
import hddm # 0.6.1
import sys # 3.5.5
import os
import numpy as np # 1.11.3
import matplotlib.pyplot as plt
import pandas # 0.23.4

# Adding a monkey in order to save the models dbs:
import pickle
# Get around a problem with saving regression outputs in Python 3
def savePatch(self, fname):
    with open(fname, 'wb') as f:
        pickle.dump(self, f)
hddm.HDDM.savePatch = savePatch

saved_model_name = 'miceq'



# Load data:

data = hddm.load_csv('goodAllMice_id_stim_rtNorm_choiceQ_0.csv')
data = hddm.utils.flip_errors(data)
docdir = 'DDM'
os.chdir(docdir)



# DDM fit:

# Note: Because we were interested in the effect of backgrounds on both the 
# signal and noise and as we used a custom code, in which the diffusion 
# variance is set to 1 (ref [2] above), the fitted parameters were scaled to 
# quantify the backgrounds effect on the drift and diffusion variance (see 
# Method and Matlab code: Figure4FigureS1FigureS2.m).  

# Running separate models for each individual mouse:

subj_params_chi = []
subj_params_ML = []
for subj_idx, subj_data in data.groupby('subj_idx'):
    # chi-square optimization:
   m_subj_chi = hddm.HDDM(subj_data, include=('z'), 
                      depends_on={'v': ['stim'], 
                                  'a': ['stim']}, 
                                  p_outlier=.05)
   subj_params_chi.append(m_subj_chi.optimize('chisquare', quantiles=(.1, .3, .5, .7, .9 )))
   # ML optimization:
   m_subj_ML = hddm.HDDM(subj_data, include=('z'), 
                      depends_on={'v': ['stim'], 
                                  'a': ['stim']}, 
                                  p_outlier=.05)
   subj_params_ML.append(m_subj_ML.optimize('ML'))
   
# Saving subject parameters: 
params_by_subject_chi = pandas.DataFrame(subj_params_chi)
params_by_subject_ML = pandas.DataFrame(subj_params_ML)
params_by_subject_chi.to_csv('subj_params_'+saved_model_name+'_chi.csv')
params_by_subject_ML.to_csv('subj_params_'+saved_model_name+'_ML.csv')



# simulate data, based on point estimates:

stimNames = pandas.unique(data.stim)

df_chi = pandas.read_csv('subj_params_'+saved_model_name+'_chi.csv')
df_ML = pandas.read_csv('subj_params_'+saved_model_name+'_ML.csv')

dict_chi = {}
dict_ML = {}

for i in range(6):
    # Chi-square fit:
    df_chi_sub = df_chi.loc[i]
    t_chi_sub = df_chi_sub['t']
    z_chi_sub = df_chi_sub['z']    
    # ML fit:
    df_ML_sub = df_ML.loc[i]
    t_ML_sub = df_ML_sub['t']
    z_ML_sub = df_ML_sub['z']
    
    for s in range(len(stimNames)):
        # Chi-square fit:
        a_chi_sub = df_chi_sub['a('+stimNames[s]+')']
        v_chi_sub = df_chi_sub['v('+stimNames[s]+')']
        dict_chi['subj_'+str(i)+'__stim_'+stimNames[s]] = {'v':v_chi_sub, 
                'a':a_chi_sub, 't':t_chi_sub, 'z':z_chi_sub}
        # ML fit:
        a_ML_sub = df_ML_sub['a('+stimNames[s]+')']
        v_ML_sub = df_ML_sub['v('+stimNames[s]+')']
        dict_ML['subj_'+str(i)+'__stim_'+stimNames[s]] = {'v':v_ML_sub, 
                'a':a_ML_sub, 't':t_ML_sub, 'z':z_ML_sub}

# Chi-square fit:
simData_chi, params = hddm.generate.gen_rand_data(dict_chi,size=1000)
sorted_simData_chi = simData_chi.sort_values(by=['condition'])
del sorted_simData_chi['subj_idx']
sorted_simData_chi.to_csv('simData_'+saved_model_name+'_chi.csv')
# ML fit:
simData_ML, params = hddm.generate.gen_rand_data(dict_ML,size=1000)
sorted_simData_ML = simData_ML.sort_values(by=['condition'])
del sorted_simData_ML['subj_idx']
sorted_simData_ML.to_csv('simData_'+saved_model_name+'_ML.csv')


