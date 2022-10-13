#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on 29/09/21
@author: Rebecca Varney, University of Exeter (rmv203@exeter.ac.uk)

"""

#%%

# Analysis imports
import numpy as np

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
from matplotlib.lines import Line2D


#%%
## Loading in correlation arrays
rmse_cSoil_cmip5 = np.load('saved_variables/rmse_cSoil_cmip5_array.npy')
rmse_cSoil_cmip6 = np.load('saved_variables/rmse_cSoil_cmip6_array.npy')

rmse_NL_cSoil_cmip5 = np.load('saved_variables/rmse_NL_cSoil_cmip5_array.npy')
rmse_NL_cSoil_cmip6 = np.load('saved_variables/rmse_NL_cSoil_cmip6_array.npy')

rmse_npp_cmip5 = np.load('saved_variables/rmse_npp_cmip5_array.npy')
rmse_npp_cmip6 = np.load('saved_variables/rmse_npp_cmip6_array.npy')

rmse_tau_cmip5 = np.load('saved_variables/rmse_tau_cmip5_array.npy')
rmse_tau_cmip6 = np.load('saved_variables/rmse_tau_cmip6_array.npy')
rmse_NL_tau_cmip6 = np.load('saved_variables/rmse_NL_tau_cmip6_array.npy')

#
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
n_models = len(cmip6_models)
cmip6_dictionary_cSoil = {}
cmip6_dictionary_cSoil_NL = {}
cmip6_dictionary_npp = {}
cmip6_dictionary_tau = {}
cmip6_dictionary_tau_NL = {}
for model_i in range(n_models):
    cmip6_dictionary_cSoil[cmip6_models[model_i]] = rmse_cSoil_cmip6[model_i]
    cmip6_dictionary_cSoil_NL[cmip6_models[model_i]] = rmse_NL_cSoil_cmip6[model_i]
    cmip6_dictionary_npp[cmip6_models[model_i]] = rmse_npp_cmip6[model_i]
    cmip6_dictionary_tau[cmip6_models[model_i]] = rmse_tau_cmip6[model_i]
    cmip6_dictionary_tau_NL[cmip6_models[model_i]] = rmse_NL_tau_cmip6[model_i]

#
cmip5_models = ['BNU-ESM', 'CCSM4', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M']
n_models2 = len(cmip5_models)
cmip5_dictionary_cSoil = {}
cmip5_dictionary_cSoil_NL = {}
cmip5_dictionary_npp = {}
cmip5_dictionary_tau = {}
for model_j in range(n_models2):
    cmip5_dictionary_cSoil[cmip5_models[model_j]] = rmse_cSoil_cmip5[model_j]
    cmip5_dictionary_cSoil_NL[cmip5_models[model_j]] = rmse_NL_cSoil_cmip5[model_j]
    cmip5_dictionary_npp[cmip5_models[model_j]] = rmse_npp_cmip5[model_j]
    cmip5_dictionary_tau[cmip5_models[model_j]] = rmse_tau_cmip5[model_j]


#%%
# Setting up the figure

fig_figure1 = plt.figure(1, figsize=(106,24))
gs = gspec.GridSpec(1, 3, figure=fig_figure1, hspace=0.2, wspace=0.2)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':3,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':90,
    'xtick.labelsize':90,
    'ytick.labelsize':90,
    'font.size':90,
}
plt.rcParams.update(params)

label_list = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'BNU-ESM', 'CanESM2 / CanESM5', 'CCSM4 / CESM2', 'CNRM-ESM2-1', 'GFDL-ESM2G / GFDL-ESM4', 'GISS-E2-R', 'IPSL-CM5A-LR / IPSL-CM6A-LR', 'MIROC-ESM / MIROC-ES2L', 'MPI-ESM-LR / MPI-ESM1-2-LR', 'NorESM1-M / NorESM2-LM', 'HadGEM2-ES / UKESM1-0-LL', 'Ensemble means']

# cSoil
column = 0
row = 0
ax = fig_figure1.add_subplot(gs[row, column])

x = np.arange(14)  # the label locations
width = 0.45  # the width of the bars

cmip5_list = [0, 0, cmip5_dictionary_cSoil['BNU-ESM'], cmip5_dictionary_cSoil['CanESM2'], cmip5_dictionary_cSoil['CCSM4'], 0, cmip5_dictionary_cSoil['GFDL-ESM2G'], cmip5_dictionary_cSoil['GISS-E2-R'], cmip5_dictionary_cSoil['IPSL-CM5A-LR'], cmip5_dictionary_cSoil['MIROC-ESM'], cmip5_dictionary_cSoil['MPI-ESM-LR'], cmip5_dictionary_cSoil['NorESM1-M'], cmip5_dictionary_cSoil['HadGEM2-ES'], 8.93]
cmip6_list = [cmip6_dictionary_cSoil['ACCESS-ESM1-5'], cmip6_dictionary_cSoil['BCC-CSM2-MR'], 0, cmip6_dictionary_cSoil['CanESM5'], cmip6_dictionary_cSoil['CESM2'], cmip6_dictionary_cSoil['CNRM-ESM2-1'], cmip6_dictionary_cSoil['GFDL-ESM4'], 0, cmip6_dictionary_cSoil['IPSL-CM6A-LR'], cmip6_dictionary_cSoil['MIROC-ES2L'], cmip6_dictionary_cSoil['MPI-ESM1-2-LR'], cmip6_dictionary_cSoil['NorESM2-LM'], cmip6_dictionary_cSoil['UKESM1-0-LL'], 8.26]

rects1 = ax.bar(x - width/2, cmip5_list, width, color='b', label='CMIP5')
rects2 = ax.bar(x + width/2, cmip6_list, width, color='g', label='CMIP6')

cmip5_list_nitrogen = [0, 0, 0, 0, cmip5_dictionary_cSoil['CCSM4'], 0, 0, 0, 0, 0, 0, cmip5_dictionary_cSoil['NorESM1-M'], 0, 0]
cmip6_list_nitrogen = [cmip6_dictionary_cSoil['ACCESS-ESM1-5'], 0, 0, 0, cmip6_dictionary_cSoil['CESM2'], 0, 0, 0, 0, cmip6_dictionary_cSoil['MIROC-ES2L'], cmip6_dictionary_cSoil['MPI-ESM1-2-LR'], cmip6_dictionary_cSoil['NorESM2-LM'], cmip6_dictionary_cSoil['UKESM1-0-LL'], 0]

rects1 = ax.bar(x - width/2, cmip5_list_nitrogen, width, color='b',  hatch='/', label='nitrogen cycle')
rects2 = ax.bar(x + width/2, cmip6_list_nitrogen, width, color='g',  hatch='/', label='nitrogen cycle')

ax.text(0.5,1.09, '(a) Soil carbon', ha='center', transform=ax.transAxes, fontweight = 'bold')
ax.set_ylabel(r'$C_{s}$ (kg m$^{-2}$)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
#ax.legend()
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


# NPP
column = 1
row = 0
ax = fig_figure1.add_subplot(gs[row, column])

x = np.arange(14)  # the label locations
width = 0.45  # the width of the bars

cmip5_list = [0, 0, cmip5_dictionary_npp['BNU-ESM'], cmip5_dictionary_npp['CanESM2'], cmip5_dictionary_npp['CCSM4'], 0, cmip5_dictionary_npp['GFDL-ESM2G'], cmip5_dictionary_npp['GISS-E2-R'], cmip5_dictionary_npp['IPSL-CM5A-LR'], cmip5_dictionary_npp['MIROC-ESM'], cmip5_dictionary_npp['MPI-ESM-LR'], cmip5_dictionary_npp['NorESM1-M'], cmip5_dictionary_npp['HadGEM2-ES'], 0.277]
cmip6_list = [cmip6_dictionary_npp['ACCESS-ESM1-5'], cmip6_dictionary_npp['BCC-CSM2-MR'], 0, cmip6_dictionary_npp['CanESM5'], cmip6_dictionary_npp['CESM2'], cmip6_dictionary_npp['CNRM-ESM2-1'], cmip6_dictionary_npp['GFDL-ESM4'], 0, cmip6_dictionary_npp['IPSL-CM6A-LR'], cmip6_dictionary_npp['MIROC-ES2L'], cmip6_dictionary_npp['MPI-ESM1-2-LR'], cmip6_dictionary_npp['NorESM2-LM'], cmip6_dictionary_npp['UKESM1-0-LL'], 0.200]

rects1 = ax.bar(x - width/2, cmip5_list, width, color='b', label='CMIP5')
rects2 = ax.bar(x + width/2, cmip6_list, width, color='g', label='CMIP6')

cmip5_list_nitrogen = [0, 0, 0, 0, cmip5_dictionary_npp['CCSM4'], 0, 0, 0, 0, 0, 0, cmip5_dictionary_npp['NorESM1-M'], 0, 0]
cmip6_list_nitrogen = [cmip6_dictionary_npp['ACCESS-ESM1-5'], 0, 0, 0, cmip6_dictionary_npp['CESM2'], 0, 0, 0, 0, cmip6_dictionary_npp['MIROC-ES2L'], cmip6_dictionary_npp['MPI-ESM1-2-LR'], cmip6_dictionary_npp['NorESM2-LM'], cmip6_dictionary_npp['UKESM1-0-LL'], 0]

rects1 = ax.bar(x - width/2, cmip5_list_nitrogen, width, color='b',  hatch='/', label='nitrogen cycle')
rects2 = ax.bar(x + width/2, cmip6_list_nitrogen, width, color='g',  hatch='/', label='nitrogen cycle')

ax.text(0.5,1.09, '(b) NPP', ha='center', transform=ax.transAxes, fontweight = 'bold')
ax.set_ylabel(r'NPP (kg m$^{-2}$ yr$^{-1}$)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
#ax.legend()
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


# legend
handels = []
handels.extend([Line2D([0,0],[0,0], linewidth=20, color='g', label='CMIP6')])
handels.extend([Line2D([0,0],[0,0], linewidth=20, color='b', label='CMIP5')])
#handels.extend([Line2D([0,0],[0,0], linewidth=20, color='white', hatch='/', label='nitrogen cycle')])
labels = ['CMIP6', 'CMIP5']#, 'nitrogen cycle']
leg1 = ax.legend(handels, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -1.025), fontsize=90)#title='Model colors',
plt.gca().add_artist(leg1)



# tau
column = 2
row = 0
ax = fig_figure1.add_subplot(gs[row, column])

x = np.arange(14)  # the label locations
width = 0.45  # the width of the bars

cmip5_list = [0, 0, cmip5_dictionary_tau['BNU-ESM'], cmip5_dictionary_tau['CanESM2'], cmip5_dictionary_tau['CCSM4'], 0, cmip5_dictionary_tau['GFDL-ESM2G'], cmip5_dictionary_tau['GISS-E2-R'], cmip5_dictionary_tau['IPSL-CM5A-LR'], cmip5_dictionary_tau['MIROC-ESM'], cmip5_dictionary_tau['MPI-ESM-LR'], cmip5_dictionary_tau['NorESM1-M'], cmip5_dictionary_tau['HadGEM2-ES'], 228]
cmip6_list = [cmip6_dictionary_tau['ACCESS-ESM1-5'], cmip6_dictionary_tau['BCC-CSM2-MR'], 0, cmip6_dictionary_tau['CanESM5'], cmip6_dictionary_tau['CESM2'], cmip6_dictionary_tau['CNRM-ESM2-1'], cmip6_dictionary_tau['GFDL-ESM4'], 0, cmip6_dictionary_tau['IPSL-CM6A-LR'], cmip6_dictionary_tau['MIROC-ES2L'], cmip6_dictionary_tau['MPI-ESM1-2-LR'], cmip6_dictionary_tau['NorESM2-LM'], cmip6_dictionary_tau['UKESM1-0-LL'], 223]

rects1 = ax.bar(x - width/2, cmip5_list, width, color='b', label='CMIP5')
rects2 = ax.bar(x + width/2, cmip6_list, width, color='g', alpha=1, label='CMIP6')

cmip5_list_nitrogen = [0, 0, 0, 0, cmip5_dictionary_tau['CCSM4'], 0, 0, 0, 0, 0, 0, cmip5_dictionary_tau['NorESM1-M'], 0, 0]
cmip6_list_nitrogen = [cmip6_dictionary_tau['ACCESS-ESM1-5'], 0, 0, 0, cmip6_dictionary_tau['CESM2'], 0, 0, 0, 0, cmip6_dictionary_tau['MIROC-ES2L'], cmip6_dictionary_tau['MPI-ESM1-2-LR'], cmip6_dictionary_tau['NorESM2-LM'], cmip6_dictionary_tau['UKESM1-0-LL'], 0]

rects1 = ax.bar(x - width/2, cmip5_list_nitrogen, width, color='b',  hatch='/', label='nitrogen cycle')
rects2 = ax.bar(x + width/2, cmip6_list_nitrogen, width, color='g',  hatch='/', label='nitrogen cycle')

ax.text(0.5, 1.09, '(c) Soil carbon turnover time', ha='center', transform=ax.transAxes, fontweight = 'bold')
ax.set_ylabel(r'$\tau_{s}$ (yr)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
#ax.legend()
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')


#%%
fig_figure1.savefig('figures/fig05', bbox_inches='tight')
plt.close()