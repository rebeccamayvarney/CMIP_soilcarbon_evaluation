#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wednesday 4th August 2021
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


## Comparing the spatial correlations between models and obs to see if there are any relationships
# i.e. if a model simulates something well does it do something else badly? etc.

## Loading in correlation arrays
#CsMs_correlations_cmip6 = np.load('saved_variables/cSoilmoisture_correlations_cmip6.npy')
CsT_correlations_cmip6 = np.load('saved_variables/cSoiltau_correlations_cmip6.npy')
CsNPP_correlations_cmip6 = np.load('saved_variables/cSoilNPP_correlations_cmip6.npy')

NPPMs_correlations_cmip6 = np.load('saved_variables/NPPmrso_correlations_cmip6.npy')
NPPT_correlations_cmip6 = np.load('saved_variables/NPPtas_correlations_cmip6.npy')

tausMs_correlations_cmip6 = np.load('saved_variables/tausmrsos_correlations_cmip6.npy')
tausT_correlations_cmip6 = np.load('saved_variables/tautas_correlations_cmip6.npy')


#%% cSoilAbove1m
#CsTau_Above1m_C = -0.00142
#CsNPP_Above1m_C = 0.134
#CsTau_Above1m_N = 0.153
#CsNPP_Above1m_N = 0.261
#
#tauMs_Above1m_C = 0.211
#tausT_Above1m_C = -0.2542
#tauMs_Above1m_N = 0.335
#tausT_Above1m_N = -0.434


#%% Ensemble means

cmip6_CsT = -0.0605
cmip6_CsNPP = 0.424
cmip6_NPPMs = 0.494
cmip6_NPPT = 0.257
cmip6_tausMs = 0.198
cmip6_tausT = -0.398


#%%
#obs_CsMs = CsMs_correlations_cmip6[0]
obs_CsT = CsT_correlations_cmip6[0]
obs_CsNPP = CsNPP_correlations_cmip6[0]
obs_NPPMs = NPPMs_correlations_cmip6[0]
obs_NPPT = NPPT_correlations_cmip6[0]
obs_tausMs = tausMs_correlations_cmip6[0]
obs_tausT = tausT_correlations_cmip6[0]


cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
n_models = len(cmip6_models)
model_shapes = ['o', '^', 'v', '1', 's', '*', 'x', '+', 'd', '<', 'h']
colors_cmip6 = ['peachpuff', '#fb8072', '#80b1d3', 'dodgerblue', 'red', 'darkcyan', 'darkgreen', 'olive', 'gold', 'orange', 'darkseagreen']
                
#%%
# Setting up the figure

fig_figure1 = plt.figure(1, figsize=(78,20))
gs = gspec.GridSpec(1, 3, figure=fig_figure1, hspace=0.35, wspace=0.35)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':3,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':72,
    'xtick.labelsize':72,
    'ytick.labelsize':72,
    'font.size':72,
}
plt.rcParams.update(params)


#%%
# plotting

min_axis_value = 1
max_axis_value = -0.5

for model_i in range(n_models):
    model = cmip6_models[model_i] # seleting the models

    model_j = model_i+1

    # 1
    column = 0
    row = 0
    ax = fig_figure1.add_subplot(gs[row, column])

    x1 = CsNPP_correlations_cmip6[model_j]
    y1 = CsT_correlations_cmip6[model_j]
    ax.axhline(y=0, color='k', alpha=0.01)
    ax.axvline(x=0, color='k', alpha=0.01)
    plt.plot(x1, y1, marker='o', color=colors_cmip6[model_i], markersize=50)
    plt.plot(obs_CsNPP, obs_CsT, marker='*', color='k', markersize=60)
    plt.plot(cmip6_CsNPP, cmip6_CsT, marker='o', color='grey', markersize=50)
    #plt.plot(CsNPP_Above1m_C, CsTau_Above1m_C, marker='o', color='darkseagreen', markersize=50)
    #plt.plot(CsNPP_Above1m_N, CsTau_Above1m_N, marker='o', color='darkseagreen', markersize=50)

    # one to one line         
    one_to_one_line = np.linspace(min_axis_value, max_axis_value, 100)
    #plt.plot(one_to_one_line, one_to_one_line, 'grey', linewidth=1)
    
    ax.text(-0.25, 0.95, '(a)',transform=ax.transAxes, ha='center', fontweight = 'bold', fontsize=72)
    
    plt.xlabel(r'$C_{s}$-NPP')
    plt.ylabel(r'$C_{s}$-$\tau_{s}$')
    #plt.title(r'$C_{s}$-$\tau_{s}$ Vs $C_{s}$-NPP')
    plt.xlim((-1,1))
    plt.ylim((-1,1))
    # plt.xticks([0,0.5,1,1.5,2])
    # plt.yticks([0,0.5,1,1.5,2])

    # if model_i==0:
    #     # Calculating r-squared value
    #     correlation_matrix = ma.corrcoef(CsNPP_correlations_cmip6, CsT_correlations_cmip6)
    #     correlation_xy = correlation_matrix[0,1]
    #     r_squared = correlation_xy**2
    #     print('1:', r_squared, correlation_matrix)
    #     plt.text(-0.4, 0.85, r'r=%0.3f' % correlation_xy, fontsize=20)


    # 2
    column = 1
    row = 0
    ax = fig_figure1.add_subplot(gs[row, column])

    x2 = NPPMs_correlations_cmip6[model_j]
    y2 = NPPT_correlations_cmip6[model_j]
    ax.axhline(y=0, color='k', alpha=0.01)
    ax.axvline(x=0, color='k', alpha=0.01)
    
    plt.plot(x2, y2, marker='o', color=colors_cmip6[model_i], markersize=50)
    plt.plot(obs_NPPMs, obs_NPPT, marker='*', color='k', markersize=60)
    plt.plot(cmip6_NPPMs, cmip6_NPPT, marker='o', color='grey', markersize=50)

    # one to one line         
    one_to_one_line = np.linspace(min_axis_value, max_axis_value, 100)
    #plt.plot(one_to_one_line, one_to_one_line, 'grey', linewidth=1)

    ax.text(-0.25, 0.95, '(b)',transform=ax.transAxes, ha='center', fontweight = 'bold', fontsize=72)

    plt.xlabel(r'NPP-$\theta$')
    plt.ylabel(r'NPP-T')
    #plt.title(r'NPP-T Vs NPP-$\theta$')
    plt.xlim((-1,1))
    plt.ylim((-1,1))

    # if model_i==0:
    #     # Calculating r-squared value
    #     correlation_matrix = ma.corrcoef(NPPMs_correlations_cmip6, NPPT_correlations_cmip6)
    #     correlation_xy = correlation_matrix[0,1]
    #     r_squared = correlation_xy**2
    #     print('1:', r_squared, correlation_matrix)
    #     plt.text(-0.4, 0.85, r'r=%0.3f' % correlation_xy, fontsize=20)


    # 3
    column = 2
    row = 0
    ax = fig_figure1.add_subplot(gs[row, column])

    x3 = tausMs_correlations_cmip6[model_j]
    y3 = tausT_correlations_cmip6[model_j]
    ax.axhline(y=0, color='k', alpha=0.01)
    ax.axvline(x=0, color='k', alpha=0.01)
    plt.plot(x3, y3, marker='o', color=colors_cmip6[model_i], markersize=50)
    plt.plot(obs_tausMs, obs_tausT, marker='*', color='k', markersize=60)
    plt.plot(cmip6_tausMs, cmip6_tausT, marker='o', color='grey', markersize=50)
    #plt.plot(tauMs_Above1m_C, tausT_Above1m_C, marker='1', color='darkseagreen', markersize=50)
    #plt.plot(tauMs_Above1m_N, tausT_Above1m_N, marker='<', color='darkseagreen', markersize=50)

    # one to one line         
    one_to_one_line = np.linspace(min_axis_value, max_axis_value, 100)
    #plt.plot(one_to_one_line, one_to_one_line, 'grey', linewidth=1)

    ax.text(-0.25, 0.95, '(c)',transform=ax.transAxes, ha='center', fontweight = 'bold', fontsize=72)

    plt.xlabel(r'$\tau_{s}$-$\theta$')
    plt.ylabel(r'$\tau_{s}$-T')
    #plt.title(r'$\tau_{s}$-T Vs $\tau_{s}$-$\theta$')
    plt.xlim((-1,1))
    plt.ylim((-1,1))

    # if model_i==0:
    #     # Calculating r-squared value
    #     correlation_matrix = ma.corrcoef(tausMs_correlations_cmip6, tausT_correlations_cmip6)
    #     correlation_xy = correlation_matrix[0,1]
    #     r_squared = correlation_xy**2
    #     print('3:', r_squared, correlation_matrix)
    #     plt.text(-0.4, 0.85, r'r=%0.3f' % correlation_xy, fontsize=20)

    # # 4
    # column = 1
    # row = 0
    # ax = fig_figure1.add_subplot(gs[row, column])

    # x4 = CsNPP_correlations_cmip6[model_j]
    # y4 = CsMs_correlations_cmip6[model_j]
    # plt.plot(x4, y4, marker=model_shapes[model_i], color='b', markersize=50)
    # plt.plot(obs_CsNPP, obs_CsMs, marker='X', color='k', markersize=60)

    # # one to one line         
    # one_to_one_line = np.linspace(min_axis_value, max_axis_value, 100)
    # #plt.plot(one_to_one_line, one_to_one_line, 'grey', linewidth=1)

    # plt.xlabel(r'$C_{s}$-NPP')
    # plt.ylabel(r'$C_{s}$-$\theta$')
    # #plt.title(r'$C_{s}$-$\theta$ Vs $C_{s}$-NPP')
    # plt.xlim((-1,1))
    # plt.ylim((-1,1))

    # if model_i==0:
    #     # Calculating r-squared value
    #     correlation_matrix = ma.corrcoef(CsNPP_correlations_cmip6, CsMs_correlations_cmip6)
    #     correlation_xy = correlation_matrix[0,1]
    #     r_squared = correlation_xy**2
    #     print('3:', r_squared, correlation_matrix)
    #     plt.text(-0.4, 0.85, r'r=%0.3f' % correlation_xy, fontsize=20)


handels = []
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='peachpuff', label='ACCESS-ESM1-5')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='#fb8072', label='BCC-CSM2-MR')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='#80b1d3', label='CanESM5')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='dodgerblue', label='CESM2')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='red', label='CNRM-ESM2-1')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='darkcyan', label='GFDL-ESM4')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='darkgreen', label='IPSL-CM6A-LR')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='olive', label='MIROC-ES2L')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='gold', label='MPI-ESM1-2-LR')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='orange', label='NorESM2-LM')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='darkseagreen', label='UKESM1-0-LL')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='o', markersize=50, color='grey', label='CMIP6 ensemble mean')])
#handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='*', markersize=50, color='k', label='OBSERVATIONS')])
label = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL', 'CMIP6 ensemble mean']
#leg = ax.legend(handels, label, loc='lower center', ncol=4, bbox_to_anchor=(-0.7, -0.8), labelspacing=2, title='CMIP6', fontsize=72)
leg = ax.legend(handels, label, loc='upper center', ncol=6, bbox_to_anchor=(-0.85, -0.65), labelspacing=0.75, title='CMIP6', fontsize=72)
plt.gca().add_artist(leg)




#########################################################################################
#%% CMIP5

## Loading in correlation arrays
CsMs_correlations_cmip5 = np.load('saved_variables/cSoilmoisture_correlations_cmip5.npy')
print(len(CsMs_correlations_cmip5))
CsT_correlations_cmip5 = np.load('saved_variables/cSoiltau_correlations_cmip5.npy')
CsNPP_correlations_cmip5 = np.load('saved_variables/cSoilNPP_correlations_cmip5.npy')

NPPMs_correlations_cmip5 = np.load('saved_variables/NPPmrso_correlations_cmip5.npy')
NPPT_correlations_cmip5 = np.load('saved_variables/NPPtas_correlations_cmip5.npy')

tausMs_correlations_cmip5 = np.load('saved_variables/tausmrsos_correlations_cmip5.npy')
tausT_correlations_cmip5 = np.load('saved_variables/tautas_correlations_cmip5.npy')


# Ensemble means

cmip5_CsT = 0.108
cmip5_CsNPP = 0.222
cmip5_NPPMs = 0.258
cmip5_NPPT = 0.230
cmip5_tausMs = 0.210
cmip5_tausT = -0.314


#
#obs_CsMs = CsMs_correlations_cmip5[0]
#obs_CsT = CsT_correlations_cmip5[0]
#obs_CsNPP = CsNPP_correlations_cmip5[0]
#obs_NPPMs = NPPMs_correlations_cmip5[0]
#obs_NPPT = NPPT_correlations_cmip5[0]
#obs_tausMs = tausMs_correlations_cmip5[0]
#obs_tausT = tausT_correlations_cmip5[0]


cmip5_models = ['BNU-ESM', 'CCSM4', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M']
n_models = len(cmip5_models)
model_shapes = ['o', '^', 'v', '1', 's', '*', 'x', '+', 'd', '<']
colors_cmip5 = ['darkblue', 'dodgerblue', '#80b1d3', 'darkcyan', '#8dd3c7', 'darkseagreen', 'darkgreen', 'olive', 'gold', 'orange']


#%%
# plotting

for model_i in range(n_models):
    model = cmip5_models[model_i] # seleting the models

    model_j = model_i+1

    # 1
    column = 0
    row = 0
    ax = fig_figure1.add_subplot(gs[row, column])

    x1 = CsNPP_correlations_cmip5[model_j]
    y1 = CsT_correlations_cmip5[model_j]
    plt.plot(x1, y1, marker='X', color=colors_cmip5[model_i], markersize=50)
    plt.plot(cmip5_CsNPP, cmip5_CsT, marker='X', color='grey', markersize=50)

    # one to one line         
    one_to_one_line = np.linspace(min_axis_value, max_axis_value, 100)
    #plt.plot(one_to_one_line, one_to_one_line, 'grey', linewidth=1)

    plt.xlabel(r'$C_{s}$-NPP')
    plt.ylabel(r'$C_{s}$-$\tau_{s}$')
    #plt.title(r'$C_{s}$-$\tau_{s}$ Vs $C_{s}$-NPP')
    plt.xlim((-1,1))
    plt.ylim((-1,1))
    # plt.xticks([0,0.5,1,1.5,2])
    # plt.yticks([0,0.5,1,1.5,2])

    # if model_i==0:
    #     # Calculating r-squared value
    #     correlation_matrix = ma.corrcoef(CsNPP_correlations_cmip5, CsT_correlations_cmip5)
    #     correlation_xy = correlation_matrix[0,1]
    #     r_squared = correlation_xy**2
    #     print('1:', r_squared, correlation_matrix)
    #     plt.text(-0.6, 0.85, r'r=%0.3f' % correlation_xy, fontsize=20)


    # 2
    column = 1
    row = 0
    ax = fig_figure1.add_subplot(gs[row, column])

    x2 = NPPMs_correlations_cmip5[model_j]
    y2 = NPPT_correlations_cmip5[model_j]
    plt.plot(x2, y2, marker='X', color=colors_cmip5[model_i], markersize=50)
    plt.plot(cmip5_NPPMs, cmip5_NPPT, marker='X', color='grey', markersize=50)

    # one to one line         
    one_to_one_line = np.linspace(min_axis_value, max_axis_value, 100)
    #plt.plot(one_to_one_line, one_to_one_line, 'grey', linewidth=1)

    plt.xlabel(r'NPP-$\theta$')
    plt.ylabel(r'NPP-T')
    #plt.title(r'NPP-T Vs NPP-$\theta$')
    plt.xlim((-1,1))
    plt.ylim((-1,1))

    # if model_i==0:
    #     # Calculating r-squared value
    #     correlation_matrix = ma.corrcoef(NPPMs_correlations_cmip5, NPPT_correlations_cmip5)
    #     correlation_xy = correlation_matrix[0,1]
    #     r_squared = correlation_xy**2
    #     print('1:', r_squared, correlation_matrix)
    #     plt.text(-0.6, 0.85, r'r=%0.3f' % correlation_xy, fontsize=20)


    # 3
    column = 2
    row = 0
    ax = fig_figure1.add_subplot(gs[row, column])

    x3 = tausMs_correlations_cmip5[model_j]
    y3 = tausT_correlations_cmip5[model_j]
    plt.plot(x3, y3, marker='X', color=colors_cmip5[model_i], markersize=50)
    plt.plot(cmip5_tausMs, cmip5_tausT, marker='X', color='grey', markersize=50)

    # one to one line         
    one_to_one_line = np.linspace(min_axis_value, max_axis_value, 100)
    #plt.plot(one_to_one_line, one_to_one_line, 'grey', linewidth=1)

    plt.xlabel(r'$\tau_{s}$-$\theta$')
    plt.ylabel(r'$\tau_{s}$-T')
    ###plt.title(r'$\tau_{s}$-T Vs $\tau_{s}$-$\theta$')
    plt.xlim((-1,1))
    plt.ylim((-1,1))

    # if model_i==0:
    #     # Calculating r-squared value
    #     correlation_matrix = ma.corrcoef(tausMs_correlations_cmip5, tausT_correlations_cmip5)
    #     correlation_xy = correlation_matrix[0,1]
    #     r_squared = correlation_xy**2
    #     print('3:', r_squared, correlation_matrix)
    #     plt.text(-0.6, 0.85, r'r=%0.3f' % correlation_xy, fontsize=20)

    # # 4
    # column = 1
    # row = 0
    # ax = fig_figure1.add_subplot(gs[row, column])

    # x4 = CsNPP_correlations_cmip5[model_j]
    # y4 = CsMs_correlations_cmip5[model_j]
    # plt.plot(x4, y4, marker=model_shapes[model_i], color='r', markersize=50)

    # # one to one line         
    # one_to_one_line = np.linspace(min_axis_value, max_axis_value, 100)
    # #plt.plot(one_to_one_line, one_to_one_line, 'grey', linewidth=1)

    # plt.xlabel(r'$C_{s}$-NPP')
    # plt.ylabel(r'$C_{s}$-$\theta$')
    # #plt.title(r'$C_{s}$-$\theta$ Vs $C_{s}$-NPP')
    # plt.xlim((-1,1))
    # plt.ylim((-1,1))

    # if model_i==0:
    #     # Calculating r-squared value
    #     correlation_matrix = ma.corrcoef(CsNPP_correlations_cmip5, CsMs_correlations_cmip5)
    #     correlation_xy = correlation_matrix[0,1]
    #     r_squared = correlation_xy**2
    #     print('3:', r_squared, correlation_matrix)
    #     plt.text(-0.6, 0.85, r'r=%0.3f' % correlation_xy, fontsize=20)

handels = []
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='darkblue', label='BNU-ESM')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='dodgerblue', label='CCSM4')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='#80b1d3', label='CanESM2')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='darkcyan', label='GFDL-ESM2G')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='#8dd3c7', label='GISS-E2-R')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='darkseagreen', label='HadGEM2-ES')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='darkgreen', label='IPSL-CM5A-LR')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='olive', label='MIROC-ESM')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='gold', label='MPI-ESM-LR')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='orange', label='NorESM1-LM')])
handels.extend([Line2D([0,0],[0,0], linestyle='None', marker='X', markersize=50, color='grey', label='CMIP5 ensemble mean')])
label = ['BNU-ESM', 'CCSM4', 'CanESM2', 'GFDL-ESM2G', 'GISS-E2-R', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MPI-ESM-LR', 'NorESM1-M', 'CMIP5 ensemble mean']
leg = ax.legend(handels, label, loc='upper center', ncol=6, bbox_to_anchor=(-0.85, -0.25), labelspacing=0.75, title='CMIP5', fontsize=72)
plt.gca().add_artist(leg)
    
handels3 = []
handels3.extend([Line2D([0,0],[0,0], linestyle='None', marker='*', markersize=50, color='k', label='Benchmark datasets')])
label3 = ['Benchmark datasets']
leg = ax.legend(handels3, label3, loc='upper center', ncol=1, bbox_to_anchor=(-0.86, -1.05), labelspacing=0.75, title='Observations', fontsize=72)
plt.gca().add_artist(leg)




#%%
fig_figure1.savefig('correlation_figures/fig08', bbox_inches='tight')
plt.close()