#!/usr/bin/env python

#Title:         Module to facilitate basic reaction calculations for Sub-Chandra research
#Author:        Adam Jacobs
#Creation Date: Nov. 1, 2013

#Usage: load as a module

#Description: This module contains code to facilitate carrying out some basic reaction calculations
#             for analysis in the Sub-Chandra II paper.
#
#Revision History
#Programmer            Date                    Change
# ------------------------------------------------------------------------------
#Adam Jacobs         11/01/2013              Code created

#TODO: 
# 1) 

#Notes on conventions, programming style:
#  1) All classes are in CamelCase with the first letter capitalized.  Class methods
#     are also in CamelCase with the first letter lower-case, e.g. myMethod().
#  2) Non-class functions and all variables are lowercase with underscores (_) 
#     acting as delimiters when needed/wanted.
#  3) An underscore prefix, e.g. _variable, means the named item is intended to be
#     private (Python doesn't provide for easily enforcing data hiding, it's all by 
#     convention).
#  4) Names that are in ALL_CAPS are intended as constants.

######################
### Global Imports ###
######################
import sys

#################################
### Global Data and Constants ###
#################################
#_EXAMPLE = 'string'

###############
### Classes ###
###############

#################
### Functions ###
#################
def plot_logenuc():
  import numpy as np
  import matplotlib
  #If 'inline' backend, we're running within IPython so don't change anything.
  #Otherwise, we're running as script to generate eps, so use PS backend
  if 'inline' not in matplotlib.get_backend():
    matplotlib.use('PS')
  import matplotlib.pyplot as plt

  matplotlib.rcParams['text.usetex'] = 'True'
  # Make data based on MESA calcs
  rho_exp = np.linspace(5.0, 7.0, 11)
  rho     = np.array([10**rexp for rexp in rho_exp])
  rho_nco = np.array([10**rexp for rexp in rho_exp[-3:]])

  # 200 MK rates
  enuc_cagob2  = np.array(
      [10**17.581540508, 
       10**17.800549914,
       10**18.024481339,
       10**18.254609218,
       10**18.492973480,
       10**18.738377513,
       10**18.986558564,
       10**19.234718720,
       10**19.477401243,
       10**19.704819234,
       10**19.897451724])
  enuc_tralf2  = np.array(
      [10**11.188457831,  #rho = 10**5.0
       10**11.626476644,  #rho = 10**5.2
       10**12.074446783,  #rho = 10**5.4
       10**12.531401511,  #rho = 10**5.6
       10**12.996679600,  #rho = 10**5.8
       10**13.469768554,  #rho = 10**6.0
       10**13.948445742,  #rho = 10**6.2
       10**14.432637076,  #rho = 10**6.4
       10**14.910185702,  #rho = 10**6.6
       10**15.450434567,  #rho = 10**6.8
       10**16.013534830]) #rho = 10**7.0
  enuc_nagf2   = np.array(
      [10**10.356064284, 
       10**10.596664766,
       10**10.840074663,
       10**11.083643491,
       10**11.322076922,
       10**11.545046637,
       10**11.796568378,
       10**12.066216636,
       10**12.348913239,
       10**12.646625087,
       10**12.961628265])
  enuc_cago2   = np.array(
      [10**7.066187025, 
       10**7.303537281,
       10**7.544452034,
       10**7.787414219,
       10**8.029050316,
       10**8.262882790,
       10**8.482066479,
       10**8.743164789,
       10**9.015982418,
       10**9.302270992,
       10**9.604056810])
  enuc_nco2    = np.array(
      [10**7.621642798,
       10**7.907935225,
       10**8.209727138])
  
  # 300 MK rates
  enuc_cagob3  = np.array(
      [10**19.030792008, 
       10**19.241106991,
       10**19.454092787,
       10**19.670440935,
       10**19.891022033,
       10**20.116932101,
       10**20.349550944,
       10**20.591053956,
       10**20.837613055,
       10**21.086204225,
       10**21.333254625])
  enuc_tralf3  = np.array(
      [10**13.785826231,  #rho = 10**5.0
       10**14.206456199,  #rho = 10**5.2
       10**14.632427790,  #rho = 10**5.4
       10**15.065124085,  #rho = 10**5.6
       10**15.506286282,  #rho = 10**5.8
       10**15.957694147,  #rho = 10**6.0
       10**16.417507200,  #rho = 10**6.2
       10**16.885601880,  #rho = 10**6.4
       10**17.360983029,  #rho = 10**6.6
       10**17.842218167,  #rho = 10**6.8
       10**18.324423763]) #rho = 10**7.0
  enuc_nagf3   = np.array(
      [10**13.685591364, 
       10**13.909659660,
       10**14.141720810,
       10**14.379825361,
       10**14.621578864,
       10**14.865447555,
       10**15.108015882,
       10**15.342720865,
       10**15.567331414,
       10**15.829352083,
       10**16.103271007])
  enuc_cago3   = np.array(
      [10**9.738896065, 
       10**9.959526033,
       10**10.185497624,
       10**10.419619014,
       10**10.658290572,
       10**10.900155663,
       10**11.143164162,
       10**11.383040770,
       10**11.611826212,
       10**11.842872776,
       10**12.107808383])
  enuc_nco3    = np.array(
      [10**10.973769648,
       10**11.204818056,
       10**11.469756036])
  
  #make the plots
  #fig, ax = plt.subplots(nrows=2, ncols=1)
  fig, ax = plt.subplots(nrows=1, ncols=1)
  
  ###plot the 200MK data
  ax.plot(rho,     enuc_cagob2, label=r'CagO-by')
  ax.plot(rho,     enuc_tralf2, label=r'3$\alpha$')
  ax.plot(rho,     enuc_nagf2,  label=r'NagF')
  ax.plot(rho,     enuc_cago2,  label=r'CagO')
  ax.plot(rho_nco, enuc_nco2,   label=r'NCO')

  #ax[0].plot(rho,     enuc_tralf3, 'b--')
  #ax[0].plot(rho,     enuc_cagob3, 'g--')
  #ax[0].plot(rho,     enuc_cago3,  'r--')
  #ax[0].plot(rho_nco, enuc_nco3,   'c--')

  #tune plot settings
  ax.set_yscale('log')
  ax.set_xscale('log')
  ax.set_ylabel(r'$\dot{\epsilon}~[\rm{erg}~\rm{g}^{-1}~\rm{s}^{-1}]$', fontsize='x-large')
  ax.set_ylim(1.e6, 1.e22)

  ax.set_xlabel(r'$\rho~[\rm{g}~\rm{cm}^{-3}]$', fontsize='x-large')
  
  ax.tick_params(labelsize='x-large')
  ax.legend(loc=2, fontsize='medium')
  ax.set_title(r'$T = 200~\rm{MK}$', fontsize='x-large')

  ###plot the 300MK data
  #ax[1].plot(rho,     enuc_cagob3, label=r'CagO-by')
  #ax[1].plot(rho,     enuc_tralf3, label=r'3$\alpha$')
  #ax[1].plot(rho,     enuc_nagf3,  label=r'NagF')
  #ax[1].plot(rho,     enuc_cago3,  label=r'CagO')
  #ax[1].plot(rho_nco, enuc_nco3,   label=r'NCO')

  ##tune plot settings
  #ax[1].set_yscale('log')
  #ax[1].set_xscale('log')
  #ax[1].set_ylabel(r'$\dot{\epsilon}~[\rm{erg}~\rm{g}^{-1}~\rm{s}^{-1}]$', fontsize='xx-large')
  #ax[1].set_ylim(1.e6, 1.e22)

  #ax[1].set_xlabel(r'$\rho~[\rm{g}~\rm{cm}^{-3}]$', fontsize='xx-large')
  #
  #ax[1].tick_params(labelsize='x-large')
  #ax[1].legend(loc=2, fontsize='large')
  #ax[1].set_title(r'$T = 300 ~\rm{MK}$', fontsize='xx-large')

  fig.set_size_inches(5.0*1.25, 5.0)
  #fig.tight_layout()
  #This fixes a bug in mpl's eps generation
  matplotlib.rc('ps', usedistiller='xpdf')
  fig.savefig('rxns.eps', bbox_inches='tight')

  return



#################
### Execution ###
#################
#This is only for testing. rxncalcs.py is intended to be used as a module.
if __name__== "__main__":
  if(len(sys.argv) <= 1):
    #TODO: Add arg checks
    pass

  #Put tests here
  plot_logenuc()
