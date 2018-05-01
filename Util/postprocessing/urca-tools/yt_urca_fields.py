"""
Module for setting up some commonly used fields for MAESTRO urca datasets with yt.

Donald E. Willcox
"""

import yt
import numpy as np

class PhysicalConstants:
    N_AVO = 6.02214129e23

class UrcaShellFields(object):
    def __init__(self):
        return

    def setup(self, ds):
        # ds should be a MAESTRO dataset in yt corresponding to the urca shell variables

        try:
            ds.add_field(('boxlib','urca23_shell_unscaled'), units='',
                         function=UrcaShellFields._urca23_shell_unscaled)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: urca23_shell_unscaled field could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','urca23_shell'), units='',
                         function=UrcaShellFields._urca23_shell)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: urca23_shell field could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','weak_xrate_na23'), units='',
                         function=UrcaShellFields._weak_xrate_na23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: weak_xrate_na23 field could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','weak_xrate_ne23'), units='',
                         function=UrcaShellFields._weak_xrate_ne23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: weak_xrate_ne23 field could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','reduced_x23'), units='',
                         function=UrcaShellFields._reduced_x23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: reduced_x23 could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','enucloss_epart_ecap23'), units='',
                         function=UrcaShellFields._enucloss_epart_ecap23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: enucloss_epart_ecap23 could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','enucloss_epart_beta23'), units='',
                         function=UrcaShellFields._enucloss_epart_beta23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: enucloss_epart_beta23 could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','enucloss_epart_urca23'), units='',
                         function=UrcaShellFields._enucloss_epart_urca23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: enucloss_epart_urca23 could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','enucloss_eprho_urca23'), units='g',
                         function=UrcaShellFields._enucloss_eprho_urca23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: enucloss_eprho_urca23 could not be added because it relies on a field not in the dataset.')
            pass
        
    @staticmethod
    def _urca23_shell_unscaled(field, data):
        return data['boxlib','ecap23']*data['boxlib','beta23']*data['boxlib','X(na23)']*data['boxlib','X(ne23)']

    @staticmethod
    def _urca23_shell(field, data):
        return data['boxlib','urca23_shell_unscaled']/np.amax(data['boxlib','urca23_shell_unscaled'])

    @staticmethod
    def _weak_xrate_na23(field, data):
        return data['boxlib','beta23']*data['boxlib','X(ne23)'] - data['boxlib', 'X(na23)']*data['boxlib','ecap23']

    @staticmethod
    def _weak_xrate_ne23(field, data):
        return data['boxlib', 'X(na23)']*data['boxlib','ecap23'] - data['boxlib','beta23']*data['boxlib','X(ne23)']

    @staticmethod
    def _reduced_x23(field, data):
        return data['boxlib','X(na23)']*data['boxlib', 'X(ne23)']/(data['boxlib','X(na23)']+data['boxlib', 'X(ne23)'])

    @staticmethod
    def _enucloss_epart_ecap23(field, data):
        return -data['boxlib','X(na23)']*data['boxlib', 'epart_ecap23']*PhysicalConstants.N_AVO/23.0

    @staticmethod
    def _enucloss_epart_beta23(field, data):
        return -data['boxlib','X(ne23)']*data['boxlib', 'epart_beta23']*PhysicalConstants.N_AVO/23.0

    @staticmethod
    def _enucloss_epart_urca23(field, data):
        return data['boxlib', 'enucloss_epart_ecap23'] + data['boxlib', 'enucloss_epart_beta23']

    @staticmethod
    def _enucloss_eprho_urca23(field, data):
        # Energy loss rate due to A=23 urca reactions in erg/g/s * g/cm^3 * cm^3 = erg/s
        return data['boxlib', 'enucloss_epart_urca23'] * data['boxlib', 'density'] * data['boxlib', 'cell_volume']
