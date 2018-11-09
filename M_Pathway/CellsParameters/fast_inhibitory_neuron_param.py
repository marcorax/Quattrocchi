# -*- coding: utf-8 -*-
# @Author: mmilde
# @Date:   2018-01-10 15:35:29
# @Last Modified by:   mmilde
# @Last Modified time: 2018-01-17 15:50:45

"""This file contains default parameter for dpi neuron. For more details on model see
models/equations/dpi_neuron.py

Attributes:
    parameters (dict): Neuron parameters
"""

from brian2 import pF, nS, mV, ms, pA, nA
from NCSBrian2Lib.models.parameters import constants

fast_in_parameters = {"kn": constants.KAPPA_N,
              "kp": constants.KAPPA_P,
              "Ut": constants.UT,
              "Io": constants.I0,
              "Cmem": 1.5 * pF,
              "Ispkthr": 1. * nA,
              "refP": 1. * ms,
              "Ireset": 0.6 * pA,
              "Iconst": constants.I0,
              ##################
              "Itau": 20. * pA,  #              "Itau": 20. * pA,  #
              "Ishunt": constants.I0,  #
              "Ith": 0.9 * pA,  #              "Ith": 0.6 * pA,  #
              #  ADAPTATION  #################
              "Ica": 2. * pA,
              "Itauahp": 1 * pA,
              "Ithahp": 1 * pA,
              "Cahp": 1 * pF,
              "Iahp": constants.I0,
              #  POSTIVE FEEDBACK #################
              "Iath": 0.5 * nA,
              "Iagain": 50. * pA,
              "Ianorm": 10. * pA,
              }