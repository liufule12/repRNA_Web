__author__ = 'aleeee'

import os

ARGS = ['rec_data', 'mode', 'k', 'd', 'n', 'lamada', 'lag', 'w']

METHODS_PHYCHE_INDEX = ['DAC', 'DCC', 'DACC', 'PC-PseDNC-General', 'SC-PseDNC-General']
METHODS_LAG = ['DAC', 'DCC', 'DACC']
METHODS_LAMADA_W = ['PseSSC', 'PseDPC', 'PC-PseDNC-General', 'SC-PseDNC-General']
METHODS_PSE = ['PC-PseDNC-General', 'SC-PseDNC-General']
METHODS_ACC = ['DAC', 'DCC', 'DACC']

RNA_ALP_NUM = 4

ALPHABET_RNA = "ACGU"
MAX_SEQ_NUM = 50

MODEL_PSESSC_PATH = os.getcwd() + '/webserver/model/PseNC_RL-single_0.4.jar'
MODEL_PSEDPC_PATH = os.getcwd() + '/webserver/model/Stru_PseFV_sealed.jar'