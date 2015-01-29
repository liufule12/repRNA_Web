__author__ = 'aleeee'

import os

ARGS = ['rec_data', 'mode', 'k', 'd', 'n', 'lamada', 'lag', 'w']

METHODS_KMER = ['Nucleotide composition', 'Dinucleotide composition', 'Trinucleotide composition',
                'Tetranucleotide composition', 'Pentanucleotide composition', 'Hexanucleotide composition']
METHODS_PHYCHE_INDEX = ['Auto covariance', 'Cross covariance', 'Auto-cross covariance', 'PC-PseDNC', 'SC-PseDNC']
METHODS_LAG = ['Auto covariance', 'Cross covariance', 'Auto-cross covariance']
METHODS_LAMADA_W = ['PseSSC', 'PseDPC', 'PC-PseDNC', 'SC-PseDNC']
METHODS_PSE = ['PC-PseDNC', 'SC-PseDNC']
METHODS_ACC = ['Auto covariance', 'Cross covariance', 'Auto-cross covariance']

RNA_ALP_NUM = 4

ALPHABET_RNA = "ACGU"
MAX_SEQ_NUM = 50

MODEL_PSESSC_PATH = os.getcwd() + '/webserver/model/PseNC_RL-single_0.4.jar'
MODEL_PSEDPC_PATH = os.getcwd() + '/webserver/model/Stru_PseFV_sealed_0.4.jar'