#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pydna.utils import rc

"""
Copied from main-3d92a74abb.js

Processed with https://www.codeconvert.ai/javascript-to-python-converter

r.TmCalc.prototype.nnBr={aa:{dh:-9.1,ds:-22.2},
                         tt:{dh:-9.1,ds:-22.2},
                         at:{dh:-8.6,ds:-20.4},
                         ta:{dh:-6,ds:-21.3},
                         ca:{dh:-5.8,ds:-22.7},
                         tg:{dh:-5.8,ds:-22.7},
                         gt:{dh:-6.5,ds:-22.4},
                         ac:{dh:-6.5,ds:-22.4},
                         ct:{dh:-7.8,ds:-21},
                         ag:{dh:-7.8,ds:-21},
                         ga:{dh:-5.6,ds:-22.2},
                         tc:{dh:-5.6,ds:-22.2},
                         cg:{dh:-11.9,ds:-27.2},
                         gc:{dh:-11.1,ds:-24.4},
                         gg:{dh:-11,ds:-19.9},
                         cc:{dh:-11,ds:-19.9}}

"""

r_TmCalc_prototype_nnBr = {
    'aa': {'dh': -9.1, 'ds': -22.2},
    'tt': {'dh': -9.1, 'ds': -22.2},
    'at': {'dh': -8.6, 'ds': -20.4},
    'ta': {'dh': -6.0, 'ds': -21.3},
    'ca': {'dh': -5.8, 'ds': -22.7},
    'tg': {'dh': -5.8, 'ds': -22.7},
    'gt': {'dh': -6.5, 'ds': -22.4},
    'ac': {'dh': -6.5, 'ds': -22.4},
    'ct': {'dh': -7.8, 'ds': -21.0},
    'ag': {'dh': -7.8, 'ds': -21.0},
    'ga': {'dh': -5.6, 'ds': -22.2},
    'tc': {'dh': -5.6, 'ds': -22.2},
    'cg': {'dh': -11.9, 'ds': -27.2},
    'gc': {'dh': -11.1, 'ds': -24.4},
    'gg': {'dh': -11.0, 'ds': -19.9},
    'cc': {'dh': -11.0, 'ds': -19.9}}

prop_seqs = {}

for key, td_dict in r_TmCalc_prototype_nnBr.items():
    prop_seqs[f"{key}/{rc(key)[::-1]}"] = td_dict["dh"], td_dict["ds"]

tablekeys = """\
AA/TT
AT/TA
TA/AT
CA/GT
GT/CA
CT/GA
GA/CT
CG/GC
GC/CG
GG/CC""".splitlines()

for key in tablekeys:
    print(key, prop_seqs[key.lower()][0], prop_seqs[key.lower()][1])


# https://github.com/biopython/biopython/blob/af00cf6c79887d80df8673e5cacde0786415ce34/Bio/SeqUtils/MeltingTemp.py#L176

"""
DNA_NN3 = {
    "init": (0, 0), "init_A/T": (2.3, 4.1), "init_G/C": (0.1, -2.8),
    "init_oneG/C": (0, 0), "init_allA/T": (0, 0), "init_5T/A": (0, 0),
    "sym": (0, -1.4),
    "AA/TT": (-7.9, -22.2), "AT/TA": (-7.2, -20.4), "TA/AT": (-7.2, -21.3),
    "CA/GT": (-8.5, -22.7), "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2), "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9)}

"""
