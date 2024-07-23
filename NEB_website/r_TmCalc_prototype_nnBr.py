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

newtable = {}
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
    newtable[key] = prop_seqs[key.lower()][0], prop_seqs[key.lower()][1]


# https://github.com/biopython/biopython/blob/af00cf6c79887d80df8673e5cacde0786415ce34/Bio/SeqUtils/MeltingTemp.py#L176

from Bio.SeqUtils import MeltingTemp as mt

nn_tables = {"DNA_NN1": mt.DNA_NN1,
             "DNA_NN2": mt.DNA_NN2,
             "DNA_NN3": mt.DNA_NN3,
             "DNA_NN4": mt.DNA_NN4}

start = len(nn_tables) + 1
new_tables = {}
for i, (nn_table_name, nn_table) in enumerate(nn_tables.items()):
    nt = nn_table.copy()
    nt.update(newtable)
    newname = f"DNA_NN{start+i}"
    new_tables[newname] = nt

    print(newname, "=", nt)
