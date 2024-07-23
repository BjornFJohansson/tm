#!/usr/bin/env python
# coding: utf-8

from Bio.SeqUtils import MeltingTemp as mt
import csv
from scipy.optimize import minimize

primerlist = []
with open('tmcalc_batch.txt', newline='') as csvfile:
    neb_file_reader = csv.reader(csvfile, delimiter='\t')
    for row in neb_file_reader:
        name, sequence, nebtm, *rest = row
        primerlist.append([name, sequence, int(nebtm)])


def fun(salt, nn_table, saltcorr):
    sqs = 0
    for name, primer, nebtm in primerlist:
        bptm = mt.Tm_NN(
            primer,
            nn_table=nn_table,
            Na=salt,  # mM
            K=0,  # mM
            Tris=0,  # mM
            Mg=2.0,
            dnac1=500,  # nM
            dnac2=0,  # nM
            dNTPs=0.8,  # mM
            saltcorr=saltcorr)
        sqs+= (nebtm - bptm) ** 2
    return sqs

DNA_NN5 = {'init': (0, 0), 'init_A/T': (0, 0), 'init_G/C': (0, 0), 'init_oneG/C': (0, -16.8), 'init_allA/T': (0, -20.1), 'init_5T/A': (0, 0), 'sym': (0, -1.3), 'AA/TT': (-9.1, -22.2), 'AT/TA': (-8.6, -20.4), 'TA/AT': (-6.0, -21.3), 'CA/GT': (-5.8, -22.7), 'GT/CA': (-6.5, -22.4), 'CT/GA': (-7.8, -21.0), 'GA/CT': (-5.6, -22.2), 'CG/GC': (-11.9, -27.2), 'GC/CG': (-11.1, -24.4), 'GG/CC': (-11.0, -19.9)}
DNA_NN6 = {'init': (0.6, -9.0), 'init_A/T': (0, 0), 'init_G/C': (0, 0), 'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0), 'sym': (0, -1.4), 'AA/TT': (-9.1, -22.2), 'AT/TA': (-8.6, -20.4), 'TA/AT': (-6.0, -21.3), 'CA/GT': (-5.8, -22.7), 'GT/CA': (-6.5, -22.4), 'CT/GA': (-7.8, -21.0), 'GA/CT': (-5.6, -22.2), 'CG/GC': (-11.9, -27.2), 'GC/CG': (-11.1, -24.4), 'GG/CC': (-11.0, -19.9)}
DNA_NN7 = {'init': (0, 0), 'init_A/T': (2.3, 4.1), 'init_G/C': (0.1, -2.8), 'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0), 'sym': (0, -1.4), 'AA/TT': (-9.1, -22.2), 'AT/TA': (-8.6, -20.4), 'TA/AT': (-6.0, -21.3), 'CA/GT': (-5.8, -22.7), 'GT/CA': (-6.5, -22.4), 'CT/GA': (-7.8, -21.0), 'GA/CT': (-5.6, -22.2), 'CG/GC': (-11.9, -27.2), 'GC/CG': (-11.1, -24.4), 'GG/CC': (-11.0, -19.9)}
DNA_NN8 = {'init': (0.2, -5.7), 'init_A/T': (2.2, 6.9), 'init_G/C': (0, 0), 'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0), 'sym': (0, -1.4), 'AA/TT': (-9.1, -22.2), 'AT/TA': (-8.6, -20.4), 'TA/AT': (-6.0, -21.3), 'CA/GT': (-5.8, -22.7), 'GT/CA': (-6.5, -22.4), 'CT/GA': (-7.8, -21.0), 'GA/CT': (-5.6, -22.2), 'CG/GC': (-11.9, -27.2), 'GC/CG': (-11.1, -24.4), 'GG/CC': (-11.0, -19.9)}

nn_tables = {"DNA_NN1": mt.DNA_NN1,
             "DNA_NN2": mt.DNA_NN2,
             "DNA_NN3": mt.DNA_NN3,
             "DNA_NN4": mt.DNA_NN4,
             "DNA_NN5": DNA_NN5,
             "DNA_NN6": DNA_NN6,
             "DNA_NN7": DNA_NN7,
             "DNA_NN8": DNA_NN8}


funcvalues = {}
bnds = ((0.1, 500),)
for nn_table_name, nn_table in nn_tables.items():
    for saltcorr in (1,2,3,4,5,6,7):
        result = minimize(fun, 50, args=(nn_table, saltcorr), bounds=bnds)
        funcvalues[nn_table_name, saltcorr] = result.fun

funcvalues = dict((k, funcvalues[k]) for k in sorted(funcvalues, key=funcvalues.get))

(nn_table, saltcorr), salt = list(funcvalues.items())[0]

print(nn_table, saltcorr, salt)
print()

for name, primer, nebtm in primerlist:
    bptm = mt.Tm_NN(
        primer,
        nn_table=nn_tables[nn_table],
        Na=salt,  # mM
        K=0,  # mM
        Tris=0,  # mM
        Mg=2.0,
        dnac1=500,  # nM
        dnac2=0,  # nM
        dNTPs=0.8,  # mM
        saltcorr=saltcorr,
    )
    if round(bptm) != nebtm:
        print(name, primer, nebtm, round(bptm), nebtm - round(bptm))

print(nn_table, saltcorr, salt)
