#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the pydna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

import collections as _collections
from Bio.SeqUtils.MeltingTemp import Tm_NN
from Bio.SeqUtils.MeltingTemp import DNA_NN1
import math as _math
from pydna.tm import tmbresluc

dHBr = _collections.defaultdict(dict)
dSBr = _collections.defaultdict(dict)

dHBr[0][0] = -9100
dSBr[0][0] = -24
dHBr[0][1] = -7633.3
dSBr[0][1] = -20.7
dHBr[0][2] = -6500
dSBr[0][2] = -17.3
dHBr[0][3] = -8500
dSBr[0][3] = -22.9
dHBr[0][6] = -7800
dSBr[0][6] = -20.8
dHBr[0][7] = -8066.7
dSBr[0][7] = -21.7
dHBr[0][10] = -8200
dSBr[0][10] = -22.4
dHBr[0][12] = -7800
dSBr[0][12] = -20.7
dHBr[0][13] = -8000
dSBr[0][13] = -21.5
dHBr[0][17] = -8450
dSBr[0][17] = -22.4
dHBr[0][18] = -7150
dSBr[0][18] = -19.1
dHBr[0][19] = -8600
dSBr[0][19] = -23.9
dHBr[0][21] = -7800
dSBr[0][21] = -20.7
dHBr[0][22] = -8850
dSBr[0][22] = -24
dHBr[0][23] = -8000
dSBr[0][23] = -21.5
dHBr[0][24] = -7550
dSBr[0][24] = -20.6
dHBr[1][0] = -5800
dSBr[1][0] = -14.4
dHBr[1][1] = -8866.7
dSBr[1][1] = -21.8
dHBr[1][2] = -9233.3
dSBr[1][2] = -22.3
dHBr[1][3] = -7722.2
dSBr[1][3] = -19.2
dHBr[1][6] = -9566.7
dSBr[1][6] = -22.4
dHBr[1][7] = -7611.1
dSBr[1][7] = -19.1
dHBr[1][10] = -8683.3
dSBr[1][10] = -21.6
dHBr[1][12] = -7516.7
dSBr[1][12] = -18.4
dHBr[1][13] = -8100
dSBr[1][13] = -20
dHBr[1][17] = -7683.3
dSBr[1][17] = -18.4
dHBr[1][18] = -9400
dSBr[1][18] = -22.4
dHBr[1][19] = -7800
dSBr[1][19] = -20.7
dHBr[1][21] = -8200
dSBr[1][21] = -19.7
dHBr[1][22] = -6800
dSBr[1][22] = -17.6
dHBr[1][23] = -8100
dSBr[1][23] = -20
dHBr[1][24] = -8516.7
dSBr[1][24] = -21.5
dHBr[2][0] = -5800
dSBr[2][0] = -12.9
dHBr[2][1] = -10233.3
dSBr[2][1] = -25.1
dHBr[2][2] = -11000
dSBr[2][2] = -26.6
dHBr[2][3] = -8500
dSBr[2][3] = -20.5
dHBr[2][6] = -11900
dSBr[2][6] = -27.8
dHBr[2][7] = -8200
dSBr[2][7] = -20.1
dHBr[2][10] = -9850
dSBr[2][10] = -24.3
dHBr[2][12] = -8400
dSBr[2][12] = -19.8
dHBr[2][13] = -9125
dSBr[2][13] = -22
dHBr[2][17] = -8850
dSBr[2][17] = -20.4
dHBr[2][18] = -11450
dSBr[2][18] = -27.2
dHBr[2][19] = -7800
dSBr[2][19] = -20.8
dHBr[2][21] = -9566.7
dSBr[2][21] = -22.4
dHBr[2][22] = -6800
dSBr[2][22] = -16.9
dHBr[2][23] = -9125
dSBr[2][23] = -22
dHBr[2][24] = -9400
dSBr[2][24] = -23.7
dHBr[3][0] = -6900
dSBr[3][0] = -18.1
dHBr[3][1] = -8000
dSBr[3][1] = -20.3
dHBr[3][2] = -7733.3
dSBr[3][2] = -19.2
dHBr[3][3] = -7722.2
dSBr[3][3] = -20
dHBr[3][6] = -8200
dSBr[3][6] = -20.1
dHBr[3][7] = -7566.7
dSBr[3][7] = -19.7
dHBr[3][10] = -8133.3
dSBr[3][10] = -20.9
dHBr[3][12] = -7316.7
dSBr[3][12] = -18.7
dHBr[3][13] = -7725
dSBr[3][13] = -19.8
dHBr[3][17] = -7550
dSBr[3][17] = -19.1
dHBr[3][18] = -7966.7
dSBr[3][18] = -19.6
dHBr[3][19] = -8066.7
dSBr[3][19] = -21.7
dHBr[3][21] = -7611.1
dSBr[3][21] = -19.1
dHBr[3][22] = -7483.3
dSBr[3][22] = -19.9
dHBr[3][23] = -7725
dSBr[3][23] = -19.8
dHBr[3][24] = -7900
dSBr[3][24] = -20.5
dHBr[6][0] = -5600
dSBr[6][0] = -13.5
dHBr[6][1] = -9533.3
dSBr[6][1] = -23.5
dHBr[6][2] = -11100
dSBr[6][2] = -26.7
dHBr[6][3] = -7700
dSBr[6][3] = -19.1
dHBr[6][6] = -11000
dSBr[6][6] = -26.6
dHBr[6][7] = -7733.3
dSBr[6][7] = -19.2
dHBr[6][10] = -8750
dSBr[6][10] = -22
dHBr[6][12] = -8350
dSBr[6][12] = -20.1
dHBr[6][13] = -8550
dSBr[6][13] = -21
dHBr[6][17] = -8300
dSBr[6][17] = -20.1
dHBr[6][18] = -11050
dSBr[6][18] = -26.7
dHBr[6][19] = -6500
dSBr[6][19] = -17.3
dHBr[6][21] = -9233.3
dSBr[6][21] = -22.3
dHBr[6][22] = -6050
dSBr[6][22] = -15.4
dHBr[6][23] = -8550
dSBr[6][23] = -21
dHBr[6][24] = -8800
dSBr[6][24] = -22
dHBr[7][0] = -6966.7
dSBr[7][0] = -17.9
dHBr[7][1] = -8233.3
dSBr[7][1] = -20.8
dHBr[7][2] = -7700
dSBr[7][2] = -19.1
dHBr[7][3] = -7988.9
dSBr[7][3] = -20.4
dHBr[7][6] = -8500
dSBr[7][6] = -20.5
dHBr[7][7] = -7722.2
dSBr[7][7] = -20
dHBr[7][10] = -8500
dSBr[7][10] = -21.7
dHBr[7][12] = -7333.3
dSBr[7][12] = -18.5
dHBr[7][13] = -7916.7
dSBr[7][13] = -20.1
dHBr[7][17] = -7733.3
dSBr[7][17] = -19.2
dHBr[7][18] = -8100
dSBr[7][18] = -19.8
dHBr[7][19] = -8500
dSBr[7][19] = -22.9
dHBr[7][21] = -7722.2
dSBr[7][21] = -19.2
dHBr[7][22] = -7733.3
dSBr[7][22] = -20.4
dHBr[7][23] = -7916.7
dSBr[7][23] = -20.1
dHBr[7][24] = -8100
dSBr[7][24] = -21
dHBr[10][0] = -5800
dSBr[10][0] = -15.2
dHBr[10][1] = -8183.3
dSBr[10][1] = -20.2
dHBr[10][2] = -8350
dSBr[10][2] = -20.1
dHBr[10][3] = -7333.3
dSBr[10][3] = -18.5
dHBr[10][6] = -8400
dSBr[10][6] = -19.8
dHBr[10][7] = -7316.7
dSBr[10][7] = -18.7
dHBr[10][10] = -8100
dSBr[10][10] = -20.2
dHBr[10][12] = -7075
dSBr[10][12] = -17.7
dHBr[10][13] = -7587.5
dSBr[10][13] = -18.9
dHBr[10][17] = -7100
dSBr[10][17] = -17.5
dHBr[10][18] = -8375
dSBr[10][18] = -19.9
dHBr[10][19] = -7800
dSBr[10][19] = -20.7
dHBr[10][21] = -7516.7
dSBr[10][21] = -18.4
dHBr[10][22] = -6800
dSBr[10][22] = -17.9
dHBr[10][23] = -7587.5
dSBr[10][23] = -18.9
dHBr[10][24] = -8075
dSBr[10][24] = -20.4
dHBr[12][0] = -7450
dSBr[12][0] = -18.5
dHBr[12][1] = -8933.3
dSBr[12][1] = -22.9
dHBr[12][2] = -8750
dSBr[12][2] = -22
dHBr[12][3] = -8500
dSBr[12][3] = -21.7
dHBr[12][6] = -9850
dSBr[12][6] = -24.3
dHBr[12][7] = -8133.3
dSBr[12][7] = -20.9
dHBr[12][10] = -9025
dSBr[12][10] = -23.3
dHBr[12][12] = -8100
dSBr[12][12] = -20.2
dHBr[12][13] = -8562.5
dSBr[12][13] = -21.8
dHBr[12][17] = -8650
dSBr[12][17] = -21.4
dHBr[12][18] = -9300
dSBr[12][18] = -23.1
dHBr[12][19] = -8200
dSBr[12][19] = -22.4
dHBr[12][21] = -8683.3
dSBr[12][21] = -21.6
dHBr[12][22] = -7825
dSBr[12][22] = -20.4
dHBr[12][23] = -8562.5
dSBr[12][23] = -21.8
dHBr[12][24] = -8475
dSBr[12][24] = -22.2
dHBr[13][0] = -6625
dSBr[13][0] = -16.8
dHBr[13][1] = -8558.3
dSBr[13][1] = -21.5
dHBr[13][2] = -8550
dSBr[13][2] = -21
dHBr[13][3] = -7916.7
dSBr[13][3] = -20.1
dHBr[13][6] = -9125
dSBr[13][6] = -22
dHBr[13][7] = -7725
dSBr[13][7] = -19.8
dHBr[13][10] = -8562.5
dSBr[13][10] = -21.8
dHBr[13][12] = -7587.5
dSBr[13][12] = -18.9
dHBr[13][13] = -8075
dSBr[13][13] = -20.3
dHBr[13][17] = -7875
dSBr[13][17] = -19.4
dHBr[13][18] = -8837.5
dSBr[13][18] = -21.5
dHBr[13][19] = -8000
dSBr[13][19] = -21.5
dHBr[13][21] = -8100
dSBr[13][21] = -20
dHBr[13][22] = -7312.5
dSBr[13][22] = -19.2
dHBr[13][23] = -8075
dSBr[13][23] = -20.3
dHBr[13][24] = -8275
dSBr[13][24] = -21.3
dHBr[17][0] = -7350
dSBr[17][0] = -18.8
dHBr[17][1] = -8583.3
dSBr[17][1] = -22.1
dHBr[17][2] = -8800
dSBr[17][2] = -22
dHBr[17][3] = -8100
dSBr[17][3] = -21
dHBr[17][6] = -9400
dSBr[17][6] = -23.7
dHBr[17][7] = -7900
dSBr[17][7] = -20.5
dHBr[17][10] = -8475
dSBr[17][10] = -22.2
dHBr[17][12] = -8075
dSBr[17][12] = -20.4
dHBr[17][13] = -8275
dSBr[17][13] = -21.3
dHBr[17][17] = -8375
dSBr[17][17] = -21.2
dHBr[17][18] = -9100
dSBr[17][18] = -22.9
dHBr[17][19] = -7550
dSBr[17][19] = -20.6
dHBr[17][21] = -8516.7
dSBr[17][21] = -21.5
dHBr[17][22] = -7450
dSBr[17][22] = -19.7
dHBr[17][23] = -8275
dSBr[17][23] = -21.3
dHBr[17][24] = -8175
dSBr[17][24] = -21.3
dHBr[18][0] = -5700
dSBr[18][0] = -13.2
dHBr[18][1] = -9883.3
dSBr[18][1] = -24.3
dHBr[18][2] = -11050
dSBr[18][2] = -26.7
dHBr[18][3] = -8100
dSBr[18][3] = -19.8
dHBr[18][6] = -11450
dSBr[18][6] = -27.2
dHBr[18][7] = -7966.7
dSBr[18][7] = -19.6
dHBr[18][10] = -9300
dSBr[18][10] = -23.1
dHBr[18][12] = -8375
dSBr[18][12] = -19.9
dHBr[18][13] = -8837.5
dSBr[18][13] = -21.5
dHBr[18][17] = -8575
dSBr[18][17] = -20.2
dHBr[18][18] = -11250
dSBr[18][18] = -26.9
dHBr[18][19] = -7150
dSBr[18][19] = -19.1
dHBr[18][21] = -9400
dSBr[18][21] = -22.4
dHBr[18][22] = -6425
dSBr[18][22] = -16.1
dHBr[18][23] = -8837.5
dSBr[18][23] = -21.5
dHBr[18][24] = -9100
dSBr[18][24] = -22.9
dHBr[19][0] = -6000
dSBr[19][0] = -16.9
dHBr[19][1] = -6833.3
dSBr[19][1] = -16.8
dHBr[19][2] = -5600
dSBr[19][2] = -13.5
dHBr[19][3] = -6966.7
dSBr[19][3] = -17.9
dHBr[19][6] = -5800
dSBr[19][6] = -12.9
dHBr[19][7] = -6900
dSBr[19][7] = -18.1
dHBr[19][10] = -7450
dSBr[19][10] = -18.5
dHBr[19][12] = -5800
dSBr[19][12] = -15.2
dHBr[19][13] = -6625
dSBr[19][13] = -16.8
dHBr[19][17] = -5900
dSBr[19][17] = -14.9
dHBr[19][18] = -5700
dSBr[19][18] = -13.2
dHBr[19][19] = -9100
dSBr[19][19] = -24
dHBr[19][21] = -5800
dSBr[19][21] = -14.4
dHBr[19][22] = -7550
dSBr[19][22] = -20.5
dHBr[19][23] = -6625
dSBr[19][23] = -16.8
dHBr[19][24] = -7350
dSBr[19][24] = -18.8
dHBr[21][0] = -6833.3
dSBr[21][0] = -16.8
dHBr[21][1] = -9133.3
dSBr[21][1] = -23.1
dHBr[21][2] = -9533.3
dSBr[21][2] = -23.5
dHBr[21][3] = -8233.3
dSBr[21][3] = -20.8
dHBr[21][6] = -10233.3
dSBr[21][6] = -25.1
dHBr[21][7] = -8000
dSBr[21][7] = -20.3
dHBr[21][10] = -8933.3
dSBr[21][10] = -22.9
dHBr[21][12] = -8183.3
dSBr[21][12] = -20.2
dHBr[21][13] = -8558.3
dSBr[21][13] = -21.5
dHBr[21][17] = -8533.3
dSBr[21][17] = -20.9
dHBr[21][18] = -9883.3
dSBr[21][18] = -24.3
dHBr[21][19] = -7633.3
dSBr[21][19] = -20.7
dHBr[21][21] = -8866.7
dSBr[21][21] = -21.8
dHBr[21][22] = -7233.3
dSBr[21][22] = -18.7
dHBr[21][23] = -8558.3
dSBr[21][23] = -21.5
dHBr[21][24] = -8583.3
dSBr[21][24] = -22.1
dHBr[22][0] = -7550
dSBr[22][0] = -20.5
dHBr[22][1] = -7233.3
dSBr[22][1] = -18.7
dHBr[22][2] = -6050
dSBr[22][2] = -15.4
dHBr[22][3] = -7733.3
dSBr[22][3] = -20.4
dHBr[22][6] = -6800
dSBr[22][6] = -16.9
dHBr[22][7] = -7483.3
dSBr[22][7] = -19.9
dHBr[22][10] = -7825
dSBr[22][10] = -20.4
dHBr[22][12] = -6800
dSBr[22][12] = -17.9
dHBr[22][13] = -7312.5
dSBr[22][13] = -19.2
dHBr[22][17] = -7175
dSBr[22][17] = -18.7
dHBr[22][18] = -6425
dSBr[22][18] = -16.1
dHBr[22][19] = -8850
dSBr[22][19] = -24
dHBr[22][21] = -6800
dSBr[22][21] = -17.6
dHBr[22][22] = -8200
dSBr[22][22] = -22.2
dHBr[22][23] = -7312.5
dSBr[22][23] = -19.2
dHBr[22][24] = -7450
dSBr[22][24] = -19.7
dHBr[23][0] = -6625
dSBr[23][0] = -16.8
dHBr[23][1] = -8558.3
dSBr[23][1] = -21.5
dHBr[23][2] = -8550
dSBr[23][2] = -21
dHBr[23][3] = -7916.7
dSBr[23][3] = -20.1
dHBr[23][6] = -9125
dSBr[23][6] = -22
dHBr[23][7] = -7725
dSBr[23][7] = -19.8
dHBr[23][10] = -8562.5
dSBr[23][10] = -21.8
dHBr[23][12] = -7587.5
dSBr[23][12] = -18.9
dHBr[23][13] = -8075
dSBr[23][13] = -20.3
dHBr[23][17] = -7875
dSBr[23][17] = -19.4
dHBr[23][18] = -8837.5
dSBr[23][18] = -21.5
dHBr[23][19] = -8000
dSBr[23][19] = -21.5
dHBr[23][21] = -8100
dSBr[23][21] = -20
dHBr[23][22] = -7312.5
dSBr[23][22] = -19.2
dHBr[23][23] = -8075
dSBr[23][23] = -20.3
dHBr[23][24] = -8275
dSBr[23][24] = -21.3
dHBr[24][0] = -5900
dSBr[24][0] = -14.9
dHBr[24][1] = -8533.3
dSBr[24][1] = -20.9
dHBr[24][2] = -8300
dSBr[24][2] = -20.1
dHBr[24][3] = -7733.3
dSBr[24][3] = -19.2
dHBr[24][6] = -8850
dSBr[24][6] = -20.4
dHBr[24][7] = -7550
dSBr[24][7] = -19.1
dHBr[24][10] = -8650
dSBr[24][10] = -21.4
dHBr[24][12] = -7100
dSBr[24][12] = -17.5
dHBr[24][13] = -7875
dSBr[24][13] = -19.4
dHBr[24][17] = -7375
dSBr[24][17] = -17.6
dHBr[24][18] = -8575
dSBr[24][18] = -20.2
dHBr[24][19] = -8450
dSBr[24][19] = -22.4
dHBr[24][21] = -7683.3
dSBr[24][21] = -18.4
dHBr[24][22] = -7175
dSBr[24][22] = -18.7
dHBr[24][23] = -7875
dSBr[24][23] = -19.4
dHBr[24][24] = -8375
dSBr[24][24] = -21.2

primer = "AGTCTAGTCTGTGTAGTTTCGACTAGTCTATCG"
pydna_tm = tmbresluc(primer)
assert tmbresluc(primer) == 63.38496307044147
print(pydna_tm, "pydna_tm")

primerc = 500.0
saltc = 50

saltc = float(saltc) / 1000
pri = primerc * 1e-9

dS = -12.4
dH = -3400

# dS = -16.8
# dH = 0

STR = primer.lower()

for i in range(len(STR) - 1):
    n1 = ord(STR[i])
    n2 = ord(STR[i + 1])
    dH +=  dHBr[n1 - 97][n2 - 97]
    dS +=  dSBr[n1 - 97][n2 - 97]

tm = (dH / (1.987 * _math.log(pri / 1600) + dS)) - 273.15

# sallt correction Schildkraut 1965 = biopython method 1
corr =  (16.6 * _math.log(saltc)) / _math.log(10)

tm += corr
print(tm, "fermentas/pydna")

newdata = {'init': (-3.4, -12.4)}

DNA_NN1.update(newdata)

bp = Tm_NN(primer,
           check=True,
           strict=True,
           c_seq=None,
           shift=0,
           nn_table=DNA_NN1,
           tmm_table=None,
           imm_table=None,
           de_table=None,
           dnac1=500/1600,
           dnac2=0,
           selfcomp=False,
           Na=50,
           K=0,
           Tris=0,
           Mg=0,
           dNTPs=0,
           saltcorr=1)
print(bp, "biopython")
