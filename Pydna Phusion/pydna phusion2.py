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

primer = "AGTCTAGTCTGTGTAGTTTCGACTAGTCTATCG"

tm_pydna = tmbresluc(primer)

assert tm_pydna == 63.38496307044147

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

print("biopython", bp)

print(tm_pydna-bp)
