#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copied from main-3d92a74abb.js

Processed with https://www.codeconvert.ai/javascript-to-python-converter

r.TmCalc.prototype.saltCorrect = function() {
        var e, t, n, r;
        if (0 === this.saltc && 0 === this.disaltc) throw {
          name: "Missing salt",
          message: "Unable to calculate salt correction because mono and divalent salt are missing."
        };
        this.sc_sch = 16.6 * Math.log(this.saltc) / Math.LN10, this.sc_sl = .368 * this.wseq.length * Math.log(this.saltc), t = 3.92, n = 1.42, r = 8.31,
          (e = this.saltc > 0 ? Math.sqrt(this.disaltc) / this.saltc : -1) >= 0 && e < .22 ? this.sc_ow = 1e-5 * (4.29 * this.fgc - 3.95) * Math.log(this
            .saltc) + 94e-7 * Math.log(this.saltc) * Math.log(this.saltc) : (e >= .22 && e < 6 && (t = 3.92 * (.843 - .352 * Math.sqrt(this.saltc) * Math
            .log(this.saltc)), n = 1.42 * (1.279 - .00403 * Math.log(this.saltc) - .00803 * Math.log(this.saltc) * Math.log(this.saltc)), r = 8.31 * (
            .486 - .258 * Math.log(this.saltc) + .00525 * Math.log(this.saltc) * Math.log(this.saltc) * Math.log(this.saltc))), this.sc_ow = 1e-5 * (t +
            -.911 * Math.log(this.disaltc) + this.fgc * (6.26 + n * Math.log(this.disaltc)) + .5 / (this.wseq.length - 1) * (52.5 * Math.log(this
              .disaltc) - 48.2 + r * Math.log(this.disaltc) * Math.log(this.disaltc))))
      },
"""

# https://github.com/biopython/biopython/blob/af00cf6c79887d80df8673e5cacde0786415ce34/Bio/SeqUtils/MeltingTemp.py#L560

import math
from Bio import SeqUtils


from Bio.SeqUtils.MeltingTemp import salt_correction

seq = "AGCGGATAACAATTTCACACAGGA"

bp = salt_correction(Na=0, K=50, Tris=0, Mg=1.5, dNTPs=0, method=7, seq=seq)
print(bp, "Biopython #7")

class TmCalc:
    def __init__(self, seq):
        self.saltc = 50 / 1000
        self.disaltc = 1.5/1000
        self.wseq = seq
        self.fgc = SeqUtils.gc_fraction(seq, "ignore")

    def saltCorrect(self):
        if self.saltc == 0 and self.disaltc == 0:
            raise ValueError("Missing salt: Unable to calculate salt correction because mono and divalent salt are missing.")

        self.sc_sch = 16.6 * math.log(self.saltc) / math.log(10)
        self.sc_sl = 0.368 * len(self.wseq) * math.log(self.saltc)

        t = 3.92
        n = 1.42
        r = 8.31

        e = math.sqrt(self.disaltc) / self.saltc if self.saltc > 0 else -1

        if 0 <= e < 0.22:
            print("e", e)
            self.sc_ow = 1e-5 * (4.29 * self.fgc - 3.95) * math.log(self.saltc) + 94e-7 * math.log(self.saltc) * math.log(self.saltc)
        else:
            if 0.22 <= e < 6:
                t = 3.92 * (0.843 - 0.352 * math.sqrt(self.saltc) * math.log(self.saltc))
                n = 1.42 * (1.279 - 0.00403 * math.log(self.saltc) - 0.00803 * math.log(self.saltc) * math.log(self.saltc))
                r = 8.31 * (0.486 - 0.258 * math.log(self.saltc) + 0.00525 * math.log(self.saltc) * math.log(self.saltc) * math.log(self.saltc))

            self.sc_ow = 1e-5 * (t + -0.911 * math.log(self.disaltc) + self.fgc * (6.26 + n * math.log(self.disaltc)) + 0.5 / (len(self.wseq) - 1) * (52.5 * math.log(self.disaltc) - 48.2 + r * math.log(self.disaltc) * math.log(self.disaltc)))

        return self.sc_ow

sc = TmCalc(seq).saltCorrect()
print(sc, "r.TmCalc.prototype.saltCorrect")

assert bp==sc
