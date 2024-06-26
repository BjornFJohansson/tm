from Bio.SeqUtils import MeltingTemp as mt

"""
OneTaq DNA polymerase 

                                    Tm
P1fwd	 AGCGGATAACAATTTCACACAGGA	58
P1rev	 GTAAAACGACGGCCAGT	        54
pBeta	 AGCGGATAACAATTTCAC	        48
P3fwd	 AGCGGATAAGGGCAATTTCAC	    57
P3rev	 GTAAAACGACGGCCA	        50

Tm values are calculated using thermodynamic data from Santa Lucia
salt correction of Owczarzy
For Phusion® DNA Polymerases, the salt correction of Schildkraut.

SantaLucia (1998) PNAS 95:1460-5
Owczarzy et al (2004) Biochem 43:3537-54


"""


p = P1fwd = "AGCGGATAACAATTTCACACAGGA"


"""

1X OneTaq® Standard Reaction Buffer Pack
20 mM Tris-HCl
22 mM NH4Cl
22 mM KCl
1.8 mM MgCl2
0.06% IGEPAL® CA-630
0.05% Tween® 20
(pH 8.9 @ 25°C)

"""


"""

Salt correction

1. 16.6 x log[Na+]
   (Schildkraut & Lifson (1965), Biopolymers 3: 195-208)
2. 16.6 x log([Na+]/(1.0 + 0.7*[Na+]))
   (Wetmur (1991), Crit Rev Biochem Mol Biol 126: 227-259)
3. 12.5 x log(Na+]
   (SantaLucia et al. (1996), Biochemistry 35: 3555-3562
4. 11.7 x log[Na+]
   (SantaLucia (1998), Proc Natl Acad Sci USA 95: 1460-1465
5. Correction for deltaS: 0.368 x (N-1) x ln[Na+]
   (SantaLucia (1998), Proc Natl Acad Sci USA 95: 1460-1465)
6. (4.29(%GC)-3.95)x1e-5 x ln[Na+] + 9.40e-6 x ln[Na+]^2
   (Owczarzy et al. (2004), Biochemistry 43: 3537-3554)
7. Complex formula with decision tree and 7 empirical constants.
   Mg2+ is corrected for dNTPs binding (if present)
   (Owczarzy et al. (2008), Biochemistry 47: 5336-5353)
"""


saltcorr_dict = { 1: "Schildkraut 1965 ",
                  2: "Wetmur 1991      ",
                  3: "SantaLucia 1996  ",
                  4: "SantaLucia 1998  ",
                  5: "SantaLucia 1998  ",
                  6: "Owczarzy 2004    ",
                  7: "Owczarzy 2008    ", }

table_dict = {
str(mt.DNA_NN1): "Breslauer 1986             ", 
str(mt.DNA_NN2): "Sugimoto 1996              ", 
str(mt.DNA_NN3): "Allawi and SantaLucia 1997 ", 
str(mt.DNA_NN4): "SantaLucia & Hicks 2004    ",
}



mt.Tm_NN(p, 
         nn_table=mt.DNA_NN3,
         Na=5,
         Tris=20,
         Mg=1.8, 
         dnac1=25,    # nM
         dnac2=25,    # nM
         dNTPs=0, 
         saltcorr=6)


# for nn_table in (mt.DNA_NN1, mt.DNA_NN2, mt.DNA_NN3, mt.DNA_NN4):
#     for saltcorr in range(1, 8):
#         print(f"table {table_dict[str(nn_table)]} saltcorr {saltcorr_dict[saltcorr]}", mt.Tm_NN(
#                p, 
#                nn_table=mt.DNA_NN3,
#                Na=22,
#                Tris=20.0/2,
#                Mg=1.8, 
#                dnac1=200/2, 
#                dnac2=200/2, 
#                dNTPs=0.8, 
#                saltcorr=saltcorr) )




