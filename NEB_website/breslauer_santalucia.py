
dSBr = {
    'aa': -24,
    'tt': -24,
    'at': -23.9,
    'ta': -16.9,
    'ca': -12.9,
    'tg': -12.9,
    'gt': -17.3,
    'ac': -17.3,
    'ct': -20.8,
    'ag': -20.8,
    'ga': -13.5,
    'tc': -13.5,
    'cg': -27.8,
    'gc': -26.7,
    'gg': -26.6,
    'cc': -26.6
}

dHBr = {
    'aa': -9.1,
    'tt': -9.1,
    'at': -8.6,
    'ta': -6,
    'ca': -5.8,
    'tg': -5.8,
    'gt': -6.5,
    'ac': -6.5,
    'ct': -7.8,
    'ag': -7.8,
    'ga': -5.6,
    'tc': -5.6,
    'cg': -11.9,
    'gc': -11.1,
    'gg': -11,
    'cc': -11
}

dSSa = {
    'aa': -22.2,
    'tt': -22.2,
    'at': -20.4,
    'ta': -21.3,
    'ca': -22.7,
    'tg': -22.7,
    'gt': -22.4,
    'ac': -22.4,
    'ct': -21,
    'ag': -21,
    'ga': -22.2,
    'tc': -22.2,
    'cg': -27.2,
    'gc': -24.4,
    'gg': -19.9,
    'cc': -19.9
}



dHSa = {
    'aa': -7.9,
    'tt': -7.9,
    'at': -7.2,
    'ta': -7.2,
    'ca': -8.5,
    'tg': -8.5,
    'gt': -8.4,
    'ac': -8.4,
    'ct': -7.8,
    'ag': -7.8,
    'ga': -8.2,
    'tc': -8.2,
    'cg': -10.6,
    'gc': -9.8,
    'gg': -8,
    'cc': -8
}


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

newtable = {}
for key in tablekeys:
    newtable[key] = dHBr[key.lower()[0:2]], dSBr[key.lower()[0:2]]


newtable2 = {}
for key in tablekeys:
    newtable2[key] = dHSa[key.lower()[0:2]], dSSa[key.lower()[0:2]]
