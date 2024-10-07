from pydna.tm import tm_default
print(tm_default("agatcgactatctacttatgcatctta"))  # 58.90561619071127
print(tm_default("agatcgactatctacttatgcatctt"))  # 59.1327553388262

from Bio.SeqUtils import MeltingTemp
