from prody import *
import numpy as np
import matplotlib.pyplot as plt
import random
confProDy(auto_show=False)
confProDy(auto_secondary=True)
def get_seqs(lst_to_conv):
    primary_seqs = []
    secondary_seqs = []
    phi_seqs = []
    psi_seqs = []
    phipsi_seqs = []
    for aa in lst_to_conv:
        parser = parsePDB(aa)
        if 'A' in parser:
            chain = parser['A']
        else:
            chain = parser
        primary_sequence = chain.ca.getSequence()
        try:
            secondary_structure = chain.ca.getSecstrs()
        except:
            secondary_structure = []
            continue
        Phi = []
        Psi = []
        phipsi = []
        for res in chain.iterResidues():
            try:
                phi = calcPhi(res)
                psi = calcPsi(res)
            except:
                continue
            else:
                Phi.append(phi)
                Psi.append(psi)
                phipsi.append((phi, psi))
        primary_seqs.append(primary_sequence)
        secondary_seqs.append(secondary_structure)
        phi_seqs.append(Phi)
        psi_seqs.append(Psi)
        phipsi_seqs.append(phipsi)
    return primary_seqs, secondary_seqs, phi_seqs, psi_seqs, phipsi_seqs