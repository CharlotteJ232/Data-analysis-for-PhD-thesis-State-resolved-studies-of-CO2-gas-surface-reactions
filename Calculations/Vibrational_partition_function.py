import numpy as np
from matplotlib import pyplot as plt

#https://pubs.aip.org/jcp/article/37/7/1442/77057/Vibrational-Energies-of-the-CO2-Molecule
wn_v1 = 1355
wn_v2 = 673
wn_v3 = 2396.5

#https://doi.org/10.1006/JMSP.2000.8237
wn_v1 = (1285+1388)/2 #fermi resonance/dyad/polyad
wn_v2 = 667
wn_v3 = 2349

E_v1 = wn_v1
E_v2 = wn_v2
E_v3 = wn_v3
kb = 0.695034800 #cm-1/K

max_v1 = 50
max_v2 = 50
max_v3 = 50
T = 300 #K

thresh = 1E-5


def main():
    state=[]
    energy=[]
    degeneracy=[]
    for v1 in range(max_v1):
        for v2 in range(max_v2):
            for v3 in range(max_v3):
                E_tot = v1*E_v1+v2*E_v2+v3*E_v3
                degen = calc_degen(v2)
                state.append([v1, v2, v3])
                energy.append(E_tot)
                degeneracy.append(degen)
    energy=np.array(energy)                
    degeneracy=np.array(degeneracy)

    # print(f'{energy=}')
    boltz = boltzmann(energy, degeneracy) 
    boltz /= np.sum(boltz)

    plt.plot(energy, boltz, 'o')
    plt.axis([0,2500, 0, None])
    plt.xlabel('Wavenumber (cm$-1$')
    plt.ylabel('Population')
    plt.show()
    plt.close()


    index_sort = np.flip(np.argsort(boltz))
    boltz = boltz[index_sort]
    state = np.array(state)
    state = state[index_sort]

    above_thresh = boltz > thresh
    # print(above_thresh)
    above_thresh = np.transpose(np.argwhere(above_thresh))[0]

    # print(above_thresh.dtype)
    # print(above_thresh)

    for index in above_thresh:
        print(state[index], '{:.1e}'.format(boltz[index]))

    for index in above_thresh:
        table_row = ''
        for i in range(len(state[index])):
            table_row += str(state[index][i]) + '&'
        table_row += '{:.1e}'.format(boltz[index])+'\\\\'
        print(table_row)
       
       
 # print(np.asarray(state)[above_thresh], np.transpose(boltz[above_thresh]))

def calc_degen(n_vib):
    """
    for original degeneracy of 2
    degeneracy = n_vib+1
    because possible quanta in state 1 is n_vib+1 and number in state 2 is always determined by number in state 1
    """
    return n_vib+1

def boltzmann(E, degen, A=1):
    return A*degen*np.exp(-E/kb/T)


main()