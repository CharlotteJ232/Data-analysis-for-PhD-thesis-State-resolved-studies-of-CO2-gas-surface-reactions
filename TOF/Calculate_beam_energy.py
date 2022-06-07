import numpy as np
from matplotlib import pyplot as plt

name1 = 'He' 
m1 = 4 #amu
x1 = 25 #ml/min

name2 = 'CO2'
m2 = 44 #amu
x2 = 1 #ml/min

T = 300 #K

R = 8.314
kjmol_to_mev = 10.364

def calc_E(m, factor=5/2):
    """
    Returns E in kJ/mol
    """
    return factor*R*T*m/(x1*m1+x2*m2)*(x1+x2)/1000

E1 = calc_E(m1)
E2 = calc_E(m2)

print('E_'+name1+' = '+str(E1)+' kJ/mol')
print('E_'+name2+' = '+str(E2)+' kJ/mol')

E1_mev = E1 * kjmol_to_mev
E2_mev = E2 * kjmol_to_mev

print('E_'+name1+' = '+str(E1_mev)+' meV')
print('E_'+name2+' = '+str(E2_mev)+' meV')

xmin = 0
xmax = 0.1
x = np.linspace(xmin,xmax,100,endpoint=True)

def calc_E_range(m, x, factor=5/2):
    """
    Returns E in kJ/mol
    """
    return factor*R*T*m/(x*m1+(1-x)*m2)/1000



plt.plot(x, calc_E_range(m1, x), label=name1)
plt.legend()
plt.title('T = '+str(T)+' K')
plt.ylabel(name1+' Energy (kJ/mol)')
plt.xlabel('Fraction of '+name1+' in '+name2)
plt.axis([xmin,xmax,None,None])
plt.grid()
plt.show()
plt.close()

plt.plot(x, calc_E_range(m2, 1-x), label=name2)
plt.legend()
plt.title('T = '+str(T)+' K')
plt.ylabel(name2+' Energy (kJ/mol)')
plt.xlabel('Fraction of '+name2+' in '+name1)
plt.axis([xmin,xmax,None,None])
plt.grid()
plt.show()
plt.close()
