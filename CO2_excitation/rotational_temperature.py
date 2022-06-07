
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
folderstart = "C:/Users/Werk/Surfdrive/DATA/"
savefolder = 'CO2_on_CO2/Figures/'
save=False


theta_r = 0.561 #for CO2
J = np.array([0,2,4,6,8,10])
population = np.array([0.051116978, 0.212347685, 0.242282697, 0.231112916, 0.166937449, 0.096202275])
T_guess = 5
A_guess = 0.1
vibrational_ground_state_fraction = 0.95

def boltzmann(J, A, T):
    return A*(2*J+1)*np.exp(-1*(J*(J+1)*theta_r)/(T))

popt = curve_fit(boltzmann, J,population,p0=[A_guess,T_guess])
# print(popt[0])

many_J = np.arange(0,1000,2)
Z = np.sum(boltzmann(many_J, *popt[0]))/vibrational_ground_state_fraction

fig, ax = plt.subplots()

# plt.plot(J, population, 'o', label='Measured data')
# plt.plot(J, boltzmann(J, A_guess, T_guess), label='Guess, T='+str(T_guess)+'K')
ax.plot(J, boltzmann(J, *popt[0])/Z, label='Fit, T='+str(np.round(popt[0][1],2))+'K')
ax.plot(J, population/Z, 'o', label='Measured population (normalized)', color='black')
plt.xlabel('J')
plt.ylabel('Population')
plt.legend(loc='lower right')

ax.tick_params(top=True, direction='in')  
ax.tick_params(right=True, labelleft=True, direction='in')
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='minor', top=True, direction='in')  
ax.tick_params(which='minor', right=True, direction='in')

# plt.title('Fitted T = '+str(np.round(popt[0][1],2))+' K')
plt.axis([0,10,0,0.25])

if save:
    plt.savefig(folderstart+savefolder+'rotational_state_distribution.png',dpi=500)
plt.show()
plt.close()







