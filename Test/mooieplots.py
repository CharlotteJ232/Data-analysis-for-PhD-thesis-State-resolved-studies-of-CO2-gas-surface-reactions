# -*- coding: utf-8 -*-
"""
https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html
https://stackoverflow.com/questions/26106552/matplotlib-style-library-not-updating-when-mplstyle-files-added-deleted
https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html#customizing-with-matplotlibrc-files
@author: Charlotte
"""

import numpy as np #voor rekenen met arrays
from matplotlib import pyplot as plt #voor plotten
from matplotlib import cm #voor colormaps


#maak een testdataset, deze stap is niet nodig als je zelf al data hebt
datasetnrs = [1,2,3,4,5,10]
x = np.arange(-10,10.1,1)
print('x is ', x)
for nr in datasetnrs:
    y = nr * x #y data maak ik een lijn met helling nr
    data = np.column_stack((x,y)) #stop de twee arrays met data bij elkaar om makkelijk op te slaan
    np.savetxt('data'+str(nr)+'.txt', data, header='x, y') #sla op in bestand met naam data1.txt, data2.txt etc

#lees de data weer in (dit is dus niet nodig als je de data in python ook hebt verwerkt/gemaakt, 
# maar ik ga ervan uit dat je data met igor verwerkt en dan hier inleest)
datasetnrs = [1,2,3,4,5, 10] #dit bepaalt welke datasets worden ingelezen
datadictionary = {} #dit is een dictionary, die zijn erg handig in python https://www.codevscolor.com/python-dictionary-tutorial
for nr in datasetnrs:
    x,y=np.loadtxt('data'+str(nr)+'.txt',skiprows=1,unpack=True)
    xdata = x
    datadictionary[nr]=y

print(datadictionary[3])

nrlijst = np.array(list(datadictionary.keys()))
print('nrlijst is ',nrlijst)

print('simpele plot')
fig, ax = plt.subplots()
for nr in nrlijst:
    ax.plot(xdata, datadictionary[nr])
plt.show()
plt.close()

#kies kleuren op basis van het nummer van de dataset
min_nr = np.min(nrlijst)
max_nr = np.max(nrlijst)
#kies een colormap. De perceptually uniform zijn de beste omdat ze ook werken bij printen in zwart-wit en kleurenblindheid (zie hier korte uitleg https://jakevdp.github.io/PythonDataScienceHandbook/04.07-customizing-colorbars.html)
#hier kun je ze vinden (je hoeft de tekst niet te lezen, het gaat om de plaatjes) https://matplotlib.org/stable/tutorials/colors/colormaps.html
gekozen_colormap = cm.viridis_r 
def to_color(nr): #functie van gemaakt, omdat je anders de regel code hieronder steeds overal moet gebruiken
    return gekozen_colormap((nr-min_nr)/max_nr) #aan colormaps geef je een getal tussen 0 en 1

print('plot met gekozen kleur')
fig, ax = plt.subplots()
for nr in nrlijst:
    ax.plot(xdata, datadictionary[nr], color=to_color(nr))
plt.show()
plt.close()

print('plot met legenda')
fig, ax = plt.subplots()
for nr in nrlijst:
    ax.plot(xdata, datadictionary[nr], color=to_color(nr), label='data'+str(nr))
plt.legend()
plt.show()
plt.close()

print('plot nu met de ticks ook goed, labels en de grenzen aangegeven')
fig, ax = plt.subplots()
for nr in nrlijst:
    ax.plot(xdata, datadictionary[nr], color=to_color(nr), label='data'+str(nr))
ax.tick_params(top=True, direction='in')  
ax.tick_params(right=True, direction='in')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Title')
ax.axis([-10,10,-50,50]) #definieer de grenzen van je plot: x_min, x_max, y_min, y_max
plt.legend()
plt.show()
plt.close()

print('plot met colorbar in plaats van legenda')
fig, ax = plt.subplots()
for nr in nrlijst:
    ax.plot(xdata, datadictionary[nr], color=to_color(nr), label='data'+str(nr))
ax.tick_params(top=True, direction='in')  
ax.tick_params(right=True, direction='in')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Title')
ax.axis([-10,10,-50,50]) #definieer de grenzen van je plot: x_min, x_max, y_min, y_max

#colorbar
sm = plt.cm.ScalarMappable(cmap=gekozen_colormap, norm=plt.Normalize(vmin=min_nr, vmax=max_nr))
cbar = fig.colorbar(sm, ax=ax)
cbar.ax.set_ylabel('Nummer')  

plt.show()
plt.close()

print('plot met uniform verdeelde kleurtjes. let op: je kan hier geen colorbar gebruiken') 
fig, ax = plt.subplots()
for nr,i in zip(nrlijst,range(len(nrlijst))): #ik gebruik hier de teller i, die loopt van 0 tot de lengte van nrlijst. 
    kleur = gekozen_colormap(i/(len(nrlijst)-1)) #ik deel i hier door de maximale waarde van i, en die is len(nrlijst)-1. Zo loopt het getal dat in de colormap gaat van 0 tot 1
    ax.plot(xdata, datadictionary[nr], color=kleur, label='data'+str(nr))
ax.tick_params(top=True, direction='in')  
ax.tick_params(right=True, direction='in')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend()
ax.set_title('Title')
ax.axis([-10,10,-50,50]) #definieer de grenzen van je plot: x_min, x_max, y_min, y_max
plt.show()
plt.close()

print('plot met verticaal verschoven lijnen')
fig, ax = plt.subplots()
oplopendenrlijst = np.sort(nrlijst) #door de lijst te laten oplopen, heb je niet een willekeurige volgorde van traces
spacing = 15 #probeer hier wat getallen uit om te kijken wat mooi is
for nr,i in zip(oplopendenrlijst,range(len(oplopendenrlijst))): #ik gebruik hier de teller i, die loopt van 0 tot de lengte van nrlijst. 
    offset = i * spacing #de absolute offset wordt dus steeds groter voor opeenvolgende datasets
    ax.plot(xdata, datadictionary[nr]-np.min(datadictionary[nr])+offset, color=to_color(nr), label='data'+str(nr)) 
ax.tick_params(top=True, direction='in')  
ax.tick_params(right=True, direction='in')
ax.set_xlabel('x')
ax.set_ylabel('y (arb. units)')
ax.tick_params(labelleft=False) #je kunt de getallen op de as ook verbergen als je wilt
ax.set_title('Title')
ax.axis([None,None,None,None]) #definieer de grenzen van je plot: x_min, x_max, y_min, y_max

#colorbar
sm = plt.cm.ScalarMappable(cmap=gekozen_colormap, norm=plt.Normalize(vmin=min_nr, vmax=max_nr))
cbar = fig.colorbar(sm, ax=ax)
cbar.ax.set_ylabel('Nummer')  

plt.show()
plt.close()


