import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

folderstart = 'P:/Surfdrive/DATA/'
folderstart = 'D:/Surfdrive/DATA/'
year = '2020/'
file_kw = '/KW/Images/stickingprob_vs_position.txt'
file_TOF = '/TOF/Images/4.0/energy.txt'
savefolder = folderstart+'D2 sticking probability as a function of step density and kinetic energy/TOFparams/'
save_figs = True

datadic = {0.7:{'kw_day':'03 Mar/200310', 'TOF_day':'03 Mar/200310','c':'red'},
           0.8:{'kw_day':'05 May/200512', 'TOF_day':'01 Jan/200128','c':'orange'},
           1.2:{'kw_day':'03 Mar/200305', 'TOF_day':'02 Feb/200220','c':'gold'},
           2.0:{'kw_day':'06 Jun/200602', 'TOF_day':'05 May/200511','c':'yellow'},
           2.9:{'kw_day':'03 Mar/200303', 'TOF_day':'02 Feb/200218','c':'greenyellow'},
           4.0:{'kw_day':'03 Mar/200311', 'TOF_day':'03 Mar/200306','c':'green'},
           4.77:{'kw_day':'05 May/200519', 'TOF_day':'05 May/200519','c':'mediumseagreen'},
           6.1:{'kw_day':'05 May/200520', 'TOF_day':'02 Feb/200204','c':'teal'},
           7.8:{'kw_day':'08 Aug/200817', 'TOF_day':'08 Aug/200813_2','c':'blue'},
           9.4:{'kw_day':'08 Aug/200814', 'TOF_day':'08 Aug/200813','c':'darkblue'},
           10.7:{'kw_day':'07 Jul/200727', 'TOF_day':'07 Jul/200727','c':'indigo'},
           12.9:{'kw_day':'09 Sep/200904', 'TOF_day':'09 Sep/200904','c':'purple'},
           13.1:{'kw_day':'09 Sep/200910', 'TOF_day':'09 Sep/200910','c':'black'}}


for key in list(datadic.keys()):
    fitparams = np.loadtxt(savefolder+'fitparams_'+str(key)+'kjm.txt',usecols=2,unpack=True,skiprows=1)
    v0 = fitparams[1]
    alpha = fitparams[3]