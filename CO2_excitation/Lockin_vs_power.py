import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm 
import datetime as dt
from scipy.optimize import curve_fit

################ General Parameters #################

folderstart = "C:/Users/jansenc3/surfdrive/"
folderstart = "C:/Users/Werk/Surfdrive/"
folder = folderstart+"DATA/Power and wavelength data/"
timecmap = cm.viridis

################ Experiment Parameters #################

#5 mL/min measurement for signal vs power
offresonance = 2.15
min_power = 0.008
year = '2021'
month = '05'
day = '04'
timedic = {'firstpurge':'11:51:58',
           'reversepurge':'12:06:14',
           'start1':'14:07:34',
           'stop1':'16:35:25'}
refwl = 4252.718
wl_window = 0.01

#1 mL/min measurement for signal vs power
offresonance = 7.4
min_power = 0.008
month = '05'
day = '21'
timedic = {'start0':'11:18:49',
           'start1':'11:30:00',
           'stop1':'13:34:24'}
refwl = 4252.71
wl_window = 0.015

# #5 mL/min measurement for signal vs position (with slit) 
# month = '05'
# day = '27'
# timedic = {'start':'15:14:58',
#            'stop':'17:12:32'}

# positiondic = {42.57:'15:14:58',
#                42.0:'15:23:05',
#                45.0:'15:26:17',
#                44.5:'15:30:21',
#                44.0:'15:33:26',
#                43.5:'15:38:13',
#                43.2:'15:41:15',
#                43.1:'15:47:29',
#                42.9:'15:50:28',
#                42.7:'16:01:37',
#                42.6:'16:14:49',
#                42.5:'16:35:36',
#                42.4:'16:49:51',
#                42.3:'17:02:12',
#                'stop':'17:12:32'}
# refwl = 4252.71
# wl_window = 0.015

#5 mL/min measurement for signal vs power, with new attenuation system
offresonance = 2.05
min_power = 0.0076
year = '2021'
month = '06'
day = '01'
timedic = {'start1':'13:06:57',
           'stop1':'14:53:57'}
refwl = 4252.71
wl_window = 0.015
lens = False

# #5 mL/min measurement for signal vs power, with RAP, possibly misaligned
# offresonance = 2.05
# min_power = 0.0076
# year = '2021'
# month = '06'
# day = '04'
# timedic = {'start1':'13:22:46',
#            'stop1':'14:39:50'}
# refwl = 4252.708
# wl_window = 0.009
# lens = True

# #1 mL/min measurement for signal vs power, with RAP
# offresonance = 7.4
# min_power = 0.004
# month = '06'
# day = '11'
# timedic = {'start1':'12:22:26',
#            'stop1':'13:45:00'}
# refwl = 4252.71
# wl_window = 0.015
# lens = True

# #5 mL/min measurement for signal vs power, with RAP alignment the same as the measurement on 11/6/21
# offresonance = 2.05
# min_power = 0.0039
# year = '2021'
# month = '06'
# day = '17'
# timedic = {'start0':'15:10:55',
#            'start1':'15:23:03',
#            'stop1':'17:13:49'}
# refwl = 4252.715
# wl_window = 0.006
# lens = True

months = {'01':'Jan', '02':'Feb', '03':'Mar', '04':'Apr', '05':'May', '06':'Jun', '07':'Jul', '08':'Aug', '09':'Sep', '10':'Oct', '11':'Nov', '12':'Dec'}
savefolder = (folderstart+'DATA/'+str(year)+'/'+month+' '+months[month]+'/'+year[-2:]+month+day+'/Laser/Images/')
save = False

################ Main #################

def main():
    if save:
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)

    datetimearray, timestamparray, power, wl, lockin = read_file()
    
    main_lockin_vs_power(datetimearray, timestamparray, power, wl, lockin)
    
    # main_lockin_vs_position(datetimearray, timestamparray, power, wl, lockin)

    
################ Script functions #################

def main_lockin_vs_position(datetimearray, timestamparray, power, wl, lockin, save=False):
    """
    Script that plots the signal vs beam position
    """
    datamask = select_datawindow(timestamparray, starttime=convert_time_to_timestamp(timedic['start']), stoptime=convert_time_to_timestamp(timedic['stop']))
    timestamparray=timestamparray[datamask]
    power=power[datamask]
    wl=wl[datamask]
    lockin=lockin[datamask]

    dm = remove_offresonance_wl(timestamparray,wl,plot=True, refwl=refwl, wl_window=wl_window, invert=False)
    plot_raw_data(timestamparray[dm]-np.min(timestamparray[dm]), (power[dm],wl[dm],lockin[dm]))
    dm_invert = remove_offresonance_wl(timestamparray,wl,plot=True, refwl=refwl, wl_window=wl_window, invert=True)
    plot_raw_data(timestamparray[dm_invert]-np.min(timestamparray[dm_invert]), (power[dm_invert],wl[dm_invert],lockin[dm_invert]))
    
    #make an array with the slit positions
    positions = np.full(len(timestamparray),-1.0) #fill with -1 so I can filter later
    integrationtime = 30 #s
    avg_onresonance = {}
    avg_offresonance = {}
    for i in range(len(positiondic.keys())-1):
        pos = list(positiondic.keys())[i]
        end = list(positiondic.keys())[i+1]
        datamask = select_datawindow(timestamparray, starttime=convert_time_to_timestamp(positiondic[pos])+integrationtime, stoptime=convert_time_to_timestamp(positiondic[end]))
        positions[datamask] = pos
        avg_onresonance[pos] = np.average(lockin[datamask*dm])
        avg_offresonance[pos] = np.average(lockin[datamask*dm_invert])
    
    plot_raw_data(timestamparray-timestamparray[0], (power,wl,lockin,positions))
    
    #plot signal vs positions
    datamask = positions > 0 #to remove the -1 positions that I did not want to use
    plt.plot(positions[datamask*dm], lockin[datamask*dm], '.', label='On resonance')
    plt.plot(positions[datamask*dm_invert], lockin[datamask*dm_invert], '.', label='Off resonance')
    plt.plot(list(avg_onresonance.keys()),list(avg_onresonance.values()),'o', c='black', label='Average on resonance')
    plt.plot(list(avg_offresonance.keys()),list(avg_offresonance.values()),'o', c='r', label='Average off resonance')
    plt.xlabel('Slit position (mm)')
    plt.ylabel('lock-in signal')
    plt.legend()
    plt.show()
    plt.close()
    
    #plot ratio of on/off resonance vs position
    for pos in avg_offresonance.keys():
        plt.plot(pos,avg_onresonance[pos]/avg_offresonance[pos],'o',c='black')
    plt.axis([None,None,1,None])
    plt.xlabel('Slit position (mm)')
    plt.ylabel('Ratio lock-in signal on/off resonance')
    plt.show()
    plt.close()
 
def main_lockin_vs_power(datetimearray, timestamparray, power, wl, lockin):
    """
    Script that plots and fits the signal vs laser power
    """
    #plot all data
    plot_raw_data(timestamparray-timestamparray[0], (power,wl,lockin),save=save,savenote='all_raw')
    
    #Select time window
    datamask = select_datawindow(timestamparray, starttime=convert_time_to_timestamp(timedic['start1']), stoptime=convert_time_to_timestamp(timedic['stop1']))
    timestamparray=timestamparray[datamask]
    power=power[datamask]
    wl=wl[datamask]
    lockin=lockin[datamask]
    #Plot raw data in window
    plot_raw_data(timestamparray-timestamparray[0], (power,wl,lockin),save=save,savenote='window_raw')
    
    #Remove data where laser was off resonance
    datamask = remove_offresonance_wl(timestamparray,wl,plot=True, refwl=refwl, wl_window=wl_window)
    timestamparray=timestamparray[datamask]
    power=power[datamask]
    wl=wl[datamask]
    lockin=lockin[datamask]
    #Plot raw data with off resonance removed
    plot_raw_data(timestamparray-timestamparray[0], (power,wl,lockin), save=save, savenote='onresonance')
    
    #Plot signal vs laser power
    plot_lockin_vs_power(power,lockin,timestamparray,offresonancevalue=offresonance,save=save)
    for function in ['exp', 'exp_pol']:
        max_sig = fit_lockin_vs_power(power,lockin,timestamparray,function=function,save=save)


    #Plot population vs actual laser power
    correctedpower = absolute_laser_power(power)
    pop = population(lockin, max_sig, offresonance, lens=lens)
    plot_lockin_vs_power(correctedpower, pop, timestamparray, save=False)
    


################ Other functions #################
    
def convert_time_to_timestamp(time):
    if time==None:
        return time
    datetime=year+'-'+month+'-'+day+' '+time
    timestamp = dt.datetime.timestamp(dt.datetime.strptime(datetime,'%Y-%m-%d %H:%M:%S'))
    return timestamp

def select_datawindow(timestamparray,starttime=None,stoptime=None):
    """
    Returns a datamask for the selected time window
    """
    if starttime:
        low = timestamparray > starttime
    else:
        low = True 
    
    if stoptime: 
        high = timestamparray < stoptime
    else:
        high = True
        
    return low*high  
    
def remove_offresonance_wl(timestamparray, wl, refwl=4252.71, wl_window=0.01, integrationtime=30, plot=False, invert=False, save=False): #data points up to integrationtime seconds after relocking are also removed
    """
    Removes all data where the laser was off resonance, and also up to integrationtime
    later to account for the delay from the lock-in amplifier
    If invert: returns datamask that selects all data where the laser is off resonance
    """
    #create data mask without integration time
    low = wl > refwl - wl_window
    high = wl < refwl + wl_window
    onresonance = low*high
    #if invert, it returns the datamask for off resonance instead of on resonance
    if invert: 
        onresonance = np.invert(onresonance)
    
    if plot:
        plt.plot(timestamparray,wl,'.',label='all data')
        plt.plot(timestamparray[onresonance],wl[onresonance],'.',label='off resonance removed')
    
    #include integration time
    times_to_remove = np.argwhere(np.invert(onresonance)) #start with only the times off resonance
    for timestamp in timestamparray[np.invert(onresonance)]: #for each time off resonance, also include every time up to 30 seconds afterwards
        lo = timestamparray > timestamp
        hi = timestamparray < timestamp + integrationtime
        mask = lo*hi #these are all times between the time off resonance and integrationtime afterwards
        times_to_remove = np.concatenate((times_to_remove, np.argwhere(mask)))
    #convert list of indices to data mask
    datamask = np.ones(len(timestamparray))
    datamask[times_to_remove] = 0
    datamask = datamask.astype(bool)
    if plot:
        plt.plot(timestamparray[datamask],wl[datamask],'.',label='integration time included')
        plt.legend()
        if save:
            plt.savefig(savefolder+'offresonance_removed_invert='+str(invert)+'.png',dpi=500, bbox_inches='tight')
        plt.show()
        plt.close()
        
    return datamask

def smooth(array, rollingaverage=6):
    """
    Can only smooth the "center" of the array between rollingaverage/2 and -rollingaverage/2. 
    The edges are filled with the first and last averaged value.
    """
    if rollingaverage%2:
        print('please input an even number for rollingaverage')
        return
    sm = np.cumsum(array)
    avg = (sm[rollingaverage:]-sm[:-rollingaverage])/rollingaverage
    array2=np.copy(array)
    array2[int(rollingaverage/2):-int(rollingaverage/2)]=avg
    array2[:int(rollingaverage/2)]=avg[0]
    array2[-int(rollingaverage/2):]=avg[-1]
    return array2
  
def plot_lockin_vs_power(power,lockin,timestamparray,offresonancevalue=0,save=False):
    colorvalues = (timestamparray-timestamparray[0])/np.max(timestamparray-timestamparray[0])
    # print(colormap)
    plt.scatter(power,lockin,2,marker='o',c=colorvalues,cmap=timecmap)
    plt.plot([0,np.max(power)*1.2],[offresonancevalue,offresonancevalue],'--',c='black')
    # plt.axis([, None,0,None])
    plt.xlabel('Laser power (W)')
    plt.ylabel('Lockin signal')
    if save:
        plt.savefig(savefolder+'lockin_vs_power.png',dpi=500, bbox_inches='tight')
    plt.show()
    plt.close()

def fit_lockin_vs_power(power,lockin,timestamparray, function='exp',fixzero=False,save=False):
    """
    Fits an exponential: A*(1-exp(b(x-c)))+d, function='exp'
    or multiplies exponential with polynomial: A*(1-e*(power-f)^n*exp(b(power-c)))+d, function='exp_pol'
    """     
    def exp(power, A,b,c,d):
        return A*(1-np.exp(b*(power-c)))+d
    
    def exp_pol(power, A,b,c,d,e,n):
        return A*(1-e*(power-c)**n*np.exp(b*(power-c)))+d

   
    A_guess = np.max(lockin)-np.min(lockin)
    b_guess = -100
    c_guess = np.min(power)*0.5 #to avoid issues with dividing by 0
    d_guess = np.min(lockin)
    
    if function=='exp':
        func = exp
        guess = [A_guess,b_guess,c_guess,d_guess]
    elif function=='exp_pol':
        func = exp_pol
        e_guess = 0.5
        n_guess = -0.1
        guess = [A_guess,b_guess,c_guess,d_guess, e_guess, n_guess]

    bounds_min = np.full(len(guess),-np.inf)
    bounds_max = np.full(len(guess),np.inf)

    if fixzero:
        bounds_min[2]=0.999*c_guess
        bounds_max[2]=1.001*c_guess
        bounds_min[3]=0.999*d_guess
        bounds_max[3]=1.001*d_guess
    

    popt = curve_fit(func, power,lockin,p0=guess,bounds=(bounds_min,bounds_max))
    max_signal = func(100000, *popt[0])
    print(popt[0])

    plt.plot(power,lockin,'o',markersize=2,label='Data',c='gray')
    plt.plot(power,func(power,*guess),'.',markersize=0.5,label='guess',c='b')
    plt.plot(power,func(power,*popt[0]),'o',markersize=2,label='Fit ('+function+')',c='r')
    plt.plot(power,np.full(len(power),max_signal),'--',c='black',label='max excitation')
    plt.plot(power,np.full(len(power),offresonance),'--',c='gray',label='no excitation')
    plt.xlabel('Power')
    plt.ylabel('Lock-in signal')
    plt.legend()
    plt.title('Max excitation: '+str(np.round(max_signal,2))+'V = '+str(np.round(max_signal/offresonance*100-100,2))+'% increase')
    if save:
        plt.savefig(savefolder+'lockin_vs_power_'+function+'_fit.png',dpi=500, bbox_inches='tight')
    plt.show()
    plt.close()

    return max_signal
    
def plot_raw_data(timestamparray, data, starttime=None, endtime=None,save=False,savenote=''):
    colorvalues = (timestamparray-timestamparray[0])/np.max(timestamparray-timestamparray[0])
    fig, axes = plt.subplots(np.shape(data)[0],sharex='col')
    fig.subplots_adjust(hspace=0)
 
    def find_label(array):
        label = 'Unknown'
        if np.max(array)<1 and np.min(array)>0:
            label='Power'
        elif np.max(array)<20:
            label='Lock-in'
        elif np.max(array)>4000:
            label='Wl'
        return label
    
    for i in range(len(axes)):
        # axes[i].plot(timestamparray, data[i],label=find_label(data[i]))
        axes[i].scatter(timestamparray,data[i],2,c=colorvalues,cmap=timecmap,label=find_label(data[i]))
        axes[i].set_ylabel(find_label(data[i]))

    axes[i].set_xlabel('Time (s)')
    axes[i].axis([starttime,endtime,None,None])
    if save:
        plt.savefig(savefolder+'data_'+savenote+'.png',dpi=500, bbox_inches='tight')
    plt.show()    
    plt.close()

def read_file():
    power, wl, lockin = np.loadtxt(folder+str(year)+'-'+str(month)+'-'+str(day)+".txt", usecols=(2,3,6),unpack=True)

    nans = ~np.isnan(power)
    power = power[nans]
    wl = wl[nans]
    lockin = lockin[nans]

    datearray = np.genfromtxt(folder+str(year)+'-'+str(month)+'-'+str(day)+'.txt',dtype='str')[:,0]
    datearray = datearray[nans]
    timearray = np.genfromtxt(folder+str(year)+'-'+str(month)+'-'+str(day)+'.txt',dtype='str')[:,1]
    timearray = timearray[nans]
    datetimearray = np.core.defchararray.add(datearray,' ')
    datetimearray = np.core.defchararray.add(datetimearray,timearray)
    datetimearray = [dt.datetime.strptime(x,'%m-%d-%Y %H:%M:%S') for x in datetimearray]
    timestamparray = np.array([dt.datetime.timestamp(x) for x in datetimearray])

    return datetimearray, timestamparray, power, wl, lockin

    # timestamparray = np.zeros(len(timearray))

    # for i in range(len(timearray)):
        # timestamparray[i] = 3600*int(timearray[i][:2]) + 60*int(timearray[i][3:5]) + int(timearray[i][6:8])
    
def absolute_laser_power(power, T_air=0.90, T_window=0.85): #power_before in mW, power_loss_air in %
    """
    T_air: transmission of the air pocket between the window and power meter \n
    T_window: transmission of the IR window
    """
    return (power-np.min(power))/T_air/T_window

def population(lockin, max_signal, offres, lens=False):
    if lens:
        return (lockin-offres)/(max_signal-offres)
    else:
        return (lockin-offres)/(max_signal-offres)/2



################ End functions #################

main()