import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm 
import datetime as dt
from scipy.optimize import curve_fit

################ General Parameters #################

folderstart = "C:/Users/jansenc3/surfdrive/"
folderstart = "C:/Users/Werk/Surfdrive/"
folder = folderstart+"DATA/Power and wavelength data/"
subfolder = ""
timecmap = cm.viridis

################ Experiment Parameters #################

#20 mL/min CO2, R(2)
min_power = 0.002
year = '2022'
month = '01'
day = '24'
timedic = {'start1':'12:43:50',
           'stop1':'13:37:23'}
integration_time = 10 #s
lens = True

# #20 mL/min CO2, R(2)
# min_power = 0.0004
# year = '2022'
# month = '01'
# day = '21'
# timedic = {'start1':'14:45:57',
#            'stop1':'15:21:55'}
# integration_time = 30 #s
# lens = True

# #20 mL/min CO2, R(2) with slit
# min_power = 0.0036
# year = '2022'
# month = '02'
# day = '15'
# timedic = {'start1':'11:24:30',
#            'stop1':'14:32:29'}
# integration_time = 30 #s
# subfolder = 'R2/'
# lens = True

# #20 mL/min CO2, R(0) with slit
# min_power = 0.0036
# year = '2022'
# month = '02'
# day = '15'
# timedic = {'start1':'14:46:04',
#            'stop1':'16:37:23'}
# integration_time = 30 #s
# subfolder = 'R0/'
# lens = True

#20 mL/min CO2, R(0)
min_power = 0.0036
year = '2022'
month = '02'
day = '24'
timedic = {'start1':'12:10:58',
           'stop1':'14:27:58'}
integration_time = 10 #s
subfolder = 'R0/'
lens = True

#20 mL/min CO2, R(6)
min_power = 0.005
year = '2022'
month = '02'
day = '25'
timedic = {'start1':'11:41:15',
           'stop1':'13:18:04'}
integration_time = 10 #s
subfolder = 'R6/'
lens = True

#20 mL/min CO2, R(8)
min_power = 0.005
year = '2022'
month = '02'
day = '25'
timedic = {'start1':'13:54:23',
           'stop1':'15:57:45'}
integration_time = 10 #s
subfolder = 'R8/'
lens = True

#20 mL/min CO2, R(10)
min_power = 0.005
year = '2022'
month = '03'
day = '03'
timedic = {'start1':'14:21:47',
           'stop1':'16:12:23'}
integration_time = 10 #s
subfolder = 'R10/'
lens = True

#20 mL/min CO2, R(2)
min_power = 0.0046
year = '2022'
month = '03'
day = '04'
timedic = {'start1':'11:42:27',
           'stop1':'14:09:20'}
integration_time = 10 #s
subfolder = 'R2/'
lens = True

#20 mL/min CO2, R(4)
min_power = 0.0041
year = '2022'
month = '03'
day = '04'
timedic = {'start1':'15:43:05',
           'stop1':'17:15:00'}
integration_time = 10 #s
subfolder = 'R4/'
lens = True



############### Stuff #################

months = {'01':'Jan', '02':'Feb', '03':'Mar', '04':'Apr', '05':'May', '06':'Jun', '07':'Jul', '08':'Aug', '09':'Sep', '10':'Oct', '11':'Nov', '12':'Dec'}
savefolder = (folderstart+'DATA/'+str(year)+'/'+month+' '+months[month]+'/'+year[-2:]+month+day+'/Laser/Images/'+subfolder)
save = True

################ Main #################

def main():
    if save:
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)

    datetimearray, timestamparray, power, wl, lockin, mct = read_file()
    
    main_lockin_vs_power(datetimearray, timestamparray, power, wl, lockin, mct)


    
################ Script functions #################
 
def main_lockin_vs_power(datetimearray, timestamparray, power, wl, lockin, mct):
    """
    Script that plots and fits the signal vs laser power
    """
    #plot all data
    plot_raw_data(timestamparray-timestamparray[0], (power,wl,lockin, mct),save=save,savenote='all_raw')
    
    #Select time window
    datamask = select_datawindow(timestamparray, starttime=convert_time_to_timestamp(timedic['start1']), stoptime=convert_time_to_timestamp(timedic['stop1']))
    timestamparray=timestamparray[datamask]
    power=power[datamask]
    wl=wl[datamask]
    lockin=lockin[datamask]
    mct=mct[datamask]
    #Plot raw data in window
    plot_raw_data(timestamparray-timestamparray[0], (power,wl,lockin, mct),save=save,savenote='window_raw')
    
    #Find data where laser was off resonance
    datamask = remove_offresonance_mct(timestamparray,wl,mct, plot=True, integrationtime=integration_time, remove_onresonance=True)
    lockin_offresonance = lockin[datamask]
    if save:
        np.savetxt(savefolder+'lockin_offresonance.txt', [np.average(lockin_offresonance)])
    plt.plot(power[datamask], lockin_offresonance, 'o')
    plt.show()
    plt.close()


    #Remove data where laser was off resonance
    datamask = remove_offresonance_mct(timestamparray,wl,mct, plot=True, integrationtime=integration_time)
    timestamparray=timestamparray[datamask]
    power=power[datamask]
    wl=wl[datamask]
    lockin=lockin[datamask]
    mct=mct[datamask]
    if save:
        np.savetxt(savefolder+'lockin_vs_power.txt',np.column_stack((power,lockin)))
    #Plot raw data with off resonance removed
    plot_raw_data(timestamparray-timestamparray[0], (power,wl,lockin, mct), save=save, savenote='onresonance')


    
    #Plot signal vs laser power
    plot_lockin_vs_power(power,lockin,timestamparray,offresonancevalue=np.average(lockin_offresonance),save=save)
    #Fit signal vs laser power
    for function in ['exp', 'exp_pol', 'dbl_exp']:
        max_sig = fit_lockin_vs_power(power,lockin,timestamparray,function=function,save=save, offresonance=np.average(lockin_offresonance))


    #Plot population vs actual laser power
    correctedpower = absolute_laser_power(power)
    pop = population(lockin, max_sig, np.average(lockin_offresonance), lens=lens)
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
    

def remove_offresonance_mct(timestamparray, wl, mct, integrationtime=10, plot=False, remove_onresonance=False, save=False, threshold=None):
    if not threshold:
        threshold = (2*np.max(mct) + np.min(mct))/3

    if remove_onresonance:
        removes = mct < threshold
    else:
        removes = mct > threshold

    remove_mask = np.argwhere(removes)
    for timestamp in timestamparray[removes]: #for each time off resonance, also include every time up to 30 seconds afterwards
        lo = timestamparray > timestamp
        hi = timestamparray < timestamp + integrationtime*3
        mask = lo*hi #these are all times between the time off resonance and integrationtime*3 afterwards
        remove_mask = np.concatenate((remove_mask, np.argwhere(mask)))
    #convert list of indices to data mask
    datamask = np.ones(len(timestamparray))
    datamask[remove_mask] = 0
    datamask = datamask.astype(bool)

    if plot:
        plt.plot(timestamparray,wl,'.',label='all data')
        plt.plot(timestamparray[np.invert(removes)],wl[np.invert(removes)],'.',label='off resonance removed')
        plt.plot(timestamparray[datamask],wl[datamask],'.',label='integration time included')
        plt.legend()
        if save:
            plt.savefig(savefolder+'offresonance_removed_invert='+str(invert)+'.png',dpi=500, bbox_inches='tight')
        plt.show()
        plt.close()
        

    return datamask
    
    

    #include integrationtime
        

def remove_offresonance_wl(timestamparray, wl, refwl=4252.71, wl_window=0.01, integrationtime=10, plot=False, invert=False, save=False): #data points up to integrationtime seconds after relocking are also removed
    """
    Removes all data where the laser was off resonance, and also up to integrationtime*3
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
        hi = timestamparray < timestamp + integrationtime*3
        mask = lo*hi #these are all times between the time off resonance and integrationtime*3 afterwards
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

def fit_lockin_vs_power(power,lockin,timestamparray, function='exp',fixzero=False,save=False, offresonance=0):
    """
    Fits an exponential: A*(1-exp(b(x-c)))+d, function='exp'
    or multiplies exponential with polynomial: A*(1-e*(power-f)^n*exp(b(power-c)))+d, function='exp_pol'
    or double exponential: A1*(1-exp(b1(x-c1)))+A2*(1-exp(b2(x-c2)))+d, function='dbl_exp'
    """     
    def exp(power, A,b,c,d):
        return A*(1-np.exp(b*(power-c)))+d
    
    def exp_pol(power, A,b,c,d,e,n):
        return A*(1-e*(power-c)**n*np.exp(b*(power-c)))+d

    def dbl_exp(power, A1, A2, b1, b2, c1, c2, d):
        return A1*(1-np.exp(b1*(power-c1)))+A2*(1-np.exp(b2*(power-c2)))+d
   
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
        n_guess = -0.08 #-0.1 for older data
        guess = [A_guess,b_guess,c_guess,d_guess, e_guess, n_guess]
    elif function=='dbl_exp':
        func = dbl_exp
        guess = [A_guess/2,A_guess/2,b_guess,b_guess,c_guess,c_guess,d_guess]

    bounds_min = np.full(len(guess),-np.inf)
    bounds_max = np.full(len(guess),np.inf)

    if fixzero:
        bounds_min[2]=0.999*c_guess
        bounds_max[2]=1.001*c_guess
        bounds_min[3]=0.999*d_guess
        bounds_max[3]=1.001*d_guess

    plt.plot(power,lockin,'o',markersize=2,label='Data',c='gray')
    plt.plot(power,func(power,*guess),'.',markersize=0.5,label='guess',c='b')

    popt = curve_fit(func, power,lockin,p0=guess,bounds=(bounds_min,bounds_max))
    max_signal = func(100000, *popt[0])
    print(popt[0])

    plt.plot(power,func(power,*popt[0]),'o',markersize=2,label='Fit ('+function+')',c='r')
    plt.plot(power,np.full(len(power),max_signal),'--',c='black',label='max excitation')
    plt.plot(power,np.full(len(power),offresonance),'--',c='gray',label='no excitation')
    plt.xlabel('Power')
    plt.ylabel('Lock-in signal')
    plt.legend()
    plt.title('Max excitation: '+str(np.round(max_signal,2))+'V = '+str(np.round(max_signal/offresonance*100-100,2))+'% increase')
    if save:
        plt.savefig(savefolder+'lockin_vs_power_'+function+'_fit.png',dpi=500, bbox_inches='tight')
        np.savetxt(savefolder+'fitparameters_'+function+'.txt',popt[0])
        
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
    power, wl, lockin, mct = np.loadtxt(folder+str(year)+'-'+str(month)+'-'+str(day)+".txt", usecols=(2,3,6,7),unpack=True)

    nans = ~np.isnan(power)
    power = power[nans]
    wl = wl[nans]
    lockin = lockin[nans]
    mct = mct[nans]

    datearray = np.genfromtxt(folder+str(year)+'-'+str(month)+'-'+str(day)+'.txt',dtype='str')[:,0]
    datearray = datearray[nans]
    timearray = np.genfromtxt(folder+str(year)+'-'+str(month)+'-'+str(day)+'.txt',dtype='str')[:,1]
    timearray = timearray[nans]
    datetimearray = np.core.defchararray.add(datearray,' ')
    datetimearray = np.core.defchararray.add(datetimearray,timearray)
    datetimearray = [dt.datetime.strptime(x,'%m-%d-%Y %H:%M:%S') for x in datetimearray]
    timestamparray = np.array([dt.datetime.timestamp(x) for x in datetimearray])

    return datetimearray, timestamparray, power, wl, lockin, mct

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