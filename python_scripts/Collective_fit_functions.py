import numpy as np






#############################################
#
#Data pre processing 
#
#############################################

def Single_list_cut(listToCut, cutsStart,cutsEnd):
    '''Function that cuts a list according to the given specifications
    listTocut: list that is used for cutting
    cutsStart: List that has a number of start values for the individual cuts
    cutsEnd: List that has the same number as cutsStart for the individual ends
    Returns a list where all cuts are placed in one list
    '''
    if len(cutsStart) != len(cutsEnd):
        raise ValueError('inconsistent cuts chosen')
    listTmp = []
    for (x,y) in zip(cutsStart, cutsEnd):
        listTmp = listTmp + listToCut[x:y]

    return listTmp    

def list_shrinker(toshrinkList:list, sizeLimit:int):
    '''until the list is smaller than the given sizeLimit, delete every second element of a list in a loop, for improving the speed of e.g. a fit of the data'''
    if not isinstance(toshrinkList,list):
        return
    while len(toshrinkList) > sizeLimit:
        del toshrinkList[::2]




def Find_washings_via_peaks(data:list, nrWashings:int, peakRecognition:float, peakSize) ->list:       
    '''Trying to determine the peaks by their very fast speed
    data: to be searched for peaks
    nrWashings: nr of expected peaks
    peakRecognition: difference between adjacent sites that one can recognize as a peak
    peakSize: give approximate size of a peak, will skip this for finding further peak, after the first one
    '''
    def individual_peak_finder(data, startposition, peakRecognition:float):
        '''Using data[startPosition:], here it attempts to calculate the first peak position, when it starts
        data: to be used
        startposition: startposition when it should try to find the peak in data
        peakRecognition: change needed for recogniting a peak, if the difference between two consecutive values is larger than this, found the peak
        '''
        tmpdata = data[startposition:]
        for i in range(len(tmpdata) - 1):
            if abs(tmpdata[i + 1] - tmpdata[i]) > peakRecognition:
                return i + startposition
        raise IndexError('could not find a peak')
    
    peakPositions = []
    nextStart = 0
    for i in range(nrWashings):

        peakPositions.append(individual_peak_finder(data,nextStart, peakRecognition))
        nextStart = peakPositions[-1] + peakSize

    return peakPositions





def cut_peaks_and_shrink_data(data:dict,timeKey,usedKeys, nrWashings:int,endRunTime = np.inf,offsetAfterPeak=[10,10], offsetBeforePeak = [50], nrOfDataPoints = 300,peakrecognition = 0.7, peakSize = 100):
    '''Here trying to cut out the peaks in the data
    data: dictionary of time and concentration data
    timeKey: designate which one of the data elements is the time
    usedKeys: designate all the other used keys for the data (maybe not all elements in data should be used)
    nrWashings: the number of washings used in the data, for peak finding
    endRUnTime: optional, give the time in s until the data should be used(cuts out all later times from data)
    offsetAfterPeak: optional, defines how many points are removed after the registered peak positions 
    offsetBeforePeak: optional, defines how many data points are removed before the program detected peaks (to remove washing related irrelevant behavior)
    nrOfDataPoints: optional, data is cleaned at the end and reduced to this number of points, such that fitting is easier (don't want to fit 10k+ points)
    peakrecognition: optional, value that decides when a peak is registered. If concurrent data points change by that much, a peak is registered (the last data in data is used, further assumed to be the same for all sets in data)
    peakSize: optional, designate the width of a peak, after finding the first peak, this number of data points is skipped and the next peak is search for data points ~last_peak_posi + peakSize
   
     returns datavalsdict, timevals, t0_estimates
    '''
    peakstartPosis = Find_washings_via_peaks(data[list(data.keys())[-1]],nrWashings,peakrecognition,peakSize)

    #will break with default values if different than 2 peaks is considered, as right now that is kind of hardcoded here...
    startcuts =  [x + y for (x,y) in zip(peakstartPosis,offsetAfterPeak)] 
    endcuts = [x - y for (x,y) in zip(peakstartPosis[1:],offsetBeforePeak)]

    if endRunTime  == np.inf:
        endcuts.append(-1)
    else:
        for i in range(len(data[timeKey])):
            if data[timeKey][i] > endRunTime:  
                endcuts.append(i)
                break
        #add final time if larger time chosen    
        if len(endcuts) != len(startcuts):
            endcuts.append(-1)     


    print('using the start and endcuts:')
    print(startcuts)
    print(endcuts)
    print('for total length' + str(len(data[timeKey])))

    datavalsdict = dict() #contains all cutted data sets
    t0_estimates = [data[timeKey][x] for x in peakstartPosis] #get the times of the cuts as t0 estimates
    timevals = Single_list_cut(data[timeKey], startcuts,endcuts)
    list_shrinker(timevals,nrOfDataPoints)
    dataKeys = []
    for x in usedKeys:
        if x !=timeKey:
            dataKeys.append(x)
    for x in dataKeys:
        yvals = Single_list_cut(data[x], startcuts,endcuts)
        list_shrinker(yvals,nrOfDataPoints)
        datavalsdict[x] = yvals

    return datavalsdict, timevals, t0_estimates


def Cut_timeslot_from_dataDict(dataDict:dict, timeInterval:list, timeKey:str)->dict:
    '''returns a new dict that consiers only a certain time interval'''
    if len(timeInterval) != 2:
        raise ValueError('interval to be cut needs exactly 2 values')
    if not timeKey in dataDict.keys():
        raise ValueError('wrong timeKey given')
    
    tvals = dataDict[timeKey]
    cutPosis = []
    for y in timeInterval:
        for i in range(len(tvals)):
            if tvals[i] > y:
                cutPosis.append(i)
                break
            if i + 1 == len(tvals):
                cutPosis.append(i)

    returndict = dict()
    for x in dataDict.keys():
        returndict[x] = dataDict[x][cutPosis[0]:cutPosis[1]]
    return returndict

def lin_fit(vals, a,b):
    return a*np.array(vals) + b




#############################################
#############################################
#
#
# The fitting process
#
#
#############################################
#############################################

import matplotlib.pyplot as plt
import python_scripts.Rate_equation_solution as rateEqs
import lmfit


def collective_spr_fit_linear_k23(t,C1,C2,C3,a_k12_1,a_k12_2,k21,a_k23_1,a_k23_2,k32,t0_1,t0_2,motor_concentrations,weight_C3=2):
    '''linear_k23_fit'''
    listAllData = []
    for x in motor_concentrations:
        listAllData.append(rateEqs.Fit_3_states_discrete_rate_changes(t,C1_start=C1,C2_start=C2,C3_start=C3,k12=[a_k12_1*x,a_k12_2*x],k21=[k21,k21],k13=[0,0],k23=[a_k23_1*x,a_k23_2*x],k32=[k32,k32],k31=[0,0],t_0=[t0_1,t0_2]))
   
    result=np.array([])
    for x in listAllData:
        result=np.hstack([result, np.array(x[1]) + weight_C3*np.array(x[2])])

    if np.NaN in np.log(result):
        print('found nan in collective somehow')    
    return np.log(result)



def collective_spr_fit_constant_k23(t,C1,C2,C3,a_k12_1,a_k12_2,k21:float,a_k23_1,a_k23_2,k32:float,t0_1,t0_2,motor_concentrations,weight_C3=1):
    '''constant_k23_fit'''
    listAllData = []
    for x in motor_concentrations:
        listAllData.append(rateEqs.Fit_3_states_discrete_rate_changes(t,C1_start=C1,C2_start=C2,C3_start=C3,k12=[a_k12_1*x,a_k12_2*x],k21=[k21,k21],k13=[0,0],k23=[a_k23_1,a_k23_2],k32=[k32,k32],k31=[0,0],t_0=[t0_1,t0_2]))
   
    result=np.array([])
    for x in listAllData:
        result=np.hstack([result, np.array(x[1]) + weight_C3*np.array(x[2])])
        # result=np.hstack([result, np.array(x[1]) ])

    if np.NaN in np.log(result):
        print('found nan in collective somehow')    
    return np.log(result)



def fit_spr_data_to_model(fitparams,times,fitdata,fitfunction,concentrations,only_show_start_fit = False):
    '''Fit function that fits to the spr data using a supplied function'''
    #bring parameters into parameter form for lmfit
    parametersdict = dict()
    for x in fitparams:      
            parametersdict[x] = lmfit.Parameter(x,**fitparams[x])

    #show for referene the plot using the start parameter, then leave function
    if only_show_start_fit == True:
        datavals = np.reshape(fitdata,(len(concentrations),len(times)))
        for y in datavals:
            plt.plot(times,y, 'g')
        for (y,z) in zip(np.reshape(fitfunction(times,**parametersdict),(len(concentrations),len(times))),concentrations):
            plt.plot(times,y, '--', label=z)
        plt.title('start values')
        plt.legend()
        plt.show()
        return

    #define the model and fit to the data
    model = lmfit.Model(fitfunction, independent_vars=['t'])
    result = model.fit(fitdata, t=times, **parametersdict)

    return (result,fitfunction)


def plot_and_save_result(result,data,times,title,concentrations):
    '''Using the result of fit_spr_data_model, it is plotted and the result saved'''
    fitvals = result[1](times,**result[0].params)
    datavals = np.reshape(np.log(data),(len(concentrations),len(times)))



    for y in datavals:
        plt.plot(times,y, 'g')
    for (y,z) in zip(np.reshape(fitvals,(len(concentrations),len(times))),concentrations):
        plt.plot(times,y, '--', label=z)
    plt.legend(title='C motor in nM')
    plt.xlabel('t in s')
    plt.ylabel('response')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(title + '.pdf')
