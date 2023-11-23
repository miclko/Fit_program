import numpy as np
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





def cut_peaks_and_shrink_data(data:dict,timeKey,usedKeys, nrWashings:int,endRunTime = np.inf,cutStartOffset=[10,10], cutEndOffset = 50, nrOfDataPoints = 300,peakrecognition = 0.7, peakSize = 100):
    '''Here trying to cut out the peaks in the data
    data: dictionary of time and concentration data
    timeKey: designate which one of the data elements is the time
    usedKeys: designate all the other used keys for the data (maybe not all elements in data should be used)
    nrWashings: the number of washings used in the data, for peak finding
    endRUnTime: optional, give the time in s until the data should be used(cuts out all later times from data)
    cutStartOffset: optional, defines how many data points are removed before the program detected peaks (to remove washing related irrelevant behavior)
    cutEndOffset: optional, defines how many points are removed after the registered peak positions
    nrOfDataPoints: optional, data is cleaned at the end and reduced to this number of points, such that fitting is easier (don't want to fit 10k+ points)
    peakrecognition: optional, value that decides when a peak is registered. If concurrent data points change by that much, a peak is registered (the last data in data is used, further assumed to be the same for all sets in data)
    peakSize: optional, designate the width of a peak, after finding the first peak, this number of data points is skipped and the next peak is search for data points ~last_peak_posi + peakSize
    returns datavalsdict, timevals, t0_estimates
    '''
    peakstartPosis = Find_washings_via_peaks(data[list(data.keys())[-1]],nrWashings,peakrecognition,peakSize)

    startcuts =  [x + y for (x,y) in zip(peakstartPosis,cutStartOffset)] 
    endcuts = [x - cutEndOffset for x in peakstartPosis[1:]]

    if endRunTime  == np.inf:
        endcuts.append(-1)
    else:
        for i in range(len(data[timeKey])):
            if data[timeKey][i] > endRunTime:  
                endcuts.append(i)
                break

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