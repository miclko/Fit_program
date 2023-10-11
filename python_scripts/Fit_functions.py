from python_scripts.Rate_equation_solution import *
from scipy.optimize import curve_fit
import warnings

class Storage_of_fit_choice:
    '''This class generates a way to call the fit function with the chosen fit parameters
    for that, it tries to translate in the init phase the startfitdict with its true/false elements into a consistent function call that can be called by curve fit
    It is assumed that startfitvaldict has already gone through Order_fitvalDict() function and is well defined, all true elements correspond to the fit parameters
    The chronological elements in the kwargs dict correspond to the chronological fit values
    It also has the capability to get the resulting fitted curves and stores the results from the fit, so that one can get
    everything by using functions from this class.
    It also generates valid boundaries from the given ones'''
    # nrFitArgs = 0
    # argumentStorageDict = dict() #stores all elements in a dict, also has a tuple that says where to put the variable elements upon fitting
    # listOfFitPositions = [] #stores the path to each arg used for fitting

    def __init__(self,startfitvaldict:dict):
        self.argumentStorageDict = dict() #contains the kwargs for the function call
        self.listOfFitPositions = []   #contains the positions where the fit params are being placed, in the format [ ['kvalname', listPosition:int], ...]

        for x in startfitvaldict:
            if isinstance(startfitvaldict[x][0], list):
                self.argumentStorageDict[x] = startfitvaldict[x][0]
            else:
                self.argumentStorageDict[x] = [startfitvaldict[x][0]]

        for x in  self.argumentStorageDict:
            if not isinstance(startfitvaldict[x][1], list):
                if startfitvaldict[x][1] == True:
                    self.listOfFitPositions.append([x, 0]) #assumes that always lists are used for true/false statements, not iplemented the case of shortcut notation
            else:
                for i in range(len(startfitvaldict[x][1])):
                    if startfitvaldict[x][1][i] == True:
                        self.listOfFitPositions.append([x, i])
       
        self.nrFitArgs = len(self.listOfFitPositions)

        # print('nrFitargs: ' + str(self.nrFitArgs))  
        # print('listofposis: ' )               
        # print(self.listOfFitPositions)

    def __modify_parameter_dict(self, *args):
        '''uses the arg list to change all parameters for the fit function'''
        for (x,y) in zip(self.listOfFitPositions, args):
            self.argumentStorageDict[x[0]][x[1]] =  y

    def fit_positions_and_names(self):
        '''return the names of all variable param positions used in the calculation'''
        return self.listOfFitPositions

    def current_variable_params_as_list(self):
        #get the parameters currently stored in the object, used for variable fitting
        paramsList = []
        for x in self.listOfFitPositions:
            paramsList.append(self.argumentStorageDict[x[0]][x[1]])
        return  paramsList
    
    def generate_boundaries_for_fit(self, boundaryDict:dict):
        '''uses a given boundary dict to generate consistent boundaries that can be used for the fit function
        will use default arguments if no boundaries given ([0],[np.inf] typically) in total the boundary list is of the form [[lower bounds], [upper bounds]] with len(lower bounds)=nrFitArgs'''
        # boundValDict['C1_start'] = [5,10] #if no fit, then no keyword stored of it, see C2,C3
        # boundValDict['t_0'] = [[1.7,5],[1.8, 5.3]]
        listOfBounds = [[],[]]
        for x in self.listOfFitPositions:
            if x[0] in boundaryDict:
                if not isinstance(boundaryDict[x[0]][0],list):
                    listOfBounds[0].append(boundaryDict[x[0]][0]) #append lower bound value
                    listOfBounds[1].append(boundaryDict[x[0]][1]) #append upper bound value
                else:
                    listOfBounds[0].append(boundaryDict[x[0]][0][x[1]]) #append lower bound value
                    listOfBounds[1].append(boundaryDict[x[0]][1][x[1]]) #append upper bound value
            else:
                listOfBounds[0].append(0) #append default lower bound value
                listOfBounds[1].append(np.inf) #append default upper bound value
        return listOfBounds            
                
                


    def call_full_function(self, times, *args):
        if len(args) != self.nrFitArgs:
            raise ValueError('wrong number of arguments given for the expected fit function')
        self.__modify_parameter_dict(*args)
        return Fit_3_states_discrete_rate_changes(times, **self.argumentStorageDict) 

    def call_fit_log_function(self, times, *args):
        if len(args) != self.nrFitArgs:
            raise ValueError('wrong number of arguments given for the expected fit function')
        self.__modify_parameter_dict(*args)
        result = Fit_3_states_discrete_rate_changes(times, **self.argumentStorageDict) 
        return np.log(np.array(result[1]) + np.array(result[2]) )




#complete fit of a full protocol for one data set
def Order_fitvalDict(startfitvaldict:dict,boundvaldict:dict, nrWashings:int):
    '''Order the fit vals according to specific order defined here, i.e. how k_ij etc should be used as args, for curve fit
       will only take true values as fit params
       returns ordered fit dictionary with all relevant parameters and boundaries (dicstart,dicboundaries)
    '''
    startFitDict = dict() #contain params and boundaries for which should be fitted
    boundariesDict = dict()
    tupleOfOrder = ('C1_start','C2_start','C3_start', 'k12', 'k21', 'k13', 'k31', 'k23', 'k32', 't_0' ) #Order of keywords used for the args later
    concentrationTuple = ( 'C1_start','C2_start','C3_start')
    for x in tupleOfOrder:
        if x in startfitvaldict:
            startFitDict[x] = startfitvaldict[x]
            if x in boundvaldict:
                boundariesDict[x] = boundvaldict[x]
        else:
            if not x in concentrationTuple:
                startFitDict[x] = [[0 for _ in range(nrWashings)] , [False for _ in range(nrWashings)] ]
            else:
                startFitDict[x] = [0,False]    

            
    return (startFitDict,boundariesDict)

            # startValDict['k12'] = [[1,1.5],[True,False]]


# def Generate_args_fit_function_from_dicts(startfitvaldict:dict,boundvaldict:dict, nrWashings:int):
#     '''Generate a function that takes as argument exactly the  number of args that we want'''
#     # dstart, dbounds = Order_fitvalDict(startfitvaldict,boudnvaldict, nrWashings)
    

#     return lambda x, *args: Fit_3_states_discrete_rate_changes(t_vals, A1_start, A2_start,A3_start, k12,k13,k23,k21,k31,k32,t_0)




def Fit_to_function_multiple_washes(timedata,posidata,startfitvaldict:dict,boundvaldict:dict,nrWashings:int):
    '''Fit to single data posidata.
        Takes dictionaries for startfitvaldict and boundvaldict, 
        startfitvaldict takes a dict that contains all parameters in the form 
            'k12': [[val1,val2], [true/false, true/false]], or 'k12': [[val1,val2], true/false] if true, then takes it as a fit variable, false no fit, default is using all for fitting
        nr Washings is the number of washings used for fitting, will give size of parameters to fit
     '''
     
    # print('used time start' + str(timedata[0]))
    if len(timedata) != len(posidata):
        raise ValueError('inconsistent position and time data given')
    if len(timedata) > 300:
        warnings.warn("data is quite large, best to reduce data set, still trying to fit...")
        
    #Fit_3_states_discrete_rate_changes(t_vals:List[float], A1_start:float, A2_start:float,A3_start:float, k12:List,k13:List,k23:List,k21:List,k31:List,k32:List,t_0:List)
    d1,d2 = Order_fitvalDict(startfitvaldict,boundvaldict,nrWashings)
    usedFitFunction = Storage_of_fit_choice(d1)
    startfitVals = usedFitFunction.current_variable_params_as_list()
    boundsOfFit = usedFitFunction.generate_boundaries_for_fit(boundvaldict)
    # print(startfitVals)
    # print(boundsOfFit)
    popt, pcov = curve_fit(lambda x, *args: usedFitFunction.call_fit_log_function(x,*args), timedata, np.log(posidata),p0=startfitVals,bounds=boundsOfFit)

    # legendTxt = ''
    # for x in tuple(popt):
    #     legendTxt = legendTxt + f'{x:3.3f}' + ','
    # plt.plot(timevals, yvalslist[2], 'r-',label='exp curve')
    # plt.plot(timevals, Wash_2_two_state_function(timevals, *popt), 'r-', label=legendTxt  )

    return popt,pcov, usedFitFunction