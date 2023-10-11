#Numeric solution of arbitrary first order rate equations, with discrete changes of rate eq.

import numpy as np
from bisect import bisect
from typing import List
from scipy.linalg import expm
def First_order_rate_eq_solution(t_vals:List[float],t_0:float, M_mat:np.ndarray, A_start:np.ndarray ) ->List[List[np.ndarray]]:
    '''Solution of the rate equation for n concentrations
        solves the linear rate eq. eq. dA/dt = M * A, with M the rate matrix, A=(a_1,a_2,...,a_n) the concentrations
        returns a list of n lists that contain the corresponding concentrations over time'''
    if not M_mat.ndim == 2:
        raise TypeError('matrix not the right dimension')
    if not M_mat.shape[0] == M_mat.shape[1]:
        raise TypeError('need symmetric rate matrix')
    if not A_start.ndim == 1:
        raise TypeError('Only 1 axis allowed for A values')
    if not A_start.shape[0] == M_mat.shape[0]:
        raise TypeError('rates have different size to start values')
    finallist = [[] for x in A_start]
    for x in t_vals:
        if x <= t_0:
            for (y,z) in zip(finallist, A_start):
                y.append(z)
        else:
            nextvals = expm(M_mat * (x-t_0)).dot(A_start)
            for (y,z) in zip(finallist,nextvals):
                y.append(z)

    return np.array(finallist)


def Discrete_rate_changes_first_order_rate_eq(t_vals:List[float],t_changes:List[float],M_mats:List[np.ndarray],A_start:np.ndarray) ->List[List[np.ndarray]]:
    '''Solves the first order rate equations for discrete rate changes and arbitrary, fixed, number of concentrations/states
    Input: t_vals: list of calc. times, t_changes: list of time changes, M_mats: list of rate matrices, A_start: start concentration list
    returns: list of lists containing the concentrations over time.
    '''
    if not all(x.ndim == 2 for x in M_mats):
        raise TypeError('rate matrices are not all of the same dimension')
    if not all(x.shape[0]==x.shape[1] for x in M_mats):
        raise TypeError('rate matrices not symmetric')
    if not A_start.ndim == 1:
        raise TypeError('Only 1 axis allowed for A values')
    if not all(A_start.shape[0] == x.shape[0] for x in M_mats):
        raise ValueError('rates have different size to start values')
    if not len(t_changes) == len(M_mats):
        raise TypeError('inconsistent lists of changes and rate matrices')
    # if not all(t_vals[i+1] >= t_vals[i] for i in range(len(t_vals)-1)):
    #     raise ValueError('ascending time list needed')

    i = 0
    finallist = [[] for x in A_start]

    A_intermediary = A_start.astype('float')
    for time in t_vals:
        stepResult = []
        
        if time < t_changes[0]:
            for j in range(len(finallist)):
                finallist[j].append(A_intermediary[j])
            continue

        if  i < len(t_changes) and time > t_changes[i]:
            i = i + 1
            for j in range(len(finallist)):
                A_intermediary[j] = finallist[j][-1]  


        stepResult = First_order_rate_eq_solution(t_vals=[time],t_0 = t_changes[i - 1], M_mat = M_mats[i-1], A_start = A_intermediary)
        for x,y in zip(finallist,stepResult):
            x.append(y[0]) 
    return finallist
        
def rate_matrix_creator_3_states(k12:float,k21:float,k13:float,k31:float,k23:float,k32:float) -> np.array:  
        ''''Rate matrix for three different states'''
        return  np.array([[- k12 - k13, k21, k31]
                                ,[k12, - k21 - k23, k32]
                                ,[k13, k23, - k31 - k32]])




def Fit_3_states_fix_rates_function(t_vals:List[float], C1_start, C2_start, C3_start, k12,k13,k23,k21,k31,k32, t_0):
    '''General solution of the three state rate equation
        Takes a list of time values t_vals and evolves the rate eq. with the parameters A1,A2,C3_start and te=he rates k_ij. 
        The run starts at time t_0, if there are t_vals that are smaller than t_0, then there is assumed to be only vanishing rates
        returns the full solution of all three concentrations
    '''
    M_matrix = rate_matrix_creator_3_states(k12=k12,k21=k21,k13=k13,k31=k31, k23=k23, k32=k32)
    if isinstance(C1_start, list) and isinstance(C2_start, list) and isinstance(C3_start, list):
        A_vec = np.array([C1_start[0],C2_start[0],C3_start[0]],dtype=np.float64)
    elif isinstance(C1_start, float) and isinstance(C2_start, float) and isinstance(C3_start, float):
        A_vec = np.array([C1_start,C2_start,C3_start],dtype=np.float64)
    else: 
       raise TypeError('Inconsistent definition of the concentration start values')

    return First_order_rate_eq_solution(t_vals, t_0, M_matrix, A_vec)


def Fit_3_states_discrete_rate_changes(t_vals:List[float], C1_start:float, C2_start:float,C3_start:float, k12:List,k13:List,k23:List,k21:List,k31:List,k32:List,t_0:List) ->List[List]:
    ''''Receive all concentrations evolved over time given a vector of start times t_0
        It is assumed, that there is a discrete number of rate changes over time, i.e. len(t_0) = len(k_ij) for all i, j
        t_vals: ordered list of times in the evolution of the system
        C1_start: Start concentration of the first state
        C2_start: Start concentration of the second state
        C3_start: Start concentration of the third state
        k_ij: list of changes of the rates going from state i to j
        returns a list of lists that contain each the concentration over the times, for the three different concentrations [[c1],[c2],[c3]]
        The size of the first list is given by the size of the t_0/ k_ij lists
    '''
    if len(t_0) != len(k12) or len(t_0) != len(k13) or len(t_0) != len(k23) or len(t_0) != len(k32) or len(t_0) != len(k31) or len(t_0) != len(k21):
        raise ValueError('inconsistently defined list of parameters')
    
    #store all solutions for the concentrations here
    resultslist = []

    #position of of when the time list is in the next interval of t_0
    intervalSwitches = [] 
    for x in t_0:
        intervalSwitches.append(bisect(t_vals,x))
    intervalSwitches.append(len(t_vals))    #end of the data
    #storage of the start positions needed for the next washing, contains solution of the corresponding endtime using  values of i-1 and ending at t_0[i]
    washingstart = [] #need to be calculated extra, as t_vals does not necessarily contain the washing times (t_0's)

    for i in range(len(t_0)):
      t_interval = []
      if i == 0: 
         t_interval = t_vals[0:intervalSwitches[i+1]]
         results = Fit_3_states_fix_rates_function(t_interval, C1_start, C2_start,C3_start,k12=k12[i],k21=k21[i],k13=k13[i], k31=k31[i],k23=k23[i],k32=k32[i],t_0=t_0[i])
         resultslist.append(results)
         #Storage of the start concentration for each washing in dependence of the chosen washing times t_0, only relevant if more than 1 washing
         if len(t_0) > 1:
            washingstart.append(Fit_3_states_fix_rates_function( [t_0[i+1]], C1_start, C2_start,C3_start,k12=k12[i],k21=k21[i],k13=k13[i], k31=k31[i],k23=k23[i],k32=k32[i],t_0=t_0[i]))
      else:
         t_interval=t_vals[intervalSwitches[i]:intervalSwitches[i+1]]
         results = Fit_3_states_fix_rates_function(t_interval, washingstart[-1][0][0], washingstart[-1][1][0],washingstart[-1][2][0],k12=k12[i],k21=k21[i],k13=k13[i], k31=k31[i],k23=k23[i],k32=k32[i],t_0=t_0[i])
         resultslist.append(results)
         if len(t_0) > i + 1:
            washingstart.append(Fit_3_states_fix_rates_function( [t_0[i+1]], washingstart[-1][0], washingstart[-1][1],washingstart[-1][2],k12=k12[i],k21=k21[i],k13=k13[i], k31=k31[i],k23=k23[i],k32=k32[i],t_0=t_0[i]))  

    c1 = []
    c2 = []
    c3 = []
    firstskip = 0
    for x in resultslist:
       c1.extend( x[0])
       c2.extend( x[1])
       c3.extend( x[2])
        
    return [c1,c2,c3]


def Fit_3_states_discrete_rate_changes_log(t_vals:List[float], C1_start:float, C2_start:float,C3_start:float, k12:List,k13:List,k23:List,k21:List,k31:List,k32:List,t_0:List) ->List[List]:
   '''Calls Fit_3_states_discrete_rate_changes_log and outputs the results with log applied to them'''
   return np.log(Fit_3_states_discrete_rate_changes(t_vals,C1_start,C2_start,C3_start, k12,k13,k23,k21,k31,k32,t_0))


