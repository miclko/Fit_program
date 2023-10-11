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