from xlsx_file_read_extractor import *
import matplotlib.pyplot as plt


####################################################
####################################################
###get the data used from xlsx file
####################################################
####################################################



filename = "../230911_data_for_Michael" #xlsx file that contains all data in the different sheets
#extract the data to analyse from the xlsx file
sheetNr = 12
#TrpR: [2,10,12]
#CoCl:  [4,6,8,14,16]

used_workbook = Get_workbook_from_xlsx_file(filename, True)
usedSheet = used_workbook[used_workbook.sheetnames[sheetNr]]   #sheet decides the experiment, use 2+ for experiments, 2/3, 4/5 etc are the same, just different channels
dataDict = Get_data_from_sheet(usedSheet) #contains all data in form of dictionary, splitting it up by the used concentrations (and also time extra)
allkeys = list(dataDict.keys()) #keys for the data dict
print('name of the used sheet:')
print(usedSheet.title)
print('found key names:')
for x in allkeys:
    print(x)









