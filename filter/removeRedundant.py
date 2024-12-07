#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 2/1/2018
# compareDatabases.py

import sys
import os
import datetime
import numpy
from numpy import random
import gzip
import math
import re

##########################
##### MAIN FUNCTIONS #####
##########################

def removeRedundantTrios():
    # nodeToLeaves = {}
    # with gzip.open('optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.leaves.txt.gz') as f:
    #     for line in f:
    #         splitLine = (line.decode('utf8').strip()).split('\t')
    #         if splitLine[0].isdigit():
    #             nodeToLeaves[int(splitLine[0])] = int(splitLine[1])

    myTrios = []
    trioToPVal = {}
    trioToSites = {}
    trioToLeaves = {}
    trioToLine = {}
    lc = 0
    with open('combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClustersNewTiebreak3seqP02RussPval005.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            myTrios.append([int(splitLine[0]),int(splitLine[3]),int(splitLine[6])])
            trioToLine[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])] = splitLine
            if splitLine[13].startswith('0/'):
                splitLine[13] = (1.0/float(splitLine[13][2:]))
            elif splitLine[13].startswith('NA'):
                continue
            trioToPVal[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])] = float(splitLine[13])
            trioToSites[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])] = len(splitLine[16])
            # trioToLeaves[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])] = nodeToLeaves[int(splitLine[3])]+nodeToLeaves[int(splitLine[6])]
            lc += 1

    toRemove = {}
    for i in range(0,len(myTrios)):
        for j in range(i+1,len(myTrios)):
            if checkTwo(myTrios[i],myTrios[j]) == True: # if circular logic:
                if trioToPVal[joinerU(myTrios[i])] < trioToPVal[joinerU(myTrios[j])]: # remove case with lower pval
                    toRemove[joinerU(myTrios[j])] = True
                elif trioToPVal[joinerU(myTrios[i])] > trioToPVal[joinerU(myTrios[j])]:
                    toRemove[joinerU(myTrios[i])] = True
                elif trioToPVal[joinerU(myTrios[i])] == trioToPVal[joinerU(myTrios[j])]:

                    if trioToSites[joinerU(myTrios[i])] > trioToSites[joinerU(myTrios[j])]: # remove case with more informative sites
                        toRemove[joinerU(myTrios[i])] = True
                    elif trioToSites[joinerU(myTrios[i])] < trioToSites[joinerU(myTrios[j])]:
                        toRemove[joinerU(myTrios[j])] = True
                    elif trioToSites[joinerU(myTrios[i])] == trioToSites[joinerU(myTrios[j])]:

                        if trioToLeaves[joinerU(myTrios[i])] < trioToLeaves[joinerU(myTrios[j])]:
                             toRemove[joinerU(myTrios[i])] = True
                        if trioToLeaves[joinerU(myTrios[i])] > trioToLeaves[joinerU(myTrios[j])]:
                             toRemove[joinerU(myTrios[j])] = True
                        elif trioToLeaves[joinerU(myTrios[i])] == trioToLeaves[joinerU(myTrios[j])]:
                            print(myTrios[i], myTrios[j], trioToPVal[joinerU(myTrios[i])])

    myOutString = ''
    for t in trioToLine:
        if not t in toRemove:
            myOutString += joiner(trioToLine[t])+'\n'
    open('combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClustersNewTiebreak3seqP02RussPval005RemoveCircular.txt','w').write(myOutString)



##########################
#### HELPER FUNCTIONS ####
##########################

def checkTwo(list1, list2):
    mySame = 0
    for k in list1:
        if k in list2:
            mySame += 1
    if mySame == 3:
        return(True)
    elif mySame == 2:
        if list1[0] not in list2 and list2[0] not in list1:
            return(False)
        else:
            return(True)
    else:
        return(False)


def toInt(myList):
    myReturn = []
    for k in myList:
        myReturn.append(int(k))
    return(myReturn)

def joinerU(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return('_'.join(newList))

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return('\t'.join(newList))

#########################
##### FUNCTION CALL #####
#########################

def main():
    removeRedundantTrios()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit















