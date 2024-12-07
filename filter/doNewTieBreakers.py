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

"""
Strategy:
- take only the greatest parsimony score improvement for each putative recombinant node but don't break ties (keep all)
- merge all breakpoint predictions that are adjacent/overlapping and have identical parents
- then take fewest breakpoints first
- if still tied, take smallest 3seq string
- if still tied, take smaller 3seq p-value
- if they have the same parents, take the larger breakpoint interval
- if still tied, larger total number of descendants
"""

##########################
##### MAIN FUNCTIONS #####
##########################

def applyPval():
    myOutString = ''
    with open('combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters3SeqP02.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if float(splitLine[20]) <= 0.2:
                if splitLine[13].startswith('0/') or splitLine[13].startswith('NA') or float(splitLine[13]) < 0.05:
                    myOutString += joiner(splitLine)+'\n'
    open('combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters3seqP02RussPval005.txt','w').write(myOutString)


def doNewTiebreakers():
    # nodeToLeaves = {}
    # with gzip.open('optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.leaves.txt.gz') as f:
    #     for line in f:
    #         splitLine = (line.decode('utf8').strip()).split('\t')
    #         if splitLine[0].isdigit():
    #             nodeToLeaves[int(splitLine[0])] = int(splitLine[1])

    bp1 = {}
    bp2 = {}
    recombToBPs = {}
    recombToStringSize = {}
    recombToLines = {}
    with open('combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters3seqP02RussPval005.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if not int(splitLine[0]) in recombToBPs:
                    recombToBPs[int(splitLine[0])] = []
                    recombToStringSize[int(splitLine[0])] = []
                    recombToLines[int(splitLine[0])] = []
                tempPreStart = 0
                tempStart = 0
                tempMid = 0
                tempEnd = 0
                tempPostEnd = 0
                myInfSites = toInt(splitLine[15].split(','))
                myStart1 = int(splitLine[1].split(',')[0][1:])
                myStart2 = int(splitLine[1].split(',')[1][:-1])
                myEnd1 = int(splitLine[2].split(',')[1][:-1])
                myEnd2 = int(splitLine[2].split(',')[0][1:])
                for k in myInfSites:
                    if k <= myStart1:
                        tempPreStart += 1
                    elif k >= myStart1 and k <= myStart2:
                        tempStart += 1
                    elif k > myStart2 and k < myEnd1:
                        tempMid += 1
                    elif k >= myEnd1 and k <= myEnd2:
                        tempEnd += 1
                    elif k > myEnd2:
                        tempPostEnd += 1
                if tempPreStart == 0 or tempPostEnd == 0:
                    recombToBPs[int(splitLine[0])].append(1)
                    recombToStringSize[int(splitLine[0])].append(len(myInfSites))
                    recombToLines[int(splitLine[0])].append(splitLine)
                else:
                    recombToBPs[int(splitLine[0])].append(2)
                    recombToStringSize[int(splitLine[0])].append(len(myInfSites))
                    recombToLines[int(splitLine[0])].append(splitLine)

    myOutString = ''
    for k in recombToBPs:
        tempB = []
        tempS = []
        tempL = []
        for i in range(0,len(recombToBPs[k])):
            tempB.append(recombToBPs[k][i])
            tempS.append(recombToStringSize[k][i])
            tempL.append(recombToLines[k][i])

        if len(tempB) == 1: # if only one, just print it
            myOutString += joiner(tempL[0])+'\n'
        else:
            if tempB.count(1) == 1:
                myOutString += joiner(tempL[tempB.index(1)])+'\n' # if one best, print that

            else: # else, get all tied and go to next tiebreaker: smallest 3seq string
                if tempB.count(1) == 0:
                    newB = tempB
                    newS = tempS
                    newL = tempL
                elif tempB.count(1) > 1:
                    newB = []
                    newS = []
                    newL = []
                    for i in range(0,len(tempB)):
                        if tempB[i] == 1:
                            newB.append(tempB[i])
                            newS.append(tempS[i])
                            newL.append(tempL[i])

                minS = min(newS)
                if newS.count(minS) == 1:
                    myOutString += joiner(newL[newS.index(minS)])+'\n'
                elif newS.count(minS) > 1:
                    tempB = []
                    tempS = []
                    tempL = []
                    tempLeaves = []
                    tempP = []
                    tempPval = []
                    for i in range(0,len(newS)):
                        if newS[i] == minS:
                            tempB.append(newB[i])
                            tempS.append(newS[i])
                            tempL.append(newL[i])
                            # tempLeaves.append(nodeToLeaves[int(newL[i][3])]+nodeToLeaves[int(newL[i][6])])
                            tempP.append(set([int(newL[i][3]),int(newL[i][6])]))
                            tempPval.append(float(newL[i][-1]))

                    minPval = min(tempPval)
                    if tempPval.count(minPval) == 1:
                        myOutString += joiner(tempL[tempPval.index(minPval)])+'\n'
                    elif tempPval.count(minPval) > 1:
                        if getNumUnique(tempP) == 1: # if all have the same parents, print biggest breakpoint interval
                            myOutString += joiner(getBiggestBreakpointInterval(tempL))+'\n'
                        else:
                            print(tempP)
                            newB = []
                            newS = []
                            newL = []
                            newLeaves = []
                            newP = []
                            newPval = []
                            for i in range(0,len(tempPval)):
                                if tempPval[i] == minPval:
                                    newB.append(tempB[i])
                                    newS.append(tempS[i])
                                    newL.append(tempL[i])
                                    # newLeaves.append(nodeToLeaves[int(tempL[i][3])]+nodeToLeaves[int(tempL[i][6])])
                                    newP.append(set([int(tempL[i][3]),int(tempL[i][6])]))

                            minLeaves = min(newLeaves)
                            if newLeaves.count(minLeaves) == 1:
                                myOutString += joiner(newL[newLeaves.index(minLeaves)])+'\n'
                            elif newLeaves.count(minLeaves) > 1:
                                tempB = []
                                tempS = []
                                tempL = []
                                tempLeaves = []
                                tempP = []
                                tempPval = []
                                for i in range(0,len(newLeaves)):
                                    if newLeaves[i] == minLeaves:
                                        tempB.append(newB[i])
                                        tempS.append(newS[i])
                                        tempL.append(newL[i])
                                        # tempLeaves.append(nodeToLeaves[int(newL[i][3])]+nodeToLeaves[int(newL[i][6])])
                                        tempP.append(set([int(newL[i][3]),int(newL[i][6])]))
                                myOutString += joiner(getBiggestBreakpointInterval(tempL))
    open('combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClustersNewTiebreak3seqP02RussPval005.txt','w').write(myOutString)


##########################
#### HELPER FUNCTIONS ####
##########################

def toInt(myList):
    myReturn = []
    for k in myList:
        myReturn.append(int(k))
    return(myReturn)

def getBiggestBreakpointInterval(myList):
    temp = []
    for k in myList:
        temp.append(numpy.sum(toInt(k[1][1:-1].split(',')))+numpy.sum(toInt(k[2][1:-1].split(','))))
    return(myList[temp.index(max(temp))])


def getNumUnique(myList):
    myReturn = 1
    for i in range(0,len(myList)):
        for j in range(i+1,len(myList)):
            if myList[i] != myList[j]:
                #print(myList[i], myList[j])
                myReturn += 1
    return(myReturn)


def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return('\t'.join(newList))

#########################
##### FUNCTION CALL #####
#########################

def main():
    applyPval()
    doNewTiebreakers()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit
















