#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 16:19:16 2016

@author: tom
"""

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

with open('interesting/interesting.txt') as f:
    primers = f.read().splitlines()


oryzr = importr("oryzr")

bof = robjects.Vector(primers)

foo = oryzr.LocToRefSeq(bof)

# print(foo.rx2(1)[0][i])
# print(foo.rx2(1)[1][i])

print(foo.rx2(2))

print(foo)

# make a dict of successful refSeqIds
myDict = {}
for i in range(len(foo.rx2(1).rx2(1))):
    # print(i)
    myDict[foo.rx2(1)[0][i]] = foo.rx2(1)[1][i]

print(myDict)

# list of locs that didn't get matched to refSeq
list(foo.rx2(2)[0])
