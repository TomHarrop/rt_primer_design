# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 16:19:16 2016

@author: tom
"""

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

oryzr = importr("oryzr")

bof = robjects.Vector(['LOC_Os07g46670', 'LOC_Os02g47660', 'LOC_Os02g40410'])

foo = oryzr.LocToRefSeq(bof)

# print(foo.rx2(1)[0][i])
# print(foo.rx2(1)[1][i])

myDict = {}
for i in range(len(foo.rx2(1)) - 1):
    # print(i)
    myDict[foo.rx2(1)[0][i]] = foo.rx2(1)[1][i]

print(myDict)
