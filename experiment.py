#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Nov 14 2016
@Author: Dawei Wang

'''
import numpy as np
from sympy import *
import math
from cable_project import *

if __name__ == "__main__":
    '''
    #############Benchmark 1 Calculation#####################
    ########Homework Benchmark. Lambda squre = 644, p' = 0.27
    '''
    EA = 10**7
    L = 51.05
    mass = 15.0
    l = 50                  #Horizontal Displacement
    h = 0.0
    C1 = Cable(L,l,h,mass,EA,boundaryCondition1 = (0,0),boundaryCondition2 = (l,h))

    load1 = singlePoint_Load(25, 2000)
    getCableResponse(C1, load1)
    '''
    #############Benchmark 2 Calculation#####################
    ########Heavy Cable, p' = 0.097 < 0.1####################
    '''
    EA = 10**7
    L = 51.05
    mass = 15.0
    l = 50                  #Horizontal Displacement
    h = 0.0
    C1 = Cable(L,l,h,mass,EA,boundaryCondition1 = (0,0),boundaryCondition2 = (l,h))

    load1 = singlePoint_Load(25, 2000/0.28*0.1)
    getCableResponse(C1, load1)
    '''
    #############Benchmark 3 Calculation#####################
    ########Taut Flat Cable, p' = 0.1 < 100####################
    '''
    #print str(3)+'dw'
    EA = 10**5
    L = 100
    mass = 20
    l = 100  # Horizontal Displacement
    h = 0.0
    C1 = Cable(L, l, h, mass, EA, boundaryCondition1=(0, 0), boundaryCondition2=(l, h))

    load1 = singlePoint_Load(50, 2000)
    getCableResponse(C1, load1)
    '''
    #############Benchmark 4 Calculation#####################
    ########Taut Flat Cable, p' = 0.1 > 100##################
    '''
    #print str(4) + 'dw'
    EA = 10 ** 5
    L = 1
    mass = 20
    l = 1  # Horizontal Displacement
    h = 0.0
    C1 = Cable(L, l, h, mass, EA, boundaryCondition1=(0, 0), boundaryCondition2=(l, h))

    load1 = singlePoint_Load(0.5, 20000)
    getCableResponse(C1, load1)
