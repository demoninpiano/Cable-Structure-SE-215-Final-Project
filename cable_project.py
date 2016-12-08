#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Nov 14 2016
@Author: Dawei Wang

'''

import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from pylab import *
#plt.interactive(True)

import math

class Cable(object):

    def __init__(self,
                 L, l, h, mass, stiffness, **kwargs):
        self.L = L
        self.H_distance = l
        self.V_distance = h
        self.m = mass * 9.81
        self.EA = stiffness
        self.bc1 = kwargs.get('boundaryCondition1',(0,0))
        self.bc2 = kwargs.get('boundaryCondition2',(l,h))
        self.bc1_type = kwargs.get('type of bc1','fix end')
        self.bc2_type = kwargs.get('type of bc2','fix end')
        #self.g = 9.81


    def solveHV(self):
        H, V = symbols('H,V')
        eq1 = (H*self.L/self.EA) + (H/self.m)*(asinh(V/H)\
            - asinh((V - self.m*self.L)/H))\
            -self.H_distance

        eq2 = (self.m*self.L/self.EA) * (V/self.m- self.L/2.0)\
            + (H/self.m) * (((V/H)**2+1)**0.5 - (((V-self.m*self.L)/H)**2+1)**0.5)\
            - self.V_distance
        #nsolve([x + y ** 2 - 4, exp(x) + x * y - 3], [x, y], [1, 1])[1]
        [Hforce, Vforce] = nsolve([eq1,eq2],[H,V],[1,1])
        return Hforce, Vforce

    def addForce(self,H,V):

        self.H = H
        self.V = V

        pass

    def addSegment(self, segment):
        self.d = segment
        self.Le = self.H_distance*(1+8*(self.d/float(self.H_distance))**2)

        pass

    def getCableType(self):
        self.Lambda = (self.m*self.H_distance/self.H)**2*self.H_distance/float(self.H*self.Le/float(self.EA))


        if self.Lambda < 24:
            self.type = 'taut flat cable'

        elif self.Lambda >24:
            self.type = 'heavy cable'



class Load(object):

    '''

    Inherent Parent Load class. Function overwrited in sub-class

    '''

    def solveHV(self,cable):

        pass

    def solveTension(self, cable):

        pass

    def solveDisplacement_X(self, cable):

        pass

    def solveDisplacement_Z(self, cable):

        pass

    def solveHVLinear(self,cable):

        pass



    pass

class singlePoint_Load(Load):

    def __init__(self,
                 location, load):
        self.location = location
        self.magnitude = load

    def solveHV(self, cable):
        s = np.arange(0, cable.L + 0.05, 0.05)
        H, V = symbols('H,V')
        eq1 = (H*s[0]/cable.EA) + (H/cable.m)\
                     * (asinh(V/H) - asinh((V - cable.m*s[0])/H))

        eq2 = (cable.m*s[0]/cable.EA) * (V/cable.m - s[0]/2)\
                     + (H/cable.m)*((1+(V/H)**2)**0.5-
                    (1+((V-cable.m*s[0])/H)**2)**0.5)

        eq3 = (H * s[-1] / cable.EA) + (H / cable.m) \
                     * (asinh(V / H) - asinh((V - cable.m * s[-1] - self.magnitude) / H))\
                     + (H / cable.m)*(asinh((V-self.magnitude-cable.m*self.location)/H )\
                     - asinh((V-cable.m*self.location)/H)) - cable.H_distance

        eq4 = (cable.m*s[-1]/cable.EA) * (V/cable.m - s[-1]/2) \
                     + (H / cable.m) * ((1 + (V / H) ** 2) ** 0.5 -
                            (1 + ((V - cable.m * s[-1]-self.magnitude) / H) ** 2) ** 0.5)\
                     + (H/cable.m) * ( self.magnitude/H*cable.m/cable.EA*(self.location-s[-1])
                     + (1+((V-cable.m*self.location-self.magnitude)/H)**2)**0.5
                     - (1+((V-cable.m*self.location)/H)**2)**0.5) - cable.V_distance
        [Hforce, Vforce] = nsolve([eq1, eq2, eq3, eq4], [H, V], [500, 500])

        cable.loadH = Hforce
        cable.loadV = Vforce

        return Hforce, Vforce

        pass

    def solveTension(self, cable):
        s = np.arange(0, cable.L+0.05, 0.05)
        T = np.zeros_like(s)

        for i in range(len(s)):
            if s[i] < self.location:
                T[i] = (self.H**2 + (self.V - cable.m*s[i])**2)**0.5
            else:
                T[i] = (self.H**2 + (self.V - cable.m*s[i] - self.magnitude)**2)**0.5

        return T

    def solveDisplacement_X(self,cable):
        s = np.arange(0, cable.L+0.05, 0.05)
        X = np.zeros_like(s)

        for i in range(len(s)):
            if s[i] < self.location:
                X[i] = (cable.loadH*s[i]/cable.EA) + (cable.loadH/cable.m)\
                     * (asinh(cable.loadV/cable.loadH) - asinh((cable.loadV - cable.m*s[i])/cable.loadH))
            else:
                X[i] = (cable.loadH * s[i] / cable.EA) + (cable.loadH / cable.m) \
                     * (asinh(cable.loadV / cable.loadH) - asinh((cable.loadV - cable.m * s[i] - self.magnitude) / cable.loadH))\
                     + (cable.loadH / cable.m)*(asinh((cable.loadV-self.magnitude-cable.m*self.location)/cable.loadH )\
                     - asinh((cable.loadV-cable.m*self.location)/cable.loadH))

        return X


    def solveDisplacement_Z(self, cable):
        s = np.arange(0, cable.L + 0.05, 0.05)
        Z = np.zeros_like(s)

        for i in range(len(s)):
            if s[i] < self.location:
                Z[i] = (cable.m*s[i]/cable.EA) * (cable.loadV/cable.m - s[i]/2)\
                     + (cable.loadH/cable.m)*((1+(cable.loadV/cable.loadH)**2)**0.5-
                                          (1+((cable.loadV-cable.m*s[i])/cable.loadH)**2)**0.5)
            else:
                Z[i] = (cable.m*s[i]/cable.EA) * (cable.loadV/cable.m - s[i]/2) \
                     + (cable.loadH / cable.m) * ((1 + (cable.loadV / cable.loadH) ** 2) ** 0.5 -
                            (1 + ((cable.loadV - cable.m * s[i]-self.magnitude) / cable.loadH) ** 2) ** 0.5)\
                     + (cable.loadH/cable.m) * ( self.magnitude/cable.loadH*cable.m/cable.EA*(self.location-s[i])
                     + (1+((cable.loadV-cable.m*self.location-self.magnitude)/cable.loadH)**2)**0.5
                     - (1+((cable.loadV-cable.m*self.location)/cable.loadH)**2)**0.5)

        cable.addSegment(Z[cable.H_distance*10])
        return Z

    def solveHVLinear(self, cable):
        self.p_prime = self.magnitude / float(cable.m * cable.H_distance)

        cable.linearH = cable.m*cable.H_distance**2/8.0/float(cable.d)

        cable.getCableType()

        print 'The cable type is: ' +cable.type
        print '\n'
        print "The value of p' is:" +str(self.p_prime)

        print "The value of lambda squre is:" +str(cable.Lambda)

        x1 = self.location/float(cable.H_distance)

        self.h_prime = 6 * self.p_prime * x1*(1-x1)/(1+(12/float(cable.Lambda)))

        self.linearh = self.h_prime * cable.linearH

        print 'Extra Linear solution of Horizontal force due to external load is : '+str(self.linearh)
        print 'Linear solution of Horizontal force due to external load is : '+str(self.linearh+cable.linearH)


class multiplePointLoad(Load):

    '''

    subclass of inherent parent Load class
    Furture work of response to multiple points load implement here

    '''
    pass


def getCableResponse(cable, load):

    H, V = cable.solveHV()

    print 'The self-weight Horizontal Force is: ' + str(H)
    print 'The self-weight Vertical Force is: ' + str(V)
    print "\n"
    cable.addForce(H, V)

    H, V = load.solveHV(cable)
    print 'The response Non-linear Horizontal Force is: ' + str(H)
    print 'The response Non-linear Vertical Force is: ' + str(V)

    xDisp = load.solveDisplacement_X(cable)
    zDisp = load.solveDisplacement_Z(cable)


    plt.figure(1)
    plt.grid()
    plt.plot(xDisp, -zDisp)
    plt.axis([0,cable.H_distance,-20,0])
    ax = plt.gca()
    ax.set_autoscale_on(False)
    plt.xlabel('Horizontal Length')
    plt.ylabel('Vertical Displacement of Cable')
    plt.annotate('Max Segment of Cable', xy=(load.location, -5), xytext=(30, -10),
                 arrowprops=dict(facecolor='black', shrink=0.05),
                 )
    plt.show()

    print 'The response cable segment is: '+str(cable.d)+'\n'


    load.solveHVLinear(cable)

    pass








