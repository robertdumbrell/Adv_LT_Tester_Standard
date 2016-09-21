
#-------------------------------------------------------------------------------
# Name:        Numerical Differentiation
# Purpose:     To Differentiate
#
# Author:      Mattias
#
# Created:     26/06/2013
# Copyright:   (c) Mattias 2012
# Licence:     <my licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env pytho
#import numpy as np
from numpy import linspace, matrix,interp,zeros,float
import scipy.sparse as SS
from scipy.sparse.linalg import spsolve
#from matplotlib import  *


class Finite_Difference():

    def TwoPointCentral(self,x,y):
#2-point formula
        dyf = zeros(x.shape[0])
        for i in range(y.shape[0]-1):
            dyf[i] = (y[i+1] - y[i])/(x[i+1]-x[i])
#set last element by backwards difference
        dyf[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])
        #plot(x,dyf,'r-',label='2pt-forward diff')
        #semilogy()
        #show()
        return dyf

    def FourPointCentral(self,x,y):

        '''
        Assumes evenly spaced points!!

        calculate dy by 4-point center differencing using array slices

        \frac{y[i-2] - 8y[i-1] + 8[i+1] - y[i+2]}{12h}

        y[0] and y[1] must be defined by lower order methods
        and y[-1] and y[-2] must be defined by lower order methods
        '''

        dy = zeros(y.shape,float) #we know it will be this size
        h = x[1]-x[0] #this assumes the points are evenely spaced!
        dy[2:-2] = (y[0:-4] - 8*y[1:-3] + 8*y[3:-1] - y[4:])/(12.*h)

        dy[0] = (y[1]-y[0])/(x[1]-x[0])
        dy[1] = (y[2]-y[1])/(x[2]-x[1])
        dy[-2] = (y[-2] - y[-3])/(x[-2] - x[-3])
        dy[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])
        return dy

class Regularisation():

    """Provides the first derivitave of something, WITH EQUALLY SPACED POINTS"""
    def FirstDerivative(self,X,Y,lam):

        Xreal = matrix(X).T
        Yreal =self.DerivitiveMatrix(1,Xreal[1]-Xreal[0],matrix(Y).T.shape[0])* matrix(Y).T

        #obtaining Differential matrix
        D=self.DerivitiveMatrix(2,Xreal[1]-Xreal[0],Yreal.shape[0])

        #print X.shape,self.Smoothed(lam,Yreal,D).shape,(X[:-1]+X[1:]).shape
        return interp(X,(X[:-1]+X[1:])/2,self.Smoothed(lam,Yreal,D)) #usually would lose a point, as evaultes inbetween points, so conver back to total number, with error at the end areas

    def DerivitiveMatrix(self,Order,Deltax,n):

        #Have to use SS as numpy can't deal with these big arrays.
        #Also its a smaller demand on the system as it only records/does opperations of values with number
        D= -SS.eye(n-1,n)+SS.eye(n-1,n,1)

        #Computing the correct Order (the order is the differential to make smooth) Matrix
        #It should be 2 orders higher than the higest differential you are going to use
        for i in range(1,Order):
            D= D[:-1,:-i]*D
            #print D


        D=(Deltax[0,0]**(-Order))*D

        #print D

        return D

    def Smoothed(self,lam,Y,D):

        return spsolve((SS.identity(D.T.shape[0]) + lam*D.T*D),Y).T[:,0]



#Comparing
'''
X = linspace(-10,10,10000)

Y =     -X**5   +4*X**4    +5e1*X**3   -1e2*X**2   +5e2*X   + 1e4
Yd =    -5*X**4 +16*X**3   +1.5e2*X**2   -2e2*X      +5e2

#Y = exp(X/0.02585)
#Yd = exp(X/0.02585)/0.02585


plot(X,Regularisation().FirstDerivative(X,Y,0),'b--',label='Regularisation Derivative')
ax.plot(X,Finite_Difference().TwoPointCentral(X,Y),'r-.',label='Two Point Finite Difference')
#ax.plot(X,Finite_Difference().FourPointCentral(X,Y),'r-.',label='Four Point Finite Difference')
#ax.plot(X,Yd,'g--',label='Real derivative')
#grid(True,'major')
#legend(loc=2)
#semilogy()
#twinx()
#plot(X,Y,'k,',label='Original function')
#legend(loc=1)
#semilogy()

show()
'''


