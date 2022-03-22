###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        Define the class to be associated with an invariancy region.
#
################################################################################

from matrix_manipulation import *
from read_problem import printMatrix
from collections import Counter
from math import floor, ceil
import time
import copy


# Define the Invariancy Region Class

class InvRgn:
    def __init__(   self, 
                    pari, 
                    gMatrix, 
                    basis, 
                    xVar, 
                    startingPoint, 
                    epsilon, 
                    paramSpace,
                    endPoints):
        self.grads      = [ [] for _ in range(len(basis) + len(paramSpace)) ]
        self.defIneq    = [ pari.zero() for _ in range(len(basis) + len(paramSpace)) ]
        self.eps        = epsilon
        self.tableau    = gMatrix
        self.basis      = basis
        self.xVar       = xVar
        self.startPnt   = startingPoint
        self.rhs        = []
        self.point      = []
        self.endPoints  = copy.deepcopy(endPoints)
        self.GetIneqAndGradients(pari, paramSpace, True)


    # Processing the invariancy region is finished. Only store the basis and the
    # defining constraints. No need to continue storing the entire tableau
    def Finalize(self):
        del self.tableau
    
    # Getters
    def Tableau(self):
        return self.tableau
        
    def Basis(self):
        return self.basis
    
    def RHS(self):
        if len(self.rhs) == 0:
            for row in self.tableau:
                self.rhs.append(row[-1])
        return self.rhs
        
    def Grads(self):
        return self.grads
        
    def DefIneq(self):
        return self.defIneq
        
    def EndPoints(self):
        return self.endPoints
    
    # Use Polynomial Roots to Compute the Endpoints of an Interval
    def GetExtremes(self, pari):
        retVal = 0
        leftVal = self.endPoints[0]
        rightVal = self.endPoints[1]
        startVal = self.startPnt[0]
#        print("midpoint: ",startVal)

        for k in range(len(self.defIneq)):
            if pari.poldegree(self.defIneq[k]) > 0:
                roots = pari.polrootsreal(self.defIneq[k], [floor(leftVal), ceil(rightVal)])
#                print("roots: ",roots)
                roots = Counter(roots)
                for r, mult in roots.items():
                    if mult % 2 != 0:
                        if r > leftVal and r < startVal:
                            leftVal = r
                        elif r < rightVal and r > startVal:
                            rightVal = r
                        elif r == startVal:
                            if pari.subst(pari.deriv(self.defIneq[k]), self.xVar[0], startVal) < 0:
                                leftVal = r
                            else:
                                rightVal = r
                            
        self.endPoints[0] = leftVal
        self.endPoints[1] = rightVal
        
        return leftVal, rightVal
        
    
    # Use pari to compute the defining inequalities of the invariancy region 
    # (stored in less-than-or-equal-to form). Then, compute the gradient of each
    # These are used for the ellipsoid method.
    def GetIneqAndGradients(self, pari, paramSpace, storeGrads):
        zeros = [0]*len(self.xVar)
        for i in range(len(self.basis)):
            val = pari.substvec(pari.denominator(self.tableau[i][-1]), self.xVar[0:-1], self.startPnt)
            if val > 0.0:
                self.defIneq[i] = -1*pari.numerator(self.tableau[i][-1])
            else:
                self.defIneq[i] = pari.numerator(self.tableau[i][-1])
            if storeGrads:
                for v in self.xVar:
                    self.grads[i].append(pari.deriv(self.defIneq[i],v))
        for i in range(len(paramSpace)):
            self.defIneq[i + len(self.basis)] = paramSpace[i][0] - paramSpace[i][1]
            if storeGrads:
                for v in self.xVar:
                    self.grads[i + len(self.basis)].append(pari.deriv(self.defIneq[i + len(self.basis)],v))
