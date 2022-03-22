###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        To read in and solve an instance of the uni-parametric
#                   Linear Complementarity Problem (upLCP) using a binary search
#                   style methodology.
#
################################################################################

# Getting Started
from cypari2 import Pari
import sys
import re
import logging
import time
import multiprocessing
from read_flags import *
from read_problem import *
from matrix_manipulation import *
from crisscross import *
from up_inv_region import *
import random
import os

# Initialize pari
pari = Pari()

# Declare Variables, etc
numVar          = 0
numParam        = 0
gxInitialized   = False
mIsNumeric      = True
parallelStart   = False
showProgress    = True
gMatrix         = 0
xVar            = 0
paramSpace      = []
startingPoint   = []
epsilon         = 0.000001
parallel        = True
feasible        = False
numThreads      = multiprocessing.cpu_count()
originalGmatrix = []
originalBasis   = []
outputFilename  = "Solution.txt"

# Define function for parallel processing
def ProcessQ(q, finalPartition, created, finished, lock, numThreads):
    if showProgress:
        print("Activating thread", os.getpid())
    while True:
        interval, curBasis, curMat = q.get(block=True) #block=True means make a blocking call to wait for items in queue
        if interval is None:
            break

        if showProgress:
            print("Thread", os.getpid(), "is processing interval", interval)
        
        mult = 0.5
        point = [mult*interval[0] + (1.0 - mult)*interval[1], 0]
        basis, mat, feasible = CrissCross(pari, logging, numVar, curMat, xVar, point, epsilon, curBasis)

        if not feasible:
            sys.exit("Criss Cross failed. Exiting.")

        rgn = InvRgn(pari, mat, basis, xVar, point, epsilon, paramSpace, interval)
        lval, rval = rgn.GetExtremes(pari)
        finalPartition.put(rgn)

        if lval - interval[0] > epsilon:
            q.put( ([interval[0], lval], copy.deepcopy(basis), copy.deepcopy(mat)) )
            with lock:
                created.value += 1
        if rval - interval[1] < -epsilon:
            q.put( ([rval, interval[1]], copy.deepcopy(basis), copy.deepcopy(mat)) )
            with lock:
                created.value += 1
        with lock:
            finished.value += 1
        if created.value == finished.value:
            for i in range(numThreads):
                q.put((None,None,None))

# Set parameters using command line flags
if len(sys.argv) > 2:
    numThreads,     \
    parallelStart,  \
    showProgress  = ReadFlags(  sys, 
                                logging, 
                                numThreads,
                                parallelStart,
                                showProgress)

if numThreads <= 1:
    parallelStart = False
elif numThreads > multiprocessing.cpu_count():
    numThreads = multiprocessing.cpu_count()

# Read in the problem instance
t = time.time()
numVar, numParam, gMatrix, xVar, paramSpace, mIsNumeric, probType, numRow, numCol = ReadFile(   pari, 
                                                                                                sys, 
                                                                                                logging, 
                                                                                                re, 
                                                                                                numVar, 
                                                                                                numParam, 
                                                                                                gMatrix, 
                                                                                                xVar, 
                                                                                                paramSpace, 
                                                                                                gxInitialized, 
                                                                                                mIsNumeric)
originalGmatrix = [row[:] for row in gMatrix] #deep copy
originalBasis = list(range(numVar))

totalTime = time.time() - t
if showProgress:
    print("Time to read problem: " + str(round(totalTime, 2)) + "s")

if mIsNumeric:
    logging.warning("Warning: The data entered consists of an M matrix containing no parameters. While the method implemented here is applicable for this problem, a more efficient procedure exists. See Adelgren and Wiecek's 'A two phase algorithm for the multiparametric linear complementarity problem' (2016). This method may implemented here in a future release, but is not as of now. Continuing ... ")

# Initialize
startingPoint = [float((paramSpace[1][1] - paramSpace[0][1])/2.0), 0]
ones = [1 for i in xVar]
endPoints = []
for i in range(len(paramSpace[0])):
    if pari.substvec(paramSpace[i][0], xVar, ones) > 0:
        endPoints.append(paramSpace[i][1])
    else:
        endPoints.append(-paramSpace[i][1])
leftVal = startingPoint[0]
rightVal = startingPoint[0]
endPoints.sort()

# Set and initialize tools for parallelization
m = multiprocessing.Manager()
q = m.Queue()
finalPartition = m.Queue()
created = m.Value('i', 0)
finished = m.Value('i', 0)
lock = m.Lock()
if parallelStart:
    leftEnd = endPoints[0]
    n = numThreads - 1
    for i in range(n):
        rightEnd = (i+1)*endPoints[1]/(n*1.0)
        newInterval = [leftEnd, rightEnd]
        q.put( (copy.deepcopy(newInterval), copy.deepcopy(originalBasis), copy.deepcopy(gMatrix)) )
        with lock:
            created.value += 1
        leftEnd = rightEnd
else:
    q.put( (endPoints, copy.deepcopy(originalBasis), copy.deepcopy(gMatrix)) )
    with lock:
        created.value += 1
        
if __name__ == '__main__':
    pool = multiprocessing.Pool(numThreads, ProcessQ, (q, finalPartition, created, finished, lock, numThreads, ))

    # prevent adding anything more to the queue and wait for queue to empty
    pool.close()
    pool.join()
    pool.terminate()


totalTime = time.time() - t

print("Solution Computed. Elapsed Time: " + str(round(totalTime, 2)) + "s")


# Write the solution

outputFile = open(outputFilename, 'w')

k = 1
if probType == "LCP":
    print("The problem entered was an instance of spLCP having the form\n", file = outputFile)
    print("\tw - M(x)z = q(x)\n\tw'z = 0\n\tw,z >= 0\n", file = outputFile)

    mx = max((len(str(ele.Str())) for row in originalGmatrix for ele in row[numVar:-1]))
    print("with M(x) =\n", file = outputFile)
    for row in originalGmatrix:
        print("\t[ " + "  ".join(["{:<{mx}}".format(str((-1*ele).Str()),mx=mx) for ele in row[numVar:-1]]) + " ]", file = outputFile)
        
    mx = max((len(str(ele.Str())) for row in originalGmatrix for ele in row[2*numVar:]))
    print("\nand q(x) =\n", file = outputFile)
    for row in originalGmatrix:
        print("\t[ " + "  ".join(["{:<{mx}}".format(str(ele.Str()),mx=mx) for ele in row[2*numVar:]]) + " ]", file = outputFile)
        
    mx = max((len(str(ele.Str())) for row in paramSpace for ele in row))
    print("\nsubject to the additional restriction that 'x' must satisfy:\n", file = outputFile)
    for row in paramSpace:
        print("\t" + " <= ".join(["{:<{mx}}".format(str(ele.Str()),mx=mx) for ele in row]), file = outputFile)
        
    print("\n\n\n**************************************************************************************************\n\nThe solution was computed in " + str(round(totalTime, 2)) + " seconds and consists of the following regions.\n\n**************************************************************************************************\n\n", file = outputFile)
        
    k = 1
    while not finalPartition.empty():
        rgn = finalPartition.get()
        point = rgn.EndPoints()
#        print(point)
        if point[0] != point[1]:
            print("\n\nRegion " + str(k) + ":\n", file = outputFile)
            rhs = rgn.RHS()
            basis = rgn.Basis()
            mx = max((len(str(row.Str())) for row in rhs))
            for i in range(len(rhs)):
                var = ""
                if basis[i] < numVar:
                    var = "w_" + str(i + 1)
                else:
                    var = "z_" + str(i + 1)
                print("\t" + var + " = " + " ".join(["{:<{mx}}".format(str(rhs[i].Str()),mx=mx)]) + " >= 0 ", file = outputFile)
            print("\n\tValid over:\t" + str('%.15g'%point[0]) + " <= " + str(xVar[0].Str()) + " <= " + str('%.15g'%point[1]), file = outputFile)
            k = k + 1
            
    print("\n\n\n\nNote: The region descriptions above do not include the 'additional restrictions' listed at the top of this document, although these restrictions do, of course, apply to all regions. Additionally, all omitted variables should be assumed to be zero.", file = outputFile)
else:
    print("The problem entered was an instance of sp" + probType + " having the form\n", file = outputFile)
    if probType == "LP":
        print("\tmin \tc(x)'y\n\ts.t.\tA(x)y <= b(x)\n\t    \ty >= 0\n", file = outputFile)
    else:
        print("\tmin \tc(x)'y + (1/2)y'Q(x)y\n\ts.t.\tA(x)y <= b(x)\n\t    \ty >= 0\n", file = outputFile)
        
    mx = max((len(str(ele.Str())) for row in originalGmatrix[numRow:] for ele in row[2*numVar:])) + 1
    print("with c(x) =\n", file = outputFile)
    for row in originalGmatrix[numRow:]:
        print("\t[ " + "  ".join(["{:<{mx}}".format(str(ele.Str()),mx=mx) for ele in row[2*numVar:]]) + " ]", file = outputFile)
        
    if probType == "QP":
        mx = max((len(str(ele.Str())) for row in originalGmatrix[numRow:] for ele in row[(numVar+numRow):-1])) + 1
        print("\nand Q(x) =\n", file = outputFile)
        for row in originalGmatrix[numRow:]:
            print("\t[ " + "  ".join(["{:<{mx}}".format(str((-1*ele).Str()),mx=mx) for ele in row[(numVar+numRow):-1]]) + " ]", file = outputFile)
            
    mx = max((len(str(ele.Str())) for row in originalGmatrix[0:(numCol-1)] for ele in row[(numVar+numRow):-1])) + 1
    print("\nand A(x) =\n", file = outputFile)
    for row in originalGmatrix[0:(numCol-1)]:
        print("\t[ " + "  ".join(["{:<{mx}}".format(str((ele).Str()),mx=mx) for ele in row[(numVar+numRow):-1]]) + " ]", file = outputFile)
        
    mx = max((len(str(ele.Str())) for row in originalGmatrix[0:(numCol-1)] for ele in row[2*numVar:])) + 1
    print("\nand b(x) =\n", file = outputFile)
    for row in originalGmatrix[0:(numCol-1)]:
        print("\t[ " + "  ".join(["{:<{mx}}".format(str(ele.Str()),mx=mx) for ele in row[2*numVar:]]) + " ]", file = outputFile)
        
    mx = max((len(str(ele.Str())) for row in paramSpace for ele in row)) + 1
    print("\nsubject to the additional restriction that 'x' must satisfy:\n", file = outputFile)
    for row in paramSpace:
        print("\t" + " <= ".join(["{:<{mx}}".format(str(ele.Str()),mx=mx) for ele in row]), file = outputFile)
        
    print("\n\n\n**************************************************************************************************\n\nThe solution was computed in " + str(round(totalTime, 2)) + " seconds and consists of the following regions.\n\n**************************************************************************************************\n\n", file = outputFile)
        
    k = 1
    while not finalPartition.empty():
        rgn = finalPartition.get()
        point = rgn.EndPoints()
        if point[0] != point[1]:
            print("\n\nRegion " + str(k) + ":\n", file = outputFile)
            rhs = rgn.RHS()
            basis = rgn.Basis()
            mx = max((len(str(row.Str())) for row in rhs)) + 1
            for i in range(len(rhs)):
                var = ""
                if basis[i] < numVar:
                    if i >= numCol - 1:
                        var = "v_" + str(i + 2 - numCol)
                    else:
                        var = "s_" + str(i + 1)
                else:
                    if i >= numCol - 1: 
                        var = "y_" + str(i + 2 - numCol)
                    else:
                        var = "u_" + str(i + 1)
                print("\t" + var + " = " + " ".join(["{:<{mx}}".format(str(rhs[i].Str()),mx=mx)]) + " >= 0 ", file = outputFile)
            print("\n\tValid over:\t" + str('%.15g'%point[0]) + " <= " + str(xVar[0].Str()) + " <= " + str('%.15g'%point[1]), file = outputFile)
            k = k + 1
            
    print("\n\n\n\nNote 1: Above, 'y' variables represent the original variables, whereas 's' variables are slack variables on the inequality constraints, 'v' variables are duals for the non-negativity restrictions on the 'y' variables, and 'u' variables are duals for the inequality constraints. Additionally, all omitted variables should be assumed to be zero.", file = outputFile)
            
    print("\n\nNote 2: The region descriptions above do not include the 'additional restrictions' listed at the top of this document, although these restrictions do, of course, apply to all regions.", file = outputFile)
        
outputFile.close()

print("Number of intervals in the final partition: " + str(k - 1))


