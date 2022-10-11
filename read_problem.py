###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        Read in and solve an instance of the multiparametric
#                   Linear Complementarity Problem (mpLCP).
#
################################################################################


# Define Functions

# Print a matrix of Pari elements in a user-readable format
#
# Inputs:   matrix -- the matrix of Pari elements to be printed
def printMatrix(matrix):
    print('\n'.join([''.join(['{:20}'.format(str(item.Str())) for item in row]) 
      for row in matrix]))

# Print a matrix of Pari elements that represents a tableau as a LaTeX table
#
# Inputs:   matrix -- the matrix of Pari elements to be printed
def printLatexTableau(basis, matrix):
    varNames = []
    for i in range((len(matrix[0]) - 1)//2):
        varNames.append('w_' + str(i+1))
    for i in range((len(matrix[0]) - 1)//2):
        varNames.append('z_' + str(i+1))

    print("\\renewcommand{\\arraystretch}{1.5}")
    print("$\\begin{array}{c|" + 'r'*(len(matrix[0]) - 1) + '|r|}')
    print("\\cline{2-" + str(len(matrix[0]) + 1) + "}")
    print(' & ' + ' & '.join([str(n) for n in varNames]) + r' & \text{RHS}\\')
    print("\\cline{2-" + str(len(matrix[0]) + 1) + "}")
    print('\n'.join([varNames[basis[i]] + ' & ' + '&'.join(['{}'.format(str(item.Strtex())) for item in matrix[i]]) + r'\\'
          for i in range(len(matrix))])) 
    print("\\cline{2-" + str(len(matrix[0]) + 1) + "}")
    print("\\end{array}$")

# Initialize the G matrix and vector of parameters (treated as variables)
#
# Inputs:   pari    --  the Pari environment
#           numVar  --  an integer specifying the number of variables in the 
#                       problem. This becomes the number of rows of the G 
#                       matrix.
#           numParam--  an integer specifying the number of parameters in the 
#                       problem. This value is 2 less than the number of entries
#                       needed in the x vector.
#           gMatrix --  the matrix used to store the problem data
#           xVar    --  the vector used to store the Pari variables that 
#                       represent the problem's parameters
#
# Outputs:  gMatrix
#           xVar
def InitializeGandX(pari, numVar, numParam, gMatrix, xVar):
    #initialize with zeros
    gMatrix = [ [pari.zero() for j in range(2*numVar + 1)] for i in range(numVar)]
    #the first numVar by numVar submatrix is an identity matrix
    for i in range(numVar):
        gMatrix[i][i] = pari.one()
    
    #use Pari variables to represent the parameters
    xVar = []
    for i in range(numParam + 2):
        xVar.append(pari('x_' + str(i + 1)))

    return gMatrix, xVar

# Convert an input decimal value to a rational string representation
#
# Input:    val -- a string representation of a decimal value
#
# Output:   val -- converted to a string representing a rational value
def ConvertToFraction(val):
#    print("Converting a decimal to a rational. Value passed in: " + str(val))
    splitVal = val.strip().split('.')
    for i in range(len(splitVal)):
        if splitVal[i].strip() == "":
            splitVal[i] = "0"
    sign = '+'
    if splitVal[0][0] == "-":
        sign = '-'
    val = splitVal[0] + sign + splitVal[1] + '/' + str(10**len(splitVal[1]))
    return val

# Parse the input file
#
# Input:    Pari    --  the Pari environment
#           sys     --  the variable containing any information passed at the
#                       command line
#           logging --  the logging environment
#           re      --  the re environment (for parsing real expressions)
#           numVar  --  an integer specifying the number of variables in the 
#                       problem. This becomes the number of rows of the G 
#                       matrix.
#           numParam--  an integer specifying the number of parameters in the 
#                       problem. This value is 2 less than the number of entries
#                       needed in the x vector.
#           gMatrix --  the matrix used to store the problem data
#           xVar    --  the vector used to store the Pari variables that 
#                       represent the problem's parameters
#           paramSpace  --  the matrix used to store the constraints 
#                           representing portion of the parameter space to be 
#                           partitioned.
#           gxInitialized   --  a boolean representing whether or not the G 
#                               matrix and x vector have been initialized.
#           mIsNumeric      --  a boolean representing whether or not the M
#                               matrix contains any parameters.
#
# Outputs:  numVar
#           numParam
#           gMatrix
#           xVar
#           paramSpace
#           mIsNumeric
#           probType
def ReadFile(Pari, sys, logging, re, numVar, numParam, gMatrix, xVar, paramSpace, gxInitialized, mIsNumeric):
    # Read the file
    inputFile = open(sys.argv[1])

    lines = inputFile.readlines()

    i = 0
    numRow = 0
    numCol = 0
    while i < len(lines) and numVar <= 0:
        if not lines[i].strip() == "":
            probType = lines[i].strip().upper()
            if lines[i].strip().upper() == "LCP":
                numVar, numParam, gMatrix, xVar, paramSpace, mIsNumeric = ReadLCP(  lines, 
                                                                                    i + 1, 
                                                                                    Pari, 
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
            elif lines[i].strip().upper() == "LP" or lines[i].strip().upper() == "QP":
                numVar, numParam, gMatrix, xVar, paramSpace, mIsNumeric, numRow, numCol = ReadQP(   lines, 
                                                                                                    i + 1, 
                                                                                                    Pari, 
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
            else:
                logging.warning("Data file should start with a specification of the type of problem you wish to solve. Valid values are 'LCP', 'LP', and 'QP'. Proceeding as if problem type is LCP.")
                probType = "LCP"
                numVar, numParam, gMatrix, xVar, paramSpace, mIsNumeric = ReadLCP(  lines, 
                                                                                    0, 
                                                                                    Pari, 
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
        i += 1
            
    return numVar, numParam, gMatrix, xVar, paramSpace, mIsNumeric, probType, numRow, numCol
    

# Read the input file if it is associated with an mpLCP instance
#
# Input:    lines   --  remaining lines of the input file
#           start   --  index of first line to read from lines
#           Pari    --  the Pari environment
#           sys     --  the variable containing any information passed at the
#                       command line
#           logging --  the logging environment
#           re      --  the re environment (for parsing real expressions)
#           numVar  --  an integer specifying the number of variables in the 
#                       problem. This becomes the number of rows of the G 
#                       matrix.
#           numParam--  an integer specifying the number of parameters in the 
#                       problem. This value is 2 less than the number of entries
#                       needed in the x vector.
#           gMatrix --  the matrix used to store the problem data
#           xVar    --  the vector used to store the Pari variables that 
#                       represent the problem's parameters
#           paramSpace  --  the matrix used to store the constraints 
#                           representing portion of the parameter space to be 
#                           partitioned.
#           gxInitialized   --  a boolean representing whether or not the G 
#                               matrix and x vector have been initialized.
#           mIsNumeric      --  a boolean representing whether or not the M
#                               matrix contains any parameters.
#
# Outputs:  numVar
#           numParam
#           gMatrix
#           xVar
#           paramSpace
#           mIsNumeric
def ReadLCP(lines, start, Pari, sys, logging, re, numVar, numParam, gMatrix, xVar, paramSpace, gxInitialized, mIsNumeric):
    i = start
    while i < len(lines):
        if not lines[i].strip() == "":
            if lines[i].strip().upper() == "H" or lines[i].strip().upper() == "NUM_VAR":
                numVar = int(lines[i+1].strip())
                i += 1
            elif lines[i].strip().upper() == "K" or lines[i].strip().upper() == "NUM_PARAM":
                numParam = int(lines[i+1].strip())
                i += 1
            else:
                if not numVar or not numParam:
                    sys.exit("After problem type specification, the data file must next include a specification of the number of variables (via keywords 'h' or 'num_var') and the number of parameters (via keywords 'k' or 'num_param')! Please check file format and try again. Exiting!")
                if not gxInitialized:
                    gxInitialized = True
                    gMatrix, xVar = InitializeGandX(Pari, numVar, numParam, gMatrix, xVar)
                if lines[i].strip().upper() == "M_DATA":
                    i += 1
                    while not lines[i].strip() == "" and lines[i].strip()[0].isdigit():
                        vals = re.findall('[-+]?\d*\.?\d+', lines[i])
                        if len(vals) < 4:
                            sys.exit("Data for the M matrix must contain four comma delimited values: (1) the row index, (2) the column index, (3) the parameter index -- with 0 indicating a constant -- and (4) the coeficient. Please reformat your data and retry. Exiting!")
                        if '.' in vals[3]:
                            vals[3] = ConvertToFraction(vals[3])
                        if vals[2] == '0':
                            gMatrix[int(vals[0]) - 1][numVar + int(vals[1]) - 1] -= Pari(vals[3])
                        else:
                            gMatrix[int(vals[0]) - 1][numVar + int(vals[1]) - 1] -= Pari(vals[3])*xVar[int(vals[2]) - 1]
                            mIsNumeric = False
                        i += 1
                elif lines[i].strip().upper() == "Q_DATA":
                    i += 1
                    while not lines[i].strip() == "" and lines[i].strip()[0].isdigit():
                        vals = re.findall('[-+]?\d*\.?\d+', lines[i])
                        if len(vals) < 3:
                            sys.exit("Data for the 'q' vector must contain three comma delimited values: (1) the row index, (2) the parameter index -- with 0 indicating a constant -- and (3) the coeficient. Please reformat your data and retry. Exiting!")
                        if '.' in vals[2]:
                            vals[2] = ConvertToFraction(vals[2])
                        if vals[1] == '0':
                            gMatrix[int(vals[0]) - 1][2*numVar] += Pari(vals[2])
                        else:
                            gMatrix[int(vals[0]) - 1][2*numVar] += Pari(vals[2])*xVar[int(vals[1]) - 1]
                        i += 1
                elif lines[i].strip().upper() == "PARAM_SPACE":
                    i += 1
                    while not lines[i].strip() == "" and lines[i].strip()[0].isdigit():
                        vals = re.findall('[-+]?\d*\.?\d+', lines[i])
                        if len(vals) < 3:
                            sys.exit("Data for the parameter space constraints must contain three comma delimited values: (1) the row index, (2) the parameter index -- with 0 indicating a constant -- and (3) the coeficient. Please reformat your data and retry. Exiting!")
                        if '.' in vals[2]:
                            vals[2] = ConvertToFraction(vals[2])
                        if int(vals[0]) > len(paramSpace):
                            paramSpace.append([Pari.zero(), Pari.zero()])
                        paramSpace[int(vals[0]) - 1][0] += Pari(vals[2])*xVar[int(vals[1]) - 1]
                        i += 1
                elif lines[i].strip().upper() == "PARAM_SPACE_RHS":
                    j = 0
                    i += 1
                    while not lines[i].strip() == "" and ( (lines[i].strip()[0].isdigit() or lines[i].strip()[0] == ".") or (lines[i].strip()[0] == "-" and (lines[i].strip()[1].isdigit() or lines[i].strip()[1] == "."))):
                        if '.' in lines[i]:
                            lines[i] = ConvertToFraction(lines[i])
                        paramSpace[j][1] += Pari(lines[i])
                        j += 1
                        i += 1
                elif lines[i].strip().upper() == "END":
                    break
                else:
                    sys.exit("Unrecognized symbol " + lines[i] + ", please adjust and retry. Exiting!")
        i += 1
            
    return numVar, numParam, gMatrix, xVar, paramSpace, mIsNumeric
    
    
# Read the input file if it is associated with an mpLP or mpQP instance
#
# Input:    lines   --  remaining lines of the input file
#           start   --  index of first line to read from lines
#           Pari    --  the Pari environment
#           sys     --  the variable containing any information passed at the
#                       command line
#           logging --  the logging environment
#           re      --  the re environment (for parsing real expressions)
#           numVar  --  an integer specifying the number of variables in the 
#                       problem. This becomes the number of rows of the G 
#                       matrix.
#           numParam--  an integer specifying the number of parameters in the 
#                       problem. This value is 2 less than the number of entries
#                       needed in the x vector.
#           gMatrix --  the matrix used to store the problem data
#           xVar    --  the vector used to store the Pari variables that 
#                       represent the problem's parameters
#           paramSpace  --  the matrix used to store the constraints 
#                           representing portion of the parameter space to be 
#                           partitioned.
#           gxInitialized   --  a boolean representing whether or not the G 
#                               matrix and x vector have been initialized.
#           mIsNumeric      --  a boolean representing whether or not the M
#                               matrix contains any parameters.
#
# Outputs:  numVar
#           numParam
#           gMatrix
#           xVar
#           paramSpace
#           mIsNumeric
#           numRow
#           numCol
def ReadQP(lines, start, Pari, sys, logging, re, numVar, numParam, gMatrix, xVar, paramSpace, gxInitialized, mIsNumeric):
    i = start
    numCol = -1
    numRow = -1
    numParam = -1
    while i < len(lines):
        if not lines[i].strip() == "":
            if lines[i].strip().upper() == "NUM_COL":
                numCol = int(lines[i+1].strip())
                i += 1
            elif lines[i].strip().upper() == "NUM_ROW":
                numRow = int(lines[i+1].strip())
                i += 1
            elif lines[i].strip().upper() == "NUM_PARAM":
                numParam = int(lines[i+1].strip())
                i += 1
            else:
                if numCol < 0 or numRow < 0 or numParam < 0:
                    sys.exit("After problem type specification, the data file must next include a specification of the number of columns (via keyword 'num_col'), the number of rows (via keyword 'num_row'), and the number of parameters (via keyword 'num_param')! Please check file format and try again. Exiting!")
                if not gxInitialized:
                    numVar = numCol + numRow
                    gxInitialized = True
                    gMatrix, xVar = InitializeGandX(Pari, numVar, numParam, gMatrix, xVar)
                if lines[i].strip().upper() == "A_DATA":
                    i += 1
                    while not lines[i].strip() == "" and lines[i].strip()[0].isdigit():
#                        print(lines[i])
#                        print(repr(lines[i]))
#                        print(repr(lines[i].strip()))
                        vals = re.findall('[-+]?\d*\.?\d+', lines[i])
                        if len(vals) < 4:
                            sys.exit("Data for the A matrix must contain four comma delimited values: (1) the row index, (2) the column index, (3) the parameter index -- with 0 indicating a constant -- and (4) the coeficient. Please reformat your data and retry. Exiting!")
                        if '.' in vals[3]:
                            vals[3] = ConvertToFraction(vals[3])
                        if vals[2] == '0':
                            gMatrix[int(vals[0]) - 1][numVar + numRow + int(vals[1]) - 1] += Pari(vals[3])
                            gMatrix[numRow + int(vals[1]) - 1][numVar + int(vals[0]) - 1] -= Pari(vals[3])
                        else:
                            gMatrix[int(vals[0]) - 1][numVar + numRow + int(vals[1]) - 1] += Pari(vals[3])*xVar[int(vals[2]) - 1]
                            gMatrix[numRow + int(vals[1]) - 1][numVar + int(vals[0]) - 1] -= Pari(vals[3])*xVar[int(vals[2]) - 1]
                            mIsNumeric = False
                        i += 1
                elif lines[i].strip().upper() == "B_DATA":
                    i += 1
                    while not lines[i].strip() == "" and lines[i].strip()[0].isdigit():
                        vals = re.findall('[-+]?\d*\.?\d+', lines[i])
                        if len(vals) < 3:
                            sys.exit("Data for the 'b' vector must contain three comma delimited values: (1) the row index, (2) the parameter index -- with 0 indicating a constant -- and (3) the coeficient. Please reformat your data and retry. Exiting!")
                        if '.' in vals[2]:
                            vals[2] = ConvertToFraction(vals[2])
                        if vals[1] == '0':
                            gMatrix[int(vals[0]) - 1][2*numVar] += Pari(vals[2])
                        else:
                            gMatrix[int(vals[0]) - 1][2*numVar] += Pari(vals[2])*xVar[int(vals[1]) - 1]
                        i += 1
                elif lines[i].strip().upper() == "C_DATA":
                    i += 1
                    while not lines[i].strip() == "" and lines[i].strip()[0].isdigit():
                        vals = re.findall('[-+]?\d*\.?\d+', lines[i])
                        if len(vals) < 3:
                            sys.exit("Data for the 'c' vector must contain three comma delimited values: (1) the column index, (2) the parameter index -- with 0 indicating a constant -- and (3) the coeficient. Please reformat your data and retry. Exiting!")
                        if '.' in vals[2]:
                            vals[2] = ConvertToFraction(vals[2])
                        if vals[1] == '0':
                            gMatrix[int(vals[0]) + numRow - 1][2*numVar] += Pari(vals[2])
                        else:
                            gMatrix[int(vals[0]) + numRow - 1][2*numVar] += Pari(vals[2])*xVar[int(vals[1]) - 1]
                        i += 1
                elif lines[i].strip().upper() == "Q_DATA":
                    i += 1
                    while not lines[i].strip() == "" and lines[i].strip()[0].isdigit():
                        vals = re.findall('[-+]?\d*\.?\d+', lines[i])
                        if len(vals) < 4:
                            sys.exit("Data for the Q matrix must contain four comma delimited values: (1) the row index, (2) the column index, (3) the parameter index -- with 0 indicating a constant -- and (4) the coeficient. Please reformat your data and retry. Exiting!")
                        if '.' in vals[3]:
                            vals[3] = ConvertToFraction(vals[3])
                        if vals[2] == '0':
                            gMatrix[numRow + int(vals[0]) - 1][numVar + numRow + int(vals[1]) - 1] -= Pari(vals[3])
                        else:
                            gMatrix[numRow + int(vals[0]) - 1][numVar + numRow + int(vals[1]) - 1] -= Pari(vals[3])*xVar[int(vals[2]) - 1]
                            mIsNumeric = False
                        i += 1
                elif lines[i].strip().upper() == "PARAM_SPACE":
                    i += 1
                    while not lines[i].strip() == "" and lines[i].strip()[0].isdigit():
                        vals = re.findall('[-+]?\d*\.?\d+', lines[i])
                        if len(vals) < 3:
                            sys.exit("Data for the parameter space constraints must contain three comma delimited values: (1) the row index, (2) the parameter index -- with 0 indicating a constant -- and (3) the coeficient. Please reformat your data and retry. Exiting!")
                        if '.' in vals[2]:
                            vals[2] = ConvertToFraction(vals[2])
                        if int(vals[0]) > len(paramSpace):
                            paramSpace.append([Pari.zero(), Pari.zero()])
                        paramSpace[int(vals[0]) - 1][0] += Pari(vals[2])*xVar[int(vals[1]) - 1]
                        i += 1
                elif lines[i].strip().upper() == "PARAM_SPACE_RHS":
                    j = 0
                    i += 1
#                    print(lines[i])
                    while not lines[i].strip() == "" and ( (lines[i].strip()[0].isdigit() or lines[i].strip()[0] == ".") or (lines[i].strip()[0] == "-" and (lines[i].strip()[1].isdigit() or lines[i].strip()[1] == "."))):
                        if '.' in lines[i]:
                            lines[i] = ConvertToFraction(lines[i])
                        paramSpace[j][1] += Pari(lines[i])
                        j += 1
                        i += 1
                elif lines[i].strip().upper() == "END":
                    break
                else:
                    sys.exit("Unrecognized symbol " + repr(lines[i]) + ", please adjust and retry. Exiting!")
        i += 1
            
    return numVar, numParam, gMatrix, xVar, paramSpace, mIsNumeric, numRow, numCol
