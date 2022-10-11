###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        Perform matrix manipulations -- primarily principal pivots.
#                   Original intention is for use as part of a solver for multi-
#                   parametric Linear Complementarity Problems (mpLCP's).
#
################################################################################


# Define Functions

# Perform a principal pivot on the given matrix, i.e., perform elementary row 
# operations so that column j of M is an identity vector with a 1 in the i-th 
# position.
#
# Input:    M   --  the matrix to manipulate
#           i   --  the row index
#           j   --  the column index
#
# Output:   M   --  the updated matrix
def matrixPivot(M, i, j):
    temp = M[i][j]
    
    for k in range(len(M[i])):
        if k == j:
            M[i][k] = 1
        else:
            M[i][k] = M[i][k]/temp
    for row in range(len(M)):
        if row != i:
            temp2 = M[row][j]
            if temp2 != 0:
                for col in range(len(M[row])):
                    M[row][col] -= temp2*M[i][col]
        
    return(M)
    


# Perform an exchange pivot on the given matrix, i.e., perform the necessary
# principal pivots in order to simultaneously swap two basis elements with their
# complements
#
# Input:    M       --  the matrix to manipulate
#           i       --  the row index of the first basis element
#           iComp   --  the column index of i's complement
#           j       --  the row index of the second basis element
#           jComp   --  the column index of j's complement
#
# Output:   M   --  the updated matrix
def ExchangePivot(M, i, iComp, j, jComp, parallel):
    M = matrixPivot(M, i, jComp, parallel)
    M = matrixPivot(M, j, iComp, parallel)
    M[i], M[j] = M[j], M[i]
        
    return(M)



