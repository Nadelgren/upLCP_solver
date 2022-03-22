###############################################################################
#
#   Author:         Nathan Adelgren
#   Affiliation:    Andlinger Center For Energy and the Environment
#                   Princeton University
#
#   Purpose:        Read in and set parameters passed from the command line in 
#                   order to solve an instance of the multiparametric Linear 
#                   Complementarity Problem (mpLCP).
#
################################################################################


# Define Functions

# Print a warning message if a commandline flag is passed with an unrecognized
# value
#
# Input:    flag -- the command line flag
#           message --  the warning message to print
#           logging --  the logging environment
def PrintInvalidParameterMessage(flag, curValue, validValues, logging):
    logging.warning("Invalid value for flag '" + str(flag) + "', ignoring and leaving at default value of " + str(curValue) + ". For future reference, valid values are " + str(validValues) + ".")

# Parse the input flags
#
# Input:    sys --  the variable containing any information passed at the
#                   command line
#           logging --  the logging environment
#           numThreads      --  the number of threads with which the program 
#                               should be run
#           parallelStart   --  a boolean indicating whether or not the original
#                               search region should be split into 'numThread'
#                               subregions at the start, and each immediately
#                               passed to the processing queue
#           showProgress    --  a boolean indicating whether or not information
#                               about the intervals being processed should be 
#                               displayed throughout execution
#
# Outputs:  numThreads
#           parallelStart
def ReadFlags(  sys, 
                logging, 
                numThreads,
                parallelStart,
                showProgress):
    # Read the flags
    
    i = 2
    while i < len(sys.argv):
        if sys.argv[i][0] != '-':
            print("Invalid Command Line Argument (missing '-'). Exiting.\n")
            sys.exit("Got " + str(sys.argv[i]))
        else:
            if sys.argv[i] == "-parStart":
                i += 1
                if sys.argv[i].upper() == "T":
                    parallelStart = True
                elif sys.argv[i].upper() == "F":
                    parallelStart = False
                else:
                    PrintInvalidParameterMessage("-parStart", parallelPivot, "T and F", logging);
            elif sys.argv[i] == "-showProgress":
                i += 1
                if sys.argv[i].upper() == "T":
                    showProgress = True
                elif sys.argv[i].upper() == "F":
                    showProgress = False
                else:
                    PrintInvalidParameterMessage("-showProgress", showProgress, "T and F", logging);
            elif sys.argv[i] == "-numThreads":
                i += 1
                try:
                    val = int(sys.argv[i])
                    if val > 0:
                        numThreads = val
                    else:
                        PrintInvalidParameterMessage("-numThreads", numThreads, "positive integers", logging);
                except ValueError:
                    PrintInvalidParameterMessage("-numThreads", numThreads, "positive integers", logging);
            else:
                print("Invalid Command Line Argument. Exiting.\n")
                sys.exit("Got " + str(sys.argv[i]))
        i += 1
    
    return numThreads, parallelStart, showProgress
