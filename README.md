
# upLCPsolver


### Introduction

upLCPsolver is a software package designed to solve the uni-parametric Linear Complementarity Problem (upLCP):

$$
\begin{align*}
w - M(\theta)z &= q(\theta)\\\\[1mm]
w^\top z &= 0\\\\[1mm]
w,z &\geq 0
\end{align*}
$$

upLCPsolver offers a simplification of the methodology developed for multi-parametric Linear Complementarity Problems (mpLCP) by Nathan Adelgren.[^fn1] Specifics of the methodology employed within upLCPsolver are outlined in a paper submitted for publication in March 2022. This document will be updated with a reference if/when that work is published.

Additonally, as both linear programs (LPs) and quadratic programs (QPs) can be formulated as LCPs, upLCPsolver is capable of solving uni-parametric LPs (upLPs) having the form 

$$
\begin{align*}
\min & c(\theta)^\top x\\\\[1mm]
\text{s.t.} & A(\theta) x \leq b(\theta)\\\\[1mm]
& x \geq 0
\end{align*}
$$

or uni-parametric QPs (upQPs) having the form

$$
\begin{array}{lr}
\min & \frac{1}{2} x^\top Q(\theta) x + c(\theta)^\top x\\\\[1mm]
\text{s.t.} & A(\theta) x \leq b(\theta)\\\\[1mm]
& x \geq 0.
\end{array}
$$

Note that correct functionality of this version of upLCP relies on the following assumptions:

> *Assumption 1* -- The parameter $\theta$ lies within an interval $\Theta := [\alpha, \beta] \subset \mathbb{R}$ and, moreover, the associated upLCP, upLP, and/or upQP is feasible for every $\theta \in \Theta$.

> *Assumption 2* -- When the given uni-parametric problem is formulated as an instance of upLCP, the matrix $M(\theta)$ is sufficient for all $\theta \in \Theta$ (with $I$ defined as above). 

Refer to Cottle et al.[^fn2] for the definition of sufficiency.

### Background

upLCPsolver is written in Python and is released as open source code under the (enter license information here).
The code has been written by Nathan Adelgren.

upLCPsolver has only been tested on Linux systems, but should be compatible with other operating systems as long as all dependencies listed below can be met.

### Dependencies

upLCPsolver depends on:

- Python 3 -- Download from python.org or install with your favorite package manager
- [PARI](https://pari.math.u-bordeaux.fr/) and [CyPari2](https://cypari2.readthedocs.io/en/latest/) -- Can be installed using apt (or similar) and pip, respectively. **Note**, however, that testing of upLCPsolver with PARI version 2.11 (the version available via the apt repository at the time of this writing) *was not successful*. Successful testing was conducted using PARI version 2.14, compiled from source. Instructions for compiling PARI from source can be found in Section 3 of [this document](https://pari.math.u-bordeaux.fr/PDF/PARIwithWindows.pdf).

Additionally, the following Python libraries are employed by upLCPsolver:

- sys
- re
- logging
- time
- multiprocessing
- random
- os
- collections
- math
- copy

### Using upLCPsolver

upLCPsolver is called from the command line as follows:

    > python3 upLCP_solver.py /path/to/data/file (options)
  
#### Data File

The data file may have any extension, but must be a text file in one of the following formats:

##### upLP Format

This format is used to specify a uni-parametric linear program (upLP). The file must containing the following (in the specified order):

- lp -- the file should literally begin with the string "lp" to indicate that the problem to be specified is a upLP.
- num_row -- an integer -- the number of rows in the constraint matrix $A(\theta)$.
- num_col -- an integer -- the number of columns in the constraint matric $A(\theta)$.
- num_param -- an integer -- the number of parameters present in the instance of upLP. (Value must be 1)
- A_data -- a matrix -- describes the nonzero contents of $A(\theta)$ in the following format:
  - Each row of A_data must consist of four entries (comma delimited): 
    1. row index
    2. column index
    3. parameter index (0 indicates the constant term)
    4. coefficient
    
  See below for an example.
- c_data -- a matrix   -- describes the nonzero contents of $c(\theta)$ in the same format used above for A_data.
- b_data -- a matrix   -- describes the nonzero contents of $b(\theta)$ in the same format used above for A_data.
- Param_Space -- a matrix   -- it is assumed in this implementation that the set $\Theta$ of attainable parameter values can be represented as a system of linear inequalities in the form $Hv \leq r$. Hence, "Param_Space" describes the nonzero entries of matrix $H$ in the following format:
  - Each row of Param_Space consists of three entries (comma delimited):
    1. row index
    2. column index
    3. coefficient
- Param_Space_RHS -- a vector   -- provides the elements of the right-hand-side vector $r$, as described above. Note that Param_Space_RHS should be given in column format and even zero elements must be included.
- END -- specifies the end of the file.

As an example, consider the following instance of upLP:
\[
\begin{array}{lc}
\min & \left[
  \begin{array}{rrrr}
  1 & 1 & 1 & 1\\
  \end{array}
\right]^\top x\\[1mm]
\text{s.t} & \left[
  \begin{array}{rrrr}
  -2 & -1         & -6 & 1\\
  -2 & 3          & -1 & \theta - 2\\
  3  & \theta - 4 & 5  & -1\\
  \end{array}
\right]x \leq \left[
  \begin{array}{r}
  -2\\
  7\\
  -5
  \end{array}
\right]\\[1mm]
& x \geq 0
\end{array}
\]
with $\theta \in \Theta = [-2, 2]$.

The data file for this instance would be as follows:
        
        lp
        
        num_row
        3
         
        num_col
        4
        
        num_param
        1
         
        A_data 
         1,1,0,-2
         1,2,0,-1
         1,3,0,-6
         1,4,0,1
         2,1,0,-2
         2,2,0,3
         2,3,0,-1
         2,4,0,-2
         2,4,1,1
         3,1,0,3
         3,2,0,-4
         3,2,1,1
         3,3,0,5
         3,4,0,-1
         
        c_data
         1,0,1
         2,0,1
         3,0,1
         4,0,1
                 
        b_data 
         1,0,-2
         2,0,7
         3,0,-5
              
        Param_Space 
         1,1,-1
         2,1,1
               
        Param_Space_RHS 
         2
         2
        
        END  

##### upQP Format

This format is used to specify a uni-parametric quadratic program (upQP). The file must containing the following (in the specified order):

- qp -- the file should literally begin with the string "qp" to indicate that the problem to be specified is a upQP.
- num_row -- an integer -- the number of rows in the constraint matrix $A(\theta)$.
- num_col -- an integer -- the number of columns in the constraint matric $A(\theta)$.
- num_param -- an integer -- the number of parameters present in the instance of upLP. (Value must be 1)
- A_data -- a matrix -- describes the nonzero contents of $A(\theta)$ in the following format:
  - Each row of A_data must consist of four entries (comma delimited): 
    1. row index
    2. column index
    3. parameter index (0 indicates the constant term)
    4. coefficient
    
  See below for an example.
- Q_data -- a matrix -- describes the nonzero contents of $Q(\theta)$ in the same format used above for A_data.
- c_data -- a matrix   -- describes the nonzero contents of $c(\theta)$ in the same format used above for A_data.
- b_data -- a matrix   -- describes the nonzero contents of $b(\theta)$ in the same format used above for A_data.
- Param_Space -- a matrix   -- it is assumed in this implementation that the set $\Theta$ of attainable parameter values can be represented as a system of linear inequalities in the form $Hv \leq r$. Hence, "Param_Space" describes the nonzero entries of matrix $H$ in the following format:
  - Each row of Param_Space consists of three entries (comma delimited):
    1. row index
    2. column index
    3. coefficient
- Param_Space_RHS -- a vector   -- provides the elements of the right-hand-side vector $r$, as described above. Note that Param_Space_RHS should be given in column format and even zero elements must be included.
- END -- specifies the end of the file.

As an example, consider the following instance of upQP:
\[
\begin{array}{lc}
\min & \frac{1}{2}x^\top\left[
  \begin{array}{rrrr}
  -9\theta + 22 & -11\theta+6  & -24\theta+16  & -25\theta+18\\
  -11\theta+6   & -14\theta+20 & 4\theta - 2   & -6\theta+15\\
  -24\theta+16  & 4\theta - 2  & -8\theta + 18 & -5\theta+10\\
  -25\theta+18  & -6\theta+15  & -5\theta+10   & -3\theta+21
  \end{array}
\right]x + \left[
  \begin{array}{rrrr}
  1 & 1 & 1 & 1\\
  \end{array}
\right]^\top x\\[1mm]
\text{s.t} & \left[
  \begin{array}{rrrr}
  -2 & -1 & -6 & 1\\
  -2 & 3  & -1 & -2\\
  3  & -4 & 5  & -1\\
  \end{array}
\right]x \leq \left[
  \begin{array}{r}
  -2\\
  7\\
  -5
  \end{array}
\right]\\[1mm]
& x \geq 0
\end{array}
\]
with $\theta \in \Theta = [0, 1]$.

The data file for this instance would be as follows:
        
        qp
        
        num_row
        3
         
        num_col
        4
        
        num_param
        1
         
        A_data 
         1,1,0,-2
         1,2,0,-1
         1,3,0,-6
         1,4,0,1
         2,1,0,-2
         2,2,0,3
         2,3,0,-1
         2,4,0,-2
         3,1,0,3
         3,2,0,-4
         3,3,0,5
         3,4,0,-1
         
         Q_data
         1,1,0,22
         1,1,1,-9
         1,2,0,6
         1,2,1,-11
         1,3,0,16
         1,3,1,-24
         1,4,0,18
         1,4,1,-25
         2,1,0,6
         2,1,1,-11
         2,2,0,20
         2,2,1,-14
         2,3,0,-2
         2,3,1,4
         2,4,0,15
         2,4,1,-6
         3,1,0,16
         3,1,1,-24
         3,2,0,-2
         3,2,1,4
         3,3,0,18
         3,3,1,-8
         3,4,0,10
         3,4,1,-5
         4,1,0,18
         4,1,1,-25
         4,2,0,15
         4,2,1,-6
         4,3,0,10
         4,3,1,-5
         4,4,0,21
         4,4,1,-3
         
         c_data
         1,0,1
         2,0,1
         3,0,1
         4,0,1
                 
        b_data 
         1,0,-2
         2,0,7
         3,0,-5
              
        Param_Space 
         1,1,-1
         2,1,1
               
        Param_Space_RHS 
         0
         1
        
        END
        
        
##### upLCP Format

- lcp -- the file should literally begin with the string "lcp" to indicate that the problem to be specified is a upLCP.
- h -- an integer -- the dimension of upLCP decision variable vectors $w$ and $z$.
- k -- an integer -- the number of parameters present in the instance of upLCP. (Value must be 1)
- M_data -- a matrix   -- describes the nonzero contents of $M(\theta)$ in the following format:
  - Each row of M_data must consist of four entries (comma delimited): 
    1. row index
    2. column index
    3. parameter index (0 indicates the constant term)
    4. coefficient
    
  See below for an example.
- q_data -- a matrix   -- describes the nonzero contents of $q(\theta)$ in the same format used above for M_data.
- Param_Space -- a matrix   -- it is assumed in this implementation that the set $\Theta$ of attainable parameter values can be represented as a system of linear inequalities in the form $Hv \leq r$. Hence, "Param_Space" describes the nonzero entries of matrix $H$ in the following format:
  - Each row of Param_Space consists of three entries (comma delimited):
    1. row index
    2. column index
    3. coefficient
- Param_Space_RHS -- a vector   -- provides the elements of the right-hand-side vector $r$, as described above. Note that Param_Space_RHS should be given in column format and even zero elements must be included.
- END -- specifies the end of the file.
                       
As an example, consider the following instance of upLCP:
\[
\begin{array}{c}
w - \left[
  \begin{array}{rrrr}
  0 & 0             & -2  & -1\\
  0 & 0             & -5  & \theta + 7\\
  1 & 3             & 0   & 0\\
  1 & -\theta - 5 &     & 0\\
  \end{array}
\right]z = \left[
  \begin{array}{r}
  2\\
  \theta + 2\\
  20\\
  10\\
  \end{array}
\right]\\[1mm]
w^\top z = 0\\[1mm]
w,z \geq 0
\end{array}
\]
with $\theta \in \Theta = [-3, 1]$.

The data file for this instance would be as follows:
        
        lcp
        
        h 
        4
        
        k
        1
        
        M_data
        1,3,0,-2
        1,4,0,-1
        2,3,0,-5
        2,4,0,7
        2,4,1,1
        3,1,0,1
        3,2,0,3
        4,1,0,1
        4,2,0,-5
        4,2,1,-1
        
        q_data
        1,0,2
        2,0,2
        2,1,1
        3,0,20
        4,0,10
        
        Param_Space
        1,1,-1
        2,1,1
        
        Param_Space_RHS 
        3
        1
        
        END  
  
#### Options

Options can be passed to upLCPsolver as command line flags in the form "-Flag Value". Most options are used to set the value of an individual parameter. Available options are:

- -numThreads -- A positive integer, less than or equal to the number of available threads, specifying the number of threads to use when partitioning the parameter space. (Default: the number of available threads)
- -parStart -- A boolean indicating whether or not the initial search interval should be divided into "numThreads" subintervals and have the partitioning search begin by processing each subinterval in parallel. We note that during testing we observed that finding initial bases was quite time consuming, and thus, using all available threads to find initial bases caused poor performance. (Default: False)
- -showProgress -- A boolean indicating whether or not information about the intervals being processed should be displayed throughout execution. (Default: True)


**Note**: At the command line, appropriate values for booleans are assumed to be "T" and "F".


#### Full Example of Calling upLCPsolver from the Command Line:

    > python3 upLCP_solver.py /path/to/data/file -numThreads 4 -parStart F -showProgress T

### Licensing

upLCPsolver is free software. You are welcome to redistribute it and/or modify it under the
terms of the GNU General Public License version 3 (or later) as published by the Free Software
Foundation. upLCPsolver is distributed as a resource for the research community, but HAS NO WARRANTY WHATSOEVER.

Please refer to the License for details. It should be contained in the file 'LICENSE'. If not, it can be obtained [here](https://www.gnu.org/licenses/gpl-3.0.en.html).


[^fn1]: Adelgren N. Advancing Parametric Optimization: On Multiparametric Linear Complementarity Problems with Parameters in General Locations. Springer, 2021. ([DOI](https://doi.org/10.1007/978-3-030-61821-6))

[^fn2]: Richard W Cottle, Jong-Shi Pang, and Richard E Stone. The Linear Complementarity
Problem. SIAM, 2009. ([DOI](https://doi.org/10.1137/1.9780898719000))
