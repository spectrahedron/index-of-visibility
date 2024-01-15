# index-of-visibility
Numerical experiments approximating the index of visibility

Run **LP-approximation.py** to reproduce the calculations. 

- The optimisation part of the code relies on the gurobi package.
- Scroll down past the function definitions to change parameters.
Set dimension d to 2 or 3.
- Choose the number of steps: this is the number of values of lambda for which 
the approximate solution is calculated.
- Choose the number of **subdivision**s. This is the size of the grid that
approximates the measure. Note that subdivisions = 1000 results in a grid of 
the size 1000x1000, which is the number of variables in the linear program.
- The code outputs the results in the current directory. To change the directory,
modify the **path** variable.
