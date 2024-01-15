# -*- coding: utf-8 -*-
"""
15 Jamuary 2024
@author: Vera Roshchina https://github.com/spectrahedron/

Numerical calculation of the optimal measure that minimises the integral of 
1 - κ(φ,ψ).

Note that the code is written with clarity in mind, and where code optimisation
won't bring substantial benefits clarity is chosen over efficiency,
to minimise errors and improve readability.

The optimisation part of the code relies on the gurobi package.

Scroll down past the function definitions to change parameters.

Set dimension d to 2 or 3.

Choose the number of steps: this is the number of values of lambda for which 
the approximate solution is calculated.

Choose the number of subdivisions. This is the size of the grid that
approximates the measure. Note that subdivisions = 1000 results in a grid of 
the size 1000x1000, which is the number of variables in the linear program.

The code outputs the results in the current directory. To change the directory,
modify the path variable.

"""

import math
import gurobipy as gp
from gurobipy import GRB
import datetime
import os


# The function sin(θ+2 η)/cos(η)
def rootFunction(theta,eta):
    return math.sin(theta+ 2*eta)/math.cos(eta) 

# Calculate the root η of the equation Λ=sin(θ+2 η)/cos(η) that
# lies in the interval [(π/2-θ)+, (π-θ)/2]. 
# The function is using bisection to calculate the root.
# This is reliable while not time consuming, as we only need to precalculate 
# a small number of values before running the linear programming routine.
def getRoot(Lambda, theta,accuracy = 0.00000000001):
    # minimal and maximal values on the interval
    left = max([0,math.pi/2 - theta])
    right = (math.pi-theta)/2
    # initial guess
    midpoint = (left+right)/2
    
    delta = (right-left)/2
    while (delta>accuracy):
        leftval = rootFunction(theta, left)-Lambda
        midval = rootFunction(theta, midpoint)-Lambda

        if (leftval*midval<0):
            # the root is on the left
            right = midpoint
        else:
            left = midpoint
        
        midpoint = (right+left)/2
        delta = (right-left)/2
    return midpoint


# the integral coefficient
def kappa(phi,psi,Lambda):
    theta = phi+psi
    if (Lambda>1):
        return 1+ math.cos(theta)
    else:
        if (theta<math.pi-math.asin(Lambda)):
            eta = getRoot(Lambda,theta)
        else:
            eta = 0
        return 1+ math.cos(theta+ 2*eta)+2*Lambda*math.sin(eta) 


# Solves an LP approximation to the problem for a given Lambda.
# Divisions is the number of subdivisions of psi and phi axes.
# divisions = 10 gives 10 subdivisions, so it is a 10x10 grid.
def minIntegral(Lambda,dimension = 2,divisions = 10):
    
    # length of each division segment
    # along each side: pi/2 divided by the number of divisions
    delta = (math.pi/2)/divisions
    
    model = gp.Model('visibility')
    
    # add variables: each variable correspond to a cell, indexed as x[i,j],
    # with i and j from 0 to divisions-1.
    x = model.addVars(divisions,divisions,vtype=GRB.CONTINUOUS,name = "x")
    
    
    # add constraints: the sum of x[k,*] should be equal to the projection of the 
    # measure on the relevant interval. 
    for k in range(divisions):

        angle1 = k*delta
        angle2 = (k+1)*delta
        rhs = (math.sin(angle2))**(dimension-1)-(math.sin(angle1))**(dimension-1)
        model.addConstr(gp.quicksum(x[k,j] for j in range(divisions)) ==rhs)
        model.addConstr(gp.quicksum(x[i,k] for i in range(divisions)) ==rhs)
    
    
    indices = [(i,j) for j in range(divisions) for i in range(divisions)]

    for (i,j) in indices:
        model.addConstr(x[i,j]>=0)

    
    objective = [[kappa(math.pi*(i+0.5)/(2 *divisions),math.pi*(j+0.5)/(2*divisions),Lambda) for j in range(divisions)] for i in range(divisions)]
    model.setObjective(gp.quicksum(x[i,j]*objective[i][j] for (i, j) in indices),GRB.MINIMIZE)
    
    model.optimize()
    
    obj = model.getObjective().getValue()
   
    solution = list()
    for i in range(divisions):
        solution.append(list())
        for j in range(divisions):
            solution[i].append(x[i,j].X)

    return obj,solution



# output files
now = datetime.datetime.now()
path = now.strftime("%Y-%m-%d-%H-%M-%S")
os.mkdir(path)
path=path+"/"
f = open(path+"lambda-integral-values.txt", "a")


# number of values of lambda for which the problem is solved
steps = 100
# side subdivisions for the linear programming problem, for both angles
subdivisions = 1000
# dimension of the problem (note code only works for d=2 and d=3)
d = 2

if (d==2):
    maxlambda = 3*math.pi/8 
else:
    maxlambda = maxlambda = 4.0/3
    
step = maxlambda/steps

delta = (math.pi/2)/subdivisions


# generate the list of λ's from step to maxlambda
lambdas = [step*(m+1) for m in range(steps)]


it = 0

for l in lambdas:
     
    it = it+1
    
    print('Step '+str(it)+' out of '+str(len(lambdas))+', lambda='+"{:.12f}".format(l)+'.')
    
    if (d==2):
        Lambda = (8*l)/(3*math.pi)
    else:
        Lambda = 3*l/4.0
    
    val,solution = minIntegral(Lambda,d,subdivisions)
    
    # write the value of lambda and the value of integral in the format useful
    # for Mathematica plotting
    f.write('{'+"{:.12f}".format(l)+', '+"{:.12f}".format(val)+'},\n')
    
    g = open(path+"measure-"+str(it)+'-lambda-'+"{:.12f}".format(l).replace('.','-')+".txt","a")
    g.write("diimension="+str(d)+'\n')
    g.write("lambda="+"{:.12f}".format(l)+'\n')
    g.write("integral="+"{:.12f}".format(val)+'\n')
            
    # only write the values that are nonzero, otherwise it is impossible 
    # to plot with Mathematica
    for i in range(len(solution)):
        line = solution[i]
        for j in range(len(line)):
            zval = line[j]
            if (zval>0.00000001):
                x = "{:.12f}".format(i*delta+delta/2)
                y = "{:.12f}".format(j*delta+delta/2)
                z = "{:.12f}".format(line[j])
                g.write('{'+x+','+y+','+z+'},')            
    
    g.close()

f.close()