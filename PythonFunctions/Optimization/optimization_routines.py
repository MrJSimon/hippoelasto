# -*- coding: utf-8 -*-
"""
Created on Sun Jul 27 19:13:42 2025

@author: jeg_e
"""
##############################################################################
##
## Author:      Jamie E. Simon
##
## Description: 
##              
##############################################################################

## Import modulues
import numpy as np
from scipy.optimize import minimize

def PredictionStatementTension(params,StressFunction,Xi,nu = 0.50):
    ## Compute the stretches in the x- and z-directions, considering incompressibility
    lam2 = 1.0 + Xi
    lam1 = lam2**(-nu)
    lam3 = lam2**(-nu)
    ## Collect all inputs: stretches and params
    input_args = [lam1,lam2,lam3] + list(params)
    ## Compute prediction statement based on input
    Ypred = StressFunction(*input_args)
    return Ypred
    
def ObjectiveFunctionSSD(params,PredictionStatement_i,StressFunction_i,Xi,Yi):
    ## Compute prediction statement
    Yj = PredictionStatement_i(params,StressFunction_i,Xi)
    ## Compute sum of squared differences
    SSD = (1/len(Yi))*np.sum((Yj-Yi)**2)
    return SSD

def EnergyConstraintTension(params,EnergyFunction,Xi,nu = 0.50):
    ## Compute the stretches in the x- and z-directions, considering incompressibility
    lam2 = 1.0 + Xi
    lam1 = lam2**(-nu)
    lam3 = lam2**(-nu)
    ## Collect all inputs: stretches and params
    input_args = [lam1,lam2,lam3] + list(params)
    ## Compute the energy function
    Wvals = EnergyFunction(*input_args)   
    return Wvals


def OptimizationSLSQP(ObjectiveFunction, coefs, args, constraints = False,
                      method = 'SLSQP', 
                      options = {'ftol': 10e-30, 'disp': True, 'maxiter': 3000}):

    ## Call minimization/optimization 
    solution = minimize(ObjectiveFunction, coefs, args=args, 
                        constraints=constraints, 
                        method='SLSQP', 
                        options=options)
    
    ## Return fitting parameters
    return solution.x
