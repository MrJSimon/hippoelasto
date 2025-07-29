# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 21:48:01 2025

@author: jeg_e
"""
##############################################################################
##
## Author:      Jamie E. Simon
##
## Description: 
##              
##############################################################################

## Import packages
import numpy as np


def NumericalElasticModulus(Xi,PredictionFunction,StressFunction, params, nu):
    """
    Computes the numerical tangent modulus (derivative of stress with respect to strain)
    from a predicted stress-strain response using logarithmically spaced strain values.
    
    Parameters
    ----------
    Xi : array_like
        Original strain data (1D array). Used to define bounds for log-spaced sampling.
    
    PredictionFunction : callable
        A function that predicts the stress given material parameters, a stress function,
        strain values, and Poisson's ratio. Should have the signature:
        `PredictionFunction(params, StressFunction, strain_array, nu) -> stress_array`.
    
    StressFunction : callable
        A material model or constitutive function used inside `PredictionFunction`.
    
    params : array_like
        Parameters to be passed into the `PredictionFunction` (e.g., material constants).
    
    nu : float
        Poisson’s ratio used in the prediction model.
    
    Returns
    -------
    tangent_modulus : ndarray
        Array of computed tangent modulus values (dσ/dε) over the artificial strain range.
    
    Notes
    -----
    - Strain values (`Xi_p`) are generated using `np.logspace` for higher resolution
      near the origin, improving numerical accuracy in the elastic region.
    - A small offset is applied to `Xi_min` to avoid issues with `log10(0)`.
    - Useful for estimating the elastic modulus from the initial slope of the stress-strain curve.
    """
    
    # Add a small positive offset to avoid log10(0)
    Xi_min = np.min(Xi[np.nonzero(Xi)])  # get smallest non-zero strain
    Xi_max = np.max(Xi)

    # Define log bounds safely
    Xlog1 = np.log10(Xi_min * 0.1)
    Xlog2 = np.log10(Xi_max * 1.1)

    # Make artifical X-range
    Xi_p = np.logspace(Xlog1,Xlog2/10,num=500,endpoint = True)
    
    # Generate prediction statement
    Yi_p = PredictionFunction(params, StressFunction, Xi_p, nu)

    # Compute tangent modulus
    tangent_modulus = np.gradient(Yi_p, Xi_p)
    
    return tangent_modulus[0]