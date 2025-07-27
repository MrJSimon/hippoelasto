##############################################################################
##
## Author:      Jamie E. Simon
##
## Description: 
##              
##############################################################################


## Import packages
import sympy as sp
from ..StretchDescription.stretches_invariants import Invariants

def FirstPiolaKirschoffStress(Wi,
                              I1b, I2b, Jac,
                              L11, L22, L33):
    """ This module computes the deviatoric and volumetric
        stress using the first Piola-Kirschoff definition.
        Under the assumption of uniaxial tension.
        
    input
    ---------
    Wi:  sympy, strain energy density [J]
    
    I1b: sympy, 1. modified invariant [-]
    
    I2b: sympy, 2. modified invariant [-]
    
    Jac: sympy, jacobian determinant of the deformation gradient tensor [-]
    
    L11: sympy, 1. stretch [-]
    
    L22: sympy, 2. stretch [-]
    
    L33: sympy, 3. stretch [-]
        
    output
    ---------
    P22: sympy, first Piola Kirshcoff stress [MPa]
    
    """
    
    ## Derivative of W with respect to the modified invariants
    dWdI1b = sp.diff(Wi, I1b) 
    dWdI2b = sp.diff(Wi, I2b)    
    
    ## Derivative of W with respect to the jacobian J
    dWdJ = sp.diff(Wi, Jac)
    
    ## Get the invariants
    I1b_expr, I2b_expr, Jac_expr = Invariants(L11, L22, L33)
    
    ## Derivative of the modified invariants with respect to the stretches
    dI1bdL22 = sp.diff(I1b_expr, L22)
    dI2bdL22 = sp.diff(I2b_expr, L22)
    
    ## Derivative of the jacobian with respect to the stretches
    dJdL22 = sp.diff(Jac_expr, L22)
    
    ## Compute the deviatoric stress using the chain rule
    P22_dev = dWdI1b * dI1bdL22 + dWdI2b * dI2bdL22
    
    ## Compute the volumetric stress
    P22_vol = dWdJ * dJdL22
    
    ## Compute the total stress
    P22 = P22_dev + P22_vol
    
    ## Substitute the definition of I1bar, I1bar and Jac into the total stress
    P22 = P22.subs({I1b: I1b_expr, I2b:I2b_expr, Jac:Jac_expr})   
    
    return P22


def FirstPiolaKirschoffStressIncompressible(Wi,
                                            I1b, I2b, Jac,
                                            L11, L22, L33):
    """ This module computes the deviatoric and volumetric
        stress using the first Piola-Kirschoff definition.
        Under the assumption of uniaxial tension.
        
    input
    ---------
    Wi:  sympy, strain energy density [J]
    
    I1b: sympy, 1. modified invariant [-]
    
    I2b: sympy, 2. modified invariant [-]
    
    Jac: sympy, jacobian determinant of the deformation gradient tensor [-]
    
    L11: sympy, 1. stretch [-]
    
    L22: sympy, 2. stretch [-]
    
    L33: sympy, 3. stretch [-]
        
    output
    ---------
    P22: sympy, first Piola Kirshcoff stress [MPa]
    
    """
    
    ## Derivative of W with respect to the modified invariants
    dWdI1b = sp.diff(Wi, I1b) 
    dWdI2b = sp.diff(Wi, I2b)    
    
    ## Get the invariants
    I1b_expr, I2b_expr, _ = Invariants(L11, L22, L33)
    
    ## Derivative of the modified invariants with respect to the stretches
    dI1bdL22 = sp.diff(I1b_expr, L22)
    dI2bdL22 = sp.diff(I2b_expr, L22)
    dI1bdL33 = sp.diff(I1b_expr, L33)
    dI2bdL33 = sp.diff(I2b_expr, L33)
    
    ## Lagrangian multiplier/external pressure 2D-plane-strain definition
    ## derived from P33_vol = 0.0
    p_lagrange  = dWdI1b * dI1bdL33 + dWdI2b * dI2bdL33
     
    ## Compute the deviatoric stress using the chain rule
    P22_dev = dWdI1b * dI1bdL22 + dWdI2b * dI2bdL22
    
    ## Compute the total stress
    P22 = P22_dev - p_lagrange
    
    ## Substitute the definition of I1bar, I1bar
    P22 = P22.subs({I1b: I1b_expr, I2b:I2b_expr})   
    
    return P22
