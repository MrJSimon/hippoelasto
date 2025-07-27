##############################################################################
##
## Author:      Jamie E. Simon
##
## Description: 
##              
##############################################################################

## Import packages
from ..StretchDescription.stretches_invariants import Invariants

def EnergyInvariantModified(Wi,
                            I1b, I2b, Jac,
                            L11, L22, L33):
    """ This module substitutes the stress description of the modified
        invariants into the defined energy formulation
        
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
    Wi_modified: sympy, first Piola Kirshcoff stress [MPa]
    
    """
    
    ## Get the invariants
    I1b_expr, I2b_expr, Jac_expr = Invariants(L11, L22, L33)
    
    ## Substitute the definition of I1bar, I1bar and Jac
    Wi = Wi.subs({I1b: I1b_expr, I2b:I2b_expr, Jac:Jac_expr})
    
    return Wi