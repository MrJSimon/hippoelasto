##############################################################################
##
## Author:      Jamie E. Simon
##
## Description: 
##              
##############################################################################


def Invariants(L11,L22,L33):
    """ This module computes the first and second invariants and modified
        invariants and the jacobian of the deformation gradient using
        the stretches.
        
    input
    ---------
    L11:  sympy, 1. stretch [-]
    
    L22:  sympy, 2. stretch [-]
    
    L33:  sympy, 3. stretch [-]
        
    output
    ---------
    I1b_expr:  sympy, 1. modified invariant [-]
    
    I2b_expr:  sympy, 2. modified invariant [-]
    
    jac_expr:  sympy, jacobian determinant of the deformation gradient tensor [-]
    
    """
    
    ## Compute the determinant J (volume change)
    Jac_expr = L11 * L22 * L33

    ## Compute invariants
    I1_expr = L11**2 + L22**2 + L33**2
    I2_expr = L11**2 * L22**2 + L22**2 * L33**2 + L33**2 * L11**2
    
    ## Compute modified invariants
    I1b_expr = Jac_expr**(-2/3) * I1_expr
    I2b_expr = Jac_expr**(-4/3) * I2_expr
    
    return I1b_expr, I2b_expr, Jac_expr