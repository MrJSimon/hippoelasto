##############################################################################
##
## Author:      Jamie E. Simon
##
## Description: 
##              
##############################################################################

## -------- ABAQUS INPUT FORMATTING --------------- ##
from sympy.printing.fortran import fcode
from sympy import Integer, Float

def fortran_d0_lines(names, exprs, value_range=range(-100, 101), indent=6):
    """
    Convert symbolic expressions to Fortran assignment lines with .D0 double precision literals.
    
    Parameters:
        names: list of variable names (strings)
        exprs: list of sympy expressions
        value_range: range of integers to promote (default -100 to 100)
        indent: number of spaces to prefix each line
        
    Returns:
        A multiline Fortran code string
    """
    # Promote integers to double precision floats
    subs_map = {Integer(n): Float(n, 64) for n in value_range}
    promoted_exprs = [expr.xreplace(subs_map) for expr in exprs]
    
    # Generate clean Fortran code lines
    d_line = "\n".join(
        f"		 {name:<5} = {fcode(expr, assign_to=None, source_format='free', standard=95)}"
        for name, expr in zip(names, promoted_exprs)
    )

    return d_line