##############################################################################
##
## Author:      Jamie E. Simon
##
## Description:
##
##############################################################################

## Load in modules
import numpy as np
import sympy as sp
from pathlib import Path
from sympy import lambdify
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from PythonFunctions.StressDescription.piola_kirschoff_stress import FirstPiolaKirschoffStress
from PythonFunctions.EnergyDescription.energy_substitution import EnergyInvariantModified
from PythonFunctions.Abaqus.generate_vumat import GenerateVumatHyperelasticity

# -------------------- LOAD IN DATA -------------------- #

## Specimen tollerances
L0, T0, W0 = 60.0, 0.75, 5.0 # Specimen tollerances

## Create temperature list
temp_list = [-30,-20,-10,0,10,20]

## Create strain rate array
eps_rate = np.array([10,100,1000,5000])*(1.0/60.0)*(1/L0)

## Set epsilon thresshold
eps_thress = 0.1

## create repeated array for output purposes.
eps_rate_repeated  = np.tile(eps_rate,len(temp_list))
temp_list_repeated = np.repeat(temp_list,len(eps_rate))

## Set fontsize of labels - plotter
fontsize_plot = 13
markersize_plot = 12

## Load in variable names
variable_names = np.loadtxt('var_names.txt',dtype=str,delimiter=None)

## Load in stress strain array
nominal_strain_array = np.loadtxt('nominal_strain.txt')
nominal_stress_array = np.loadtxt('nominal_stress.txt')

## Run through all variables in variable names
for i in range(0,1): #len(variable_names)
    ## Initiate stress and strain
    sig_n = nominal_stress_array[:,i]
    eps_n = nominal_strain_array[:,i]

    ## Remove values equaivalent to zero
    eps_n = eps_n[sig_n!=0]
    sig_n = sig_n[sig_n!=0]
    
    ## Zero out stress and strain
    sig_n = sig_n - sig_n[0]
    eps_n = eps_n - eps_n[0]
    
    ## Only get acting stress below 0.02 strain
    sig_n = sig_n[eps_n<=eps_thress]
    eps_n = eps_n[eps_n<=eps_thress]

# ---------------- Define symbolic variables ---------------- #
# Material parameters
C10, C01, C20, D1 = sp.symbols('C10 C01 C20 D1')

# ---------------- Define placeholders for invariants ---------------- #
I1b, I2b, J_sym = sp.symbols('I1b I2b J')

# Principal stretches
lambda_11, lambda_22, lambda_33 = sp.symbols('lambda_11 lambda_22 lambda_33', positive=True)

# ---------------- Define strain energy function ---------------- #
W = (
    C10 * (I1b - 3) +
    C01 * (I2b - 3) +
    C20 * (I1b - 3)**2 +
    (1 / D1) * (J_sym - 1)**2
)

# ---------- Create combined list for symbolic deriviation -----  #
symbolic_combi_list = (lambda_11,lambda_22,lambda_33,
                       C10, C01, C20, D1)

# ---------- Create symbolic list for the VUMAT.f file ---------- #
symbolic_param_list = [C10,C01,C20,D1]
symbolic_deriv_list = [sp.diff(W, I1b), sp.diff(W, I2b), sp.diff(W, J_sym)]
symbolic_namin_list = ['dWdI1','dWdI2','dWdJ1']
  
def PredictionStatement(params,StressFunction,Xi,nu = 0.495):
    ## Set material parameters/coefficients
    C10, C01, C20 = params[0], params[1], params[2]
    ## Compute D1 based on the poissions ratio and C10 and C01 parameters
    D1 = 3.0*(1.0 - 2.0*nu)/(2.0*(1.0+nu)*(C01+C10))
    ## Compute the stretches in the x- and z-directions, considering volumetric changes
    lam2 = 1.0 + Xi
    lam1 = lam2**(-nu)
    lam3 = lam2**(-nu)
    ## Compute prediction statement based on input
    Ypred = StressFunction(lam1,lam2,lam3,C10,C01,C20,D1)
    ## Return output
    return Ypred
    
def ObjectiveFunction(params,PredictionStatement_i,StressFunction_i,Xi,Yi):
    ## Compute prediction statement
    Yj = PredictionStatement_i(params,StressFunction_i,Xi)
    ## Compute sum of squared differences
    SSD = (1/len(Yi))*np.sum((Yj-Yi)**2)
    return SSD

def EnergyConstraint(params,EnergyFunction,Xi,nu = 0.495):
    ## Set material parameters/coefficients
    C10, C01, C20 = params[0], params[1], params[2]
    ## Compute D1 based on the poissions ratio and C10 and C01 parameters
    D1 = 3.0*(1.0 - 2.0*nu)/(2.0*(1.0+nu)*(C01+C10))
    ## Compute the stretches in the x- and z-directions, considering volumetric changes
    lam2 = 1.0 + Xi
    lam1 = lam2**(-nu)
    lam3 = lam2**(-nu)
    ## Compute the energy function
    Wvals = EnergyFunction(lam1,lam2,lam3,C10,C01,C20,D1)   
    return Wvals

       
## Get the stress function
P22_total = FirstPiolaKirschoffStress(W,I1b,I2b,J_sym,lambda_11,lambda_22,lambda_33)

## Get the modified energy-formulation
W22_total = EnergyInvariantModified(W,I1b,I2b,J_sym,lambda_11,lambda_22,lambda_33)

## Create function statement for evaluation
P22_func = lambdify(symbolic_combi_list,P22_total,modules='numpy')

## Get the modified energy-formulation
W22_func = lambdify(symbolic_combi_list,W22_total,modules='numpy')
 
## Initial guess
coefs = np.array([1,1,1])

## Update constraints
constraints = ({'type': 'ineq', 'fun': lambda params: EnergyConstraint(params, W22_func, eps_n)},
               {'type': 'ineq', 'fun': lambda params: -params[0]},
               {'type': 'ineq', 'fun': lambda params: params[1]},
               {'type': 'ineq', 'fun': lambda params: params[2]})

#Call minimize with SLSQP and constraints
solution = minimize(ObjectiveFunction, coefs, args=(PredictionStatement, P22_func, eps_n, sig_n), constraints=constraints, method='SLSQP', options={'disp': True, 'maxiter': 3000})

## Set fitting parameters to the objective function
model_coef_opt = solution.x

print('The optimization parameters are: ')
print(model_coef_opt)

plt.figure()
plt.plot(eps_n,sig_n)
plt.plot(eps_n,PredictionStatement(model_coef_opt, P22_func, eps_n))
plt.show()
    
# %% Generate VUMAT fortran file
GenerateVumatHyperelasticity(W,
                             symbolic_param_list,
                             symbolic_deriv_list,
                             symbolic_namin_list,
                             template_name = 'VUMAT_2D_planestrain_incompressible_template.f')
