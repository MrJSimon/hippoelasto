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
from sympy import lambdify
from PythonFunctions.StressDescription.piola_kirschoff_stress import FirstPiolaKirschoffStress
from PythonFunctions.EnergyDescription.energy_substitution import EnergyInvariantModified
from PythonFunctions.Abaqus.generate_vumat import GenerateVumatHyperelasticity
from PythonFunctions.Optimization.optimization_routines import PredictionStatementTension
from PythonFunctions.Optimization.optimization_routines import ObjectiveFunctionSSD
from PythonFunctions.Optimization.optimization_routines import EnergyConstraintTension
from PythonFunctions.Optimization.optimization_routines import OptimizationSLSQP
from PythonFunctions.PlottingFunctions.plotting_functions import plotStressStrainCurve
from PythonFunctions.TangentModulus.numerical_elastic_modulus import NumericalElasticModulus

## ------------------------------ DATA INPUT ------------------------------ ##

# Load in data
data = np.loadtxt('data\\nominal_stress_strain_data.txt',delimiter=',')

# Set strain- and stress data
eps_n, sig_n = data[:,0], data[:,1]

## ----------------------- STRAIN ENERGY DEFINITION ----------------------- ##

# Assume poisons ration
nu = 0.1

# Define symbolic variables
C10, C01, C20, D = sp.symbols('C10 C01 C20 D')

# Define placeholders for modified invariants
I1b, I2b, J_sym = sp.symbols('I1b I2b detJ')

# Define stretches for present in a uniaxial stretch state
lambda_11, lambda_22, lambda_33 = sp.symbols('lambda_11 lambda_22 lambda_33', positive=True)

# Define strain energy function
W = C10 * (I1b - 3) + C01 * (I2b - 3) + C20 * (I1b - 3)**2 + (1 / D) * (J_sym - 1)**2

# Create combined list for symbolic deriviation
symbolic_combi_list = (lambda_11,lambda_22,lambda_33,
                       C10, C01, C20, D)

# Create symbolic list for the VUMAT.f file
symbolic_param_list = [C10,C01,C20,D]
symbolic_deriv_list = [sp.diff(W, I1b), sp.diff(W, I2b), sp.diff(W, J_sym)]
symbolic_namin_list = ['dWdI1','dWdI2', 'dWdJ']
symbolic_mater_list = [str(param) for param in symbolic_param_list] + ['E', 'nu']
        
# Get the stress function
P22_total = FirstPiolaKirschoffStress(W,I1b,I2b,J_sym,lambda_11,lambda_22,lambda_33)

# Get the modified energy-formulation
W_modified = EnergyInvariantModified(W,I1b,I2b,J_sym,lambda_11,lambda_22,lambda_33)

# Create function statement for evaluation
P22_func = lambdify(symbolic_combi_list,P22_total,modules='numpy')

# Get the modified energy-formulation
W_func = lambdify(symbolic_combi_list,W_modified,modules='numpy')
 
# Initial guess
coefs = np.ones(len(symbolic_combi_list[3:]))

# Construct constraints
constraints = ({'type': 'ineq', 'fun': lambda params: EnergyConstraintTension(params, W_func, eps_n, nu = nu)})

# Construct args
args = (PredictionStatementTension, P22_func, eps_n, sig_n, nu)

# Conduct optimization and get best parameters
model_coef_opt = OptimizationSLSQP(ObjectiveFunctionSSD, coefs, args, constraints = constraints)

print('The optimization parameters are: ')
print(model_coef_opt)

## --------------------------- GENERATE VUMAT ------------------------------ ##

# % Generate VUMAT fortran file
GenerateVumatHyperelasticity(W,
                          symbolic_mater_list,
                          symbolic_deriv_list,
                          symbolic_namin_list,
                          template_name = 'VUMAT_2D_planestrain_template.f')


## ------------------- PLOT AND SAVE FIGURE TO OUTPUT----------------------- ##

# Make artifical X-range
eps_p = np.linspace(np.min(eps_n)-np.min(eps_n)/10,np.max(eps_n)+np.min(eps_n)/10,num=50,endpoint = True)

# Compute prediction statement
sig_p = PredictionStatementTension(model_coef_opt, P22_func, eps_p, nu)

# Plot stress strain curve and save to output
plotStressStrainCurve(eps_n,sig_n,eps_p,sig_p)

##
tangent_modulus = np.gradient(sig_p, eps_p)

import matplotlib.pyplot as plt

plt.figure()
plt.plot(eps_p,tangent_modulus)
plt.show()

# %%

## Tangent modulus estimation

# Add a small positive offset to avoid log10(0)
eps_min = np.min(eps_n[np.nonzero(eps_n)])  # get smallest non-zero strain
eps_max = np.max(eps_n)

# Define log bounds safely
Xlog1 = np.log10(eps_min * 0.1)
Xlog2 = np.log10(eps_max * 1.1)

# Make artifical X-range
eps_p = np.logspace(Xlog1,Xlog2,num=500,endpoint = True)

sig_p = PredictionStatementTension(model_coef_opt, P22_func, eps_p, nu)

tangent_modulus = np.gradient(sig_p, eps_p)

plt.figure()
plt.plot(eps_p,tangent_modulus)
plt.xlabel('eps_{nominal} [-]')
plt.ylabel('tangent moduli [MPa]')
plt.show()


oink = NumericalElasticModulus(eps_n,PredictionStatementTension, P22_func, model_coef_opt, nu)

print(oink)