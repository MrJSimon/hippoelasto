# Quickstart: Custom Hyperelastic Model Calibration & VUMAT Export

This guide walks you through the full workflow in **hippoelasto**, from defining a custom hyperelastic model to exporting a VUMAT file for use in Abaqus.

## 🔧 Step 3: Compute Symbolic Quantities

We start by computing the symbolic stress and energy terms needed for calibration and VUMAT generation.

```python
# Compute First Piola–Kirchhoff stress
P22_total = FirstPiolaKirschoffStressIncompressible(
    W, I1b, I2b, J_sym,
    lambda_11, lambda_22, lambda_33
)

# Compute modified strain energy for stability
W_modified = EnergyInvariantModified(
    W, I1b, I2b, J_sym,
    lambda_11, lambda_22, lambda_33
)
```

## 🧠 Step 4: Create Callable Functions

Lambdify the symbolic expressions to evaluate them numerically during optimization.

```python
# Create callable stress function
P22_func = lambdify(symbolic_combi_list, P22_total, modules='numpy')

# Create callable energy function
W_func = lambdify(symbolic_combi_list, W_modified, modules='numpy')
```

## 🎯 Step 5: Optimize Model Parameters

Define the objective function and constraints, then perform parameter optimization using SLSQP.

```python
# Initial parameter guess
coefs = np.ones(len(symbolic_combi_list[3:]))

# Define energy-based constraint
constraints = ({
    'type': 'ineq',
    'fun': lambda params: EnergyConstraintTension(params, W_func, eps_n)
})

# Define optimization arguments
args = (PredictionStatementTension, P22_func, eps_n, sig_n)

# Run optimization
model_coef_opt = OptimizationSLSQP(
    ObjectiveFunctionSSD,
    coefs,
    args,
    constraints=constraints
)
```

##  📝 Step 6: Generate VUMAT

Automatically export a VUMAT subroutine based on the symbolic model and fitted parameters.

```python
GenerateVumatHyperelasticity(
    W,
    symbolic_param_list,
    symbolic_deriv_list,
    symbolic_namin_list,
    template_name='VUMAT_2D_planestrain_incompressible_template.f'
)
```

## 📉 Step 7: Visualize Fit

Plot the experimental and predicted stress-strain curves.

```python
# Create artificial prediction range
eps_p = np.linspace(
    np.min(eps_n) - np.min(eps_n)/10,
    np.max(eps_n) + np.min(eps_n)/10,
    num=50,
    endpoint=True
)

# Evaluate prediction
sig_p = PredictionStatementTension(model_coef_opt, P22_func, eps_p)

# Plot and save
plotStressStrainCurve(eps_n, sig_n, eps_p, sig_p)
```
