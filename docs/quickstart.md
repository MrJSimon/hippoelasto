# Quickstart: Custom Hyperelastic Model Calibration & VUMAT Export

This guide walks you through the full workflow in **hippoelasto**, from defining a custom hyperelastic model to exporting a VUMAT file for use in Abaqus.

## 📂 Step 1: Load data
The data file must be comma-separated and contain the following columns:

| strain [-] | stress [MPa] |
|---------------:|--------------:|
| 0.0          | 0.069           |
| 0.011          | 4.109         |
| 0.022          | 6.654         |
| ⋮              | ⋮              |
| 0.066          | 10.246        |
| 0.077          | 10.483        |
| 0.088          | 10.643        |
| ⋮              | ⋮              |

Place your tensile data file (e.g., yourfile.txt) in the data/ folder of the repository.

```python
# Load in data
data = np.loadtxt('data\\nominal_stress_strain_data.txt',delimiter=',')

# Set strain- and stress data
eps_n, sig_n = data[:,0], data[:,1]
```

## ⚙️ Step 2: Choose Material Model
**Specify a custom hyperelastic model** (e.g., Neo-Hookean, Mooney-Rivlin, Ogden) in terms of strain energy potential

```python
# Define symbolic variables
C10, C01, C20 = sp.symbols('C10 C01 C20')

# Define placeholders for modified invariants
I1b, I2b, J_sym = sp.symbols('I1b I2b J')

# Define stretches uniaxial state
lambda_11, lambda_22, lambda_33 = sp.symbols('lambda_11 lambda_22 lambda_33', positive=True)

# Define strain energy function
W = C10 * (I1b - 3) + C01 * (I2b - 3) + C20 * (I1b - 3)**2

# Create combined list for symbolic derivation
symbolic_combi_list = (lambda_11,lambda_22,lambda_33,
                       C10, C01, C20)

# Create symbolic list for the VUMAT.f file
symbolic_param_list = [C10,C01,C20]
symbolic_deriv_list = [sp.diff(W, I1b), sp.diff(W, I2b)]
symbolic_namin_list = ['dWdI1','dWdI2']
```

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
