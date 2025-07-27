# hippoelasto
**hippoelasto** is a Python package for hyperelastic material characterization and calibration, and VUMAT generation for finite element modeling in Abaqus.

## Purpose  
**hippoelasto** is built for researchers and engineers working with soft material modeling and simulation. It provides a full pipeline—from experimental data fitting to VUMAT generation—for hyperelastic material models.

1. Load and manage experimental stress-strain data 
1. Select and configure hyperelastic material models 
3. Calibrate model parameters using optimization routines
4. Visualize fitted results and assess model performance
5. Automatically generate a VUMAT subroutine for Abaqus

# Installation
Install **hippoelasto** by cloning the repository onto your local machine using the following command

    git clone https://github.com/MrJSimon/hippoelasto.git

### Requirements
This program was tested using

    python 3.10.9

Install the required Python packages, listed in requirements.txt, with:  

    pip install -r requirements.txt

# Getting started
This package is intended to be run through a Python IDE (PyCharm, VSCode, Spyder, etc.). 

    python asm\asm-package\start_GUI.py

# Basic workflow

### 📂 Step 1: Load data
The data file must be comma-separated and contain the following columns:

| strain [-] | stress [MPa] |
|---------------:|--------------:|
| 0.0          | 0.069           |
| 0.011          | 4.109         |
| 0.022          | 6.654         |
| ⋮              | ⋮              |
| 0.033          | 8.379         |
| 0.044          | 9.44          |
| 0.055          | 9.977         |
| ⋮              | ⋮              |
| 0.066          | 10.246        |
| 0.077          | 10.483        |
| 0.088          | 10.643        |
| ⋮              | ⋮              |

Place your tensile data file (e.g., yourfile.txt) in the data/ folder of the repository.

⚙️ Step 2: Choose Material Model
Specify custom hyperlastic model (e.g., Neo-Hookean, Mooney-Rivlin, Ogden) in terms of strain energy potential



TBD: Select from built-in hyperelastic models (e.g., Neo-Hookean, Mooney-Rivlin, Ogden).

You can also specify custom strain energy functions.




### 🔖 Step 2: Label Manager
Use the plus/+ or minus/- sign to add or remove labels. The program currently supports 9 labels.

### 🎨 Step 3: Paint tools
To mask an image start by:

&nbsp;&nbsp;&nbsp;&nbsp; **3.1** Activate the label of interest by pressing the added label-button  
&nbsp;&nbsp;&nbsp;&nbsp; **3.2** Press the brush tool button  
&nbsp;&nbsp;&nbsp;&nbsp; **3.3** Guide the mouse to the main image window and paint on top of the image 

Please note that the brush and eraser tool will only work on the image, not on the entire canvas.

### 📊 Step 3: Select Features
In the **Random Forest Classifier** panel, choose the features to include in training (e.g., Sobel, Canny Edge, Gaussian filters).

### 🌲 Step 4: Random Forest Classifier
Click **Train** to build a Random Forest model using the selected features and masks.

### 🔮 Step 5: Predict Segmentation
Use the **Predict** button to apply the trained model to the image shown in the main GUI window. 

Predictions will appear in a **separate** window next to the original image.

### 🤖 Step 5: Automation Manager
To apply the trained model to multiple images:

&nbsp;&nbsp;&nbsp;&nbsp; 5.1 Use the **Automation Manager** to select or deselect the images to be predicted  
&nbsp;&nbsp;&nbsp;&nbsp; 5.2 Press the **Predict** button to apply the trained model to the selected images  
&nbsp;&nbsp;&nbsp;&nbsp; 5.3 Lean back and relax while **asm** processes the batch

### 📝 Step 6: Save Results
Click the **save json** button to save the current session information in a 'config.json' file. The file contains all necessary metadata to restore the session later.
