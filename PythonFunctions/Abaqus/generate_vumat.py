##############################################################################
##
## Author:      Jamie E. Simon
##
## Description: 
##              
##############################################################################


## -------- ABAQUS INPUT FORMATTING --------------- ##
from pathlib import Path
from PythonFunctions.Abaqus.fortran_formatting import fortran_d0_lines

def GenerateVumatHyperelasticity(StrainEnergyDensity,
                                 MaterialPropsParam,
                                 StrainEnergyDerivativeExprs,
                                 StrainEnergyDerivativeNames,
                                 template_name = 'VUMAT_2D_planestrain_template.f'):
  
    ## Set main path
    main_path = "PythonFunctions\\Abaqus\\templates\\"
    
    ## Set modified template name
    output_name = template_name.split('_template.f')[0] + "_modified.f"
    
    ## Load template
    template_path = Path(main_path + template_name)
      
    ## Set output path
    output_path = Path("output\\"+output_name)
    
    ## Set identifiers for input text
    ID_1 = "*** INPUT FROM PYTHON PROGRAM *** STRAIN ENERGY DEFINITION"
    ID_2 = "C 	  *** INPUT FROM PYTHON PROGRAM *** MATERIAL PARAMETERS"
    ID_3 = "C		 *** INPUT FROM PYTHON PROGRAM *** DERIVATIVE OF STRAIN-ENERGY FUNCTION"
    ID_4 = "C 	  *** INPUT FROM PYTHON PROGRAM *** MATERIAL INITIATION"
    
    ## Open template and read it in
    with template_path.open("r", encoding="utf-8") as f:
        output = f.read()
    
    ## Set strain energy density for input
    W_line = 'W = ' + str(StrainEnergyDensity)
    
    ## Set strain energy parameter properties for input
    M_line = "\n".join(
        f"      {param:<3} = props({i+1})"
        for i, param in enumerate(MaterialPropsParam)
    )
    
    ## Set elastic material parameter properties for input
    P_line = ",".join(
        f" {param}"
        for param in MaterialPropsParam
    )
    
    ## Add The badass real 8 baby !!!
    P_line = '	  ' + 'Real*8' + P_line
    
    ## Set partial derivative of the strain energy w. respect to the invariatns
    D_line = fortran_d0_lines(StrainEnergyDerivativeNames,
                              StrainEnergyDerivativeExprs)
        
    ## Replace Identifiers with the appropriate definitions
    output = output.replace(ID_1,W_line,1)
    output = output.replace(ID_2,M_line,1)
    output = output.replace(ID_3,D_line,1)
    output = output.replace(ID_4,P_line,1)
      
    ## Write to output
    with output_path.open("w", encoding="utf-8") as f:
        f.write(output)

    return print('Fortran file generated and saved to output')
