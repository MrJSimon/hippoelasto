##############################################################################
##
## Author:      Jamie E. Simon
##
## Description: 
##
##############################################################################

## Import modulues
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#import matplotlib.lines as mlines


def initiate_font_settings():
    """
    This function initiates a latex style fonts to be used for plotting

    Returns
    -------
    None.

    """
    
    # Global settings for nicer fonts
    mpl.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'STIX', 'Georgia'],
        'font.size': 14,                     # base font size
        'axes.labelsize': 16,               # axis label font size
        'axes.titlesize': 16,               # title font size
        'legend.fontsize': 14,              # legend font size
        'xtick.labelsize': 12,              # x-axis tick font size
        'ytick.labelsize': 12,              # y-axis tick font size
        'mathtext.fontset': 'stix',         # for LaTeX-style math
    })


def plotStressStrainCurve(Xi, Yi, Xp, Yp,
                          picture_name = 'output//predictionvsdata.pdf',
                          fontsize_plot   = 16,
                          markersize_plot = 12,
                          labelsize_plot  = 15,
                          linewidth_plot  = 5):
    
    ## Initiate font settings
    initiate_font_settings()

    ## Set x- and ylabels
    xlabel, ylabel=r"$\epsilon_{\mathrm{nominal}}$ [-]", r"$\sigma_{\mathrm{nominal}}$ [MPa]"

    ## Set figure
    fig,axes = plt.subplots(1,1,figsize = (6,6), tight_layout = True)

    ## Plot prediction and data
    axes.plot(Xp,Yp,linestyle='-',linewidth=linewidth_plot,color='black',alpha=1.0,label='Prediction')
    axes.plot(Xi,Yi,linestyle='--',linewidth=linewidth_plot-2,color='red',alpha=1.0,label='Data')

    ## Set legend
    axes.legend(fontsize = fontsize_plot,edgecolor='gray',facecolor='white',framealpha=1.0)

    ## Set grid
    axes.grid(color='gray', linestyle='-', linewidth=1) 

    ## Set x and y labels
    axes.set_xlabel(xlabel,fontsize=fontsize_plot)
    axes.set_ylabel(ylabel,fontsize=fontsize_plot) 

    ## Set tick sizes
    axes.tick_params(axis='both', which='major', labelsize=labelsize_plot)

    ## Set Xscale
    Xmin, Xmax = np.min(Xp), np.max(Xp)
    Ymin, Ymax = np.min(Yp), np.max(Yp)
    axes.set_xlim([Xmin,Xmax])
    axes.set_ylim([Ymin,Ymax])

    ## Save figure !!!
    plt.savefig(picture_name, bbox_inches='tight')
    
    ## Print saving output
    print('Figure ' + picture_name + ' saved')
    
    plt.close()
    return



def plotTangenmodulus(Xp, Yp, Yt,
                          picture_name = 'output//tangentmodulus.pdf',
                          fontsize_plot   = 16,
                          markersize_plot = 12,
                          labelsize_plot  = 15,
                          linewidth_plot  = 5):  
    ## Compute gradient
    dXpdYp = np.gradient(Yp, Xp)
    
    ## Set Xscale
    Xmin, Xmax = np.min(Xp), np.max(Xp)
    Ymin, Ymax = np.min(Yp), np.max(Yp)
    
    ## Initiate font settings
    initiate_font_settings()

    ## Set x- and ylabels
    xlabel, ylabel=r"$\epsilon_{\mathrm{nominal}}$ [-]", r"$E_{\mathrm{tangent}}$ [MPa]"

    ## Set figure
    fig,axes = plt.subplots(1,1,figsize = (6,6), tight_layout = True)

    ## Plot prediction and data
    axes.plot(Xp,dXpdYp,linestyle='-',linewidth=linewidth_plot,color='black',alpha=1.0,label='Prediction')
    axes.plot([Xmin,Xmax],[Yt, Yt],linestyle='--',linewidth=linewidth_plot-2,color='red',alpha=1.0,label='Elastic modulus')

    ## Set legend
    axes.legend(fontsize = fontsize_plot,edgecolor='gray',facecolor='white',framealpha=1.0)

    ## Set grid
    axes.grid(color='gray', linestyle='-', linewidth=1) 

    ## Set x and y labels
    axes.set_xlabel(xlabel,fontsize=fontsize_plot)
    axes.set_ylabel(ylabel,fontsize=fontsize_plot) 

    ## Set tick sizes
    axes.tick_params(axis='both', which='major', labelsize=labelsize_plot)

    
    axes.set_xlim([Xmin,Xmax])
    #axes.set_ylim([Ymin,Ymax])

    ## Save figure !!!
    plt.savefig(picture_name, bbox_inches='tight')
    
    ## Print saving output
    print('Figure ' + picture_name + ' saved')
    
    plt.close()
    return


def saveMaterialParameters(param_names,param_values):
    
    # Combine into rows
    table_data = list(zip(param_names, param_values))

    # Create figure
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axis('off')

    # Table
    table = ax.table(cellText=table_data,
                     colLabels=["Parameter", "Value"],
                     loc='center',
                     cellLoc='center',
                     colLoc='center')
    table.scale(1, 1.5)


    # Save to PDF
    with PdfPages("output/model_parameters.pdf") as pdf:
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
    plt.close()
    return
