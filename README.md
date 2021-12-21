Slower organisms exhibit sudden population disappearances in a reddened world
=========

#### Authors
Christopher J. Greyson-Gaito<sup>*1</sup>, Gabriel Gellner<sup>1</sup>, Kevin S. McCann<sup>1</sup>
----------

### Affiliations
*Corresponding Author - christopher@greyson-gaito.com

1. Department of Integrative Biology, University of Guelph, Guelph, ON, Canada

## ORCID
* CJGG &ndash; 0000-0001-8716-0290
* GG &ndash; 0000-0001-8170-1463
* KSM &ndash; 0000-0001-6031-7913

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5796525.svg)](https://doi.org/10.5281/zenodo.5796525)

## Julia scripts and datasets

### Folder and file structure
* data &ndash; empty folder for data files to be placed (created in slowfast_whitenoise.jl and slowfast_rednoise.jl)
* figs &ndash; empty folder for figures to be placed (created in slowfast_figures.jl)
* scripts
    * packages.jl &ndash; list of packages required (file used in other scripts)
    * quasipotential.R &ndash; R script file to produce quasipotential figures in Supporting Information
    * slowfast_canardfinder.jl &ndash; julia script file containing code for the canard finder algorithm (file used in other scripts)
    * slowfast_commoncode.jl &ndash; julia script file containing common code used in other scripts
    * slowfast_metabolicmodel.jl &ndash; julia script file containing code producing the Yodzis and Innes metabolic model (file used in other scripts)
    * slowfast_figures.jl &ndash; julia script to produce the figures in the manuscript
    * slowfast_rednoise.jl &ndash; julia script to run the canard finder algorithm with red noise
    * slowfast_whitenoise.jl &ndash; julia script to run the canard finder algorithm with white noise
* .gitignore &ndash; file containing files and folders that git should ignore
* LICENSE.txt &ndash; CC by 4.0 License for this repository
* README.md &ndash; this file

### Instructions

1. Download the GitHub/Zenodo repo
2. Open the repo in Visual Studio Code (if you haven't already done so, set up [Julia in Visual Studio Code](https://www.julia-vscode.org/))
3. Most of the analysis here requires multiple cores. Thus to set up multiple cores on your computer, in Visual Studio Code find the Julia: Num Threads setting in the Extension Settings of the Visual Studio Code Julia Extension. Change this setting to at most the number of logical cores in your computer. Restart Julia. All parallel computing will run automatically regardless of the number of cores selected.
4. Run slowfast_whitenoise.jl and slowfast_rednoise.jl to create the data required for slowfast_figures.jl. Note, depending on the number of cores in your computer, this will take a long time.
5. Run slowfast_figures.jl to produce the figures in the manuscript.