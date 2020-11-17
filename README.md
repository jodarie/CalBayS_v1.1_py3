# CalBayS_v1.1_py3
### Calorimetry data analysis using Bayesian Statistics for python 3

@author Jose David Rivera

## Introduction

Script to perform statistical analysis via Bayesian approach of data obtained from ITC experiments.
CalBayS_v1 allows using more than one curve with different Wiseman C-parameter [DOI:10.1016/0003-2697(89)90213-3](https://www.sciencedirect.com/science/article/pii/0003269789902133?via%3Dihub) in a global fit, in order to compute the affinity and other thermodynamics parameters with high precision.

## Dependencies

CalBayS_v1.1 was run in Ubuntu 20.04 using the nohup option. The script is python 3.7.6 compatible and it requires the following packages:

- numpy
- matplotlib  
- scipy
- corner
- emcee
- random

## Configuration Parameters

The parameters should be set in params_all.py file.

#### Parallelization parameter

- run_nohup --> (True/False) Run in background the process, parallel when more than one model is used.

#### Input data   

To run the script the user must create files with each corresponding column to each experiment with different Wiseman C-parameter.

- Conc_protein_data --> Path + file name of protein concetration (&mu;M).
- Conc_ligand_data --> Path + file name of ligand concetration (&mu;M).
- DH_data --> Path + file name of enthalpy/molar.
- Volume_injection_data --> Path + file name of injection volume (ml).
- Volume_cell --> Vector of cell volume (ml).
- Conc_injection --> Vector of injection concentration (&mu;M).

#### Output options 

- color --> Vector of the colors for each curve profile and data scatter figures.
- File_out --> Initial name of out file.

#### MCMC parameters

- Number_of_walkers --> Number of chains.
- Number_of_points --> Number of points in each chain.
- burn_in --> Number of points of burn in.

Note that the run time and reliability of calculation depend on number of chains and  points set here. For more details, the user is suggested to read [emcee documentation](https://emcee.readthedocs.io/en/stable/).  

#### Hill parameters 

- Hill --> (True/False) Set if the Hill model will be used.
- Boundary_conditions_Hill --> Tuple of boundary conditions (minnHill,maxnHill,minKd,maxKd,minDH,maxDH,minq0,maxq0).
- ePrior_Hill --> (True/False) Set Gaussian prior information for &epsilon; parameters, see the section below.
- eP_Hill --> Tuple of center values for &epsilon;.
- deP_Hill --> Tuple of sigma values for &epsilon;.

#### _n_-Independent parameters

- Indp --> (True/False) Set if _n_-independent sites model will be used.
- Boundary_conditions_Indp --> Tuple of boundary conditions (minn,maxn,minK,maxK,minDH,maxDH,minq0,maxq0).
- nPrior_Indp --> (True/False) Set Gaussian prior information for _n_ parameter, see the section below.
- nP_Indp --> Center value for _n_.
- dnP_Indp --> Sigma value for _n_.
- ePrior_Indp --> (True/False) Set Gaussian prior information for &epsilon; parameters, see the section below.
- eP_Indp --> Tuple of center values for &epsilon;.
- deP_Indp --> Tuple of sigma values for &epsilon;.

#### Stepwise parameters

- Stepwise --> (True/False) Set if Stepwise model will be used.
- Boundary_conditions_sw --> Tuple of boundary conditions (minlogK,maxlogK,minDH,maxDH,minq0,maxq0). 
- n_model_sw --> Define the number of Kd in Step-Wise model.
- ePrior_sw --> (True/False) Set Gaussian prior information for &epsilon; parameters, see the section below.
- eP_sw --> Tuple of center values for &epsilon;.
- deP_sw --> Tuple of sigma values for &epsilon;.

#### Stepwise parameters with equal &Delta;H<sub>bind</sub>

- Stepwise_eqDH --> (True/False) Set if stepwise model with equal DH for all binding sites will be used.
- Boundary_conditions_sweqDH --> Tuple of boundary conditions (minlogK,maxlogK,minDH,maxDH,minq0,maxq0). 
- n_model_sweqDH --> Define the number of Kd in Step-Wise model.
- ePrior_sweqDH --> (True/False) Set Gaussian prior information for &epsilon; parameters, see the section below.
- eP_sweqDH --> Tuple of center values for &epsilon;.
- deP_sweqDH --> Tuple of sigma values for &epsilon;.

## Parameters to fit

#### Parameters of interest

- _n_<sub>Hill</sub>: Hill coefficient.
- _K<sup>(i)</sup><sub>d</sub>_: Dissociation constant (_i_ for each event in stepwise model).
- _n_: Stoichiometry value (_n_-independent sites model). 
- &Delta;_H_<sup>(_i_)</sup><sub>bind</sub>: Molar binding enthalpy (_i_ for each event in stepwise model). For stepwise model with equal &Delta;_H_<sub>bind</sub>, the script computes a mean molar binding enthalpy for every binding site.

#### Nuisance parameters 

Here _i_ correponds to each measurement.

- &sigma;<sub>i</sub> : Measurement error.
- [_L_]<sup>(i)</sup><sub>0</sub>: Concentration of ligand in the syringe.
- &epsilon;<sub>i</sub>: Calibration factor of the initial protein concentration in the cell.
- &Delta;_H_<sup>(_i_)</sup><sub>0</sub>: Molar enthalpy offset. 


## Running the code

After download the code and setting the file params_all.py, the user can run it through the following command line:

```
$ python MAIN_RUN.py
```

The code has an example by using 4 titration curves. 


## Attribution

Please cite Cardoso, Rivera et al. (2020) if you find this code useful in your research.
The BibTeX entry for the paper is:

```
@article{CARDOSO2020,
title = "CALX-CBD1 Ca2+-binding cooperativity studied by NMR spectroscopy and ITC with Bayesian statistics",
journal = "Biophysical Journal",
year = "2020",
issn = "0006-3495",
doi = "https://doi.org/10.1016/j.bpj.2020.05.031",
url = "http://www.sciencedirect.com/science/article/pii/S0006349520304549",
author = "M.V.C. Cardoso and J.D. Rivera and P.M. Vitale and M.F.S. Degenhardt and L.A. Abiko and C.L.P. Oliveira and R.K. Salinas",
abstract = "he Na+/Ca2+ exchanger of Drosophila melanogaster, CALX, is the main Ca2+-extrusion mechanism in olfactory sensory neurons and photoreceptor cells. Na+/Ca2+ exchangers have two Ca2+ sensor domains, CBD1 and CBD2. In contrast to the mammalian homologues, CALX is inhibited by Ca2+-binding to CALX-CBD1, while CALX-CBD2 does not bind Ca2+ at physiological concentrations. CALX-CBD1 consists of a β-sandwich and displays four Ca2+ binding sites at the tip of the domain. In this study, we used NMR spectroscopy and isothermal titration calorimetry (ITC) to investigate the cooperativity of Ca2+-binding to CALX-CBD1. We observed that this domain binds Ca2+ in the slow exchange regime at the NMR chemical shift time scale. Ca2+-binding restricts the dynamics in the Ca2+-binding region. Experiments of 15N CEST and 15N R2 dispersion allowed the determination of Ca2+ dissociation rates (≈30s-1). NMR titration curves of residues in the Ca2+-binding region were sigmoidal due to the contribution of chemical exchange to transverse magnetization relaxation rates, R2. Hence, a novel approach to analyze NMR titration curves was proposed. Ca2+-binding cooperativity was examined assuming two different stoichiometric binding models and using a Bayesian approach for data analysis. Fittings of NMR and ITC binding curves to the Hill model yielded nHill∼2.9, near maximum cooperativity (nHill=4). By assuming a stepwise model to interpret the ITC data, we found that the probability of binding from 2 up to 4 Ca2+ is approximately three orders of magnitude higher than that of binding a single Ca2+. Hence, four Ca2+ ions bind almost simultaneously to CALX-CBD1. Cooperative Ca2+-binding is key to enable this exchanger to efficiently respond to changes in the intracellular Ca2+-concentration in sensory neuronal cells."
}
```
