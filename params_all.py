#################################
#          CalBayS_v1           #
#################################
#      Jose D Rivera 2019       #
#################################

##################################
#   Write your parameters here   #
##################################

run_nohup = True # Run in background the process, parallel when more than one model is used

##################################
#          Input data            #
##################################

Volume_cell = [1.4397,1.4397,1.4397,1.4397] # ml Volume of cell
Conc_injection = [500,500,500,500] # uM
Conc_protein_data = "Data/Protein.dat" # Path + file name of protein concentration data 
Conc_ligand_data = "Data/Ligand.dat" # Path + file name of ligand concentration data
DH_data = "Data/DH.dat" # Path + file name of enthalpy/molar data
Volume_injection_data = "Data/Volume_inj.dat" # Path + file name of injection volume data

##################################
#         Output options         #
##################################

color = ["black","red","blue","green","gold","pink","cyan","peru","limegreen"] # Define the colors of the curve profile and data scatter figures
File_out = "Results" # Initial name of out file

##################################
#        MCMC parameters         #
##################################

Number_of_walkers = 100 # Number of chains
Number_of_points = 1500 # Number of points in each chain
burn_in = 1000 # Burn in

##################################
#        Hill parameters         #
##################################

Hill = True # True --> Hill model is used. False --> Hill model is not used
Boundary_conditions_Hill = (0,8,0.01,100,-30,2,-5,5) # (minnHill,maxnHill,minKd,maxKd,minDH,maxDH,minq0,maxq0)
ns = 4.0 # Set the stoichiometry

##################################
#   Set Gaussian priors (Hill)   #
##################################

ePrior_Hill = True # Set Gaussian prior information for eps_i parameters
eP_Hill = (0.8988835016857637,0.9155636199093159,1.1499269245534458,1.0389285760370357) # Center values
deP_Hill = (0.0067897586725015335,0.005619155810149712,0.004799042776802542,0.0038638395087877697) # Sigma values

##################################
#    nIndependent parameters     #
##################################

Indp = True # True --> Independent model is used. False --> Independent model is not used
Boundary_conditions_Indp = (0,8,0.01,100,-30,2,-5,5) # (minn,maxn,minK,maxK,minDH,maxDH,minq0,maxq0)

##################################
#   Set Gaussian priors (Indep)  #
##################################

nPrior_Indp = True # Set Gaussian prior information for n parameter
nP_Indp = 4.0 # Center value 
dnP_Indp = 0.01 # Sigma value
ePrior_Indp = False # Set Gaussian prior information for eps_i parameters
eP_Indp = (0.9158764649722917,1.0539921125160618,1.394476252895587,1.2692705895627583) # Center values
deP_Indp = (0.11328679607774439,0.10832811807268983,0.1320580864208163,0.11818704512176903) # Sigma values

##################################
#      Stepwise parameters       #
##################################

Stepwise = True # True --> Stepwise model is used. False --> Stepwise model is not used
Boundary_conditions_sw = (-6,6,-20,2,-5,5) #(minlogK,maxlogK,minDH,maxDH,minq0,maxq0) 
n_model_sw = 4 # Define the number of Kd in Step-Wise model

##################################
#    Set Gaussian priors (Sw)    #
##################################

ePrior_sw = True # Set Gaussian prior information for eps_i parameters
eP_sw = (0.8988835016857637,0.9155636199093159,1.1499269245534458,1.0389285760370357) # Center values
deP_sw = (0.0067897586725015335,0.005619155810149712,0.004799042776802542,0.0038638395087877697) # Sigma values


##################################
# Stepwise (equal DH) parameters #
##################################

Stepwise_eqDH = True # True --> Stepwise model with equal DH for all binding sites is used. False --> Stepwise model with equal DH for all binding sites is not used
Boundary_conditions_sweqDH = (-6,6,-20,2,-5,5) #(minlogK,maxlogK,minDH,maxDH,minq0,maxq0) 
n_model_sweqDH = 4 # Define the number of Kd in Step-Wise model

##################################
#  Set Gaussian priors (Sw eqDH) #
##################################

ePrior_sweqDH = True # Set Gaussian prior information for eps_i parameters
eP_sweqDH = (0.8988835016857637,0.9155636199093159,1.1499269245534458,1.0389285760370357) # Center values
deP_sweqDH = (0.0067897586725015335,0.005619155810149712,0.004799042776802542,0.0038638395087877697) # Sigma values

