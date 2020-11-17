import params_all as pall
import os

rnohup = pall.run_nohup
vcell = pall.Volume_cell
inj = pall.Conc_injection
prot = pall.Conc_protein_data
lig = pall.Conc_ligand_data
DH = pall.DH_data
vol_inj = pall.Volume_injection_data
color = pall.color
out_name = pall.File_out
num_w = pall.Number_of_walkers
num_p = pall.Number_of_points
burnin = pall.burn_in

if pall.Hill:
  f = open('params.py','w')
  f.write('Boundary_conditions = '+str(pall.Boundary_conditions_Hill)+'\n')
  f.write('Volume_cell = '+str(vcell)+'\n')
  f.write('File_out = "Hill_'+str(out_name)+'"\n')
  f.write('Number_of_walkers = '+str(num_w)+'\n')
  f.write('Number_of_points = '+str(num_p)+'\n')
  f.write('burn_in = '+str(burnin)+'\n')
  f.write('Conc_injection = '+str(inj)+'\n')
  f.write('Conc_protein_data = "'+str(prot)+'"\n')
  f.write('Conc_ligand_data = "'+str(lig)+'"\n')
  f.write('DH_data = "'+str(DH)+'"\n')
  f.write('Volume_injection_data = "'+str(vol_inj)+'"\n')
  f.write('color = '+str(color)+'\n')
  f.write('ns = '+str(pall.ns)+'\n')
  f.write('ePrior = '+str(pall.ePrior_Hill)+'\n')
  f.write('eP = '+str(pall.eP_Hill)+'\n')
  f.write('deP = '+str(pall.deP_Hill)+'\n')
  f.close()
  os.system('mv params.py Models/Hill/')
  if rnohup == True:
    os.system('nohup python Models/Hill/Fit_Hill.py > out_Hill 2> err_Hill &')
  else:
    os.system('python Models/Hill/Fit_Hill.py')

if pall.Indp:
  f = open('params.py','w')
  f.write('Boundary_conditions = '+str(pall.Boundary_conditions_Indp)+'\n')
  f.write('Volume_cell = '+str(vcell)+'\n')
  f.write('File_out = "Indp_'+str(out_name)+'"\n')
  f.write('Number_of_walkers = '+str(num_w)+'\n')
  f.write('Number_of_points = '+str(num_p)+'\n')
  f.write('burn_in = '+str(burnin)+'\n')
  f.write('Conc_injection = '+str(inj)+'\n')
  f.write('Conc_protein_data = "'+str(prot)+'"\n')
  f.write('Conc_ligand_data = "'+str(lig)+'"\n')
  f.write('DH_data = "'+str(DH)+'"\n')
  f.write('Volume_injection_data = "'+str(vol_inj)+'"\n')
  f.write('color = '+str(color)+'\n')
  f.write('nPrior = '+str(pall.nPrior_Indp)+'\n')
  f.write('nP = '+str(pall.nP_Indp)+'\n')
  f.write('dnP = '+str(pall.dnP_Indp)+'\n')
  f.write('ePrior = '+str(pall.ePrior_Indp)+'\n')
  f.write('eP = '+str(pall.eP_Indp)+'\n')
  f.write('deP = '+str(pall.deP_Indp)+'\n')
  f.close()
  os.system('mv params.py Models/Indep/')
  if rnohup == True:
    os.system('nohup python Models/Indep/Fit_independent.py > out_Indep 2> err_Indep &')
  else:
    os.system('python Models/Indep/Fit_independent.py')

if pall.Stepwise:
  f = open('params.py','w')
  f.write('Boundary_conditions = '+str(pall.Boundary_conditions_sw)+'\n')
  f.write('Volume_cell = '+str(vcell)+'\n')
  f.write('File_out = "SW_'+str(out_name)+'"\n')
  f.write('Number_of_walkers = '+str(num_w)+'\n')
  f.write('Number_of_points = '+str(num_p)+'\n')
  f.write('burn_in = '+str(burnin)+'\n')
  f.write('Conc_injection = '+str(inj)+'\n')
  f.write('Conc_protein_data = "'+str(prot)+'"\n')
  f.write('Conc_ligand_data = "'+str(lig)+'"\n')
  f.write('DH_data = "'+str(DH)+'"\n')
  f.write('Volume_injection_data = "'+str(vol_inj)+'"\n')
  f.write('color = '+str(color)+'\n')
  f.write('n_model = '+str(pall.n_model_sw)+'\n')
  f.write('ePrior = '+str(pall.ePrior_sw)+'\n')
  f.write('eP = '+str(pall.eP_sw)+'\n')
  f.write('deP = '+str(pall.deP_sw)+'\n')
  f.close()
  os.system('mv params.py Models/SW/')
  if rnohup == True:
    os.system('nohup python Models/SW/Fit_StepWise.py > out_SW 2> err_SW &')
  else:
    os.system('python Models/SW/Fit_StepWise.py')

if pall.Stepwise_eqDH:
  f = open('params.py','w')
  f.write('Boundary_conditions = '+str(pall.Boundary_conditions_sweqDH)+'\n')
  f.write('Volume_cell = '+str(vcell)+'\n')
  f.write('File_out = "SWeqDH_'+str(out_name)+'"\n')
  f.write('Number_of_walkers = '+str(num_w)+'\n')
  f.write('Number_of_points = '+str(num_p)+'\n')
  f.write('burn_in = '+str(burnin)+'\n')
  f.write('Conc_injection = '+str(inj)+'\n')
  f.write('Conc_protein_data = "'+str(prot)+'"\n')
  f.write('Conc_ligand_data = "'+str(lig)+'"\n')
  f.write('DH_data = "'+str(DH)+'"\n')
  f.write('Volume_injection_data = "'+str(vol_inj)+'"\n')
  f.write('color = '+str(color)+'\n')
  f.write('n_model = '+str(pall.n_model_sweqDH)+'\n')
  f.write('ePrior = '+str(pall.ePrior_sweqDH)+'\n')
  f.write('eP = '+str(pall.eP_sweqDH)+'\n')
  f.write('deP = '+str(pall.deP_sweqDH)+'\n')
  f.close()
  os.system('mv params.py Models/SW_eqDH/')
  if rnohup == True:
    os.system('nohup python Models/SW_eqDH/Fit_StepWise.py > out_SW_eqDH 2> err_SW_eqDH &')
  else:
    os.system('python Models/SW_eqDH/Fit_StepWise.py')
