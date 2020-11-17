import numpy as np
import bin.Functions as func
import matplotlib.pyplot as plt
#plt.rc('text', usetex=True)
import params 
import os

boCond = params.Boundary_conditions
Vc = params.Volume_cell
out = params.File_out
Nw = params.Number_of_walkers
Nc = params.Number_of_points
burnin = params.burn_in
os.system('mkdir '+out)
So = np.array(params.Conc_injection)

nconc = len(So)
if nconc==1:
  P = np.loadtxt(params.Conc_protein_data)*1.0e3 # uM
  X = np.loadtxt(params.Conc_ligand_data)*1.0e3 # uM
  DH = np.loadtxt(params.DH_data)*1.0e-3 # kcal/mol
  v = np.loadtxt(params.Volume_injection_data)*1.0e-3 # ml
  ntmp = len(P)
  Pmatrix = np.zeros((ntmp,2))
  Xmatrix = np.zeros((ntmp,2))
  DHmatrix = np.zeros((ntmp,2))
  vi = np.zeros((ntmp,2))
  Pmatrix[:,0] = P
  Xmatrix[:,0] = X
  DHmatrix[:,0] = DH
  vi[:,0] = v
else: 
  Pmatrix = np.loadtxt(params.Conc_protein_data)*1.0e3 # uM
  Xmatrix = np.loadtxt(params.Conc_ligand_data)*1.0e3 # uM
  DHmatrix = np.loadtxt(params.DH_data)*1.0e-3 # kcal/mol
  vi = np.loadtxt(params.Volume_injection_data)*1.0e-3 # ml

bo,ndim,lab = func.tools_fit(So,boCond[0],boCond[1],boCond[2],boCond[3],boCond[4],boCond[5],boCond[6],boCond[7])

eArgs = [params.ePrior,params.eP,params.deP]

pos_i,chain = func.Best_fit(Pmatrix,Xmatrix,DHmatrix,Vc,vi,So,bo,ndim,Nw,Nc,lab,eArgs,params.ns)
np.savetxt(out+'/Chains.dat',chain,fmt='%s')
func.Chain_figure(chain,out,ndim,Nw,Nc,lab)
func.Correlation_figure(chain,out,Nw,Nc,burnin,lab)
res = func.MCMC_analysis(chain,out,Nw,Nc,burnin)
P = res[:,0]
func.Figure_fit(P,out,Pmatrix,Xmatrix,DHmatrix,Vc,vi,So,params.color,params.ns)
