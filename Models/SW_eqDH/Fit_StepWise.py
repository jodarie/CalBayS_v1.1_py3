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

bo,ndim,lab = func.tools_fit(So,params.n_model,boCond[0],boCond[1],boCond[2],boCond[3],boCond[4],boCond[5])

eArgs = [params.ePrior,params.eP,params.deP]

pos_i,chain = func.Best_fit(Pmatrix,Xmatrix,DHmatrix,Vc,vi,So,bo,ndim,Nw,Nc,lab,params.n_model,eArgs)
np.savetxt(out+'/Chains.dat',chain,fmt='%s')
#chain = np.loadtxt(out+'/Chains.dat')
func.Chain_figure(chain,out,ndim,Nw,Nc,lab)
func.Correlation_figure(chain,out,Nw,Nc,burnin,lab)
res = func.MCMC_analysis(chain,out,Nw,Nc,burnin)
#res = np.loadtxt('Results_CBD1/MCMC_analysis.dat')
p = res[:,0]
Kv = 10**p[1:1+params.n_model]
Kd = np.zeros(params.n_model)
Kd[0] = 1.0/Kv[0]
print('Kd1=',Kd[0])
for i in range(1,params.n_model):
  Kd[i] = 1.0/(Kv[i]*np.prod(Kd[:i]))
  print('Kd'+str(i+1)+'=',Kd[i])
func.Figure_fit(p,out,Pmatrix,Xmatrix,DHmatrix,Vc,vi,So,params.color,params.n_model)
