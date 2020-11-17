import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import emcee
import corner
from random import randint

errSo = 0.1

def tools_fit(So,n_model,minK,maxK,minDH,maxDH,minq0,maxq0):
  nconc = len(So)
  ndim = nconc*4+2*n_model
  bo = []
  lab = []
  for i in range(n_model):
    bo.append((minDH,maxDH))
    lab.append(r'$\Delta H^{('+str(i+1)+r')}_{\rm bind}$')
  for i in range(n_model):
    bo.append((minK,maxK))
    lab.append(r'${\rm log}\beta^{('+str(i+1)+')}$')
  for i in range(nconc):  
    bo.append((0.000001,1.0))
    lab.append(r'$\sigma_{'+str(i)+'}$')
  for i in range(nconc):
    bo.append((So[i]*0.01,So[i]*3.0))
    lab.append(r'$[L]_0^{('+str(i)+')}$')
  for i in range(nconc):
    bo.append((0.0001,2.0))
    lab.append(r'$\epsilon_{'+str(i)+'}$')
  for i in range(nconc):
    bo.append((minq0,maxq0))
    lab.append(r'$\Delta H_0^{('+str(i)+')}$')
  return bo,ndim,lab

def ln_Gauss(x,mu,sigma):
  return -np.log(sigma)-0.5*np.log(2*np.pi)-((x-mu)**2/(2.0*sigma**2))

def Xfree(Xtot,Kv,Ptot):
  n = len(Kv)
  def pol(x):
    out = x*(Kv[0]*(Ptot-Xtot)+1.0)
    for i in range(1,n):
      out = out + x**(i+1)*(Kv[i]*((i+1)*Ptot-Xtot)+Kv[i-1])
    out = out + x**(n+1)*Kv[n-1] - Xtot
    return out
  try:
    res = op.brentq(pol,0,Xtot,maxiter=500)
    if 0<res<Xtot:
      return res
    else:
      return Xtot
  except:
    return Xtot

def Z(x,Kv):
  n = len(Kv)
  out = 1.0
  for i in range(n):
    out = out + x**(i+1)*Kv[i]
  return out

def DHm(Xtot,Kv,Ptot,H):
  xf =  Xfree(Xtot,Kv,Ptot)
  z = Z(xf,Kv)
  n = len(Kv)
  out = 0
  for i in range(n):
    out = out + H[i]*xf**(i+1)*Kv[i]
  return out*(1.0/z)

def Best_fit(Pmatrix,Xmatrix,DHmatrix,Vc,vi,So,bo,dim,nw,nc,lab,nmodel,eArgs):
  nconc = len(So)
  DH_exp = DHmatrix[1:,:]
  Ndata = len(DH_exp)
  def Nlogprob(p):

    DH = p[0:nmodel]
    logKv = p[nmodel:2*nmodel]
    Kv = 10**logKv
    s = p[2*nmodel:2*nmodel+nconc]
    CS = p[2*nmodel+nconc:2*nmodel+2*nconc]
    eps = p[2*nmodel+2*nconc:2*nmodel+3*nconc]
    q0 = p[2*nmodel+3*nconc:2*nmodel+4*nconc]
    Ptot = np.zeros(Pmatrix.shape)
    DHm_th = np.zeros(Xmatrix.shape)
    DH_th = np.zeros(DH_exp.shape) 

    NProb = 0
    for i,e in enumerate(eps):
      Ptot[:,i] = Pmatrix[:,i]*e
    for i,Xrow in enumerate(Xmatrix):
      for j,Xtot in enumerate(Xrow):
        DHm_th[i][j] = DHm(Xtot,Kv,Ptot[i][j],DH) 
    if ~np.isnan(np.sum(DHm_th)):
      for i in range(nconc):
        DH_th[:,i] = (1.0/CS[i])*((Ptot[1:,i]*DHm_th[1:,i]*(0.5+(Vc[i]/vi[1:,i])))+(Ptot[:-1,i]*DHm_th[:-1,i]*(0.5-(Vc[i]/vi[:-1,i]))))
        if eArgs[0]:
          NProb = NProb+(Ndata)*np.log(s[i])+(1.0/(2.0*s[i]*s[i]))*np.sum((DH_th[:,i]-DH_exp[:,i]+q0[i])**2)-ln_Gauss(CS[i],So[i],errSo*So[i])-ln_Gauss(eps[i],eArgs[1][i],eArgs[2][i])
        else:
          NProb = NProb+(Ndata)*np.log(s[i])+(1.0/(2.0*s[i]*s[i]))*np.sum((DH_th[:,i]-DH_exp[:,i]+q0[i])**2)-ln_Gauss(CS[i],So[i],errSo*So[i])
      return NProb
    else:
      return np.inf
  bounds = bo
#  result = op.differential_evolution(Nlogprob, bounds, strategy='best1bin',polish=True)
  tmp = False
  k=0
  while tmp==False:
    try:
      result = op.differential_evolution(Nlogprob, bounds, strategy='best1bin',polish=True)
      tmp = True
    except:
      print('Fail... Trying again!!')
      tmp = False
      if k==99:
        print('Verify your input files')
        break
    k+=1
  print(result)
  pos_i = result.x
  def logprob(p):
    condition = []
    for i in range(len(p)):
      condition.append(bo[i][0]<p[i]<bo[i][1])
    condition = np.array(condition)
    if np.all(condition):
      try:
        return -1.0*Nlogprob(p)
      except:
        return -np.inf
    else:
      return -np.inf
  ndim, nwalkers, npoints = dim, nw, nc
  pos = [pos_i + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
  sampler = emcee.EnsembleSampler(nwalkers,ndim,logprob)
  print('Running MCMC...')
  sampler.run_mcmc(pos,npoints)
  print('Done.')
  samples = sampler.chain[:, :, :].reshape((-1, ndim))
  return pos_i,samples

def Chain_figure(chain,out,Nparam,Nw,Nc,lab):  
  Vchain = np.arange(1,Nc+1,1)
  plt.rc('font', family='serif',size=12)
  plt.figure(figsize=(12,22))
  plt.subplots_adjust(hspace=0.2,wspace=0)
  for i in range(Nparam):
    plt.subplot(Nparam,1,i+1)
    for j in range(Nw):
      plt.plot(Vchain,chain[j*Nc:(j+1)*Nc,i],'k-',c='black',lw=0.5)
    plt.ylabel(lab[i],fontsize=20)
    if i<Nparam-1:
      plt.xticks(np.arange(0,Nc,Nc/(10)))
      plt.tick_params(labelbottom=False)
    else:
      plt.xticks(np.arange(0,Nc+1,Nc/(10)))
      plt.xlabel('Step number', fontsize=18)
  plt.savefig(out+'/Figure_chains.png',bbox_inches='tight')
  plt.close()
  return None

def MCMC_analysis(chain,out,Nw,Nc,burn_in):
  mask = np.ones(Nc*Nw,dtype=bool)
  maskf = np.zeros(burn_in,dtype=bool)
  for i in range(Nw):
    mask[i*Nc:i*Nc+burn_in] = maskf
  samples= chain[mask]
  outdata = np.array(list(map(lambda v: (v[1], v[2]-v[1], v[0]-v[1]),np.percentile(samples,[16,50,84],axis=0).T)))
  np.savetxt(out+'/MCMC_analysis.dat',outdata,fmt='%s')
  return outdata
 
def Correlation_figure(chain,out,Nw,Nc,burn_in,lab):
  mask = np.ones(Nc*Nw,dtype=bool)
  maskf = np.zeros(burn_in,dtype=bool)
  for i in range(Nw):
    mask[i*Nc:i*Nc+burn_in] = maskf
  samples = chain[mask]
  fig = corner.corner(samples, labels=lab, plot_density=False,levels=(1-np.exp(-0.5),1-np.exp(-2.0),),smooth1d=2,smooth=2,bins=20,label_kwargs={"fontsize": 16},fontsize=16)
  fig.savefig(out+'/Correlation.pdf')
  plt.close()
  return None

def Figure_fit(p,out,Pmatrix,Xmatrix,DHmatrix,Vc,vi,So,color_code,nmodel):
  nconc = len(So)
  DH_exp = DHmatrix[1:,:]
  Ndata = len(DH_exp)

  DH = p[0:nmodel]
  logKv = p[nmodel:2*nmodel]
  Kv = 10**logKv
  s = p[2*nmodel:2*nmodel+nconc]
  CS = p[2*nmodel+nconc:2*nmodel+2*nconc]
  eps = p[2*nmodel+2*nconc:2*nmodel+3*nconc]
  q0 = p[2*nmodel+3*nconc:2*nmodel+4*nconc]
  Ptot = np.zeros(Pmatrix.shape)
  DHm_th = np.zeros(Xmatrix.shape)
  DH_th = np.zeros(DH_exp.shape)
  
  for i,e in enumerate(eps):
    Ptot[:,i] = Pmatrix[:,i]*e
  for i,Xrow in enumerate(Xmatrix):
    for j,Xtot in enumerate(Xrow):
      DHm_th[i][j] = DHm(Xtot,Kv,Ptot[i][j],DH)
  for i in range(nconc):
    DH_th[:,i] = (1.0/CS[i])*((Ptot[1:,i]*DHm_th[1:,i]*(0.5+(Vc[i]/vi[1:,i])))+(Ptot[:-1,i]*DHm_th[:-1,i]*(0.5-(Vc[i]/vi[:-1,i]))))+q0[i]
  plt.rc('font', family='serif',size=24)
  fig = plt.figure(figsize=(10,8))
  for i in range(nconc):
    outres = np.zeros((len(Xmatrix[1:,i]),3))
    outres[:,0] = Xmatrix[1:,i]
    outres[:,1] = DH_exp[:,i]
    outres[:,2] = DH_th[:,i]
    np.savetxt(out+'/Profile_ITC'+str(i)+'.dat',outres,fmt='%s')
    plt.scatter(Xmatrix[1:,i],DH_exp[:,i],s=80, facecolors='none', edgecolors=(color_code[i]),alpha=None,label=r'$[P]_{0}=$'+str(round(Pmatrix[:,i][0],0)))
    plt.plot(Xmatrix[1:,i],DH_th[:,i],'-',color=(color_code[i])) 
  #plt.axis([None,None,0,10])
  plt.legend(frameon=False,loc='best',fontsize='xx-small')
  plt.xlabel(r'$[L] \; (\mu{\rm M})$')
  plt.ylabel(r'$\Delta H $ (kcal/mol of injectant)')
  plt.savefig(out+'/ITC_figure.pdf',bbox_inches='tight')
  plt.close()
  plt.clf()
