from matplotlib.colors import LogNorm
import numpy as np
import scipy, h5py
import os,sys
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d 
import seaborn as sns
import numpy as np
import pandas as pd
# initial recon and truth
recon_total = np.empty((0,6)) # recon (E,x,y,z,t0,taud)
truth_total = np.empty((0,4)) # truth (E,x,y,z)
# read file with .h5
if len(sys.argv)<2:
    filepath = r'./output' # target dir, can be changed!
    fig_path = './figure'
    log_path = '.'
else:
    filepath = sys.argv[1]
    fig_path = sys.argv[1]+'/figure'
    log_path = filepath
files = os.listdir(filepath)
for i in files:                                     # all files
    if os.path.splitext(i)[1] == ".h5":             # choose .h5 file
        print(i)                                    # print name
        h = h5py.File(os.path.join(filepath,i),'r') # full path
        result_recon = h['recon'][...]              # key 'recon'
        result_truth = h['truth'][...]              # key 'truth'
        recon_total = np.vstack((recon_total, result_recon)) # joint recon
        truth_total = np.vstack((truth_total, result_truth)) # joint truth
# Energy correct
recon_total[:,0] = 1/np.mean(recon_total[:,0])*np.mean(truth_total[:,0])*recon_total[:,0]
# truth taud
tau_real = 26*(np.log10(np.floor(10*truth_total[:,0])/10)*0.3+1)
# output recon result as recon.txt
with open(log_path + '/Recon.txt','w') as f:
    f.write('Energy bias: %f MeV\n' %(np.mean(np.abs(recon_total[:,0]-truth_total[:,0]))))
    f.write('Energy resolution: %f MeV\n' %(np.sqrt(np.var(recon_total[:,0]-truth_total[:,0]))))
    f.write('x bias: %f m\n' %(np.mean(np.abs(recon_total[:,1]-truth_total[:,1]))))
    f.write('x resolution: %f m\n' %(np.sqrt(np.var(recon_total[:,1]-truth_total[:,1]))))
    f.write('y bias: %f m\n' %(np.mean(np.abs(recon_total[:,2]-truth_total[:,2]))))
    f.write('y resolution: %f m\n' %(np.sqrt(np.var(recon_total[:,2]-truth_total[:,2]))))
    f.write('z bias: %f m\n' %(np.mean(np.abs(recon_total[:,3]-truth_total[:,3]))))
    f.write('z resolution: %f m\n' %(np.sqrt(np.var(recon_total[:,3]-truth_total[:,3]))))
    f.write('tau bias: %f ns\n' %(np.mean(np.abs(recon_total[:,5]-tau_real))))
    f.write('tau resolution: %f ns\n' %(np.sqrt(np.var(recon_total[:,3]-tau_real))))
f.close()
# plot figure, path definition move to beginning of program
# fig_path = './figure'
if not os.path.exists(fig_path): # prevent exist
    os.mkdir(fig_path)

# fig_0 E and tau_d
fig=plt.figure(num=0)
plt.hist2d(recon_total[:, 0], recon_total[:, 5],bins=100)
plt.colorbar()
plt.xlabel('Energy: [MeV]', size=20)
plt.ylabel(r'$\tau_d$: [ns]', size=20)
plt.savefig('./figure/E_vs_taud.png')
# if need to be showed
# plt.ion()
# plt.pause(5)  #显示秒数
# plt.close()

dis = np.sqrt(np.sum(truth_total[:,1:4]**2,1)) # vertex to center distance
# plot 2d 
EE = np.arange(0.2, 2.3, 0.1) # energy
xx = np.max(dis)/10**(1/3)*(np.arange(0,10,1)**(1/3)) # distance to center
# relative bias and resolution
E_bias_rlt = np.zeros((len(EE)-1, len(xx)-1))
E_res_rlt = np.zeros((len(EE)-1, len(xx)-1))
x_bias_rlt = np.zeros((len(EE)-1, len(xx)-1))
x_res_rlt = np.zeros((len(EE)-1, len(xx)-1))
tau_bias_rlt = np.zeros((len(EE)-1, len(xx)-1))
tau_res_rlt = np.zeros((len(EE)-1, len(xx)-1))
# absolute bias and resolution
E_bias_abs = np.zeros((len(EE)-1, len(xx)-1))
E_res_abs = np.zeros((len(EE)-1, len(xx)-1))
x_bias_abs = np.zeros((len(EE)-1, len(xx)-1))
x_res_abs = np.zeros((len(EE)-1, len(xx)-1))
tau_bias_abs = np.zeros((len(EE)-1, len(xx)-1))
tau_res_abs = np.zeros((len(EE)-1, len(xx)-1))
# histogram 2d count
for i in range(len(EE) - 1):
    for j in range(len(xx) - 1):
        index = (truth_total[:,0]>EE[i]) & (truth_total[:,0] < EE[i+1]) \
            & (dis > xx[j]) & (dis > xx[j+1])
        E_bias_rlt[i,j] = np.mean(np.abs(truth_total[index,0] - recon_total[index,0]))/np.sqrt((EE[i]+EE[i+1])/2)
        E_res_rlt[i,j] = np.sqrt(np.var(truth_total[index,0] - recon_total[index,0]))/np.sqrt((EE[i]+EE[i+1])/2)
        x_bias_rlt[i,j] = np.mean(np.abs(truth_total[index,1] - recon_total[index,1]))/np.sqrt((EE[i]+EE[i+1])/2)
        x_res_rlt[i,j] = np.sqrt(np.var(truth_total[index,1] - recon_total[index,1]))/np.sqrt((EE[i]+EE[i+1])/2)
        tau_bias_rlt[i,j] = np.mean(np.abs(tau_real[index] - recon_total[index,5]))/np.sqrt((EE[i]+EE[i+1])/2)
        tau_res_rlt[i,j] = np.sqrt(np.var(tau_real[index] - recon_total[index,5]))/np.sqrt((EE[i]+EE[i+1])/2)
        
        E_bias_abs[i,j] = np.mean(np.abs(truth_total[index,0] - recon_total[index,0]))
        E_res_abs[i,j] = np.sqrt(np.var(truth_total[index,0] - recon_total[index,0]))
        x_bias_abs[i,j] = np.mean(np.abs(truth_total[index,1] - recon_total[index,1]))
        x_res_abs[i,j] = np.sqrt(np.var(truth_total[index,1] - recon_total[index,1]))
        tau_bias_abs[i,j] = np.mean(np.abs(tau_real[index] - recon_total[index,5]))
        tau_res_abs[i,j] = np.sqrt(np.var(tau_real[index] - recon_total[index,5]))
# if nan exist
# *_rlt[np.isnan(*_rlt)] = 0
# *_abs[np.isnan(*_abs)] = 0

# fig_1 Energy bias
fig=plt.figure(num=1,figsize=(15,7))
plt.subplot(121)
sns.heatmap(pd.DataFrame(E_bias_rlt), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60, fontsize=15)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0, fontsize=15)
plt.xlabel('Radius: [m]', size=20)
plt.ylabel('Energy: [MeV]', size=20)
plt.title(r'Energy bias/$\sqrt{E}$', size=20)
plt.subplot(122)
sns.heatmap(pd.DataFrame(E_bias_abs), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60, fontsize=15)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0, fontsize=15)
plt.xlabel('Radius: [m]', size=20)
plt.ylabel('Energy: [MeV]', size=20)
plt.title(r'Energy bias: [MeV]', size=20)
plt.savefig(figpath + '/Energy_bias.png')

# fig_2 Energy resolution
fig=plt.figure(num=2,figsize=(15,7))
plt.subplot(121)
sns.heatmap(pd.DataFrame(E_res_rlt), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60, fontsize=15)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0, fontsize=15)
plt.xlabel('Radius: [m]', size=20)
plt.ylabel('Energy: [MeV]', size=20)
plt.title(r'Energy resolution: [MeV]', size=20)
plt.subplot(122)
sns.heatmap(pd.DataFrame(E_res_abs), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60, fontsize=15)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0, fontsize=15)
plt.xlabel('Radius: [m]', size=20)
plt.ylabel('Energy: [MeV]', size=20)
plt.title(r'Energy resolution: [MeV]', size=20)
plt.savefig(fig_path + '/Energy_res.png')

# fig_3 x bias
fig=plt.figure(num=3,figsize=(15,7))
plt.subplot(121)
sns.heatmap(pd.DataFrame(x_bias_rlt), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60, fontsize=15)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0, fontsize=15)
plt.xlabel('Radius: [m]', size=20)
plt.ylabel('Energy: [MeV]', size=20)
plt.title(r'$x$ bias: [m/MeV]', size=20)
plt.subplot(122)
sns.heatmap(pd.DataFrame(x_bias_abs), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60, fontsize=15)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0, fontsize=15)
plt.xlabel('radius: [m]', size=20)
plt.ylabel('energy: [MeV]', size=20)
plt.title(r'$x$ bias: [m]', size=20)
plt.savefig(fig_path + '/x_bias.png')

# fig_4 x resolution
fig=plt.figure(num=4,figsize=(15,7))
plt.subplot(121)
sns.heatmap(pd.DataFrame(x_res_rlt), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60, fontsize=15)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0, fontsize=15)
plt.xlabel('Radius: [m]', size=20)
plt.ylabel('Energy: [MeV]', size=20)
plt.title(r'$x$ resolution: [m/MeV]', size=20)
plt.subplot(122)
sns.heatmap(pd.DataFrame(x_res_abs), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60, fontsize=15)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0, fontsize=15)
plt.xlabel('Radius: [m]', size=20)
plt.ylabel('Energy: [MeV]', size=20)
plt.title(r'$x$ resolution: [m]', size=20)
plt.savefig(fig_path + '/x_res.png')

# fig_5 tau bias
fig=plt.figure(num=5,figsize=(15,7))
plt.subplot(121)
sns.heatmap(pd.DataFrame(tau_bias_rlt), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60, fontsize=15)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0, fontsize=15)
plt.xlabel('Radius: [m]', size=20)
plt.ylabel('Energy: [MeV]', size=20)
plt.title(r'$\tau_d$ bias: [ns/MeV]', size=20)
plt.subplot(122)
sns.heatmap(pd.DataFrame(tau_bias_abs), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60, fontsize=15)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0, fontsize=15)
plt.xlabel('Radius: [m]', size=20)
plt.ylabel('Energy: [MeV]', size=20)
plt.title(r'$\tau_d$ bias: [ns]', size=20)
plt.savefig(fig_path + '/tau_bias.png')

# fig_6 tau resolution
fig=plt.figure(num=6,figsize=(15,7))
plt.subplot(121)
sns.heatmap(pd.DataFrame(tau_res_rlt), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60, fontsize=15)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0, fontsize=15)
plt.xlabel('Radius: [m]', size=20)
plt.ylabel('Energy: [MeV]', size=20)
plt.title(r'$\tau_d$ resolution: [ns/MeV]', size=20)
plt.subplot(122)
sns.heatmap(pd.DataFrame(tau_res_abs), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60, fontsize=15)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0, fontsize=15)
plt.xlabel('Radius: [m]', size=20)
plt.ylabel('Energy: [MeV]', size=20)
plt.title(r'$\tau_d$ resolution: [ns]', size=20)
plt.savefig(fig_path + '/tau_res.png')
