from matplotlib.colors import LogNorm
from matplotlib.pyplot import *
import numpy as np
import scipy, h5py
import os,sys
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d 
import seaborn as sns
import numpy as np
import pandas as pd


recon_total = np.empty((0,6))
truth_total = np.empty((0,4))


filepath = r'./output/'
files = os.listdir(filepath)
for i in files:                             # 循环读取路径下的文件并筛选输出
    if os.path.splitext(i)[1] == ".h5":   # 筛选csv文件
        print(i)  
        h = h5py.File(os.path.join(filepath,i),'r')
        result_recon = h['recon'][...]
        result_truth = h['truth'][...]
        recon_total = np.vstack((recon_total, result_recon))
        truth_total = np.vstack((truth_total, result_truth))

# Energy correct
recon_total[:,0] = 1/np.mean(recon_total[:,0])*np.mean(truth_total[:,0])*recon_total[:,0]
tau_real = 26*(np.log10(np.floor(10*truth_total[:,0])/10)*0.3+1)

with open('Recon.txt','w') as f:
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

fig_path = './figure'
if not os.path.exists(fig_path):
    os.mkdir(fig_path)

# fig_0 E and tau_d
fig=plt.figure(num=0)
plt.hist2d(recon_total[:, 0], recon_total[:, 5],bins=100)
plt.colorbar()
plt.xlabel('Energy: [MeV]')
plt.ylabel(r'$\tau_d$: [ns]')
plt.savefig('./figure/E_vs_taud.png')
plt.ion()
plt.pause(5)  #显示秒数
plt.close()

EE = np.arange(0.2, 2.3, 0.1)
dis = np.sqrt(np.sum(truth_total[:,1:4]**2,1))
xx = np.max(dis)/10**(1/3)*(np.arange(0,10,1)**(1/3))

E_bias = np.zeros((len(EE)-1, len(xx)-1))
E_res = np.zeros((len(EE)-1, len(xx)-1))
x_bias = np.zeros((len(EE)-1, len(xx)-1))
x_res = np.zeros((len(EE)-1, len(xx)-1))
tau_bias = np.zeros((len(EE)-1, len(xx)-1))
tau_res = np.zeros((len(EE)-1, len(xx)-1))

for i in range(len(EE) - 1):
    for j in range(len(xx) - 1):
        index = (truth_total[:,0]>EE[i]) & (truth_total[:,0] < EE[i+1]) \
            & (dis > xx[j]) & (dis > xx[j+1])
        E_bias[i,j] = np.mean(np.abs(truth_total[index,0] - recon_total[index,0]))/np.sqrt((EE[i]+EE[i+1])/2)
        E_res[i,j] = np.sqrt(np.var(truth_total[index,0] - recon_total[index,0]))/np.sqrt((EE[i]+EE[i+1])/2)
        x_bias[i,j] = np.mean(np.abs(truth_total[index,1] - recon_total[index,1]))/np.sqrt((EE[i]+EE[i+1])/2)
        x_res[i,j] = np.sqrt(np.var(truth_total[index,1] - recon_total[index,1]))/np.sqrt((EE[i]+EE[i+1])/2)
        tau_bias[i,j] = np.mean(np.abs(tau_real[index] - recon_total[index,5]))/np.sqrt((EE[i]+EE[i+1])/2)
        tau_res[i,j] = np.sqrt(np.var(tau_real[index] - recon_total[index,5]))/np.sqrt((EE[i]+EE[i+1])/2)
E_bias[np.isnan(E_bias)] = 0
E_res[np.isnan(E_res)] = 0

fig=plt.figure(num=1)
sns.heatmap(pd.DataFrame(E_bias), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0)
plt.xlabel('radius: [m]')
plt.ylabel('energy: [MeV]')
plt.title(r'Energy bias/$\sqrt{E}$')
plt.savefig('./figure/Energy_bias.png')

fig=plt.figure(num=2)
sns.heatmap(pd.DataFrame(E_res), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0)
plt.xlabel('radius: [m]')
plt.ylabel('energy: [MeV]')
plt.title(r'Energy resolution/$\sqrt{E}$')
plt.savefig('./figure/Energy_res.png')

fig=plt.figure(num=3)
sns.heatmap(pd.DataFrame(x_bias), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0)
plt.xlabel('radius: [m]')
plt.ylabel('energy: [MeV]')
plt.title(r'$x$ bias: [m]')
plt.savefig('./figure/x_bias.png')

fig=plt.figure(num=4)
sns.heatmap(pd.DataFrame(x_res), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0)
plt.xlabel('radius: [m]')
plt.ylabel('energy: [MeV]')
plt.title(r'$x$ resolution: [m]')
plt.savefig('./figure/x_res.png')

fig=plt.figure(num=5)
sns.heatmap(pd.DataFrame(tau_bias), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0)
plt.xlabel('radius: [m]')
plt.ylabel('energy: [MeV]')
plt.title(r'$\tau_d$ bias: [ns]')
plt.savefig('./figure/tau_bias.png')

fig=plt.figure(num=6)
sns.heatmap(pd.DataFrame(tau_res), xticklabels= True, yticklabels= True, 
            square=True, cmap="YlGnBu")
plt.xticks(np.arange(len(xx)),np.round(xx,1),color='black',rotation=60)
plt.yticks(np.arange(len(EE)),np.round(EE,1),color='black',rotation=0)
plt.xlabel('radius: [m]')
plt.ylabel('energy: [MeV]')
plt.title(r'$\tau_d$ resolution: [ns]')
plt.savefig('./figure/tau_res.png')