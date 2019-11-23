import numpy as np
import scipy, h5py
import os,sys
import matplotlib.pyplot as plt

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
    f.write('Energy bias: %f MeV\n' %(np.mean(recon_total[:,0]-truth_total[:,0])))
    f.write('Energy resolution: %f MeV\n' %(np.sqrt(np.var(recon_total[:,0]-truth_total[:,0]))))
    f.write('x bias: %f m\n' %(np.mean(recon_total[:,1]-truth_total[:,1])))
    f.write('x resolution: %f m\n' %(np.sqrt(np.var(recon_total[:,1]-truth_total[:,1]))))
    f.write('y bias: %f m\n' %(np.mean(recon_total[:,2]-truth_total[:,2])))
    f.write('y resolution: %f m\n' %(np.sqrt(np.var(recon_total[:,2]-truth_total[:,2]))))
    f.write('z bias: %f m\n' %(np.mean(recon_total[:,3]-truth_total[:,3])))
    f.write('z resolution: %f m\n' %(np.sqrt(np.var(recon_total[:,3]-truth_total[:,3]))))
    f.write('tau bias: %f ns\n' %(np.mean(recon_total[:,5]-tau_real)))
    f.write('tau resolution: %f ns\n' %(np.sqrt(np.var(recon_total[:,3]-tau_real))))
f.close()

fig_path = './figure'
if not os.path.exists(fig_path):
    os.mkdir(fig_path)

fig=plt.figure(num=0)
plt.hist2d(recon_total[:, 0], recon_total[:, 5],bins=100)
plt.colorbar()
plt.xlabel('Energy: [MeV]')
plt.ylabel(r'$\tau_d$: [ns]')
plt.savefig('./figure/E_vs_taud.png')
plt.ion()
plt.pause(5)  #显示秒数
plt.close()
