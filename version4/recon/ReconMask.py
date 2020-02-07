from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import os
import tables

path = './version4/recon'
path_list=os.listdir(path)
E_total = []
taud_total = []
x_total = []
y_total = []
z_total = []
for filename in path_list:
    if os.path.splitext(filename)[1] == '.h5':
        print(os.path.join(path,filename))
        h = tables.open_file(os.path.join(path,filename),'r')
        recondata = h.root.Recon
        E = recondata[:]['E_sph']
        E_total = np.hstack((E_total, E))
        taud = recondata[:]['tau_d']
        taud_total = np.hstack((taud_total, taud))
        x = recondata[:]['x']
        x_total = np.hstack((x_total, x))
        y = recondata[:]['y']
        y_total = np.hstack((y_total, y))
        z = recondata[:]['z']
        z_total = np.hstack((z_total, z))
        h.close()

x = E_total
x = np.nan_to_num(x)
y = taud_total
y = np.nan_to_num(y)
r = np.sqrt(x_total**2 + y_total**2 + z_total**2)
print(np.max(x))
a = y>5
b = y<100
c = x>1.0
d = x<3.0
e = r<0.6
print(np.sum(e), np.sum(b))
index = b
print(sum(b), max(x),x.shape, r.shape)
# index = y>5 & y<200

plt.figure(1)
for i in np.arange(0.7,0.3,-0.05):
    plt.hist(x[r<i],bins=100)
plt.legend(["0.65","0.60","0.55","0.50","0.45","0.40","0.35","0.30"])
plt.xlabel('Energy/MeV')
plt.ylabel('Count')
plt.savefig('./EnergySpectrum.jpg')

x = x[index]
y = y[index]
r = r[index]
#print(x.shape, y.shape)
#print(np.min(x), np.max(x), np.min(y), np.max(y))

H, xedges, yedges = np.histogram2d(x, y, bins=50)

X, Y = np.meshgrid(xedges[1:],yedges[1:])

plt.figure(2)
plt.contourf(X,Y,np.transpose(H),cmap='Greys')
plt.xlabel('Energy/MeV')
plt.ylabel('Tau_d/ns')
plt.savefig('./EnergyVsTau.jpg')
plt.show()

'''
f = interpolate.interp2d(xedges[1:], yedges[1:], H, kind='cubic')
xnew = np.arange(np.min(xedges), np.max(xedges), 1e-3)
ynew = np.arange(np.min(yedges), np.max(yedges), 1e-2)
znew = f(xnew, ynew)

plt.subplot(122)
plt.contourf(xnew,ynew,znew)
plt.show()
'''