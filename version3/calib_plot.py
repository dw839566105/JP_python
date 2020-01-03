import numpy as np 
import h5py
import matplotlib.pyplot as plt
filepath = './calib/output/'
radius = np.linspace(-0.63, 0.63, 127)
record = np.zeros((np.size(radius), 7))
for i in np.arange(0, np.size(radius)):
    k = '%+.2f' % radius[i]
    filename = filepath + 'result' + k + '.h5'
    print(filename)
    h = h5py.File(filename,'r')
    record_tmp = h['record'][...]
    #print(record_tmp.shape)
    record[i,:] = record_tmp
fit = np.zeros((11,7))
for i in np.arange(0,7):
    plt.figure(num=i)
    pfit = np.polyfit(radius, record[:,i],10)
    fitfun = np.poly1d(pfit)
    plt.plot(radius,record[:,i], '*')
    plt.plot(radius,fitfun(radius))
    for j in np.arange(0, np.size(pfit)):
        figpar = plt.gcf()
        ax = plt.gca()
        t1 = "%.2f" % pfit[j]
        textstr = "a"+str(j)+":" + t1
        plt.text(0.7,0.95-0.05*j,textstr, size = 15, alpha = 1, transform=ax.transAxes)
        plt.legend(['Legendre coefficient', 'Polynomial fit'])
        fit[:,i] = pfit
with h5py.File('coeff.h5','w') as out:
    out.create_dataset("record", data=fit)
#plt.show()
