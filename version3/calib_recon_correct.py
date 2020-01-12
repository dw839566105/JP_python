import numpy as np 
import h5py, sys
import matplotlib.pyplot as plt
import tables
from scipy import interpolate
from numpy.polynomial import legendre as LG
from scipy.optimize import minimize

def Recon():
    fun = lambda vertex: calib(vertex)
    return fun

def calib(vertex):
    global coeff, PMT_pos, event_pe, cut
    y = event_pe
    # fixed axis
    z = np.sqrt(np.sum(vertex[1:4]**2))
    cos_theta = np.sum(vertex[1:4]*PMT_pos,axis=1)\
        /np.sqrt(np.sum(vertex[1:4]**2)*np.sum(PMT_pos**2,axis=1))
    # accurancy and nan value
    cos_theta = np.nan_to_num(cos_theta)
    cos_theta[cos_theta>1] = 1
    cos_theta[cos_theta<-1] =-1
    size = np.size(PMT_pos[:,0])
    x = np.zeros((size, cut))
    # legendre coeff
    for i in np.arange(0,cut):
        c = np.zeros(cut)
        c[i] = 1
        x[:,i] = LG.legval(cos_theta,c)

    k = np.zeros((np.size(coeff[0,:])))
    #print(np.size(coeff[0,:]))
    for i in np.arange(0,np.size(coeff[0,:])):
        # polyfit
        # fitfun = np.poly1d(coeff[:,i])
        # k[i] = fitfun(z)
        # cubic interp

        xx = np.arange(-1,1,0.01)
        yy = np.zeros(np.size(xx))
        # print(np.where(np.abs(x+0.63)<1e-5))
        # print(np.where(np.abs(x-0.64)<1e-5))
        # z = y[np.where(np.abs(x+0.63)<1e-5):np.where(np.abs(x-0.64)<1e-5)]
        # z = y[37:164]
        # print(z.shape)
        # y[np.where(np.where(np.abs(x+0.63)<1e-5)):np.where(np.where(np.abs(x-0.64)<1e-5))] = 
        yy[37:164] = coeff[:,i]
        if z>0.99:
            z = 0.99
        elif z<-0.99:
            z = -0.99
        # print(z)
        f = interpolate.interp1d(xx, yy, kind='cubic')
        k[i] = f(z)
        k[0] = k[0] + np.log(vertex[0])
    # print(k) 
    # print('haha')
    expect = np.exp(np.dot(x,k))
    L = - np.sum(np.sum(np.log((expect**y)*np.exp(-expect))))
    return L

def ReconSph(fid, fout):
    global PMT_pos, coeff, event_pe
    filepath = '../input/data/type9'
    # k = '%+.2f' % fid
    k = fid
    filename = filepath + '/calib' + k + '.h5'

    filepath = './calib/output/'
    radius = np.linspace(-0.63, 0.63, 127)
    result_total = np.empty((1,4))
    record = np.zeros((1,4))
    h = h5py.File('./calib/coeff.h5','r')
    coeff = h['coeff'][...]

    h1 = tables.open_file(filename,'r')
    truthtable = h1.root.GroundTruth
    EventID = truthtable[:]['EventID']
    ChannelID = truthtable[:]['ChannelID']

    for k in np.arange(1, max(EventID)):
        event_pe = np.zeros(np.size(PMT_pos[:,0]))
        pe = np.zeros((np.size(PMT_pos[:,0]),1))
        hit = ChannelID[EventID == k]
        tabulate = np.bincount(hit)
        event_pe[0:np.size(tabulate)] = tabulate
        pe[:,0] = event_pe
        x0 = np.sum(pe*PMT_pos,axis=0)/np.sum(pe)
        theta0 = np.array([1,0.1,0.1,0.1])
        theta0[1:4] = x0
        result = minimize(Recon(),theta0, method='SLSQP')  
        record[0,:] = np.array(result.x, dtype=float)
        result_total = np.vstack((result_total,record))
        print(record)

    with h5py.File(fout,'w') as out:
        out.create_dataset("result", data=result_total)
    

## read data from calib files
global PMT_pos, cut
f = open(r"../input/PMT/PMT1t.txt")
line = f.readline()
data_list = []
while line:
    num = list(map(float,line.split()))
    data_list.append(num)
    line = f.readline()
f.close()
PMT_pos = np.array(data_list)

cut = 7 # Legend order
ReconSph(sys.argv[1],sys.argv[2])