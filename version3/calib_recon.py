import numpy as np 
import h5py
import matplotlib.pyplot as plt
import tables
from numpy.polynomial import legendre as LG
from scipy.optimize import minimize

def Recon():
    fun = lambda vertex: calib(vertex)
    return fun

def calib(vertex):
    global coeff, PMT_pos, event_pe
    cut = 7
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
    # print((coeff.shape))
    #print(np.size(coeff[0,:]))
    for i in np.arange(0,np.size(coeff[0,:])):
        fitfun = np.poly1d(coeff[:,i])
        #print(z)
        k[i] = fitfun(z)
        #print(k)
    expect = np.exp(np.dot(x,k))
    L = - np.sum(np.sum(np.log((expect**y)*np.exp(-expect))))
    return L

global PMT_pos, coeff, event_pe

f = open(r"./calib/PMT1t.txt")
line = f.readline()
data_list = []
while line:
    num = list(map(float,line.split()))
    data_list.append(num)
    line = f.readline()
f.close()
PMT_pos = np.array(data_list)

filepath = './calib/output/'
radius = np.linspace(-0.63, 0.63, 127)
record = np.zeros((np.size(radius), 7))

h = h5py.File('./coeff.h5','r')
coeff = h['record'][...]

fid = './calib/type9/calib+0.40.h5'

h1 = tables.open_file(fid,'r')
truthtable = h1.root.GroundTruth
EventID = truthtable[:]['EventID']
ChannelID = truthtable[:]['ChannelID']

for k in np.arange(1, max(EventID)):
    event_pe = np.zeros(np.size(PMT_pos[:,0]))
    hit = ChannelID[EventID == k]
    tabulate = np.bincount(hit)
    event_pe[0:np.size(tabulate)] = tabulate
    theta0 = np.array([1,0.1,0.1,0.1])
    result = minimize(Recon(),theta0, method='SLSQP')  
    record = np.array(result.x, dtype=float)
    print(record)