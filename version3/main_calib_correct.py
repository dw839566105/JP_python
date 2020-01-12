import numpy as np 
import scipy, h5py
import tables
import sys
from scipy.optimize import minimize
from numpy.polynomial import legendre as LG
import matplotlib.pyplot as plt

def Recon():
    fun = lambda theta: calib(theta)
    return fun

def calib(theta):
    global total_pe, PMT_pos, cut, mu0
    y = total_pe
    # fixed axis
    vertex = np.array([0,0,2,10])
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
    # Poisson regression
    L = - np.sum(np.sum(np.transpose(y)*np.transpose(np.dot(x, theta) - mu0) \
        - np.transpose(np.exp(np.dot(x, theta)))))
    return L

def calib_tmp(theta):
    global total_pe, PMT_pos, cut, mu0
    y = total_pe
    # fixed axis
    vertex = np.array([0,0,2,10])
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
    # Poisson regression
    print((np.dot(x, theta)).shape)
    print(mu0.shape)
    L = - np.sum(np.sum(np.transpose(y)*np.transpose(np.dot(x, theta) - mu0) \
        - np.transpose(np.exp(np.dot(x, theta)))))
    return x

def main_calib(fid, fout):
    global total_pe, PMT_pos, cut
    filepath = '../input/data/type9'
    # k = '%+.2f' % fid
    k = fid
    filename = filepath + '/calib' + k + '.h5'

    # read files by table
    h1 = tables.open_file(filename,'r')
    print(filename)
    truthtable = h1.root.GroundTruth
    EventID = truthtable[:]['EventID']
    ChannelID = truthtable[:]['ChannelID']
    
    # read file series
    try:
        for j in np.arange(1,10,1):
            filename = filepath + '/calib' + k + '_' + str(j)+ '.h5'           
            h1 = h5py.File(filename,'r')
            print(filename)
            truthtable = h1.root.GroundTruth
            EventID_tmp = truthtable[:]['EventID']
            ChannelID_tmp = truthtable[:]['ChannelID']
            EventID = np.hstack((EventID, EventID_tmp))
            ChannelID = np.hstack((ChannelID, ChannelID_tmp))
    except:
        j = j - 1
    total_pe = np.zeros((np.size(PMT_pos[:,0]),max(EventID)))
    for k in np.arange(1, max(EventID)):
        event_pe = np.zeros(np.size(PMT_pos[:,0]))
        hit = ChannelID[EventID == k]
        tabulate = np.bincount(hit)
        event_pe[0:np.size(tabulate)] = tabulate
        total_pe[:,k-1] = event_pe
    theta0 = np.zeros(cut) # initial value
    result = minimize(Recon(),theta0, method='SLSQP')  
    record = np.array(result.x, dtype=float)

    # check
    L1 = np.exp(np.dot(calib_tmp(result.x),result.x) - mu0)        ## generalized Poisson regression
    L2 = np.mean(total_pe, axis=1)
    mu = np.log(L2)-np.log(L1)
    print(np.sum(mu0))
    print(mu)
    with h5py.File(fout,'w') as out:
        out.create_dataset("record", data=record)
    return record, mu

## read data from calib files
global total_pe, PMT_pos, cut, mu0
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
mu0 = np.array([0])
record, mu0 = main_calib('+0.00', sys.argv[2])
#print(record, mu0)
main_calib(sys.argv[1],sys.argv[2])
