import numpy as np 
import scipy, h5py
import tables
from scipy.optimize import minimize
from numpy.polynomial import legendre as LG
import matplotlib.pyplot as plt

def Recon():
    fun = lambda theta: calib(theta)
    return fun

def calib(theta):
    global total_pe, PMT_pos, cut
    y = total_pe
    vertex = np.array([0,0,2,10])
    # print(vertex[1:4]*PMT_pos)
    cos_theta = np.sum(vertex[1:4]*PMT_pos,axis=1)\
        /np.sqrt(np.sum(vertex[1:4]**2)*np.sum(PMT_pos**2,axis=1))
    # cos_theta[isnan(cos_theta)) = 0
    # if(max(abs(cos_theta))>1.001)
    #    disp('cos theta error')
    # end
    cos_theta = np.nan_to_num(cos_theta)
    cos_theta[cos_theta>1] = 1
    cos_theta[cos_theta<-1] =-1
    size = np.size(PMT_pos[:,0])
    # x = np.zeros((np.size(PMT_pos[:,0],cut)))
    # x = np.linspace(-1,1,size)
    x = np.zeros((size, cut))
    for i in np.arange(0,cut):
        c = np.zeros(cut)
        c[i] = 1
        x[:,i] = LG.legval(cos_theta,c)
    #print(cos_theta)
    #print(x)

    L1 = - np.sum(np.sum(np.transpose(y)*np.transpose(np.dot(x, theta)) \
        - np.transpose(np.exp(np.dot(x, theta)))))
    #print(L1)
    #L2 = - np.sum(np.sum(np.transpose(y)*np.transpose(np.dot(x, theta)))) \
    #    + np.size(y[0])*np.sum(np.transpose(np.exp(np.dot(x, theta))))

    #print(np.sum(np.sum(np.transpose(y)*np.transpose(np.dot(x, theta)))))
    # print(np.sum(np.transpose(np.exp(np.dot(x, theta)))))
    # print(np.sum(L1-L2))
    return L1

def calib_tmp(theta):
    global total_pe, PMT_pos, cut
    # y = total_pe
    vertex = np.array([0,0,2,10])
    cos_theta = np.sum(vertex[1:4]*PMT_pos,axis=1)/ \
        np.sqrt(np.sum(vertex[1:4]**2)*np.sum(PMT_pos**2,axis=1))
    cos_theta[cos_theta>1] = 1
    cos_theta[cos_theta<-1] =-1
    size = np.size(PMT_pos[:,0])
    x = np.zeros((size, cut))
    for i in np.arange(0,cut):
        c = np.zeros(cut)
        c[i] = 1
        x[:,i] = LG.legval(cos_theta,c)

    #print(np.sum(np.transpose(np.exp(np.dot(x, theta)))))

    return x

## read data from calib files
global total_pe, PMT_pos, cut
f = open(r"./calib/PMT1t.txt")
# f = open(r"./PMT_1t_sphere.txt")


line = f.readline()
data_list = []
while line:
    num = list(map(float,line.split()))
    data_list.append(num)
    line = f.readline()
f.close()

PMT_pos = np.array(data_list)
#PMT_pos[:,1:4] = PMT_pos[:,1:4]/1000
cut = 7
radius = np.linspace(-0.60,0.60,121)
record = np.zeros((cut, np.size(radius)))
for ii in np.arange(0,np.size(radius)):
    i = radius[ii]
    filepath = './calib/type9'
    k = '%+.2f' % i
    # h = h5py.File(os.path.join(filepath,k),'r')
    filename = filepath + '/calib' + k + '.h5'
    print(filename)
    
    h1 = tables.open_file(filename,'r')
    truthtable = h1.root.GroundTruth
    EventID = truthtable[:]['EventID']
    ChannelID = truthtable[:]['ChannelID']
    
    #h = h5py.File(filename,'r')
    #print(np.array(GroundTruth[0]).shape)
    #print((GroundTruth[0][0]))
    #GroundTruth = h['GroundTruth'][...]
    #TruthData = h['TruthData'][...] 
    try:
        for j in np.arange(1,10,1):
            filename = filepath + '/calib' + k + '_' + str(j)+ '.h5'           
            
            h1 = h5py.File(filename,'r')
            print(filename)
            truthtable = h1.root.GroundTruth
            EventID_tmp = truthtable[:]['EventID']
            ChannelID_tmp = truthtable[:]['ChannelID']
            #GroundTruth_tmp = h['GroundTruth'][...] 
            #TruthData_tmp = h['TruthData'][...] 
            EventID = np.hstack((EventID, EventID_tmp))
            ChannelID = np.hstack((ChannelID, ChannelID_tmp))
    except:
        j = j - 1

    total_pe = np.zeros((np.size(PMT_pos[:,0]),max(EventID)))
    # print(GroundTruth[np.size(GroundTruth)-1][0])
    # print((GroundTruth.shape))
    # print((GroundTruth[0].shape))
    # print(GroundTruth[:,0])
    '''
    for k in np.arange(1, min(GroundTruth[-1][0], 500)):  # event number begin from 1
        event_pe = np.zeros(np.size(PMT_pos[:,0]))
        for l in np.arange(0, np.size(GroundTruth)):
            if(GroundTruth[l][0] == k):
                event_pe[GroundTruth[l][1]] = event_pe[GroundTruth[l][1]] + 1
        #print(k)
        total_pe[:,k] = event_pe
    '''
    for k in np.arange(1, max(EventID)):
        event_pe = np.zeros(np.size(PMT_pos[:,0]))
        hit = ChannelID[EventID == k]
        tabulate = np.bincount(hit)
        event_pe[0:np.size(tabulate)] = tabulate
        total_pe[:,k-1] = event_pe

    theta = np.zeros(cut)
    args = (theta)
    print('begin fitting')
    result = minimize(Recon(),args, method='SLSQP')
    #print(result.x.shape)
    #print(type(result.x.shape))
    
    record[:,ii] = np.array(result.x, dtype=float)
    #print(result.x)
    #L1 = np.exp(np.dot(calib_tmp(result.x),result.x))        ## generalized Poisson regression
    #L2 = np.mean(total_pe, axis=1)
    #print(np.mean(L1))
    #print(np.mean(L2))
    #L2 = 1
print(record)
plt.plot(np.transpose(record))
plt.show()