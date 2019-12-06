# coding=utf-8

import numpy as np
import scipy, h5py
import scipy.stats as stats
import math
import ROOT
import os,sys
import scipy.io as scio
from scipy.optimize import minimize
from scipy import special

np.set_printoptions(suppress=True)

Light_yield = 4285*0.88 # light yield
Att_Wtr = 300 # attenuation length of water
tau_r = 1.6 # fast time constant
TTS = 5.5/2.355
QE = 0.20
PMT_radius = 0.254
c = 2.99792e8
n = 1.48
count = 0
EE = 1
xx = 0
yy = 0
zz = 0
iterNo = 0
event_count = 0
shell = 12.2 # Acrylic

def Recon():
    fun = lambda x: Likelihood(x)
    return fun

def Likelihood(args):
    Energy,\
    x,\
    y,\
    z,\
    t,\
    tau_d,\
    Att_LS\
    = args
    global shell
    distance = np.sqrt((PMT_pos[:,1] - x)**2 + (PMT_pos[:,2] - y)**2 + (PMT_pos[:,3] - z)**2)
    # LS distance
    d1 = np.tile(np.array([x,y,z]),[len(PMT_pos[:,1]),1])
    d2 = PMT_pos[:,1:4]
    d3 = d2 - d1
    # cons beyond shell ?
    lmbd = (-2*np.sum(d3*d1,1) \
        + np.sqrt(4*np.sum(d3*d1,1)**2 \
        - 4*np.sum(d3**2,1)*(-np.abs((np.sum(d1**2,1)-shell**2))))) \
        /(2*np.sum(d3**2,1))
    expect = Energy*Light_yield*np.exp(-distance*lmbd/Att_LS - distance*(1-lmbd)/Att_Wtr)*SolidAngle(x,y,z,distance)*QE
    p_pe = - expect + pe_array*np.log(expect) - np.log(special.factorial(pe_array));
    # p_pe = - np.log(stats.poisson.pmf(pe_array, expect))
    Likelihood_pe = - np.nansum(p_pe)
    p_time = TimeProfile(time_array, distance[fired_PMT], tau_d, t)
    Likelihood_time = - np.nansum(p_time)
    
    Likelihood_total = Likelihood_pe + Likelihood_time
    return Likelihood_total

def SolidAngle(x, y, z, distance):
    radius_O1 = PMT_radius # PMT bottom surface
    radius_O2 = 0.315 # PMT sphere surface
    PMT_vector = - PMT_pos[:,1:4]/np.transpose(np.tile(np.sqrt(np.sum(PMT_pos[:,1:4]**2,1)),[3,1]))
    O1 = np.tile(np.array([x,y,z]),[len(PMT_pos[:,1]),1])
    O2 = PMT_pos[:,1:4]
    d1 = np.sqrt(radius_O2**2 - radius_O1**2)
    O3 = (O2 + PMT_vector*d1)
    flight_vector = O2 - O1
    d2 = np.sqrt(np.sum(flight_vector**2,1))
    O4 = O2 - flight_vector/ \
        np.transpose(np.tile(d2*np.sqrt(radius_O2**2*(d2**2 - radius_O2**2)/d2**2),[3,1]))
    # Helen formula
    a = np.sqrt(np.sum((O4-O2)**2,1))
    b = np.sqrt(np.sum((O4-O3)**2,1))
    c = np.sqrt(np.sum((O3-O2)**2,1))
    p = (a+b+c)/2
    d = 2*a*b*c/(4*np.sqrt(p*(p-a)*(p-b)*(p-c)))
    '''
    # this part influence the speed!
    data = np.array([a,b,c])
    sorted_cols = []
    for col_no in range(data.shape[1]):
        sorted_cols.append(data[np.argsort(data[:,col_no])][:,col_no])
    sorted_data = np.column_stack(sorted_cols)

    a = sorted_data[0,:]
    b = sorted_data[1,:]
    c = sorted_data[2,:]
    d[a+b-c<=1*10**(-10)] = 1 # avoid inf
    '''
    d = np.transpose(d)
    chord = 2*np.sqrt(radius_O2**2 - d[d<radius_O2]**2)

    theta1 = np.sum(PMT_vector*flight_vector,1)/np.sqrt(np.sum(PMT_vector**2,1)*np.sum(flight_vector**2,1))
    add_area = 1/3* \
        (radius_O2 - radius_O1*np.abs(theta1[d<radius_O2]) \
        - np.sqrt(radius_O2**2 - radius_O1**2)*np.abs(np.sin(np.arccos(theta1[d<radius_O2]))))*chord
    Omega = (1-d2/np.sqrt(d2**2+radius_O1*np.abs(theta1)))/2
    Omega[d<radius_O2] = Omega[d<radius_O2] + add_area/(4*d2[d<radius_O2]**2)
    return Omega

def TimeProfile(time_array, distance, tau_d, t):
    time_correct = time_array - distance/(c/n)*1e9 - t
    time_correct[time_correct<=-8] = -8
    # if(np.sum(np.isnan(time_correct))>=1):
    #    print(t)
    # p_time
    
    p_time = TimeUncertainty(time_correct, tau_d)

    # p_time = np.exp(-time_correct/(tau_d+TTS))*(1-np.exp(-time_correct/(fast_const)))/((tau_d+TTS)**2)*((tau_d+TTS)+fast_const)
    # p_time[p_time<=0] = 1e-200
    return p_time

def TimeUncertainty(tc, tau_d):
    a1 = np.exp(((TTS**2 - tc*tau_d)**2-tc**2*tau_d**2)/(2*TTS**2*tau_d**2))
    a2 = np.exp(((TTS**2*(tau_d+tau_r) - tc*tau_d*tau_r)**2 - tc**2*tau_d**2*tau_r**2)/(2*TTS**2*tau_d**2*tau_r**2))
    a3 = np.exp(((TTS**2 - tc*tau_d)**2 - tc**2*tau_d**2)/(2*TTS**2*tau_d**2))*special.erf((tc*tau_d-TTS**2)/(np.sqrt(2)*tau_d*TTS))
    a4 = np.exp(((TTS**2*(tau_d+tau_r) - tc*tau_d*tau_r)**2 - tc**2*tau_d**2*tau_r**2)/(2*TTS**2*tau_d**2*tau_r**2))*special.erf((tc*tau_d*tau_r-TTS**2*(tau_d+tau_r))/(np.sqrt(2)*tau_d*tau_r*TTS))
    p_time  = np.log(tau_d + tau_r) - 2*np.log(tau_d) + np.log(a1-a2+a3-a4)
    
    return p_time

def con(args):
    Emin,\
    radius\
    = args

    cons = ({'type': 'ineq', 'fun': lambda x: x[0] - Emin},\
    {'type': 'ineq', 'fun': lambda x: radius**2 - (x[1]**2 + x[2]**2+x[3]**2)},\
    {'type': 'ineq', 'fun': lambda x: x[5] - 5},\
    {'type': 'ineq', 'fun': lambda x: (x[4] + 100)*(300-x[4])},\
    {'type': 'ineq', 'fun': lambda x: x[6] - 5})

    return cons

def recon_drc(time_array, fired_PMT, recon_vertex):
    time_corr = time_array - np.sum(PMT_pos[fired_PMT,1:4]-np.tile(recon_vertex[0,1:4],[len(fired_PMT),1]))/(3*10**8)
    index = np.argsort(time_corr)
    fired_PMT_sorted = fired_PMT[index]
    fired_PMT_sorted = fired_PMT_sorted[0:int(np.floor(len(fired_PMT_sorted)/10))]
    drc = np.sum(PMT_pos[fired_PMT_sorted,1:4],0)/len(fired_PMT_sorted)
    return drc

def recon(fid, fout):
    '''
    reconstruction

    fid: root reference file
    fout: output file in this step
    '''
    global event_count
    
    result_recon = np.empty((0,6))
    result_truth = np.empty((0,4))
    result_drc = np.empty((0,3))
    result_tdrc = np.empty((0,3))
    event_count = 0

    rootfile = ROOT.TFile(fid)
    TruthChain = rootfile.Get('SimTriggerInfo')
    print(fid)

    global pe_array
    global time_array
    global fired_PMT
            
    for event in TruthChain:
        pe_array = np.zeros(PMT_total)
        time_array = np.zeros(TruthChain.PEList.size())
        fired_PMT = np.zeros(TruthChain.PEList.size(), dtype = 'int')
        count = 0
        truth_vertex = np.empty((1,4))
        truth_px = np.empty((1,3))
        for truth in TruthChain.truthList:
            truth_vertex[0][0] = truth.EkMerged
            truth_vertex[0][1] = truth.x/1000
            truth_vertex[0][2] = truth.y/1000
            truth_vertex[0][3] = truth.z/1000
            for A in truth.PrimaryParticleList:
                truth_px[0][0] = A.px
                truth_px[0][1] = A.py
                truth_px[0][2] = A.pz
        for pe in TruthChain.PEList:
            if(pe.PEType != -1):
                time_array[count] = pe.PulseTime
                fired_PMT[count] = int(pe.PMTId-1); # PMTId range 1-8607
                pe_array[pe.PMTId-1] = pe_array[pe.PMTId-1] + 1
                count = count + 1

        time_array = time_array[0:count]
        fired_PMT = fired_PMT[0:count]
        
        x0 = np.zeros((1,7))

        x0[0][0] = pe_array.sum()/300
        x0[0][1] = np.sum(pe_array*PMT_pos[:,1])/np.sum(pe_array)
        x0[0][2] = np.sum(pe_array*PMT_pos[:,2])/np.sum(pe_array)
        x0[0][3] = np.sum(pe_array*PMT_pos[:,3])/np.sum(pe_array)
        x0[0][4] = 95
        x0[0][5] = 26
        x0[0][6] = 18

        Emin = 0.01
        recon_vertex = np.empty((1,7))
        args = (Emin, np.sqrt(np.sum(PMT_pos[1,:]**2)))
        cons = con(args)

        result = minimize(Recon(), x0, method='SLSQP', constraints=cons)

        recon_vertex[0,:] = result.x

        result_truth = np.vstack((result_truth, truth_vertex))
        result_recon = np.vstack((result_recon, recon_vertex))
        drc = recon_drc(time_array, fired_PMT, recon_vertex)
        result_drc = np.vstack((result_drc, drc))
        result_tdrc = np.vstack((result_tdrc, truth_px))
        # print(np.sum(result_drc*truth_px)/np.sqrt(np.sum(result_drc**2)*np.sum(truth_px**2)))
        '''
        EE = recon_vertex[0,0]
        xx = recon_vertex[0,1]
        yy = recon_vertex[0,2]
        zz = recon_vertex[0,3]
        t0 = recon_vertex[0,4]
        taud = recon_vertex[0,5]
        '''
        event_count = event_count + 1
        print(event_count)
    
    with h5py.File(fout,'w') as out:
        out.create_dataset("truth", data=result_truth)
        out.create_dataset("recon", data=result_recon)
        out.create_dataset("drc", data=result_drc)
        out.create_dataset("tdrc", data=result_tdrc)

fid = sys.argv[2]

# Load PMT positions
f = open(r"../PMT_5kt_sphere.txt")
#f = open(r"./PMT_1t_sphere.txt")

line = f.readline()
data_list = []
while line:
    num = list(map(float,line.split()))
    data_list.append(num)
    line = f.readline()
f.close()

global PMT_pos
global PMT_total
PMT_pos = np.array(data_list)
PMT_pos[:,1:4] = PMT_pos[:,1:4]/1000
PMT_total = len(PMT_pos[:,0])

event_count = 0

recon(fid, sys.argv[1])
