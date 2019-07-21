# coding=utf-8

import numpy as np
import scipy, h5py
import scipy.stats as stats
import math
import ROOT
import os,sys
import scipy.io as scio
from scipy.optimize import minimize

np.set_printoptions(suppress=True)

Light_yield = 4285*0.88 # light yield
Att_length = 18 # attenuation length
fast_const = 1.6 # fast time constant
TTS = 5.5
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
flag_break = 0

def Recon_e():
    fun = lambda x: Likelihood_e(x, tau_d)
    return fun

def Recon_m():
    fun = lambda t: Likelihood_m(t, EE, xx, yy, zz)
    return fun

def Likelihood_e(args, tau_d):
    Energy,\
    x,\
    y,\
    z,\
    t\
    = args

    distance = np.sqrt((PMT_pos[:,1] - x)**2 + (PMT_pos[:,2] - y)**2 + (PMT_pos[:,3] - z)**2)
    expect = Energy*Light_yield*np.exp(-distance/Att_length)*SolidAngle(x,y,z,distance)*QE
    p_pe = - np.log(stats.poisson.pmf(pe_array, expect))
    Likelihood_pe = np.nansum(p_pe)
    p_time = - np.log(TimeProfile(time_array, distance[fired_PMT], tau_d, t))
    Likelihood_time = np.nansum(p_time)
    
    Likelihood = Likelihood_pe + Likelihood_time
    return Likelihood

def Likelihood_m(args, Energy, x, y, z):
    tau_d,\
    t\
    = args

    distance = np.sqrt((PMT_pos[:,1] - x)**2 + (PMT_pos[:,2] - y)**2 + (PMT_pos[:,3] - z)**2)
    expect = Energy*Light_yield*np.exp(-distance/Att_length)*SolidAngle(x,y,z,distance)*QE
    p_pe = - np.log(stats.poisson.pmf(pe_array, expect))
    Likelihood_pe = np.nansum(p_pe)
    p_time = - np.log(TimeProfile(time_array, distance[fired_PMT], tau_d, t))
    Likelihood_time = np.nansum(p_time)
    Likelihood = Likelihood_pe + Likelihood_time
    return Likelihood

def SolidAngle(x, y, z, distance):
    theta = ((x - PMT_pos[:,1]) * PMT_pos[:,1] + (y - PMT_pos[:,2]) * PMT_pos[:,2] + (z - PMT_pos[:,3]) * PMT_pos[:,3]) / \
    np.sqrt((x - PMT_pos[:,1])**2 + (y - PMT_pos[:,2])**2 + (z - PMT_pos[:,3])**2) / \
    np.sqrt(PMT_pos[:,1]**2 + PMT_pos[:,2]**2 + PMT_pos[:,3]**2)
    omega = (1- distance/np.sqrt(distance**2+PMT_radius**2*abs(theta)))/2
    return omega

def TimeProfile(time_array, distance, tau_d, t):
    time_correct = time_array - distance/(c/n)*1e9 - t
    p_time = np.exp(-time_correct/(tau_d+TTS))*(1-np.exp(-time_correct/(fast_const)))/((tau_d+TTS)**2)*((tau_d+TTS)+fast_const)
    p_time[p_time<=0] = 1e-200
    return p_time

def con(args):
    Emin,\
    radius\
    = args

    cons = ({'type': 'ineq', 'fun': lambda x: x[0] - Emin},\
    {'type': 'ineq', 'fun': lambda x: radius**2 - (x[1]**1 + x[2]**2+x[3]**2)})

    return cons

def recon(fid, fout, miter=None):
    '''
    reconstruction

    fid: root reference file
    miter: intermediate files during iteration
    fout: output file in this step
    '''
    global event_count
    global flag_break
    
    result_recon = np.empty((0,5))
    result_truth = np.empty((0,4))
    result_const = np.empty((0,2))

    result_recon = np.empty((0,5))
    result_truth = np.empty((0,4))
    result_const = np.empty((0,2))

    event_count = 0

    rootfile = ROOT.TFile(fid);
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
        for truth in TruthChain.truthList:
            truth_vertex[0][0] = truth.EkMerged
            truth_vertex[0][1] = truth.x/1000
            truth_vertex[0][2] = truth.y/1000
            truth_vertex[0][3] = truth.z/1000

        for pe in TruthChain.PEList:
            if(pe.PEType != -1):
                time_array[count] = pe.PulseTime
                fired_PMT[count] = int(pe.PMTId-1); # PMTId range 1-8607
                pe_array[pe.PMTId-1] = pe_array[pe.PMTId-1] + 1
                count = count + 1

        time_array = time_array[0:count]
        fired_PMT = fired_PMT[0:count]
        global EE, xx, yy, zz, tau_d
        x0 = np.zeros((1,5))
        if(iterNo == 0):
            x0[0][0] = pe_array.sum()/300;
            x0[0][1] = np.sum(pe_array*PMT_pos[:,1])/np.sum(pe_array)
            x0[0][2] = np.sum(pe_array*PMT_pos[:,2])/np.sum(pe_array)
            x0[0][3] = np.sum(pe_array*PMT_pos[:,3])/np.sum(pe_array)
            x0[0][4] = 95
            tau_d = 25
        else:
            with h5py.File(miter) as ipt:
                result_recon_old = ipt['recon'][...]
                result_const_old = ipt["const"][...]
            x0 = result_recon_old[event_count,:]
            tau_d = result_const_old[event_count, 0]
            x0[4] = result_const_old[event_count, 1]

        Emin = 0.01
        recon_vertex = np.empty((1,5))
        recon_const = np.empty((1,2))
        args = (Emin, np.sqrt(np.sum(PMT_pos[1,:]**2)))
        cons = con(args)

        result = minimize(Recon_e(), x0, method='SLSQP',constraints=cons)

        recon_vertex[0,:] = result.x
        Likelihood_e = result.fun

        result_truth = np.vstack((result_truth, truth_vertex))
        result_recon = np.vstack((result_recon, recon_vertex))

        EE = recon_vertex[0,0]
        xx = recon_vertex[0,1]
        yy = recon_vertex[0,2]
        zz = recon_vertex[0,3]
        t0 = [tau_d, recon_vertex[0,4]];

        result = minimize(Recon_m(), t0, method='SLSQP')

        recon_const[0,:] = result.x
        vertex_tmp = recon_vertex[0,1:3]
        Likelihood_m = result.fun
        result_const = np.vstack((result_const, recon_const))
        tau_d = recon_const[0,0]

        event_count = event_count + 1
        print(event_count)

    with h5py.File(fout,'w') as out:
        out.create_dataset("truth", data=result_truth)
        out.create_dataset("recon", data=result_recon)
        out.create_dataset("const", data=result_const)

iterNo = int(sys.argv[1])
fid = sys.argv[3]

# Load PMT positions
f = open(r"./PMT_5kt_sphere.txt")
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
PMT_pos[:,1:4] = PMT_pos[:,1:4]/1000 # why is 4???
PMT_total = len(PMT_pos[:,0])

event_count = 0

if len(sys.argv) > 4:
    recon(fid, sys.argv[2], sys.argv[4])
else:
    # step 0
    recon(fid, sys.argv[2])
