import numpy as np
import scipy, h5py
import scipy.stats as stats
import math
import ROOT
import os,sys
import tables
import uproot, argparse
import scipy.io as scio
from scipy.optimize import minimize
from scipy import interpolate
from numpy.polynomial import legendre as LG
from scipy import special
from Readlog import coeff3d

# physical constant
Light_yield = 4285*0.88 # light yield
Att_LS = 18 # attenuation length of LS
Att_Wtr = 300 # attenuation length of water
tau_r = 1.6 # fast time constant
TTS = 5.5/2.355
QE = 0.20
PMT_radius = 0.254
c = 2.99792e8
n = 1.48
shell = 0.65 # Acrylic

# coeff3d
EE_tmp, radius, coeff = coeff3d()
EE = np.zeros(len(EE_tmp))
for i in np.arange(len(EE)):
    EE[i] = int(EE_tmp[i])

def Likelihood_Sph(vertex, *args):
    coeff, PMT_pos, event_pe, cut = args
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
    for i in np.arange(cut):
        # polyfit
        # fitfun = np.poly1d(coeff[:,i])
        # k[i] = fitfun(z)
        # cubic interp
        xx = radius
        yy = coeff[:,i,:]
        f = interpolate.interp2d(xx, EE, yy, kind='cubic')
        k[i] = f(vertex[0],z)
    # print(k) 
    # print('haha')
    expect = np.exp(np.dot(x,k))
    L = - np.sum(np.sum(np.log((expect**y)*np.exp(-expect))))
    return L

def Likelihood_ML(fit, *args):
    Energy,\
    x,\
    y,\
    z,\
    t,\
    tau_d\
    = fit
    PMT_pos, pe_array, time_array, fired_PMT = args
    distance, Omega = SolidAngle(x,y,z)
    lmbd = Att(x,y,z)
    # expect photons
    expect = Energy*\
        Light_yield*\
        np.exp(-distance*lmbd/Att_LS - distance*(1-lmbd)/Att_Wtr)*\
        Omega*\
        QE
    # log Poisson # p_pe = - np.log(stats.poisson.pmf(PE, expect))
    log_p_pe = - expect + pe_array*np.log(expect) 
    # this part is nonsense {- np.log(special.factorial(pe_array))}
    Likelihood_pe = - np.nansum(log_p_pe)
    # log Time profile pdf
    # log_p_time = TimeProfile(time_array, distance[fired_PMT], tau_d, t)
    # Likelihood_time = - np.nansum(log_p_time)
    # total likelihood
    Likelihood_total = Likelihood_pe
    #Likelihood_total = Likelihood_pe + Likelihood_time
    return Likelihood_total

def SolidAngle(x, y, z):
    distance = np.sqrt(np.sum((PMT_pos - np.array((x,y,z)))**2, axis=1))
    radius_O1 = PMT_radius # PMT bottom surface
    PMT_vector = - PMT_pos/np.transpose(np.tile(np.sqrt(np.sum(PMT_pos**2,1)),[3,1]))
    O1 = np.tile(np.array([x,y,z]),[len(PMT_pos[:,0]),1])
    O2 = PMT_pos
    flight_vector = O2 - O1
    d2 = np.sqrt(np.sum(flight_vector**2,1))
    theta1 = np.sum(PMT_vector*flight_vector,1)/np.sqrt(np.sum(PMT_vector**2,1)*np.sum(flight_vector**2,1))
    Omega = (1-d2/np.sqrt(d2**2+radius_O1*np.abs(theta1)))/2
    
    return distance, Omega
'''
def SolidAngle(x, y, z):
    distance = np.sqrt(np.sum((PMT_pos - np.array((x,y,z)))**2, axis=1))
    radius_O1 = PMT_radius # PMT bottom surface
    radius_O2 = 0.315 # PMT sphere surface
    PMT_vector = - PMT_pos/np.transpose(np.tile(np.sqrt(np.sum(PMT_pos**2,1)),[3,1]))
    O1 = np.tile(np.array([x,y,z]),[len(PMT_pos[:,0]),1])
    O2 = PMT_pos
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
    d = np.transpose(d)
    chord = 2*np.sqrt(radius_O2**2 - d[d<radius_O2]**2)

    theta1 = np.sum(PMT_vector*flight_vector,1)/np.sqrt(np.sum(PMT_vector**2,1)*np.sum(flight_vector**2,1))
    add_area = 1/3* \
        (radius_O2 - radius_O1*np.abs(theta1[d<radius_O2]) \
        - np.sqrt(radius_O2**2 - radius_O1**2)*np.abs(np.sin(np.arccos(theta1[d<radius_O2]))))*chord
    Omega = (1-d2/np.sqrt(d2**2+radius_O1*np.abs(theta1)))/2
    Omega[d<radius_O2] = Omega[d<radius_O2] + add_area/(4*d2[d<radius_O2]**2)
    return distance, Omega
'''

def Att(x, y, z):
    '''
    this function returns ratio in different material 
    lmbd is in the LS and 1-lmbda is the water
    '''
    # LS distance
    d1 = np.tile(np.array([x,y,z]),[len(PMT_pos[:,1]),1])
    d2 = PMT_pos
    d3 = d2 - d1
    # cons beyond shell 
    lmbd = (-2*np.sum(d3*d1,1) \
        + np.sqrt(4*np.sum(d3*d1,1)**2 \
        - 4*np.sum(d3**2,1)*(-np.abs((np.sum(d1**2,1)-shell**2))))) \
        /(2*np.sum(d3**2,1))
    lmbd[lmbd>=1] = 1
    return lmbd

def TimeProfile(time_array, distance, tau_d, t):
    time_correct = time_array - distance/(c/n)*1e9 - t
    time_correct[time_correct<=-8] = -8
    p_time = TimeUncertainty(time_correct, tau_d)
    return p_time

def TimeUncertainty(tc, tau_d):
    a1 = np.exp(((TTS**2 - tc*tau_d)**2-tc**2*tau_d**2)/(2*TTS**2*tau_d**2))
    a2 = np.exp(((TTS**2*(tau_d+tau_r) - tc*tau_d*tau_r)**2 - tc**2*tau_d**2*tau_r**2)/(2*TTS**2*tau_d**2*tau_r**2))
    a3 = np.exp(((TTS**2 - tc*tau_d)**2 - tc**2*tau_d**2)/(2*TTS**2*tau_d**2))*special.erf((tc*tau_d-TTS**2)/(np.sqrt(2)*tau_d*TTS))
    a4 = np.exp(((TTS**2*(tau_d+tau_r) - tc*tau_d*tau_r)**2 - tc**2*tau_d**2*tau_r**2)/(2*TTS**2*tau_d**2*tau_r**2))*special.erf((tc*tau_d*tau_r-TTS**2*(tau_d+tau_r))/(np.sqrt(2)*tau_d*tau_r*TTS))
    p_time  = np.log(tau_d + tau_r) - 2*np.log(tau_d) + np.log(a1-a2+a3-a4)
    
    return p_time

def con(args):
    E_min,\
    E_max,\
    tau_min,\
    tau_max,\
    t0_min,\
    t0_max\
    = args
    cons = ({'type': 'ineq', 'fun': lambda x: (x[0] - E_min)*(E_max - x[0])},\
    {'type': 'ineq', 'fun': lambda x: shell**2 - (x[1]**2 + x[2]**2 + x[3]**2)},\
    {'type': 'ineq', 'fun': lambda x: (x[5] - tau_min)*(tau_max-x[5])},\
    {'type': 'ineq', 'fun': lambda x: (x[4] - t0_min)*(t0_max-x[4])})
    return cons

def con_sph(args):
    E_min,\
    E_max,\
    tau_min,\
    tau_max,\
    t0_min,\
    t0_max\
    = args
    cons = ({'type': 'ineq', 'fun': lambda x: (x[0] - E_min)*(E_max - x[0])},\
    {'type': 'ineq', 'fun': lambda x: shell**2 - (x[1]**2 + x[2]**2 + x[3]**2)})
    return cons

def recon_drc(time_array, fired_PMT, recon_vertex):
    time_corr = time_array - np.sum(PMT_pos[fired_PMT,1:4]-np.tile(recon_vertex[0,1:4],[len(fired_PMT),1]))/(3*10**8)
    index = np.argsort(time_corr)
    fired_PMT_sorted = fired_PMT[index]
    fired_PMT_sorted = fired_PMT_sorted[0:int(np.floor(len(fired_PMT_sorted)/10))]
    drc = np.sum(PMT_pos[fired_PMT_sorted,1:4],0)/len(fired_PMT_sorted)
    return drc

def ReadPMT():
    f = open(r"./PMT1t.txt")
    line = f.readline()
    data_list = [] 
    while line:
        num = list(map(float,line.split()))
        data_list.append(num)
        line = f.readline()
    f.close()
    PMT_pos = np.array(data_list)
    return PMT_pos

def recon(fid, fout, *args):
    PMT_pos, event_count = args
    # global event_count,shell,PE,time_array,PMT_pos, fired_PMT
    '''
    reconstruction

    fid: root reference file
    fout: output file
    '''
    # Create the output file and the group
    rootfile = ROOT.TFile(fid)
    print(fid) # filename
    class ReconData(tables.IsDescription):
        EventID = tables.Int64Col(pos=0)    # EventNo
        x = tables.Float16Col(pos=1)        # x position
        y = tables.Float16Col(pos=2)        # y position
        z = tables.Float16Col(pos=3)        # z position
        t0 = tables.Float16Col(pos=4)       # time offset
        E = tables.Float16Col(pos=5)        # energy
        tau_d = tables.Float16Col(pos=6)    # decay time constant
        success = tables.Int64Col(pos=7)    # recon failure
        x_sph = tables.Float16Col(pos=8)        # x position
        y_sph = tables.Float16Col(pos=9)        # y position
        z_sph = tables.Float16Col(pos=10)        # z position
        E_sph = tables.Float16Col(pos=11)        # energy
        success_sph = tables.Int64Col(pos=12)    # recon failure
    # Create the output file and the group
    h5file = tables.open_file(fout, mode="w", title="OneTonDetector",
                            filters = tables.Filters(complevel=9))
    group = "/"
    # Create tables
    ReconTable = h5file.create_table(group, "Recon", ReconData, "Recon")
    recondata = ReconTable.row
    # Loop for event
    f = uproot.open(fid)
    a = f['SimpleAnalysis']
    for tot, chl, PEl, Pkl, nPl in zip(a.array("TotalPE"),  # total pe in an event
                    a.array("ChannelInfo.ChannelId"),       # PMT fired seq
                    a.array('ChannelInfo.PE'),              # Hit info number on PMT
                    a.array('ChannelInfo.PeakLoc'),         # Time info on PMT
                    a.array('ChannelInfo.nPeaks')):         # 
        pe_array = np.zeros(np.size(PMT_pos[:,1])) # Photons on each PMT (PMT size * 1 vector)
        fired_PMT = np.zeros(0)     # Hit PMT (PMT Seq can be repeated)
        time_array = np.zeros(0, dtype=int)    # Time info (Hit number)
        for ch, pe, pk, npk in zip(chl, PEl, Pkl, nPl):
            pe_array[ch] = pe
            time_array = np.hstack((time_array, pk))
            fired_PMT = np.hstack((fired_PMT, ch*np.ones(np.size(pk))))
        fired_PMT = fired_PMT.astype(int)
        # initial result
        result_vertex = np.empty((0,6)) # reconstructed vertex
        # initial value x[0] = [1,6]
        x0 = np.zeros((1,6))
        x0[0][0] = pe_array.sum()/300
        x0[0][1] = np.sum(pe_array*PMT_pos[:,0])/np.sum(pe_array)
        x0[0][2] = np.sum(pe_array*PMT_pos[:,1])/np.sum(pe_array)
        x0[0][3] = np.sum(pe_array*PMT_pos[:,2])/np.sum(pe_array)
        x0[0][4] = np.mean(time_array)
        x0[0][5] = 26
        # Constraints
        E_min = 0.01
        E_max = 100
        tau_min = 0.01
        tau_max = 100
        t0_min = -300
        t0_max = 300
        con_args = E_min, E_max, tau_min, tau_max, t0_min, t0_max
        cons = con(con_args)
        # reconstruction
        result = minimize(Likelihood_ML, x0, method='SLSQP', constraints=cons, \
        args = (PMT_pos, pe_array, time_array, fired_PMT))
        # result
        print(event_count, result.x, result.success)
        event_count = event_count + 1
        recondata['EventID'] = event_count
        recondata['x'] = result.x[1]
        recondata['y'] = result.x[2]
        recondata['z'] = result.x[3]
        recondata['E'] = result.x[0]
        recondata['t0'] = result.x[4]
        recondata['tau_d'] = result.x[5]
        recondata['success'] = result.success

        h = h5py.File('../calib/coeff.h5','r')
        coeff = h['coeff'][...]
        # initial value
        x0 = np.zeros((1,4))
        x0[0][0] = pe_array.sum()/300
        x0[0][1] = np.sum(pe_array*PMT_pos[:,0])/np.sum(pe_array)
        x0[0][2] = np.sum(pe_array*PMT_pos[:,1])/np.sum(pe_array)
        x0[0][3] = np.sum(pe_array*PMT_pos[:,2])/np.sum(pe_array)

        # Constraints
        # x0 = np.sum(PE*PMT_pos,axis=0)/np.sum(PE)
        theta0 = np.array([1,0.1,0.1,0.1])
        theta0[0] = x0[0][0]
        theta0[1] = x0[0][1]
        theta0[2] = x0[0][2]
        theta0[3] = x0[0][3]
        con_args = E_min, E_max, tau_min, tau_max, t0_min, t0_max
        cons_sph = con_sph(con_args)
        record = np.zeros((1,4))
        result = minimize(Likelihood_Sph, theta0, method='SLSQP',constraints=cons_sph, args = (coeff, PMT_pos, pe_array, cut))
        # record[0,:] = np.array(result.x, dtype=float)
        # result_total = np.vstack((result_total,record))

        # result
        print(event_count, result.x, result.success)
        recondata['x_sph'] = result.x[1]
        recondata['y_sph'] = result.x[2]
        recondata['z_sph'] = result.x[3]
        recondata['E_sph'] = result.x[0]
        recondata['success_sph'] = result.success
        recondata.append()

    # Flush into the output file
    ReconTable.flush()
    h5file.close()

# Automatically add multiple root files created a program with max tree size limitation.
if len(sys.argv)!=3:
    print("Wront arguments!")
    print("Usage: python Recon.py MCFileName[.root] outputFileName[.h5]")
    sys.exit(1)
# Read PMT position
PMT_pos = ReadPMT()
event_count = 0
cut = 7
ROOT.PyConfig.IgnoreCommandLineOptions = True
# Reconstruction
fid = sys.argv[1] # input file .root
fout = sys.argv[2] # output file .h5
args = PMT_pos, event_count
recon(fid, fout, *args)
