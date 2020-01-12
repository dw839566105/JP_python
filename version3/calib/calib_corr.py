import numpy as np 
import scipy, h5py
import tables
import matplotlib.pyplot as plt

radius = np.linspace(-0.63,0.63,127)
T = np.zeros((0,7))
for index in np.arange(len(radius)):
    filepath = './calib/file_corr'
    k = '%+.2f' % radius[index]
    filename = filepath + k + '.h5'
    # read files by table
    h1 = tables.open_file(filename,'r')
    #print(filename)
    a = h1.root.record
    T = np.vstack((T,a[:]))
with h5py.File('coeff_corr.h5','w') as out:
    out.create_dataset("coeff_corr", data=T)

