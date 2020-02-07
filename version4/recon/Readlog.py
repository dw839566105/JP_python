from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import os
import tables

def coeff3d():
    logs = ['0.8','1.0','1.2','1.5','1.8','2.0']
    radius = np.arange(-0.65,0.66,0.01)
    new_cell = np.zeros([len(radius), 7, len(logs)])
    print(new_cell.shape)
    for i in range(len(logs)):
        filename = './input/logs/c' + logs[i] + '.h5'
        print(filename)
        h = tables.open_file(filename,'r')
        recondata = h.root.coeff
        x_axis = h.root.x
        coeff_x = np.array(recondata[:])
        x = np.array(x_axis[:])
        for j in np.arange(len(coeff_x[0,])):
            if i == 0:
                coeff_x[0,0] = 0
            x_left1 = np.min(x)-0.01
            x_left2 = np.min(x)-0.02
            x_left3 = -1
            
            x_right1 = np.max(x)+0.01
            x_right2 = np.max(x)+0.02
            x_right3 = 1
            xx = np.hstack((x_left3, x_left2, x_left1, x[:,0], x_right1, x_right2, x_right3))
            yy = np.hstack((0, 0, 0, coeff_x[:,j], 0, 0, 0))
            f = interpolate.interp1d(xx, yy, kind='cubic')
            
            new_cell[:,j,i] = f(radius)
    return logs, radius, new_cell

if __name__ == '__main__':
    logs, radius, new_cell = coeff3d()
    print(type(logs))
    print(type(np.array(logs)[0]))
    for i in range(len(new_cell[0,:,0])):
        plt.figure(i)
        plt.plot(radius, new_cell[:,i,:])
        
        plt.xlabel('radius/m')
        plt.ylabel('coeff value')
        plt.legend(['0.8','1.0','1.2','1.5','1.8','2.0'])
        plt.title((str(i) + '-th value'))
        plt.savefig(('./fig'+str(i)+'.jpg'))