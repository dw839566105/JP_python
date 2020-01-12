import numpy as np 
import matplotlib.pyplot as plt
from numpy.polynomial import legendre as LG
cut = 7
size = 100
x = np.linspace(-1,1,size)
y = np.zeros((size, cut))
for i in np.arange(0,cut):
    c = np.zeros(cut)
    c[i] = 1
    y[:,i] = LG.legval(x,c)
plt.plot(x,y)
plt.show()
y0 = LG.legval(x, 1)
y1 = LG.legval(x, np.array((0,1)))
y2 = LG.legval(x, np.array((0,0,1)))
y3 = LG.legval(x, np.array((0,0,0,1)))
plt.plot(x,y0)
plt.plot(x,y1)
plt.plot(x,y2)
plt.plot(x,y3)
plt.show()