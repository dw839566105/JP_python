import matplotlib.pyplot as plt
import numpy as np

data = np.random.rand(10000,2)
# plt.hist2d(data[:, 0], data[:, 1],bins=100)
# plt.scatter(data[:, 0], data[:, 1], c=data[:,1])
plt.colorbar()
plt.show()