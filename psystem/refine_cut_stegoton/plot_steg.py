import numpy as np
import h5py
from matplotlib import pyplot as plt

hf = h5py.File('refn0_small_domain.h5', 'r')
#hf = h5py.File('refn0_full_domain.h5', 'r')

x = hf.get('x')
sigma = hf.get('stress')
eps = hf.get('strain')
vel = hf.get('vel')

plt.plot(x,eps,'-r',lw=2)
plt.plot(x,sigma,'-k',lw=3)
plt.savefig('steg_refn0.png')


