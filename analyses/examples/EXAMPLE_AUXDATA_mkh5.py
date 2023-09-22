#! /usr/bin/env python3

import numpy as np
jpts = np.linspace(20,90,100)
jsfs = 1.0 + 0.1*((1+np.exp(0.1*(60-jpts)))**-1 - 0.5)

import matplotlib.pyplot as plt
plt.plot(jpts, jsfs)
plt.ylim(0.0, 1.1)
plt.show()

import h5py
with h5py.File("EXAMPLE_AUXDATA.h5", "w") as f:
    dset_pts = f.create_dataset("jetpts", (100,), dtype='f', data=jpts)
    dset_sfs = f.create_dataset("jetscales", (100,), dtype='f', data=jpts)
