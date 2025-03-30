
import numpy as np
import pandas as pd
import time
from ebpca import ebpca_gaussian
from ebpca.pca import get_pca
from ebpca.preprocessing import normalize_obs

n = 50
p = 500
K = 2

X = np.genfromtxt('../data/sim1.csv', delimiter=',')
nseeds = int(X.shape[0]/n)
matr = np.zeros(shape=(nseeds * K, p))

t0 = time.time()
for seed in range(nseeds):
	nrow_from = seed * n
	nrow_to = (seed + 1) * n

	X1 = X[nrow_from:nrow_to,:]
	tau = X1[0,0]/normalize_obs(X1, K)[0,0]
	_, _, _, _, _, _, _, s, _, _ = get_pca(normalize_obs(X1, K), K)
	uest, vest = ebpca_gaussian(X1, K)

	# scale so that u has norm 1 columns
	vest[:, 0] = vest[:, 0] * np.sqrt((uest[:, 0]**2).sum()) * s[0]/n * tau
	vest[:, 1] = vest[:, 1] * np.sqrt((uest[:, 1]**2).sum()) * s[1]/n * tau

	matr[(seed * K):(seed * K + K)] = vest.transpose()
t1 = time.time()
print('average runtime: ', (t1 - t0)/nseeds)

np.savetxt('../output/ebpca_sim1.csv', matr, delimiter=",")

