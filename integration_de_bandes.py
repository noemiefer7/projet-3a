# Intégration de bandes
import numpy as np
import pandas as pd 

with open ("aphya.txt") as f:
    wl, rp = np.loadtxt(f, dtype = float, unpack = True)
    
with open ("bandes OLCI.txt") as f:
    bmin, bmax = np.loadtxt(f, dtype = int, unpack = True)


wl = wl.astype(int )
nb = len(bmin)

wlp = np.zeros(nb)
wlp = wlp.astype(int)
rpp = np.zeros(nb)



for i in range (nb):
    indmin = np.argwhere(wl==bmin[i])[0,0]
    indmax = np.argwhere(wl==bmax[i])[0,0]
    wlp[i] = round(bmin[i] + bmax[i]/2)
    rpp[i] = np.mean(rp[indmin:indmax])
    
df_sargasse = pd.DataFrame({'wlp': wlp, 'rpp': rpp})

with open('S3_aphya.txt', 'a') as f:
    dfAsString = df_sargasse.to_string(header=False, index=False)
    f.write(dfAsString)
    
        