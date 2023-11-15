# Détection des sargasses

import numpy as np
import matplotlib.pyplot as plt



# Spectre d'absorption du phytoplancton
with open ("aphya.txt") as f:
    wl, a_phy_ = np.loadtxt(f, dtype = float, unpack = True)
 
# Spectre d'absorption de l'eau pure
with open ("aw.txt") as f:
    wl, aw = np.loadtxt(f, dtype = float, unpack = True) 
    
# Réflectance des sargasses
with open ("r_Sarg.txt") as f : 
    wl, r_sarg = np.loadtxt(f, dtype = float, unpack = True)
    
    
# Modèle de Lee adapté aux sargasses
def LEE(X, wl, tetaw, tetav):
    
    ###### Paramètres d'entrée #####
    C_chl = X[0]
    C_nap = X[1]
    C_cdom = X[2]
    z = X[3]
    FC = X[4]  #fraction de couverture de sargasses
    
    ###### Constantes #####
    S_cdom = 0.0157
    S_nap = 0.0106
    a_nap = 0.0048
    b_bphy = 0.00038
    Y_phy = 0.681
    b_bnap = 0.0054
    Y_nap = 0.681
    
    
    ##### Calcul de l'absorption spectrale a #####
    a_phy = C_chl*a_phy_;
    a_nap = C_nap*a_nap*np.exp(-S_nap*(wl-440))
    a_cdom = C_cdom*np.exp(-S_cdom*(wl-440))
    
    a = aw + a_phy + a_nap + a_cdom 
    
    
    ###### Calcul de la rétro-diffusion spectrale bb #####
    b_bw = 0.00144*(wl/500)**-4.32
    b_bp = C_chl*b_bphy*(542/wl)**Y_phy + C_nap*b_bnap*(542/wl)**Y_nap
    
    b_b = b_bw + b_bp
    
    ##### Calcul de la réflectance #####
    # Coefficient d'atténuation diffuse
    K = a + b_b   
    
    u = b_b/(a + b_b)
    u_p = b_bp/(a + b_b)
    
    D_uc = 1.03*(1 + 2.4*u)**0.5
    D_ub = 1.04*(1 + 5.4*u)**0.5
    
    gp = 0.184*(1-0.602*np.exp(-3.852*u_p))
    gw = 0.115
    
    r_rsdp = gw*b_bw/(a + b_b) + gp*b_bp*(a + b_b)
    
    tetaw = tetaw*np.pi/180
    tetav = tetav*np.pi/180
    rho_f = FC*r_sarg + (1-FC)*r_rsdp
    
    r_rsC = r_rsdp*(1-np.exp(-(1/np.cos(tetaw) + D_uc/np.cos(tetav))*K*z))
    r_rsb = 1/np.pi*rho_f*np.exp(-(1/np.cos(tetaw) + D_ub/np.cos(tetav))*K*z)
    
    r_rs = r_rsC + r_rsb
    
    Rrs = 0.52*r_rs/(1 - 1.56*r_rs)
    
    return Rrs
    
    
X = [0.1, 0.1, 0.01, 0, 0.5]
tetaw = 50
tetav = 25

Rrs = LEE(X, wl, tetaw, tetav)

plt.figure()
plt.plot(wl, Rrs)
plt.title('Réflectance de surface (Modèle de Lee)')
plt.xlabel("Longueurs d'onde")
plt.ylabel('Réflectance')
plt.show()
    

    
    
