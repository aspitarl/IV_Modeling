

# Function Diodemodel(kT,I0,n1,Vs,n2)
# kT,I0,n1,Vs,n2

# wave turnon1 = root:Platueautime:turnon1
# wave turnon2 = root:Platueautime:turnon2



# Make/n=1000/o Vsd_model
# make/n=1000/o Isd_model

# wave Vsd_model = Vsd_model
# wave Isd_model = Isd_model

# j 

# for (j=0;j < dimsize(turnon1,0) ; j = j +1)

# i 
# Vrange = 1.5

# for (i=0;i<dimsize(Vsd_model,0); i=i+1)

# Vsd_model[i] = (Vrange/dimsize(Vsd_model,0))*i

# current0 = I0*(np.exp(Vsd_model[i]/(n1*kT))-1)

# if(current0 < turnon1[j])
# 	Isd_model[i] = current0
# else
# 	Isd_model[i] = turnon1[j]
# endif

# current1 = I0*(np.exp((Vsd_model[i]-Vs)/(n2*kT))-1)

# if(Vsd_model[i] > Vs)
# 	if((current1+turnon1[j]) <  turnon2[j] )
# 		Isd_model[i] += current1
# 	else
# 		Isd_model[i] += turnon2[j]
# 	endif
# endif

# duplicate/o Isd_model $("Isd_model_" + num2str(j))

# endfor

# endfor

# end



# Function diodeiterate()
# R = 1
# I0 = 1
# kT = 0.025

# i
# Niter = 3

# make/o/n=(Niter+1) Id,Vd
# make/o/n=(Niter+1) IR,VR

# Vtot,Itot

# Iavg

# Vtot = 1
# Vd[0] = 0.01
# VR[0] = 0.5

# for ( i =0; i<Niter;i=i+1)
# Id[i] = I0*np.exp(Vd[i]/kT)
# IR[i] = VR[i]/R

# Iavg = (Id[i] + IR[i])/2

# Vd[i+1] = kT*np.log(Iavg/I0)
# VR[i+1] = Iavg*R

# if(VR[i+1]>Vtot)
# VR[i+1] = Vtot
# endif

# endfor

# end


import numpy as np
import matplotlib.pyplot as plt

def diodefromI(I0pn,I0SBn,I0SBp, Vshift):
    kT = 0.025
    n = 1.3
    num = 50

    Isd = np.logspace(-13,-9,num)
    Vpn = np.array(range(num))
    Vsb = np.array(range(num))
    Vsbp = np.array(range(num))
    Vtot = np.array(range(num)) 

    for i in range(num):
        Vsbp[i] = -kT*np.log(1-(Isd[i]/I0SBp))

        if (Isd[i]< I0SBn):
            Vpn[i] = n*kT*np.log((Isd[i]/I0pn)+1)
            Vsb[i] = -kT*np.log(1-(Isd[i]/I0SBn))
        else:
            Vpn[i] = n*kT*np.log(((Isd[i]-I0SBn)/I0pn)+1)
            Vsb[i] = Vshift

    Vtot = Vpn +Vsb + Vsbp

    return Isd, Vtot


def WFmodel(phie,phih):
    taue = 5.5e-3
    tauh = 1

    e = 1.6e-19
    kb = 1.38e-23
    h = 6.626e-34
    #C = 1.686e-8
    #C = 7.311e-11
    T = 300

    for i in range(len(phie)):
        I0SBn = ((4*e*kb*taue)/h)*T*np.exp(-(e*phie[i])/(1000*kb*T))
        I0SBp = ((4*e*kb*tauh)/h)*T*np.exp(-(e*phih[i])/(1000*kb*T))
        I0pn = ((8*e*kb*1)/h)*T*np.exp(-(e*1.2*(phih[i]+phie[i]))/(1000*kb*T))
        #I0pn = 3e-14

        tpn = (h/(4*e*kb*T))*I0pn*np.exp((e*0.410)/(kb*T))

        Isd, Vtot = diodefromI(I0pn,I0SBn,I0SBp, 0.9)

    return Isd,Vtot

phie =  [ 0.3,0.25,0.2]
phih =  [ 0.1,0.15,0.2]

Isd, Vtot = WFmodel(phie,phih)

plt.plot(Vtot,Isd)
plt.show()