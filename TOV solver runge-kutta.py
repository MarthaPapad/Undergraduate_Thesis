import numpy as np
import matplotlib.pyplot as plt

def en(P,Beff):
    return(3*P+4*Beff)

def mr(r,P,Beff):
    return(11.2*10**(-6)*r**2*en(P,Beff))

def Pr(r,P,m,Beff):
    return((-1.474*en(P,Beff)*m)/r**2*(1+P/en(P,Beff))*(1+11.2*10**(-6)*r**3*P/m)*(1-2.948*m/r)**(-1))

def TOV(P0, m0, dr,Beff):
    P = [P0]
    m = [m0]
    r = [0.01]
    h = dr
    i = 0

    while P[i]>0.000000000001:
        k0 = h * Pr(r[i], P[i], m[i], Beff)
        l0 = h * mr(r[i], P[i], Beff)
        k1 = h * Pr(r[i]+h/2, P[i]+k0/2, m[i]+l0/2, Beff)
        l1 = h * mr(r[i]+h/2, P[i]+k0/2, Beff)
        k2 = h * Pr(r[i]+h/2, P[i]+k1/2, m[i]+l1/2, Beff)
        l2 = h * mr(r[i]+h/2, P[i]+k1/2, Beff)
        k3 = h * Pr(r[i]+h, P[i]+k2, m[i]+l2, Beff)
        l3 = h * mr(r[i]+h, P[i]+k2, Beff)

        P.append(P[i]+1/6*(k0+2*k1+2*k2+k3))
        m.append(m[i]+1/6*(l0+2*l1+2*l2+l3))
        r.append(r[i] + h)
        
        i += 1

    return m[-1], r[-1]


Beff_values = [13, 27, 58, 109, 209]
Mmax_list = []
all_Pep_data = []
all_Mr_data = []
condition = False

for Beff in Beff_values:
    M_list = []
    R_list = []
    ep_list = []
    P0_list = []
    M_max = -1.0
    condition = False
    
    for ep in range(4*Beff, 5001, 1):
        P0 = 1/3*(ep-4*Beff)
        final_M, final_R = TOV(P0,0.001,0.001,Beff)
        M_list.append(final_M)
        R_list.append(final_R)
        ep_list.append(ep)
        P0_list.append(P0)

        if final_M >= M_max:
            M_max = final_M
        else:
            condition = True
            break

    if not condition:
        break
     
    Mmax_list.append(M_max)

    all_Pep_data.append((ep_list, P0_list, f'$B_{{eff}} = {Beff} \, MeV \cdot fm^{{-3}}$'))
    all_Mr_data.append((R_list, M_list, f'$B_{{eff}} = {Beff} \, MeV \cdot fm^{{-3}}$'))



plt.figure()
for ep_list, P0_list, label in all_Pep_data:
    plt.plot(ep_list, P0_list, label=label)
plt.xlabel('$\epsilon \, (MeV \cdot fm^{-3})$')
plt.ylabel('$P \, (MeV \cdot fm^{-3})$')
plt.title('Pressure vs Energy Density')
plt.grid(True)
plt.legend()


plt.figure()
for R_list, M_list, label in all_Mr_data:
    plt.plot(R_list, M_list, label=label)
plt.xlabel('$R \, (km)$')
plt.ylabel('$M/M_{\odot}$')
plt.title('Mass vs Radius')
plt.grid(True)
plt.legend()


plt.figure()
for R_list, M_list, label in all_Mr_data:
    R_cubed_list = [R ** 3 for R in R_list]
    plt.plot(R_cubed_list, M_list, label=label)
plt.xlabel('$R^3 \, (km^3)$')
plt.ylabel('$M/M_{\odot}$')
plt.title('Mass vs $Radius^3$')
plt.grid(True)
plt.legend()


plt.figure()
plt.plot(Beff_values, Mmax_list, 'bo-')
plt.xlabel('$B_{eff} \, (MeV \cdot fm^{-3})$')
plt.ylabel('$M_{max}/M_{\odot}$')
plt.title('Maximum Mass vs $B_{eff}$')
plt.grid(True)

plt.show()
