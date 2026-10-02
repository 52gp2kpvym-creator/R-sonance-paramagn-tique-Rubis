#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 16:10:22 2026

@author: baptisteviu
"""

import numpy as np
import matplotlib.pyplot as plt
e=1.602e-19#C
m=9.109e-31#kg
g=2
S=3/2
D=5.73e9#Hz
pi=np.pi
hbar=1.05e-34#J.s
kb=1.380e-24#J.K
f_cavite=9.19*1e9#Hz
#Calcul de l'hamiltonien
def H(theta,B0):
    M=[[0 for k in range(4)]for k in range(4)]
    M[0][0]=2*pi*D*hbar+(3*g*e*B0/(4*m))*hbar*np.cos(theta)
    M[0][1]=g*e*B0*np.sin(theta)*np.sqrt(3)*hbar/(4*m)
    M[1][1]=-2*pi*D*hbar+(g*e*B0/(4*m))*hbar*np.cos(theta)
    M[1][0]=g*e*B0*np.sin(theta)*np.sqrt(3)*hbar/(4*m)
    M[1][2]=g*e*B0*np.sin(theta)*hbar/(2*m)
    M[2][2]=-2*pi*D*hbar-(g*e*B0/(4*m))*np.cos(theta)*hbar
    M[2][1]=g*e*B0*np.sin(theta)*hbar/(2*m)
    M[2][3]=g*e*B0*np.sin(theta)*hbar*np.sqrt(3)/(4*m)
    M[3][3]=+2*pi*D*hbar-(3*g*e*B0/(4*m))*hbar*np.cos(theta)
    M[3][2]=g*e*B0*np.sin(theta)*np.sqrt(3)*hbar/(4*m)
    M=np.array(M)
    return M
    

liste_thetadeg=np.linspace(0,190,5000)
liste_thetarad=pi*liste_thetadeg/180
B=np.linspace(0,6000e-4,5000)



#Calcul de toutes les valeurs propres et vecteurs propres possiblles pour tous les angles
# et tous les champs magnétiques
VP=[[0 for k in range(len(B))] for k in range(len(liste_thetarad))]
VecP=[[0 for k in range(len(B))] for k in range(len(liste_thetarad))]
for i in range (len(liste_thetarad)):
    for j in range(len(B)):
        VP[i][j]=np.linalg.eigh(H(liste_thetarad[i],B[j]))[0]
        VecP[i][j]=np.linalg.eigh(H(liste_thetarad[i],B[j]))[1]
VP=np.array(VP)
VecP=np.array(VecP)

#Question2
theta=40*pi/180
E=[]
for b in B:
    valeurspropres=np.linalg.eigh(H(theta,b))[0]
    E.append(valeurspropres)
E=np.array(E)
E=E/(2*pi*hbar)#Hz
E=E/1e9#GHz

# #Vecteurs propres pour B=0G
print(np.linalg.eigh(H(theta*pi/180, 0))[1])
# #Pour B=6000G
print(np.linalg.eigh(H(theta*pi/180, 6000))[1].T)






# # Représentation des niveaux d'énergie
plt.plot(B*1e4,E[:,0],label='état fondamental')
plt.plot(B*1e4,E[:,1],label='premier état excité')
plt.plot(B*1e4,E[:,2],label='deuxième état excité')
plt.plot(B*1e4,E[:,3],label='troisième état excité')
plt.xlabel("B(en G)")
plt.ylabel("fréquence(GHz)")
plt.legend()
plt.grid()



# #Transitions
Transition10=E[:,1]-E[:,0]
Transition20=E[:,2]-E[:,0]
Transition30=E[:,3]-E[:,0]
Transition21=E[:,2]-E[:,1]
Transition32=E[:,3]-E[:,2]
Transition31=E[:,3]-E[:,1]
Cavite=[f_cavite/1e9 for k in range(5000)]

#Recherche  grossière champs B de transition pour theta=40°


def intersection(L1,L2):
    n=len(L1)
    L1,L2=np.array(L1),np.array(L2)
    new_L=L1-L2
    res=[]
    for k in range(n-1):
        if new_L[k]*new_L[k+1]<=0:
            res+=[(k,k+1)]
    return res
           

indices_intersection10_min,indices_intersection10_max=intersection(Cavite,Transition10)[0]
indices_intersection21_min_1,indices_intersection21_max_1=intersection(Cavite,Transition21)[1]
indices_intersection21_min,indices_intersection21_max=intersection(Cavite,Transition21)[0]
indices_intersection32_min,indices_intersection32_max=intersection(Cavite,Transition32)[0]
B_01_1=((B[indices_intersection10_min]+B[indices_intersection10_max])/2)*1e4
B_21_2=((B[indices_intersection21_min_1]+B[indices_intersection21_max_1])/2)*1e4
B_21=((B[indices_intersection21_min]+B[indices_intersection21_max])/2)*1e4
B_32=((B[indices_intersection32_min]+B[indices_intersection32_max])/2)*1e4





print(B_01_1,B_21_2,B_21,B_32)