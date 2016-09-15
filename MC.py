# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 17:59:18 2013

This code runs Metropolis MC simulation for liquid argon. It was written
when I first started grad school, and also when I was learning to code.
This code is intended for learning purposes and has not been optimized.
(Some bugs are known to exist in this code.)
@author: byoo1
"""

import numpy as np
from pylab import *
from pylab import *
from math import *
#from mayavi import mlab
class Monte_Carlo:
    def __init__(self,sigma,epsilon,ptcl_number,ptcl_mass,box_length,
                 temperature,rcut):
        self.sigma=sigma
        self.epsilon=epsilon
        self.ptcl_number=ptcl_number
        self.box_length=box_length
        self.ptcl_mass=ptcl_mass
        self.T=temperature
        self.rcut=rcut
###############################################################################
    def LennardJones(self,r): #1-D
        #Constants
        sigma=self.sigma
        eps=self.epsilon
        #Lennard Jones Equation
        s_r=sigma/r
        W_r= 4.*eps*((s_r)**12-(s_r)**6)
        return W_r
###############################################################################       
    def Force(self,rx,ry,rz,d):
        #Constants
        sigma=self.sigma
        eps=self.epsilon
        #Force Equation
        s_d=sigma/(d)
        F_rx =24.*eps*rx*(2.*(s_d)**12-(s_d)**6)/(d**2)
        F_ry =24.*eps*ry*(2.*(s_d)**12-(s_d)**6)/(d**2)
        F_rz =24.*eps*rz*(2.*(s_d)**12-(s_d)**6)/(d**2)
        return (F_rx,F_ry,F_rz)   
############################################################################### 
    def Periodic(self,r):
        L=self.box_length
        if r<-L/2:   
            r+=L    #r=abs(r)%(L/2.)*-1.+L/2.
        elif r>=L/2:
            r-=L    #r=r%(L/2.)-L/2.
        return r
###############################################################################
    def FCC(self,rx,ry,rz):
        #Constants
        sigma=self.sigma
        N=self.ptcl_number
        L=self.box_length
        while len(rx)<N:
            for i in range(int(L/sigma)-1):
                for j in range(int(L/sigma)-1):
                    for k in range(int(L/sigma)-1):
                        #corners
                        rz=np.append(rz,i*sigma)
                        ry=np.append(ry,j*sigma)
                        rx=np.append(rx,k*sigma)
                        #center in same plane as corners
                        rz=np.append(rz,i*sigma)
                        ry=np.append(ry,(j+0.5)*sigma)
                        rx=np.append(rx,(k+0.5)*sigma)
                        #side wall faces
                        rz=np.append(rz,(i+0.5)*sigma)
                        ry=np.append(ry,(j+0.5)*sigma)
                        rx=np.append(rx,(k)*sigma)
                        #bottom and top faces
                        rz=np.append(rz,(i+0.5)*sigma)
                        ry=np.append(ry,j*sigma)
                        rx=np.append(rx,(k+0.5)*sigma)
        if len(rx)>N:
            rx=rx[:N]
            ry=ry[:N]
            rz=rz[:N]
        r=np.vstack([rx,ry,rz])
        return r
############################################################################### 
    def Initialize(self):
        #empty arrays
        rx=np.array([])
        ry=np.array([])
        rz=np.array([])
        #initialize coord with FCC lattice
        old_coord=self.FCC(rx,ry,rz)
        rx_old=old_coord[0]
        ry_old=old_coord[1]
        rz_old=old_coord[2]
        for i in range(len(rx_old)):
            rx_old[i]=self.Periodic(rx_old[i])
            ry_old[i]=self.Periodic(ry_old[i])
            rz_old[i]=self.Periodic(rz_old[i])
        #'current' coordinates after imposing initial velocities
        r_old=np.vstack([rx_old,ry_old,rz_old])
        return (r_old)
###############################################################################
    def Init_Potential(self,rx,ry,rz):
        #Constants
        N=self.ptcl_number
        rcut=self.rcut
        #Initial distance/force arrays
        rx_ab=np.array([])
        ry_ab=np.array([])
        rz_ab=np.array([])
        d_ab=np.array([])
        d_cut=np.array([])
        #Calculate interparticle distances
        for i in range(N-1):
            for j in range(i+1,N):
                #minimum image convention
                rx_ab = np.append(rx_ab,self.Periodic(rx[i]-rx[j]))
                ry_ab = np.append(ry_ab,self.Periodic(ry[i]-ry[j]))
                rz_ab = np.append(rz_ab,self.Periodic(rz[i]-rz[j]))
                d_ab=np.append(d_ab,(rx_ab[-1]**2+ry_ab[-1]**2+rz_ab[-1]**2)**0.5)
        r=np.vstack([rx,ry,rz])
        for d in d_ab:
            if d<rcut:
                d_cut=np.append(d_cut,d)
        V= sum(self.LennardJones(d_cut))
        return (r,V)
###############################################################################
    def Potential(self,rx,ry,rz,index):
        rx_ab=np.array([])
        ry_ab=np.array([])
        rz_ab=np.array([])
        d_ab=np.array([])
        d_cut=np.array([])
        for i in range(len(rx)):
            if index != i:
                rx_ab = np.append(rx_ab,self.Periodic(rx[index]-rx[i]))
                ry_ab = np.append(ry_ab,self.Periodic(ry[index]-ry[i]))
                rz_ab = np.append(rz_ab,self.Periodic(rz[index]-rz[i]))
                d_ab=np.append(d_ab,(rx_ab[-1]**2+ry_ab[-1]**2+rz_ab[-1]**2)**0.5)
        for d in d_ab:
            if d<rcut:
                d_cut=np.append(d_cut,d)
        dV= sum(self.LennardJones(d_cut))
        return dV
###############################################################################
    def Virial(self,rx,ry,rz):
        #Constants
        N=self.ptcl_number
        r_cut=self.rcut
        #Initial distance/force arrays
        rx_ab=np.array([])
        ry_ab=np.array([])
        rz_ab=np.array([])
        d_ab=np.array([])
        fx=np.zeros(N)
        fy=np.zeros(N)
        fz=np.zeros(N)
        w=0.
        #Calculate interparticle distances
        for i in range(N-1):
            for j in range(i+1,N):
                #minimum image convention
                rx_ab = np.append(rx_ab,self.Periodic(rx[i]-rx[j]))
                ry_ab = np.append(ry_ab,self.Periodic(ry[i]-ry[j]))
                rz_ab = np.append(rz_ab,self.Periodic(rz[i]-rz[j]))
                d_ab=np.append(d_ab,(rx_ab[-1]**2+ry_ab[-1]**2+rz_ab[-1]**2)**0.5)
                #Force calculations for particles within cutoff radius
                if d_ab[-1] < r_cut:
                    F=self.Force(rx_ab[-1],ry_ab[-1],rz_ab[-1],d_ab[-1])
                    fx[i]=fx[i]+F[0]
                    fy[i]=fy[i]+F[1]
                    fz[i]=fz[i]+F[2]
                    fx[j]=fx[j]-F[0]
                    fy[j]=fy[j]-F[1]     
                    fz[j]=fz[j]-F[2]
                    #virial (for pressure)
                    w+=rx_ab[-1]*F[0]+ry_ab[-1]*F[1] +rz_ab[-1]*F[2]
        w=w/3.
        return (w)
###############################################################################
    def Equil(self,max_dx,sample,nadjust):
        N=self.ptcl_number
        kb=1.
        T=self.T
        Beta=1./kb/T
        coord=self.Initialize()
        rx=coord[0]
        ry=coord[1]
        rz=coord[2]
        V=np.array([])
        Y_N=np.array([0., 0.])
        data=self.Init_Potential(rx,ry,rz)
        V=np.append(V,data[1])
        for i in range(sample):
            index=np.int(np.floor(N*np.random.uniform(0,1)))       
            dV_old=self.Potential(rx,ry,rz,index)
            dx=np.random.uniform(-1,1)*max_dx
            dy=np.random.uniform(-1,1)*max_dx
            dz=np.random.uniform(-1,1)*max_dx
            rx[index]=rx[index]+dx
            ry[index]=ry[index]+dy
            rz[index]=rz[index]+dz
            dV_new=self.Potential(rx,ry,rz,index)
            delV=dV_new-dV_old
            delBV=Beta*(delV)
            if delBV<75:
                if delBV<=0.0:
                    V=np.append(V,V[-1]+delV)
                    Y_N[0]+=1 #accept
                elif np.exp(-delBV) > np.random.uniform(0,1): #metropolis
                    V=np.append(V,V[-1]+delV)
                    Y_N[0]+=1 #accept
                else:
                    V=np.append(V,V[-1])
                    rx[index]=rx[index]-dx
                    ry[index]=ry[index]-dy
                    rz[index]=rz[index]-dz
                    Y_N[1]+=1
            else:
                V=np.append(V,V[-1])
                rx[index]=rx[index]-dx
                ry[index]=ry[index]-dy
                rz[index]=rz[index]-dz
                Y_N[1]+=1 #reject
            if i>10000 and i%nadjust==0:
                if Y_N[0]/(Y_N[0]+Y_N[1])> 0.5:
                    max_dx=max_dx*1.05
                else:
                    max_dx=max_dx*0.95                             
        return (Y_N,V,max_dx,rx,ry,rz)
###############################################################################
    def NVT(self,max_dx,sample,nadjust,P_record,init_coord=0):
        N=self.ptcl_number
        kb=1.
        T=self.T
        L=self.box_length
        rcut=self.rcut
        Beta=1./kb/T
        #coord=self.Initialize()
        coord=init_coord
        rx=coord[0]
        ry=coord[1]
        rz=coord[2]
        V=np.array([])
        P=np.array([])
        Y_N=np.array([0., 0.])
        data=self.Init_Potential(rx,ry,rz)
        V=np.append(V,data[1])
        for i in range(sample):
            index=np.int(np.floor(N*np.random.uniform(0,1)))       
            dV_old=self.Potential(rx,ry,rz,index)
            dx=np.random.uniform(-1,1)*max_dx
            dy=np.random.uniform(-1,1)*max_dx
            dz=np.random.uniform(-1,1)*max_dx
            rx[index]=rx[index]+dx
            ry[index]=ry[index]+dy
            rz[index]=rz[index]+dz
            dV_new=self.Potential(rx,ry,rz,index)
            delV=dV_new-dV_old
            #print V[-1]+delV
            delBV=Beta*(delV)
            #print delV
            if delBV<75:
                if delBV<=0.0:
                    V=np.append(V,V[-1]+delV)
                    Y_N[0]+=1 #accept
                elif np.exp(-delBV) > np.random.uniform(0,1):
                    V=np.append(V,V[-1]+delV)
                    Y_N[0]+=1 #accept
                else:
                    V=np.append(V,V[-1])
                    rx[index]=rx[index]-dx
                    ry[index]=ry[index]-dy
                    rz[index]=rz[index]-dz
                    Y_N[1]+=1
            else:
                V=np.append(V,V[-1])
                rx[index]=rx[index]-dx
                ry[index]=ry[index]-dy
                rz[index]=rz[index]-dz
                Y_N[1]+=1 #reject
            if i>10000 and i%nadjust==0:
                if Y_N[0]/(Y_N[0]+Y_N[1])> 0.5:
                    max_dx=max_dx*1.05
                else:
                    max_dx=max_dx*0.95
            if i>10000 and i%P_record==0:
                W=self.Virial(rx,ry,rz)
                inst_p=N*kb*T/(L**3)+W/(L**3)
                P=np.append(P,inst_p)
        V=V/N
        return (Y_N,V,max_dx,P,rx,ry,rz)
###############################################################################
sigma=1.#3.405e-10
kb=1.#/119.8#1.38e-23
eps=1.#119.8*kb
N=256
rho_r=0.6
L=(N/rho_r)**(1/3.)*sigma #L=10*sigma
m=1.#39.948*1.6747e-27
T_r=1.4
T=T_r/kb*eps
rcut=0.5*L
cycle=500000
nadjust=1000
P_record=1000
equil_cycle=500000

MC=Monte_Carlo(sigma,eps,N,m,L,T_r,rcut)
#equil_r=MC.Equil(0.02*sigma,cycle,nadjust)
#np.savetxt('equil_N256.txt',transpose((equil_r[3],equil_r[4],equil_r[5])))
'''
init_coord=np.loadtxt('equil_N256.txt')
data=np.vstack((init_coord[:,0],init_coord[:,1],init_coord[:,2]))
V=MC.NVT(0.02*sigma,cycle,nadjust,P_record,data)

print 'Energy: ',np.average(V[1][:])
print 'Pressure:',np.average(V[3])
print 'Yes, No: ',V[0]
print 'max displacement:',V[2]
np.savetxt('Energy_MC1.txt',transpose((V[1])))
np.savetxt('Pressure_MC1.txt',transpose((V[3])))
'''
init_coord=np.loadtxt('equil_N256.txt')
data=np.vstack((init_coord[:,0],init_coord[:,1],init_coord[:,2]))
W=MC.NVT(0.02*sigma,cycle,nadjust,P_record,data)

print 'Energy: ',np.average(W[1][:])
print 'Pressure:',np.average(W[3])
print 'Yes, No: ',W[0]
print 'max displacement:',W[2]
np.savetxt('Energy_MC2.txt',transpose((W[1])))
np.savetxt('Pressure_MC2.txt',transpose((W[3])))

init_coord=np.loadtxt('equil_N256.txt')
data=np.vstack((init_coord[:,0],init_coord[:,1],init_coord[:,2]))
Z=MC.NVT(0.02*sigma,cycle,nadjust,P_record,data)
print 'Energy: ',np.average(Z[1][:])
print 'Pressure:',np.average(Z[3])
print 'Yes, No: ',Z[0]
print 'max displacement:',Z[2]
np.savetxt('Energy_MC3.txt',transpose((Z[1])))
np.savetxt('Pressure_MC3.txt',transpose((Z[3])))


print 'FINISHED!'
        
        
