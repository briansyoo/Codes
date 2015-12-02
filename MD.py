#!/usr/bin/env python
#coding: utf8 
"""
@author: byoo1

This code runs MD simulation for liquid argon.
This code is intended for learning purposes and has not been optimized.
"""
import numpy as np #for fast array operations
import sys, os, argparse, linecache
from pylab import * #for plots
import pdb #for debugging

#from mayavi import mlab
class Molecular_Dynamics:
#stores attributes into class
    def __init__(self,sigma,epsilon,ptcl_number,ptcl_mass,box_length,
                 temperature,timestep,rcut,nout):
        self.sigma=sigma
        self.epsilon=epsilon
        self.ptcl_number=ptcl_number
        self.box_length=box_length
        self.cutoff=rcut
        self.ptcl_mass=ptcl_mass
        self.T=temperature
        self.timestep=timestep
        self.nout=nout
###############################################################################
#determines LJ Pairwise Potential
    def LennardJones(self,r): #1-D
        #Constants
        sigma=self.sigma
        eps=self.epsilon   
        #Lennard Jones Equation
        sr=sigma/r
        sr2=sr*sr
        sr6=sr2*sr2*sr2
        sr12=sr6*sr6
        W_r= 4.*eps*(sr12-sr6)
        return W_r
###############################################################################
#determines LJ Pairwise Forces
       
    def Force(self,rx,ry,rz,d):  
        #Constants
        sigma=self.sigma
        eps=self.epsilon
        #Force Equation
        sd=sigma/(d)
        sd2=sd**2
        sd6=sd2*sd2*sd2
        sd12=sd6*sd6
        d2=d*d
        F_rx =24.*eps*rx*(2.*sd12-sd6)/d2
        F_ry =24.*eps*ry*(2.*sd12-sd6)/d2
        F_rz =24.*eps*rz*(2.*sd12-sd6)/d2
        return (F_rx,F_ry,F_rz)   
###############################################################################
#Applies peridic boundary conditions provided the box length 
    def Periodic(self,r,L):      
        if r<-L/2.:   
            r+=L    #r=abs(r)%(L/2.)*-1.+L/2.
        elif r>=L/2.:
            r-=L    #r=r%(L/2.)-L/2.
        return r
###############################################################################
#Determines the kinetic energies based on the cartesian velocities 
    def kinE(self,vx,vy,vz):
        m=self.ptcl_mass
        v2x=m*sum(vx**2)
        v2y=m*sum(vy**2)
        v2z=m*sum(vz**2)
        KE=0.5*(v2x+v2y+v2z)        
        return KE
###############################################################################
#Velocity rescaling function for thermostat

    def vscale(self,T,v2):
        N=self.ptcl_number
        kb=1.38e-23
        scale=(T/(v2/(N-3)/kb/3))**0.5    
        return scale          
###############################################################################
#Returns FCC lattice coordinates provided the empty coordinate arrays
    def FCC(self,rx,ry,rz):
        #Constants
        sigma=self.sigma
        N=self.ptcl_number
        L=self.box_length
        while len(rx)<N:
            for i in range(int(L/sigma)):
                for j in range(int(L/sigma)):
                    for k in range(int(L/sigma)):
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
                        rx=np.append(rx,k*sigma)
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
#Provides initial coordinates and velocities for the FCC lattice initial configuration 
    def Initialize(self):
        #Constants
        sigma=self.sigma
        N=self.ptcl_number
        L=self.box_length
        #empty arrays
        rx=np.array([])
        ry=np.array([])
        rz=np.array([])
        #initial velocities
        vel=np.zeros((N,3))      
        for j in range(3):
            for i in range(N):
                vel[i,j]=np.random.uniform(-0.5*sigma,0.5*sigma)
        #set net momentum to 0
        sum_velx=sum(vel[:,0])
        sum_vely=sum(vel[:,1])
        sum_velz=sum(vel[:,2])
        vel[:,0]-=sum_velx/N
        vel[:,1]-=sum_vely/N
        vel[:,2]-=sum_velz/N
        #store xyz/vel coords into separate lists
        vx=np.array(vel[:,0])
        vy=np.array(vel[:,1])
        vz=np.array(vel[:,2])
        #initialize coord with FCC lattice
        old_coord=self.FCC(rx,ry,rz)
        rx_old=old_coord[0]
        ry_old=old_coord[1]
        rz_old=old_coord[2]
        for i in range(len(rx_old)):
            rx_old[i]=self.Periodic(rx_old[i],L)
            ry_old[i]=self.Periodic(ry_old[i],L)
            rz_old[i]=self.Periodic(rz_old[i],L)
        #'current' coordinates after imposing initial velocities
        v=np.vstack([vx,vy,vz])
        r_old=np.vstack([rx_old,ry_old,rz_old])
        return (r_old,v)
###############################################################################
#detemrines the pairwise distances, forces, potentials, and virials (for pressure)
# of coordinate arrays
    def Pairwise(self,rx,ry,rz):
        #Constants
        N=self.ptcl_number
        r_cutoff=self.cutoff
        L=self.box_length
        #Initial distance/force arrays
        rx_ab=np.array([])
        ry_ab=np.array([])
        rz_ab=np.array([])
        d_ab=np.array([])
        d_cut=np.array([])
        fx=np.zeros(N)
        fy=np.zeros(N)
        fz=np.zeros(N)
        w=0.
        #Wrap coordinates
        for i in range(len(rx)):
            rx[i]=self.Periodic(rx[i],L)
            ry[i]=self.Periodic(ry[i],L)
            rz[i]=self.Periodic(rz[i],L)
        #Calculate interparticle distances
        for i in range(N-1):
            for j in range(i+1,N):
                #minimum image convention
                rx_ab = np.append(rx_ab,self.Periodic((rx[i]-rx[j]),L))
                ry_ab = np.append(ry_ab,self.Periodic((ry[i]-ry[j]),L))
                rz_ab = np.append(rz_ab,self.Periodic((rz[i]-rz[j]),L))
                d_ab=np.append(d_ab,(rx_ab[-1]**2+ry_ab[-1]**2+rz_ab[-1]**2)**0.5)
                #Force calculations for particles within cutoff radius
                if d_ab[-1] < r_cutoff:
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
        r=np.vstack([rx,ry,rz])
        f=np.vstack([fx,fy,fz])
        for d in d_ab:
            if d<r_cutoff:
                d_cut=np.append(d_cut,d)
        V= np.sum(self.LennardJones(d_cut))
        return (r,f,V,w)
###############################################################################
#applies Velocity Verlet integrator (finite differences) method;
# this is the main part of the program
    def Velocity_Verlet(self,time,ensemble):
        #Constants
        sigma=self.sigma
        N=self.ptcl_number
        m=self.ptcl_mass
        L=self.box_length
        T=self.T
        rcut=self.cutoff
        dt=self.timestep
        nout=self.nout
        kb=1.38e-23
        #total number of snapshots/cycles for simulation
        cycle=int(time/dt)      
        #open empty output file for trajectories/coordinates and energy
        trjfile = open("traj.dat",'w')
        trjfile.write('             x             y             z            vx            vy            vz'+'\n') 
        edrfile = open("edr.dat",'w')
        edrfile.write('          Time          Temp      Pressure        Pot. E        Kin. E        Tot. E'+'\n')
        fmt="%14.6e"
        
        #initialize coordinates and velocities
        init=self.Initialize()
        coord=init[0]
        v_current=init[1]
        #determine initial accelerations based on initial coordinates
        data=self.Pairwise(coord[0],coord[1],coord[2])
        r_current=data[0]
        a_current=data[1]/m
        v_half=v_current+0.5*a_current*dt
        r_new=r_current+v_half*dt 
        V=data[2]
        for i in range(cycle-1):
            new_data=self.Pairwise(r_new[0],r_new[1],r_new[2])
            r_new=new_data[0]
            a_new=new_data[1]/m
            V=new_data[2]
            v_new=v_half+0.5*(a_new)*dt
            #conservation of momentum for first 200 steps
            if i%100==0 or i<100:
                sum_vx=sum(v_new[0]/sigma)
                sum_vy=sum(v_new[1]/sigma)
                sum_vz=sum(v_new[2]/sigma)        
                if sum_vx>0:
                    v_new[0]-=sum_vx/N*sigma
                    v_new[1]-=sum_vy/N*sigma
                    v_new[2]-=sum_vz/N*sigma    
            #kinetic energy
            KE=self.kinE(v_new[0],v_new[1],v_new[2])
            inst_temp=KE*2./((N-3)*kb*3.)
            inst_p=N*kb*T/(L**3)+new_data[3]/(L**3)
            #velocity scale factor
            if ensemble=='NVT' or ensemble=='nvt':
                case=cycle
            elif ensemble=='NVE' or ensemble=='nve':
                case=cycle/4. #runs nvt first for 1/4 of the nve cycle
            if i<case:
                vscale=self.vscale(T,KE*2)
                v_new=v_new*vscale
                KE=self.kinE(v_new[0],v_new[1],v_new[2])
                inst_temp=KE*2./((N-3)*kb*3.)
            #instant temperature/pressure
            #reset for next iteration
            r_current=r_new
            a_current=a_new
            v_current=v_new
            #Velocity Verlet
            v_half=v_current+0.5*a_current*dt
            r_new=r_current+dt*v_half
            #write values every nout iterations
            if i%nout ==0: 
                edrfile.write(fmt%(i*dt))
                edrfile.write(fmt%(inst_temp*kb/eps))
                edrfile.write(fmt%(inst_p*(sigma**3)/eps))
                edrfile.write(fmt%(V/eps/N))
                edrfile.write(fmt%(KE/eps/N))
                edrfile.write(fmt%((V+KE)/eps/N)+'\n')
                for coords in range(N):
                    trjfile.write(fmt%(r_new[0][coords]/sigma))
                    trjfile.write(fmt%(r_new[1][coords]/sigma))
                    trjfile.write(fmt%(r_new[2][coords]/sigma))
                    trjfile.write(fmt%(v_current[0][coords]))
                    trjfile.write(fmt%(v_current[1][coords]))
                    trjfile.write(fmt%(v_current[2][coords])+'\n')
        #PBC on last iterated coords
        for i in range(N):
            r_new[0][i]=self.Periodic((r_new[0][i]),L)/sigma
            r_new[1][i]=self.Periodic((r_new[1][i]),L)/sigma
            r_new[2][i]=self.Periodic((r_new[2][i]),L)/sigma
###############################################################################

###################### P o s t  A n a l y s i s ###############################
#code for radial distribution function
    def RDF(self, bin_number,start_cycle):
        N=self.ptcl_number
        L=self.box_length/self.sigma #reduced units
        nout=self.nout
        #loads trajectory file
        kb=1.38e-23 #kj/molecule
        eps=self.epsilon #kJ/molecule
        T=self.T #Kelvin
        T_r=T*kb/eps #reduced units
        sigma=self.sigma #meters
        xyz_trj = np.loadtxt('traj.dat')
        cycles=int(len(xyz_trj[:,0])/N)
        
        #maximum total bin length determined by half the diagonal of the cubic box
        max_bin_length = sqrt(3)/2*L #reduced units

        #create emty array for bins
        bins = np.zeros(bin_number)
        #detemines bin width
        bin_width = max_bin_length/bin_number #reduced units
        #loops over all particles to determine particle-particle distances
        k=0
        
        while k < len(xyz_trj[:,0]):
            for i in range(N-1):
                for j in range(i+1,N):
                    if k>start_cycle*N/nout:
                        rx=self.Periodic((xyz_trj[i+k,0]-xyz_trj[j+k,0]),L)
                        ry=self.Periodic((xyz_trj[i+k,1]-xyz_trj[j+k,1]),L)
                        rz=self.Periodic((xyz_trj[i+k,2]-xyz_trj[j+k,2]),L)
                        r_dist=(rx**2+ry**2+rz**2)**(0.5)
                        #pdb.set_trace()
                        bins[int(np.floor(r_dist/bin_width))]+=1                        
            k+=N
            
        normalization=np.zeros(len(bins))
        sigma_axis=np.zeros(len(bins))
        i=0
        for i in range(len(bins)):
            normalization[i]=4/3.0*3.14*0.8*(((i+1)*bin_width)**3-(i*bin_width)**3)
            sigma_axis[i]=i*bin_width

        rdf=(bins/(cycles-start_cycle/nout))/normalization/N*2
        
        
        
        figure(num=None, figsize=(8, 20), dpi=80, facecolor='w', edgecolor='k')

        subplot(5,1,1)
        plot(sigma_axis,rdf,'b',sigma_axis,np.ones(len(rdf)),'k-.')
        xlim(0,2.5)
        #xlabel('r/sigma')
        ylabel('g(r)')
        title('T*= '+str(T_r)+'$ rho$*= '+str(rho_r))
        
        #w_r=-119.8*kb*T_r*log(rdf);exp(-kb*118*T_r*((-119.8*kb*T_r*(log(rdf))-v_r)))
        
        v_r=np.zeros(len(rdf))
        for i in range(1,len(rdf)):
            v_r[i]=self.LennardJones(sigma_axis[i]*sigma) #real units

        start=45 #starting index since log(0)=-inf;guess and checked
        subplot(5,1,2)
        plot(sigma_axis,rdf*exp(v_r/T_r/119.8/kb),'b',sigma_axis,np.ones(len(rdf)),'k-.')
        #xlabel('r/sigma')
        ylabel('g(r)$_i$$_n$$_d$$_i$$_r$$_e$$_c$$_t$')
        xlim(0,2.5)
        
        subplot(5,1,3)
        plot(sigma_axis,log(rdf)*(-1/kb/T),'b',sigma_axis,np.zeros(len(rdf)),'k-.')
        xlabel('r/$\sigma$')
        ylabel('w(r)/ kJ/molecule')
        xlim(0,2.5)
        
        subplot(5,1,4)
        plot(sigma_axis[start:],v_r[start:],'b',sigma_axis,np.zeros(len(rdf)),'k-.')
        xlabel('r/$\sigma$')
        ylabel('v(r)/ kJ/molecule')
        xlim(0,2.5)
        
        subplot(5,1,5)
        plot(sigma_axis[start:],rdf[start:]-rdf[start:]*exp(v_r[start:]/T_r/119.8/kb),'b',sigma_axis,np.zeros(len(rdf)),'k-.')
        xlabel('r/$\sigma$')
        ylabel('c(r)')
        xlim(0,2.5)
        show()        

###############################################################################
###############################################################################
Ensemble='NVT'

#Parameters
sigma=3.405e-10 #LJ sigma parameter
kb=1.38e-23 #boltzmann's constant
eps=119.8*kb #LJ epsilon parameter
N=100 #number of particles
rho_r=0.8 #reduced density
L=(N/rho_r)**(1/3.)*sigma #L=10*sigma
m=39.948*1.6747e-27 #atomic mass
T_r=1.0 #reduced temperature
T=T_r/kb*eps #temperature
rcut=0.5*L #cutoff
dt=1e-14*0.86641345111/2 #timestep
time=dt*20000 #number of snapshots
nout=5 #file outputs every nout steps

#initializes class
MD=Molecular_Dynamics(sigma,eps,N,m,L,T,dt,rcut,nout)

#System description
print 'You are running %s Ensemble'%Ensemble
print 'number of particles =',N
print 'temperature set at',T,'Kelvin'
print 'set rho* =',N/(L**3)*sigma**3
print 'set L* =', L/sigma,'sigma'
print 'set T* =',T_r
print '****************************************************************'
print 'set t* =',(eps/m/sigma**2)**0.5*dt
print 'number of cycles =',int(time/dt)
print 'nout =', nout, 'steps'

############################## Run MD Code ####################################
# **comment out to run post-analysis code to avoid rerunning whole simulation**
#runs the Velocity-Verlet algorithm 

#MD.Velocity_Verlet(time,Ensemble) 


######################## Run Post-Analysis Code ###############################
#runs radial distribution function code
MD.RDF(200,18000)




print 'FINISHED!'



        
        
        
