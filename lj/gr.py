# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 17:43:27 2015

file output is in the energy section
pep8


@author: Yushan Zhang
"""
import math
import random
import copy
import matplotlib.pyplot as plot
from atom import Atom


N=864

sigma=3.4e-10
R=2.25*sigma
L=10.229*sigma
k=1.3806488e-23
Na = 6.022e23

epsilon=k*120.0
M=(39.95/Na)*(10**-3)

V=L**3

dt=1e-14

T0=90
v2=0.0

dr=sigma/10

'''Run the simulation'''
def main(total,zero):
    atoms=initialization()
    temperatures=[]
    velaclist = [] 
    energylist=[]
    PElist=[]#############################
    KElist=[]####################################
#    velacfinit=0
    originalAtoms=copy.deepcopy(atoms)    
    for t in range(0,total):
        atoms,temperature,velaclist,originalAtoms=simulation(atoms,t,temperatures, velaclist,originalAtoms,zero)
        if t>=zero:       
            energylist,PElist,KElist=total_energy(atoms,energylist,PElist,KElist)################3#########
        writeXYZ(atoms,t,'Ar.xyz')
    atoms=output(atoms,velaclist,energylist,total,PElist,KElist,zero)
    
#    return atoms,temperatures,velaclist,energylist

def simulation(atoms,step,temperatures, velaclist,originalAtoms,zero):
    atoms=update_force(atoms)
    atoms=integration(atoms)
    T,temperatures=update_temperature(atoms,temperatures)
    atoms=reset_force(atoms)
    if step==(zero+1):
        originalAtoms=copy.deepcopy(atoms)
    if step>=(zero+1):#########################################
        velaclist=velocityAutocorrelation(step,atoms,velaclist,originalAtoms,zero)
    print("Current System Temperature: " + str(T))    
    print("~~~~~~~~~~~-COMPLETED STEP " + str(step+1) + "~~~~~~~~~~~~~~")
    if step>100 and step<118:
        atoms=scale_temp(atoms,T)
    return atoms,temperatures,velaclist,originalAtoms

def initialization():
    atoms=particle_generation()
    atoms=assign_position(atoms)
    atoms=boltzmann_distribution(atoms)
    atoms=correct_momenta(atoms)    
    return atoms

'''Output functions'''
def output(atoms,velaclist,energylist,nSteps,PElist,KElist,zero):
    atom_counts=pairDistributionFunction(atoms)
    plotRDF(atom_counts)
    plotVAC(nSteps,velaclist,zero)
    plotEnergy(nSteps,energylist,zero)
    plotPE(nSteps,PElist,zero)
    plotKE(nSteps,KElist,zero)
    return atoms 

    
def pairDistributionFunction(atoms):
    atom_counts = [0]*50
        
    print("Generating rdf..."),
    for atom1 in range(0, N-1):
        for atom2 in range(1, N):
            dx = atoms[atom1].x - atoms[atom2].x
            dy = atoms[atom1].y - atoms[atom2].y
            dz = atoms[atom1].z - atoms[atom2].z

            dx -= L*round(dx/L)
            dy -= L*round(dy/L)
            dz -= L*round(dz/L)
                
            r2 = dx*dx + dy*dy + dz*dz
            r = math.sqrt(r2)

            for radius in range(0, 50):
                if (r < ((radius+1)*dr)) and (r > radius*dr):
                    atom_counts[radius] += 1
    for radius in range(1, 50):
        atom_counts[radius] *= (V/N**2)/(4*math.pi*((radius*dr)**2)*dr)
    print("done.")    
    return(atom_counts)     

def velocityAutocorrelation(step,atoms,velaclist,originalAtoms,zero):
    velacf = 0##############################
    velacfinit=0
    vx = 0
    vy = 0
    vz = 0
    if step == (zero+1):
        print('Atoms copied')
        for atom in range(0, N):
            vx += originalAtoms[atom].vx * atoms[atom].vx
            vy += originalAtoms[atom].vy * atoms[atom].vy
            vz += originalAtoms[atom].vz * atoms[atom].vz
        velacfinit += vx + vy + vz
        velacfinit /= N
        velaclist.append(velacfinit)
        print (velacfinit)
    else:   
        for atom in range(0, N):
            vx += originalAtoms[atom].vx * atoms[atom].vx
            vy += originalAtoms[atom].vy * atoms[atom].vy
            vz += originalAtoms[atom].vz * atoms[atom].vz
        velacf += vx + vy + vz
#        velacf /= N*velacfinit
        velacf /= N*velaclist[0]
        velaclist.append(velacf)
        velacf = 0
    return velaclist

def writeXYZ(atoms,step,filename):########################
    """Writes positional data to a .xyz file"""
    scale=1e-10
    with open(filename, "a") as output:
        output.write('%d \n \n' % len(atoms))
        for atom in atoms:
            output.write("Ar %s %s %s\n" % (atom.x/scale, atom.y/scale, atom.z/scale))   
    return

def total_energy(atoms,energylist,PElist,KElist):
    PE_const=4*epsilon
    KE_const=M/2
    v2=0
    for atom in range(0, N):
        KE_instant=0.0
        KE_atom=atoms[atom]
        v2=KE_atom.vx**2+KE_atom.vy**2+KE_atom.vz**2
        KE_instant+=v2
    KE=(KE_instant*KE_const)        
    for atom1 in range (0,N-1):
        for atom2 in range(atom1+1,N):
            PE_instant=0.0
            dx = atoms[atom1].x - atoms[atom2].x
            dy = atoms[atom1].y - atoms[atom2].y
            dz = atoms[atom1].z - atoms[atom2].z
            r=math.sqrt((dx)**2+(dy)**2+(dz)**2)
            
            PE_instant+=(sigma/r)**12-(sigma/r)**6 #potential energy equation      
    PE=(PE_instant*PE_const)
    energylist.append(KE+PE) 
    PElist.append(PE)
    KElist.append(KE)
    return energylist,PElist,KElist

def plotRDF(rdf):
    radiuslist=[]    
    for radius in range(0,50):
        radiuslist.append(radius*dr/sigma)
    plot.figure()
    plot.plot(radiuslist,rdf)
    plot.xlabel('Distance from certain atom   /m')
    plot.ylabel('g(r)')
    plot.show()

def plotVAC(nSteps,vac,zero):
    vac[0] = 1
    timelist=[]
    for time in range(0, len(vac)):
        timelist.append(float(time+zero) * dt)
    plot.figure()
    plot.plot(timelist, vac)
    plot.xlabel('Time   /s')
    plot.ylabel('Velocity auto-correlation')
    plot.show()

def plotEnergy(nSteps,energylist,zero):
    timelist=[]
    for time in range(0, len(energylist)):
        timelist.append(float(time+zero) * dt)
    plot.figure()
    plot.plot(timelist, energylist)
    plot.xlabel('Time   /s')
    plot.ylabel('Total energy /arbitary unit')
    plot.axis([timelist[0], timelist[-1], -1*10**(-19), 1*10**(-19)])
    plot.show() 

def plotPE(nSteps,PElist,zero):
    timelist=[]
    for time in range(0, len(PElist)):
        timelist.append(float(time+zero) * dt)
    plot.figure()
    plot.plot(timelist, PElist)
    plot.xlabel('Time   /s')
    plot.ylabel('Potential energy /arbitary unit')
    plot.axis([timelist[0], timelist[-1], -1*10**(-19), 1*10**(-19)])
    plot.show()

def plotKE(nSteps,KElist,zero):
    timelist=[]
    for time in range(0, len(KElist)):
        timelist.append(float(time+zero) * dt)
    plot.figure()
    plot.plot(timelist, KElist)
    plot.xlabel('Time   /s')
    plot.ylabel('Kinetic energy /arbitary unit')
    plot.axis([timelist[0], timelist[-1], -1*10**(-19), 1*10**(-19)])
    plot.show()






'''Simulation functions'''


def update_force(atoms):
    for atom1 in range(0, N-1):
        for atom2 in range(atom1+1, N):
#            print(atom1,atom2)
            atoms=calculateForce(atoms,atom1,atom2)
    for atom in range(0, N):
        atoms[atom].fx *= 48*epsilon
        atoms[atom].fy *= 48*epsilon
        atoms[atom].fz *= 48*epsilon       
    return atoms

def calculateForce(atoms,atom1,atom2):
    dx = atoms[atom1].x - atoms[atom2].x
    dy = atoms[atom1].y - atoms[atom2].y
    dz = atoms[atom1].z - atoms[atom2].z
    
    dx -= L*round(dx/L)
    dy -= L*round(dy/L)
    dz -= L*round(dz/L)
    
    r2 = dx*dx + dy*dy + dz*dz
    
    if r2 < (R**2):
        fr2 = (sigma**2)/r2
        fr6 = fr2**3
        force = fr6*(fr6 - 0.5)/r2
            
        forcex = force*dx
        forcey = force*dy
        forcez = force*dz
            
        atoms[atom1].fx += forcex
        atoms[atom2].fx -= forcex
        atoms[atom1].fy += forcey
        atoms[atom2].fy -= forcey
        atoms[atom1].fz += forcez
        atoms[atom2].fz -= forcez
    return atoms
    
    
def integration(atoms):
    for atom in range(0,N):
        vx0=atoms[atom].vx
        vy0=atoms[atom].vy
        vz0=atoms[atom].vz
        
        atoms[atom].vx += (atoms[atom].fx/M)*dt
        atoms[atom].vy += (atoms[atom].fy/M)*dt
        atoms[atom].vz += (atoms[atom].fz/M)*dt

            # Update positions
        newX = atoms[atom].x + (atoms[atom].vx*dt+vx0*dt)/2
        newY = atoms[atom].y + (atoms[atom].vy*dt+vy0*dt)/2
        newZ = atoms[atom].z + (atoms[atom].vz*dt+vz0*dt)/2  
        if newX < 0:
            atoms[atom].x = newX +  L
        elif newX > L:
            atoms[atom].x = newX - L
        else:
            atoms[atom].x = newX
            
        if newY < 0:
            atoms[atom].y = newY +  L
        elif newY >  L:
            atoms[atom].y = newY - L
        else:
            atoms[atom].y = newY
                
        if newZ < 0:
            atoms[atom].z = newZ +  L
        elif newZ > L:
            atoms[atom].z = newZ -  L
        else:
            atoms[atom].z = newZ    
    return atoms

def update_temperature(atoms,temperatures):
    sumv2 = 0
    for atom in atoms:
        sumv2 += atom.vx**2 + atom.vy**2 + atom.vz**2
    T = (M/(3*N*k))*sumv2
    temperatures.append(T)
    return T,temperatures

def reset_force(atoms):
    for atom in range(0, N):
        atoms[atom].fx = 0
        atoms[atom].fy = 0
        atoms[atom].fz = 0
    return atoms    
    
def scale_temp(atoms,T):
    if T > 100.0 or T < 80.0:
        print("Rescaling velocities...")
        for atom in range(0, N):
            atoms[atom].vx *= math.sqrt(T0/T)
            atoms[atom].vy *= math.sqrt(T0/T)
            atoms[atom].vz *= math.sqrt(T0/T)
    return atoms 



'''Initialization functions'''
    
def particle_generation():
    atoms=[]
    for i in range(0,N):
        atoms.append(Atom())
    return atoms
    
def assign_position(atoms):
    n=int(math.ceil(N**(1.0/3.0)))
    particle=0
    for x in range(0,n):
        for y in range(0,n):
            for z in range(0,n):
                if(particle < N):
                    atoms[particle].x=x*sigma
                    atoms[particle].y=y*sigma
                    atoms[particle].z=z*sigma    
                    particle+=1
    return atoms
    
def boltzmann_distribution(atoms):
    normDist = []
    scaling_factor = math.sqrt(k*T0/M)
    for i in range(0, 3*N):
            normDist.append(random.gauss(0,1))      
    for number in range(0, 3*N):
        normDist[number] = normDist[number]*scaling_factor          
    for atom in range(0, N):
        atoms[atom].vx = normDist[atom*3]
        atoms[atom].vy = normDist[atom*3+1]
        atoms[atom].vz = normDist[atom*3+2]        
    return atoms
    
def correct_momenta(atoms):
    sumvx = 0
    sumvy = 0
    sumvz = 0        
    for atom in range(0, N):
        sumvx += atoms[atom].vx
        sumvy += atoms[atom].vy
        sumvz += atoms[atom].vz        
    for atom in range(0, N):
        atoms[atom].vx -= sumvx/N
        atoms[atom].vy -= sumvy/N
        atoms[atom].vz -= sumvz/N
    return atoms

