#!/usr/bin/env python
# coding: utf-8

# In[22]:

#LIBRARY FOR POST-PROCESSING OF DATA OBTAINED FROM UPPASD 
#WORK IN PROGRESS

############################################################################################################################################

#functions for calculating S(T) and ΔS(T) with obtained Cv from Monte Carlo simulations

import numpy as np
import scipy 
from scipy.integrate import simps
import matplotlib.pyplot as plt
import math

def Entropy_with_simulatedCv(material, size, *Hfield):
    
    method="MC"
    
    EntropySim=[] #array of arrays for all different total entropies
    diffS=[] #array of arrays for entropy differences
    
    for field in Hfield: #does the calculations for as many fields as supplied
        
        data=open("thermal."+material+"ncell"+size+"H"+str(field)+".dat") #read data and store data
        lines=data.readlines()[1:]

        y= []
        x= []

        for line in lines:
            lines = [i for i in line.split()]
            x.append(lines[0])
            y.append(lines[4])

        temp=np.loadtxt(x)
        datay=np.loadtxt(y) 
        
        f=datay/datax #creating Cv/T

        Smag=[] #empty array for all entropies

        SmagT=scipy.integrate.trapz(f,temp) #integrating for last temperature
        Smag.insert(0,SmagT) #adding the entropy for last temperature to the array

        for i in range(1, len(temp)-1): #integrating for every temperature except last and adding to the array
            SmagT=scipy.integrate.trapz(f[:-i],temp[:-i])
            Smag.insert(0,SmagT)
    
        SmagSim=np.round(Smag,3) #rounding the values
        EntropySim.append(SmagSim)

    #calculating the difference
    
    for i in range(1,len(EntropySim)):
        diff=EntropySim[i]-EntropySim[0]
        diffS.append(diff)      
        
    return EntropySim, diffS, temp, Hfield, method


############################################################################################################################################


#functions for calculating S(T) and ΔS(T) with Cv=dU/dT 

def Entropy_with_calculatedCv(material,size, *Hfield, **kwargs):
    
    method="MC"
    
    EntropyCal=[] #an array of arrays for all the total entropies
    diffS=[] #an array of arrays for entropy differences
    
    if "rounding" not in kwargs:
        rounding=4
    
        
    for field in Hfield: #does the calculations for as many fields as supplied
        
        data=open("thermal."+material+"ncell"+size+"H"+str(field)+".dat") #read data and store data
        lines=data.readlines()[1:]

        temp= []
        Uenergy= []

        for line in lines:
            lines = [i for i in line.split()]
            temp.append(int(lines[0]))
            Uenergy.append(float(lines[5]))
        
        for i in range(len(Uenergy)): #converting energy to joule
            Uenergy[i]=(Uenergy[i]*13.606e-3)/6.242e+18
        
        Cmag=np.gradient(Uenergy,temp) #creating dU/dT

        #for i in range(len(Cmag)): #converting Cmag to kB
            #Cmag[i]=Cmag[i]/1.381e-23
        
        f=Cmag/temp #creating Cv/T

        Smag=[] #empty array for all entropies
    
        SmagT=scipy.integrate.trapz(f,temp) #integrating for last temperature 
        Smag.insert(0,SmagT*6.022e+23) #adding the entropy for last temperature to the array #converting to J/Kmol
    
        for i in range(1, len(temp)-1): #integrating for every temperature except last and adding to the array
            SmagT=scipy.integrate.trapz(f[:-i],temp[:-i])
            Smag.insert(0,SmagT*6.022e+23) #converting to J/Kmol
         
        if "rounding" not in kwargs:
            SmagCal=np.array(Smag)
            EntropyCal.append(SmagCal)
        
        else:
            SmagCal=np.round(Smag,kwargs["rounding"]) #rounding values
            EntropyCal.append(SmagCal)
        
        
    #calculating the difference
    
    for i in range(1,len(EntropyCal)):
        diff=EntropyCal[i]-EntropyCal[0]
        diffS.append(diff)  
    
    return EntropyCal, diffS, temp, Hfield, method, material


############################################################################################################################################


#function for plotting and calculating differences


def plotEntropy(*data, **kwargs):
    
    
    #plotting entropy
    
    for N in data:
        plt.figure(figsize=(10,6))
        #plotting total entropy

        plt.title("S(T) of " + N[5])
        plt.xlabel("T(K)")
        plt.ylabel("S(J/Kmol)")

        for S in N[0]:      
            
            if N[4]=="MC":    
                plt.plot(N[2][1:],S, marker="o", markersize=5)

            if N[4]=="CE" or N[4]=="CEAFM" :
                plt.plot(N[2],S,marker="o", markersize=5)

        
        if N[4]=="CEAFM":
            
            legend=[]
            legend.append("AFM 0.0H")
            
            for H in N[3]:
                legend.append(H)
            
            plt.legend(legend, title="H in Tesla", loc="lower right")
            
          
        else:
            
            plt.legend(N[3], title="H in Tesla", loc="lower right") 
        
        plt.show()


        #plotting entropy differences
        
        plt.figure(figsize=(10,6))
        plt.title("ΔS(T) of " + N[5])
        plt.xlabel("T(K)")
        plt.ylabel("ΔS(J/Kmol)")
        
        if N[4]!="CEAFM":
            
            plt.axhline(y=0, linestyle="dashed") #creating a dashed line at zero difference
            
        
        
        for S in N[1]:
            
            if N[4]=="MC":
                plt.plot(N[2][1:], S, marker="o",markersize=5)
                location="lower right"

            if N[4]=="CE" or N[4]=="CEAFM":
                plt.plot(N[2],S,marker="o", markersize=5)
                location="lower left"

        plt.legend(N[3], title="H in Tesla", loc=location) 
        
        plt.show()


############################################################################################################################################

#function for plotting magnon spectra, both adiabatic and correlatedd

def plotMagnonSpectrum(material, size, *Hfield, **kwargs):
    
    both=False
    
    
    if kwargs["spectrum"]=="both": #checking which one should be done
        
        both=True
    
    
    if kwargs["spectrum"]=="CMS" or both==True: #doing correlataed
        
        H_list=[]

        for i in range(len(Hfield)):

            H_list.append(Hfield[i])

        for field in H_list: #does the calculations for as many fields as supplied

            data=np.loadtxt("sqw."+material+"ncell"+size+"H"+str(field)+"."+str(kwargs["yscale"])+".mat") #read data and store data

            rows=data.shape[0]
            columns=data.shape[1]
            
            if "kpath" in kwargs: #setting high symmetry points if given
    
                kpoints=[]
                labels=[]
                for j in range(0,len(kwargs["kpath"]),2):
                    for i in range(0,len(kwargs["kpath"][j])):
                        kpoints.append(kwargs["kpath"][j][i])
                        labels.append(kwargs["kpath"][j+1])

                plt.xticks(kpoints,labels)
                
                for point in kpoints:
                    plt.axvline(x=point,color="w",linestyle="dashed",linewidth=0.5, alpha=1)

                                   
            if "kpath" not in kwargs:
                plt.xticks([], [])                      
                                      
            plt.xlabel("K")
            plt.ylabel("Energy (meV)")
            
            plt.imshow(data, origin="lower", extent=[0, columns, 0,rows*kwargs["yscale"]], aspect="auto", cmap="hot")
            
            if both==False:
                plt.show()
        
        
    if kwargs["spectrum"]=="AMS" or both==True: #doing adiabatic
    
        H_list=[]

        for i in range(len(Hfield)):

            H_list.append(Hfield[i])

        for field in H_list: #does the calculations for as many fields as supplied

            data=np.loadtxt("ams."+material+"ncell"+size+"H"+str(field)+".out") #read data and store data

            col=np.shape(data)[1] #getting numbers of columns
                                      
            if "kpath" in kwargs: #setting high symmetry points if given
    
                kpoints=[]
                labels=[]
                for j in range(0,len(kwargs["kpath"]),2):
                    for i in range(0,len(kwargs["kpath"][j])):
                        kpoints.append(kwargs["kpath"][j][i])
                        labels.append(kwargs["kpath"][j+1])

                plt.xticks(kpoints,labels)
                
                for point in kpoints:
                    plt.axvline(x=point,color="w",linestyle="dashed",linewidth=0.5, alpha=1)

                                   
            if "kpath" not in kwargs:
                plt.xticks([], [])                     
            
            plt.xlabel("K")
            plt.ylabel("Energy (meV)")
            
            for i in range(col-2): #plotting for every mode
                plt.plot(data[:,0],data[:,i+1], color="b", alpha=0.6)
            plt.ylim(ymin=0)
            plt.xlim(xmin=0)
            plt.show()

############################################################################################################################################


#function for plotting adiabatic magnon DOS

def plotAMSmDOS(material, size, *Hfield):
    
    H_list=[]
    
    for i in range(len(Hfield)):
        
        H_list.append(Hfield[i])
        
    for field in H_list: #does the calculations for as many fields as supplied

        magdos=open("magdos." + material + "ncell"+size+"H"+ str(field) + ".out")
        lines=magdos.readlines()

        energy= []
        DOS= [] #state occupation



        for line in lines:
            lines = [i for i in line.split()]
            energy.append(lines[0])
            DOS.append(lines[1])

            #Changing to numpy arrays
            dataEnergy=np.loadtxt(energy) 
            dataDOS=np.loadtxt(DOS)   

        for j in range(0,len(energy)): #chaning to right type
            energy[j]=float(energy[j])
            DOS[j]=float(DOS[j])


        print("Integration of "+material+"s mDOS with H=" + str(field) + "T equals:") 
        print(scipy.integrate.trapz(DOS,energy))

        plt.plot(energy, DOS)
        
        plt.show()
    


############################################################################################################################################


#function for calculating entropy with canonical ensemble

def Entropy_with_CanonicalEnsemble(material ,size, temp, *Hfield, **kwargs):
    
    method="CE"
    
    EntropyCE=[] #an array of arrays for all different total entropies
    diffS=[] #an array of arrays for all entropy differences
    
    
    
    if kwargs["transition"]=="AFM-FM": #for antiferromagnetic to ferromagnetic phase transition
            
            dataAFM=np.loadtxt("magdos." + "AFM"+material +"ncell"+size+"H0.0" + ".out") #loading data from file
            
            
            
            energy=dataAFM[:,0]
            DOS=dataAFM[:,1] #state occupation

            S=[] #empty array for entropy values

            for j in range(len(energy)): #changing a zero energy to something really small to avoid errors
                    if energy[j]==0:
                        energy[j]=1e-10

            for i in range(len(temp)): #start of the calculations


                exp=(energy/(8.617e-2*temp[i])) #the exponent in the Bose-Einstein distribution
                n=np.exp(-exp)/(1-np.exp(-exp)) #Bose-Einstein distribution created, not in the litterature way

                ln1=[] #empty arrays for the logarithms
                ln2=[]

                #the logarithms
                for i in range(len(n)): 
                    ln1.append(math.log(1+n[i]))

                for i in range(len(n)):
                    ln2.append(math.log(n[i]))

                right=(1+n)*ln1-(n*ln2) #the parantese in the integral
                inner=DOS*right #the whole integral

                integral=scipy.integrate.trapz(inner,energy)
                SCE=8.316382*integral    #multiplying the integration with kB and then with avogadros constant to get J/Kmol

                S.append(SCE) #adding to the list

            SmagH=np.round(S,kwargs["rounding"])
            EntropyCE.append(SmagH)
            
            material="FM"+material
            
            method="CEAFM"
           
    
    for field in Hfield: #does the calculations for as many fields as supplied
        
        
        data=np.loadtxt("magdos." + material +"ncell"+size+"H" + str(field) + ".out") #loading data from file
                
        energy=data[:,0]
        DOS=data[:,1] #state occupation

        S=[] #empty array for entropy values

        for j in range(len(energy)): #changing a zero energy to something really small to avoid errors
                if energy[j]==0:
                    energy[j]=1e-10

        for i in range(len(temp)): #start of the calculations


            exp=(energy/(8.617e-2*temp[i])) #the exponent in the Bose-Einstein distribution
            n=np.exp(-exp)/(1-np.exp(-exp)) #Bose-Einstein distribution created, not in the litterature way
           
            ln1=[] #empty arrays for the logarithms
            ln2=[]

            #the logarithms
            for i in range(len(n)): 
                ln1.append(math.log(1+n[i]))

            for i in range(len(n)):
                ln2.append(math.log(n[i]))

            right=(1+n)*ln1-(n*ln2) #the parantese in the integral
            inner=DOS*right #the whole integral

            integral=scipy.integrate.trapz(inner,energy)
            SCE=8.316382*integral    #multiplying the integration with kB and then with avogadros constant to get J/Kmol

            S.append(SCE) #adding to the list
        
        if "rounding" not in kwargs:
            SmagH=np.array(S)
            EntropyCE.append(SmagH)
        
        else:
            SmagH=np.round(S,kwargs["rounding"])
            EntropyCE.append(SmagH)
          
     
    #calculating differences
       
    
    for i in range(1,len(EntropyCE)):
        diff=EntropyCE[i]-EntropyCE[0]
        diffS.append(diff) 
    

    return EntropyCE, diffS, temp, Hfield, method, material


############################################################################################################################################