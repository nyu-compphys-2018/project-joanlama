# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 09:48:27 2018

@author: Joan
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 19:01:49 2018

@author: Joan
"""


import numpy as np
import matplotlib.pyplot as plt
#import time
#from random import randrange


# I download the data and store it k_mod and p_k

#First Power

power_spectrum = np.array(np.loadtxt( 'power_LD_z0.dat'), dtype=float )


k_mod = power_spectrum [:,0]

p_k_data = power_spectrum [:,1]

#I define the parameters of my grid

L = 1000

N = 150

p_k = (p_k_data)*(((2*np.pi)/L)**(-3))

# FIT: I define the function that fits my |k|^2 with the P(k) of the data


def fit(x_data,y_data,x):
    y = 0     
       
    for i in range(len(x_data)):
        if x_data[i] >= x:
            y = ((y_data[i]-y_data[i-1])/(x_data[i]-x_data[i-1]))*(x-x_data[i-1]) + y_data[i-1]
            break
    return y


# POWER SPECTRUM: I define the function that gives me the 3x3x3 tensor of sqrt(P_k/2)

def power_spec(k_mod,p_k,kx,ky,kz,phi_fourier,h):
    
    delta_k = np.zeros((len(kx),len(ky),len(kz)),complex)


    for i in range(len(kx)):
        for j in range(len(ky)):
            for k in range(len(kz)):
                
                if i==(int(len(kx)//2)) and j==(int(len(kx)//2)) and k==(int(len(kx)//2)):
                    delta_k[(i-(int(len(kx)//2)))][(j-(int(len(ky)//2)))][(k-(int(len(kz)//2)))] = 0
                else:
                    delta_k[(i-(int(len(kx)//2)))][(j-(int(len(ky)//2)))][(k-(int(len(kz)//2)))] = (-1)*phi_fourier[i][j][k]*np.sqrt(fit(k_mod,p_k,np.sqrt(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k]))/2)/(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k])
                    
                    
    return delta_k


# RANDOM: I construct the tensors A and B which are Gaussian


def random(N,kx,ky,kz):
    A = np.zeros((len(kx),len(ky),len(kz)),float)
    B = np.zeros((len(kx),len(ky),len(kz)),float)
    cos = np.zeros((len(kx),len(ky),len(kz)),float)
    sin = np.zeros((len(kx),len(ky),len(kz)),float)
    
    r = np.zeros((len(kx),len(ky),len(kz)),float)
    one = np.zeros((len(kx),len(ky),len(kz)),float)
    
    u = 0
    v = 0
    
    u_1 = 0
    v_1 = 0
    
    u_2 = 0
    v_2 = 0
    
    u_3 = 0
    v_3 = 0
    

    for i in range(int(N/2),len(kx)):
        for j in range(int(N/2),len(ky)):
            for k in range(int(N/2),len(kz)):
                             
                u = np.random.uniform()
                v = np.random.uniform()
                
                u_1 = np.random.uniform()
                v_1 = np.random.uniform()
                
                u_2 = np.random.uniform()
                v_2 = np.random.uniform()
                
                u_3 = np.random.uniform()
                v_3 = np.random.uniform()
                
                
                cos[i][j][k] = np.cos((2*np.pi)*u)
                
                cos[-i+int(N)][-j+int(N)][-k+int(N)] = np.cos((2*np.pi)*u)
                
                
                sin[i][j][k] = np.sin((2*np.pi)*u)
                
                sin[-i+int(N)][-j+int(N)][-k+int(N)] = (-1)*np.sin((2*np.pi)*u)
                
                
                r[i][j][k] = np.sqrt(-2*np.log(1-v))
                
                r[-i+int(N)][-j+int(N)][-k+int(N)] = np.sqrt(-2*np.log(1-v))
                
                               
                #Until here it works as expected!
                
                #This lines are to create the middle terms
                
                cos[i-int(N/2)][j][k] = np.cos((2*np.pi)*u_1)
                
                cos[-(i-int(N/2))+int(N)][-j+int(N)][-k+int(N)] = np.cos((2*np.pi)*u_1)
                
                
                sin[i-int(N/2)][j][k] = np.sin((2*np.pi)*u_1)
                
                sin[-(i-int(N/2))+int(N)][-j+int(N)][-k+int(N)] = - np.sin((2*np.pi)*u_1)
            
                
                r[i-int(N/2)][j][k] = np.sqrt(-2*np.log(1-v_1))
                
                r[-(i-int(N/2))+int(N)][-j+int(N)][-k+int(N)] = np.sqrt(-2*np.log(1-v_1))
                
                #########################################################################                
                
                cos[i][j-int(N/2)][k] = np.cos((2*np.pi)*u_2)
                
                cos[-i+int(N)][-(j-int(N/2))+int(N)][-k+int(N)] = np.cos((2*np.pi)*u_2)
                
                
                sin[i][j-int(N/2)][k] = np.sin((2*np.pi)*u_2)
                
                sin[-i+int(N)][-(j-int(N/2))+int(N)][-k+int(N)] = - np.sin((2*np.pi)*u_2)
                
                
                r[i][j-int(N/2)][k] = np.sqrt(-2*np.log(1-v_2))
                
                r[-i+int(N)][-(j-int(N/2))+int(N)][-k+int(N)] = np.sqrt(-2*np.log(1-v_2))
                
                #########################################################################                
                
                
                cos[i][j][k-int(N/2)] = np.cos((2*np.pi)*u_3)
                
                cos[-i+int(N)][-j+int(N)][-(k-int(N/2))+int(N)] = np.cos((2*np.pi)*u_3)
                
            
                sin[i][j][k-int(N/2)] = np.sin((2*np.pi)*u_3)
                
                sin[-i+int(N)][-j+int(N)][-(k-int(N/2))+int(N)] = - np.sin((2*np.pi)*u_3)
                
                
                r[i][j][k-int(N/2)] = np.sqrt(-2*np.log(1-v_3))
                
                r[-i+int(N)][-j+int(N)][-(k-int(N/2))+int(N)] = np.sqrt(-2*np.log(1-v_3))
                
    
                
    
    
    A = r*cos

    B = r*sin
    
    
    
    #one = cos*cos + sin*sin
    
    

    return A,B



#First Derivatives in Fourier Space

def gradient_fourier(phi,k_x,k_y,k_z):
    
    kx = np.roll(k_x,int(len(k_x))//2+1)
    
    ky = np.roll(k_y,int(len(k_y))//2+1)

    kz = np.roll(k_z,int(len(k_z))//2+1)

    
    psi_x = np.zeros((len(kx),len(ky),len(kz)),complex)
    psi_y = np.zeros((len(kx),len(ky),len(kz)),complex)
    psi_z = np.zeros((len(kx),len(ky),len(kz)),complex)

    for i in range(len(kx)):
        for j in range(len(ky)):
            for k in range(len(kz)):
                psi_x[i][j][k] = phi[i][j][k]*(1j)*kx[i]
                psi_y[i][j][k] = phi[i][j][k]*(1j)*ky[j]
                psi_z[i][j][k] = phi[i][j][k]*(1j)*kz[k]
    
    return psi_x,psi_y,psi_z



#Second Derivatives in Fourier Space


def sec_derivatives_fourier(phi,k_x,k_y,k_z):
    
    kx = np.roll(k_x,int(len(k_x))//2+1)
    
    ky = np.roll(k_y,int(len(k_y))//2+1)

    kz = np.roll(k_z,int(len(k_z))//2+1)

    partial_x_2 = np.zeros((len(kx),len(ky),len(kz)),complex)
    
    partial_y_2 = np.zeros((len(kx),len(ky),len(kz)),complex)
    
    partial_z_2 = np.zeros((len(kx),len(ky),len(kz)),complex)
    
    partial_x_y = np.zeros((len(kx),len(ky),len(kz)),complex)
    
    partial_x_z = np.zeros((len(kx),len(ky),len(kz)),complex)
    
    partial_y_z = np.zeros((len(kx),len(ky),len(kz)),complex)
    

    for i in range(len(kx)):
        for j in range(len(ky)):
            for k in range(len(kz)):
                
                partial_x_2[i][j][k] = (-1)*phi[i][j][k]*kx[i]*kx[i]
                partial_y_2[i][j][k] = (-1)*phi[i][j][k]*ky[j]*ky[j]
                partial_z_2[i][j][k] = (-1)*phi[i][j][k]*kz[k]*kz[k]
                
                partial_x_y[i][j][k] = (-1)*phi[i][j][k]*kx[i]*ky[j]
                partial_x_z[i][j][k] = (-1)*phi[i][j][k]*kx[i]*kz[k]
                partial_y_z[i][j][k] = (-1)*phi[i][j][k]*ky[j]*kz[k]
                
                
                
    
    return partial_x_2, partial_y_2, partial_z_2, partial_x_y, partial_x_z, partial_y_z



def rho_2_real(partial_x_2, partial_y_2, partial_z_2, partial_x_y, partial_x_z, partial_y_z,N):
    
    rho = np.zeros((N,N,N),complex)
    
    

    for i in range(N):
        for j in range(N):
            for k in range(N):
                
                rho[i][j][k] = (partial_x_2[i][j][k]*partial_y_2[i][j][k])\
                               +(partial_x_2[i][j][k]*partial_z_2[i][j][k])\
                               +(partial_y_2[i][j][k]*partial_z_2[i][j][k])\
                               -(partial_x_y[i][j][k]*partial_x_y[i][j][k])\
                               -(partial_x_z[i][j][k]*partial_x_z[i][j][k])\
                               -(partial_y_z[i][j][k]*partial_y_z[i][j][k])
                               
                
    
    return rho



def poisson(rho_2,k_x,k_y,k_z):
    kx = np.roll(k_x,int(len(k_x))//2+1)
    
    ky = np.roll(k_y,int(len(k_y))//2+1)

    kz = np.roll(k_z,int(len(k_z))//2+1)
    
    phi_2 = np.zeros((len(kx),len(ky),len(kz)),complex)

    
    for i in range(len(kx)):
        for j in range(len(ky)):
            for k in range(len(kz)):
                
                
                if i==(0) and j==(0) and k==(0):
                    phi_2[i][j][k] = 0
                else:
                    
                    phi_2[i][j][k] = (-1)*(rho_2[i][j][k])/(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k])
                    
    
    return phi_2



#I create the Grid for Plotting the movement of the particles

def get_grid(psi1,psi2,psi_2_x,psi_2_y,N,L,k):
    
    
    
    h = np.arange(-N/2,(N+1)/2)
    
    h_l = np.linspace(-L/2,L/2,len(h))
    
    
    mesh = [(x, y) for x in h for y in h]
    
    mesh_l = [(x, y) for x in h_l for y in h_l]
    
    
    grid = np.array(mesh)
    
    grid_l_2 = np.array(mesh_l)
    
#    grid_original = np.array(mesh_l)

    
    m = 0
    n = 0
    
    for i in range (len(h)**2):
        m = int((N/2) + grid[i][0])
        n = int((N/2) + grid[i][1])
        
    
        
        grid_l_2[i][0] = grid_l_2[i][0] - psi1[m][n][k] - psi_2_x[m][n][k]
        
        grid_l_2[i][1] = grid_l_2[i][1] - psi2[m][n][k] - psi_2_y[m][n][k]
    
        
    #I only return the final positions    
    return grid, grid_l_2


###### THE FUNCTIONS THAT COME NOW WERE USEFUL TO TEST THE CODE, I DON'T USE THEM ANYMORE ####

# POWER SPECTRUM TEST: I define the function that gives me the 3x3x3 tensor of sqrt(P_k/2)

def power_test(k_mod,p_k,kx,ky,kz):
    delta_k = np.zeros((len(kx),len(ky),len(kz)))

    for i in range(len(kx)):
        for j in range(len(ky)):
            for k in range(len(kz)):
                
                delta_k[(i-(int(len(kx)//2)))][(j-(int(len(ky)//2)))][(k-(int(len(kz)//2)))] = fit(k_mod,p_k,np.sqrt(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k]))
    
    
    return delta_k



def get_pk(phi,k_x,k_y,k_z):
    
    kx = np.roll(k_x,int(len(k_x))//2+1)
    
    ky = np.roll(k_y,int(len(k_y))//2+1)

    kz = np.roll(k_z,int(len(k_z))//2+1)
    
    
    k_module = []
    pk = []
    
    
    for i in range(len(kx)):
        for j in range(len(ky)):
            for k in range(len(kz)):
                k_module.append(np.sqrt(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k]))
                pk.append(phi[i][j][k])
                            
    
    matrix = np.zeros((len(kx)**3,2))
    for i in range(len(kx)**3):
        matrix[i][0] = k_module[i]
        matrix[i][1] = pk[i]
    matrix=matrix[np.argsort(matrix[:,0])]
    
    num = 1
    for i in range(1, len(kx)**3):
        #if matrix[i][0] != matrix[i-1][0]:
        if (matrix[i][0] >= matrix[i-1][0] - 0.00001) and (matrix[i][0] <= matrix[i-1][0] + 0.00001):
            santa = 12/25
        else:
            num += 1
    
    average = np.zeros((num,2))
    
    n = 0
    b = 1
    average[0][0] = matrix[0][0]
    average[0][1] = matrix[0][1]
    for i in range(1, len(kx)**3):
        if (matrix[i][0] >= matrix[i-1][0] - 0.00001) and (matrix[i][0] <= matrix[i-1][0] + 0.00001):
            #average[n][0] = matrix[i][0]
            average[n][1] += matrix[i][1]
            b += 1
        else:
            average[n][1] = average[n][1]/b
            n += 1
            b = 1
            average[n][0] = matrix[i][0]
            average[n][1] = matrix[i][1]

    return average


########## I run my code defining all the variables #############

# I define the -N/2< m_i <N/2
mx = np.arange(-N/2,(N+1)/2)

my = np.arange(-N/2,(N+1)/2)

mz = np.arange(-N/2,(N+1)/2)


#I define the wave vectors

kx = 2*(np.pi)/L*mx

ky = 2*(np.pi)/L*my

kz = 2*(np.pi)/L*mz


# I get phi_1 in Fourier space ---> This means solving Poisson Eqn.

A,B = random(N,kx,ky,kz)

phi_fourier = A + +1j*B


#We solve the first Poisson Eqn. in Fourier Space

phi_1_k = power_spec(k_mod,p_k,kx,ky,kz,phi_fourier,(L/N))


#DISPLACEMENT FIELDS TO FIRST ORDER IN LPT

#We calculate The Field psi_1 as the Gradient of phi_1

#I take derivatives in Fourier Space

psi_1_x_f, psi_1_y_f, psi_1_z_f = gradient_fourier(phi_1_k,kx,ky,kz)

#We go to Real Space

psi_1_x = np.fft.ifftn(psi_1_x_f)

psi_1_y = np.fft.ifftn(psi_1_y_f)

psi_1_z = np.fft.ifftn(psi_1_z_f)



#DISPLACEMENT FIELDS TO SECOND ORDER IN LPT

#Now we define "rho_2" which is a tensor composed of partial derivatives of phi_1

#This will act as a source for the Second Poisson Eqn.

#I take derivatives in Fourier Space


f_partial_x_2, f_partial_y_2, f_partial_z_2, f_partial_x_y, f_partial_x_z, f_partial_y_z = sec_derivatives_fourier(phi_1_k,kx,ky,kz)

#I get in Real Space the second derivatives

partial_x_2 = np.fft.ifftn(f_partial_x_2)

partial_y_2 = np.fft.ifftn(f_partial_y_2)

partial_z_2 = np.fft.ifftn(f_partial_z_2)


#Cross second derivatives

partial_x_y = np.fft.ifftn(f_partial_x_y)

partial_x_z = np.fft.ifftn(f_partial_x_z)

partial_y_z = np.fft.ifftn(f_partial_y_z)


#I construct the rho_2 in Real Space

rho_2_r = rho_2_real(partial_x_2, partial_y_2, partial_z_2, partial_x_y, partial_x_z, partial_y_z,len(mx))

#I get rho_2 in Fourier Space

rho_2_k = np.fft.fftn(rho_2_r) 


#I solve the Second Poisson Eqn. in Fourier Space

phi_2_k = poisson(rho_2_k,kx,ky,kz)


#We calculate The Field psi_2 as the Gradient of phi_2 in Fourier Space


psi_2_x_f, psi_2_y_f, psi_2_z_f = gradient_fourier(phi_2_k,kx,ky,kz)

#We go to Real Space

psi_2_x = (3/14)*np.fft.ifftn(psi_2_x_f)

psi_2_y = (3/14)*np.fft.ifftn(psi_2_y_f)

psi_2_z = (3/14)*np.fft.ifftn(psi_2_z_f)


#Grid showing the position of the Particles

original, grid_2 = get_grid(psi_1_x,psi_1_y,psi_2_x,psi_2_y,N,L,0)


plt.figure()
plt.hist2d(original[:,0], original[:,1], bins=(400, 400), cmap=plt.cm.jet)
plt.colorbar()
plt.title(r'$Initial \, Positions$', fontsize=30)       
plt.tick_params(axis='both', which='major', labelsize=20)
plt.show()



plt.figure()
plt.hist2d(grid_2[:,0], grid_2[:,1], bins=(500, 500), cmap=plt.cm.jet)
plt.title(r'$Final \, Positions \, 2nd \,LPT$', fontsize=30)       
plt.tick_params(axis='both', which='major', labelsize=20)
plt.colorbar()
plt.show()

