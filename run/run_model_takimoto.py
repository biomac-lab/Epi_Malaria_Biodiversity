#Import libraries
import numpy as np 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.stats import truncnorm
import pandas as pd
import os, glob

#Set the seed
np.random.seed(42)

results_path = 'results'
params_path = 'params'


#Define the parameters of the model
num_species = 4
num_scenarios = 20
mu_H = 0.01 #Número real 
lambda_H = np.random.uniform(1000/3650, 1000/1500) #Número real
lambda_V = np.random.uniform(5*10**3/3650, 5*10**3/1500,size=num_species)
gamma = 0.2 #Número real
tau_HV = np.random.rand(num_species) #Lista de números reales donde cada posición representa la especie i del vector
tau_VH = np.random.rand(num_species)
mu_V = np.random.rand(num_species)
sigma_V = np.random.rand(num_species)
sigma_H = np.random.rand(num_species)
mat_pi = np.random.rand(num_species, num_species)

dict_type = {0: 'domi', 1: 'high', 2: 'low'}

#NO SE SI ESTO VA ACA O EN OTRA CARPETA
vecTime = np.linspace(0, 500, 100)
H_total = 100
p_H_Infected = np.random.uniform(0, 0.9)
V_total = np.random.uniform(10, 50, size=num_species)
p_V_Infected = np.random.uniform(0, 0.9, size=num_species)

mat_a_test = np.random.uniform(1/(50), 1/(100), size=(num_species, num_species))
mat_a_test[0,0] = 0
mat_a_test[1,1] = 0
mat_a_test[2,2] = 0
mat_a_test[3,3] = 0

condInit = [H_total*(1-p_H_Infected), H_total*(p_H_Infected)]

for i in range(num_species):
    condInit.append(V_total[i]*(1-p_V_Infected[i]))

for i in range(num_species):
    condInit.append(V_total[i]*(p_V_Infected[i]))
    

#Define additional functions
def b_i(V, H, mat_pi):
    return (sigma_V[i]*V*sigma_H[i]*H)/(sigma_V[i]*V + sigma_H[i]*H + np.sum([mat_pi[i,j]*sigma_V[j]*V for j in range(num_species)]))

#Define the model function
def model_HostVectorNspecies(variables, t, mat_a, mat_pi):
    #The first 2 variables are the host population (humans, S-I)
    H_S = variables[0]
    H_I = variables[1]
    #The next 2*num_species variables are the vector population of different species (mosquitoes, S-I)
    V_S = [variables[i]  for i in range(2, 2+num_species)]
    V_I = [variables[i]  for i in range(2+num_species, 2+2*num_species)]
    
    derivatives = list()
    
    dH_Sdt = lambda_H - mu_H*H_S + gamma*H_I - np.sum([tau_HV[i]*(V_I[i])/(V_S[i]+V_I[i])*(b_i(V_S[i]+V_I[i], H_S+H_I, mat_pi)*H_S)/(H_S+H_I) for i in range(num_species)])
    dH_Idt = - mu_H*H_I - gamma*H_I + np.sum([tau_HV[i]*(V_I[i])/(V_S[i]+V_I[i])*(b_i(V_S[i]+V_I[i], H_S+H_I, mat_pi)*H_S)/(H_S+H_I) for i in range(num_species)])
    
    derivatives = derivatives + [dH_Sdt, dH_Idt]

    for i in range(num_species):
        dV_Sidt = -mu_V[i]*V_S[i] - tau_VH[i]*(H_I)/(H_S+H_I)*(b_i(V_S[i]+V_I[i], H_S+H_I, mat_pi)*V_S[i])/(V_S[i] + V_I[i]) + lambda_V[i]*(1-np.sum([mat_a[i,j]*(V_S[j]+V_I[j]) for j in range(num_species)]))
        derivatives.append(dV_Sidt)
        
    for i in range(num_species):
        dV_Iidt = -mu_V[i]*V_I[i] + tau_VH[i]*(H_I)/(H_S+H_I)*(b_i(V_S[i]+V_I[i], H_S+H_I, mat_pi)*V_S[i])/(V_S[i] + V_I[i]) 
        derivatives.append(dV_Iidt)
    
    return derivatives

def simulate_HostVectorNspecies(condInit, vecTime, compe_type):
    df_simu = pd.DataFrame({'Time': vecTime})

    for scenario in range(num_scenarios):
        mat_a_temp = np.load(os.path.join(params_path, str(num_species), compe_type, 'mat_a_'+str(scenario)+'.npy'))
        mat_pi_temp = np.load(os.path.join(params_path, str(num_species), compe_type, 'mat_pi_'+str(scenario)+'.npy'))
        simu_temp = odeint(model_HostVectorNspecies, condInit, vecTime, args=(mat_a_temp, mat_pi_temp))

        df_simu['H_S_'+str(scenario)] = simu_temp[:,0]
        df_simu['H_I_'+str(scenario)] = simu_temp[:,1]

        for i in range(num_species):
            df_simu['V_S_'+str(i)+'_'+str(scenario)] = simu_temp[:,2+i]
            df_simu['V_I_'+str(i)+'_'+str(scenario)] = simu_temp[:,2+num_species+i]
    
    if not os.path.exists(os.path.join(results_path, str(num_species), compe_type)):
        os.makedirs(os.path.join(results_path, str(num_species), compe_type))
    
    df_simu.to_csv(os.path.join(results_path, str(num_species), compe_type, 'simu'+str(num_scenarios)+'.csv'), index=False)
        

for key1, value1 in dict_type.items():
    for key2, value2 in dict_type.items():
        simulate_HostVectorNspecies(condInit, vecTime, value1+'_'+value2)
    
'''simu = odeint(model_HostVectorNspecies, condInit, vecTime, args=(mat_a_test, ))


#Plot the results
fig, ax = plt.subplots(2, 3, figsize=(14,8), sharex=True)

ax[0,0].plot(vecTime, simu[:,0], label='Host Susceptible')
ax[0,0].plot(vecTime, simu[:,1], label='Host Infected')
ax[0,0].grid()
ax[0,0].set_xlabel('Time')
ax[0,0].set_ylabel('Population')
ax[0,0].legend(loc='best')

ax[0,1].plot(vecTime, simu[:,2], label='Vector 1 Susceptible')
ax[0,1].plot(vecTime, simu[:,6], label='Vector 1 Infected')
ax[0,1].grid()
ax[0,1].set_xlabel('Time')
ax[0,1].set_ylabel('Population')
ax[0,1].legend(loc='best')

ax[0,2].plot(vecTime, simu[:,3], label='Vector 2 Susceptible')
ax[0,2].plot(vecTime, simu[:,7], label='Vector 2 Infected')
ax[0,2].grid()
ax[0,2].set_xlabel('Time')
ax[0,2].set_ylabel('Population')
ax[0,2].legend(loc='best')

ax[1,0].plot(vecTime, simu[:,4], label='Vector 3 Susceptible')
ax[1,0].plot(vecTime, simu[:,8], label='Vector 3 Infected')
ax[1,0].grid()
ax[1,0].set_xlabel('Time')
ax[1,0].set_ylabel('Population')
ax[1,0].legend(loc='best')

ax[1,1].plot(vecTime, simu[:,5], label='Vector 4 Susceptible')
ax[1,1].plot(vecTime, simu[:,9], label='Vector 4 Infected')
ax[1,1].grid()
ax[1,1].set_xlabel('Time')
ax[1,1].set_ylabel('Time')
ax[1,1].legend(loc='best')

plt.savefig('modelo_takimoto.png')

plt.show()'''