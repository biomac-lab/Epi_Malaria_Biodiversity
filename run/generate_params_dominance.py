#Import libraries
import numpy as np 
import os, sys
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.stats import truncnorm

#Set the seed
np.random.seed(42)

#Define the parameters of the model
num_species = 4
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

#NO SE SI ESTO VA ACA O EN OTRA CARPETA
vecTime = np.linspace(0, 500, 100)
H_total = 100
p_H_Infected = np.random.uniform(0, 0.9)
V_total = np.random.uniform(10, 50, size=num_species)
p_V_Infected = np.random.uniform(0, 0.9, size=num_species)


def generate_CompetitionMatrix(type:int, L: int):
    #Dominance scenario
    if type == 0:
        mat1 = truncnorm.rvs(a=0, b=1, loc=0.75, scale=1, size=(L,L))
        mat2 = truncnorm.rvs(a=0, b=1, loc=0.25, scale=1, size=(L,L))
        for i in range(L):
            for j in range(L):
                if i == j:
                    mat1[i,j] = 0
                    mat2[i,j] = 0
                elif i < j:
                    mat1[i,j] = 0
                elif i > j:
                    mat2[i,j] = 0
        mat_final = mat1 + mat2
    
    #High competition scenario
    elif type == 1:
        mat_final = truncnorm.rvs(a=0, b=1, loc=0.75, scale=1, size=(L,L))
        for i in range(L):
            mat_final[i,i] = 0

    #Low competition scenario
    elif type == 2:
        mat_final = truncnorm.rvs(a=0, b=1, loc=0.25, scale=1, size=(L,L))
        for i in range(L):
            mat_final[i,i] = 0
    
    return mat_final

scenarios = 20
cont = 0
dict_type = {0: 'domi', 1: 'high', 2: 'low'}

fig, axes = plt.subplots(3,3)
for type1 in range(3):
    for type2 in range(3):
        array_a = np.zeros((scenarios*num_species**2))
        array_pi = np.zeros((scenarios*num_species**2))
        while cont < scenarios:
            mat_a = generate_CompetitionMatrix(type1, num_species)
            mat_pi = generate_CompetitionMatrix(type2, num_species)

            array_a[num_species**2*cont:(num_species**2*(cont+1))] = mat_a.flatten()
            array_pi[num_species**2*cont:(num_species**2*(cont+1))] = mat_pi.flatten()

            if not os.path.exists(os.path.join('params',str(num_species),dict_type[type1]+'_'+dict_type[type2])):
                os.makedirs(os.path.join('params',str(num_species),dict_type[type1]+'_'+dict_type[type2]))
            
            np.save(os.path.join('params',str(num_species),dict_type[type1]+'_'+dict_type[type2],'mat_a_'+str(cont)), mat_a)
            np.save(os.path.join('params',str(num_species),dict_type[type1]+'_'+dict_type[type2],'mat_pi_'+str(cont)), mat_pi)
            
            cont += 1
        cont = 0

        axes[type1, type2].hist(array_a, alpha=0.5, label=r'$\alpha$')
        axes[type1, type2].hist(array_pi, alpha=0.5, label=r'$\pi$')
        
        if type1 == 2 and type2 == 2:
            axes[type1, type2].legend()

        if type2 == 0:
            axes[type1, type2].set_ylabel(dict_type[type1])

        if type1 == 0:
            axes[type1, type2].set_title(dict_type[type2])

fig.suptitle(r'Escenarios of $\pi$')
fig.supylabel(r'Escenarios of $\alpha$')
fig.show()
fig.savefig(os.path.join('params',str(num_species),'hist_experimentationScenarios.jpeg'), dpi=600)