#Import libraries
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sympy as sym
import csv

#Set the seed
np.random.seed(0)

name = 'Results_V1'

variables_population = (1000,100,0,                               # Variables iniciales de HOST (S,I,R)
                        10000,10000,10000,10000,10000,10000,      # Variables iniciales de Vectores Suceptible (V1,V2,...,V5,V6)
                        1000,1000,1000,1000,1000,1000)            # Variables inicialed de Vectores Infectados (V1,V2,...,V5,V6)

num_species = 6                                                     # Numero de especies en el modelo
r = [3/1000,3/1000,3/1000,3/1000,3/1000,3/1000]       # Variables de crecimiento de Vectores (V1,V2,...,V5,V6)
beta = [0.1/5000,0.1/5000,0.1/5000,0.1/5000,0.1/5000,0.1/5000]          # Variables de biting rate de Vectores (V1,V2,...,V5,V6)
rho = [0.02,0.02,0.02,0.02,0.02,0.02]                                   # Variables de probabilidad de infectar de los Vectores (V1,V2,...,V5,V6)
mu_S = [0.4/1000,0.4/1000,0.4/1000,0.4/1000,0.4/1000,0.4/1000]                      # Variables de muerte de Vectores Suceptibles (V1,V2,...,V5,V6)
mu_I = [1.1/1000,1.1/1000,1.1/1000,1.1/1000,1.1/1000,1.1/1000]                      # Variables de muerte de Vectores Infectados (V1,V2,...,V5,V6)

np.random.seed(0)
#gamma =np.zeros((6,6))                                          # Matriz de interaccion NULA
gamma = np.zeros((num_species, num_species))

for i in range(num_species):
    for j in range(i+1, num_species):
        random = np.random.randint(1,9)
        gamma[i][j] =  random
gamma = gamma + gamma.T
gamma = gamma/50000000

r_H = 1/1000      # Tasa de crecimiento de Host
g_S = 1/1000  # Tasa de muerte de Host Suceptible
g_I = 1/1000      # Tasa de muerte de Host Infectado
g_R = 1/1000    # Tasa de muerte de Host Recuperado
tau = 1/10000    # Tasa de suceptibilisacion ¿?
lambda_ = 1/10000 # Tasa de recuperacion de Host Infectado
time = 4000
vecTime = np.linspace(0, time,time)      # Vector de Tiempo


#Define the model function
def model(variables_population, t, num_species):

    #The first 3 variables are the host population (humans, S-I-R)
    HS = variables_population[0]
    HI = variables_population[1]
    HR = variables_population[2]

    #The next 6 variables are the vector population of different species (mosquitoes, S-I)
    VS = [variables_population[i] for i in range(3, 3+num_species)]
    VI = [variables_population[i] for i in range(3+num_species, 3+num_species+num_species)]

    #The lists to store the derivatives of the vector population
    derivadasVS = list()
    derivadasVI = list()

    infection_H = 0

    #Iteration over the total number of species
    for i in range(num_species):
        #The interaction term for each species
        interaction_V = 0
        for j in range(num_species):
            #The if removes the competition of the species with itself
            #If you want to include the intraespecific competition, remove the if
            if i != j:
                interaction_V += gamma[i][j] * (VS[i] + VI[i])  * (VS[j] + VI[j])
            
        #The derivatives of the vector population
        dtVS = r[i] * (VS[i] + VI[i]) - beta[i] * rho[i] * VS[i] * HI - mu_S[i]* VS[i] - interaction_V
        dtVI = beta[i] * rho[i] * VS[i] * HI - mu_I[i] * VI[i]

        derivadasVS.append(dtVS)
        derivadasVI.append(dtVI)

        #Calculate the infection term for the host population     
        infection_H += beta[i] * rho[i] * VI[i]

    #The derivatives of the host population
    dtHS = r_H * (HS + HI + HR) + tau * HR - infection_H * HS - g_S * HS
    dtHI = infection_H * HS - g_I * HI - lambda_ * HI
    dtHR = lambda_ * HI - tau * HR - g_R * HR

    #Store the derivatives in a list
    derivadas = list()
    derivadas.append(dtHS)
    derivadas.append(dtHI)
    derivadas.append(dtHR)

    
    for i in range(len(VS)):
        derivadas.append(derivadasVS[i])

    for i in range(len(VI)):
        derivadas.append(derivadasVI[i])

    return derivadas


solve = odeint(model,variables_population, vecTime, args=(num_species,))

fig, ax = plt.subplots(1,3, figsize = (18,3))

ax[0].plot(vecTime, solve[:,0] + solve[:,1] +solve[:,2], 'k--', label = 'Total population')
ax[0].plot(vecTime, solve[:,0], color = 'b', label = 'Suceptible population')
ax[0].plot(vecTime, solve[:,1], color = 'r',label = 'Infected population')
ax[0].plot(vecTime, solve[:,2], color = 'c', label = 'Recovered population')

ax[1].plot(vecTime, solve[:,3], 'b', label = 'Vector 1')
ax[1].plot(vecTime, solve[:,4], 'y', label = 'Vector 2')
ax[1].plot(vecTime, solve[:,5], 'r', label = 'Vector 3')
ax[1].plot(vecTime, solve[:,6], 'g', label = 'Vector 4')
ax[1].plot(vecTime, solve[:,7], 'm', label = 'Vector 5')
ax[1].plot(vecTime, solve[:,8], 'c', label = 'Vector 6')

ax[2].plot(vecTime, solve[:,9], 'b', label = 'Vector 1')
ax[2].plot(vecTime, solve[:,10], 'y', label = 'Vector 2')
ax[2].plot(vecTime, solve[:,11], 'r', label = 'Vector 3')
ax[2].plot(vecTime, solve[:,12], 'g', label = 'Vector 4')
ax[2].plot(vecTime, solve[:,13], 'm', label = 'Vector 5')
ax[2].plot(vecTime, solve[:,14], 'c', label = 'Vector 6')

ax[1].plot(vecTime, solve[:,9] + solve[:,3], 'b--')
ax[1].plot(vecTime, solve[:,10] + solve[:,4], 'y--')
ax[1].plot(vecTime, solve[:,11] + solve[:,5], 'r--')
ax[1].plot(vecTime, solve[:,12] + solve[:,6], 'g--')
ax[1].plot(vecTime, solve[:,13] + solve[:,7], 'm--')
ax[1].plot(vecTime, solve[:,14] + solve[:,8], 'c--')

ax[0].grid(1)
ax[1].grid(1)
ax[2].grid(1)
#ax[3].grid(1)

ax[0].legend(loc = 'best')
ax[1].legend(loc = 'best')
ax[2].legend(loc = 'best')
#ax[3].legend(loc = 'best')

ax[0].set_title('HOST')
ax[1].set_title('Suceptible Vectors')
ax[2].set_title('Infected Vectors')


plt.savefig('{}.jpg'.format(name))
plt.show()
print(gamma)

solve_data = pd.DataFrame(solve, columns = ['SUCEPTIBLE HOST','INFECTED HOST','RECOVERED HOST',
                                            'SUCEPTIBLE VECTOR 1','SUCEPTIBLE VECTOR 2',
                                            'SUCEPTIBLE VECTOR 3','SUCEPTIBLE VECTOR 4',
                                            'SUCEPTIBLE VECTOR 5','SUCEPTIBLE VECTOR 6',
                                            'INFECTED VECTOR 1','INFECTED VECTOR 2',
                                            'INFECTED VECTOR 3','INFECTED VECTOR 4',
                                            'INFECTED VECTOR 5','INFECTED VECTOR 6'])

Variables_total_dict = {'num_species':num_species,
                        'H': variables_population[0:3],
                        'VS': variables_population[3:],
                        'VI': variables_population[3+num_species:3+num_species*2],
                        'r': r ,
                        'beta': beta,
                        'rho': rho,
                        'mu_S': mu_S,
                        'mu_I': mu_I,
                        'gamma': gamma,
                        'r_H': r_H,
                        'g_S': g_I,
                        'g_I': g_I,
                        'g_R': g_R,
                        'tau': tau,
                        'lambda': lambda_,
                        'time': time
                        }

solve_data.to_csv('{}.csv'.format(name))

with open('Variables_{}.csv'.format(name), 'w', newline='') as archivo_csv:
    escritor_csv = csv.DictWriter(archivo_csv, fieldnames=Variables_total_dict.keys())
    # Escribir el encabezado
    escritor_csv.writeheader()
    # Escribir los datos del diccionario
    escritor_csv.writerow(Variables_total_dict)

#Puntos de equilibrio
'''HS,HI,HR = sym.symbols('H^S H^I H^R')

VS_1,VS_2,VS_3,VS_4,VS_5,VS_6 = sym.symbols('V^S_1 V^S_2 V^S_3 V^S_4 V^S_5 V^S_6')
VI_1,VI_2,VI_3,VI_4,VI_5,VI_6 = sym.symbols('V^I_1 V^I_2 V^I_3 V^I_4 V^I_5 V^I_6')

interaction_V = sym.symbols('MI_V')
infection_H = sym.symbols('I_H')

dtVS1 = r[0] * (VS_1 + VI_1) - beta[0] * rho[0] * VS_1 * HI - mu_S[0]* VS_1 + alpha * (VS_1 + VI_1)  * HI - interaction_V
dtVI1 = beta[0] * rho[0] * VS_1 * HI - mu_I[0] * VI_1

dtVS2 = r[1] * (VS_2 + VI_2) - beta[1] * rho[1] * VS_2 * HI - mu_S[1]* VS_2 + alpha * (VS_2 + VI_2)  * HI - interaction_V
dtVI2 = beta[1] * rho[1] * VS_2 * HI - mu_I[1] * VI_2

dtVS3 = r[2] * (VS_3 + VI_3) - beta[2] * rho[2] * VS_3 * HI - mu_S[2]* VS_3 + alpha * (VS_3 + VI_3)  * HI - interaction_V
dtVI3 = beta[2] * rho[2] * VS_3 * HI - mu_I[2] * VI_3

dtVS4 = r[3] * (VS_4 + VI_4) - beta[3] * rho[3] * VS_4 * HI - mu_S[3]* VS_4 + alpha * (VS_4 + VI_4)  * HI - interaction_V
dtVI4 = beta[3] * rho[3] * VS_4 * HI - mu_I[3] * VI_4

dtVS5 = r[4] * (VS_5 + VI_5) - beta[4] * rho[4] * VS_5 * HI - mu_S[4]* VS_5 + alpha * (VS_5 + VI_5)  * HI - interaction_V
dtVI5 = beta[4] * rho[4] * VS_5 * HI - mu_I[4] * VI_5

dtVS6 = r[5] * (VS_6 + VI_6) - beta[5] * rho[5] * VS_6 * HI - mu_S[5]* VS_6 + alpha * (VS_6 + VI_6)  * HI - interaction_V
dtVI6 = beta[5] * rho[5] * VS_6 * HI - mu_I[5] * VI_6

dtHS = r_H * (HS + HI) + tau * HI - infection_H * HS - g_S * HS
dtHI = infection_H * HS - g_I * HI
dtHR = landa * HI - tau * HR - g_R * HR

SOLVER = sym.solve([dtHS,dtHI,dtHR,dtVS1,dtVS2,dtVS3,dtVS4,dtVS5,dtVS6,dtVI1,dtVI2,dtVI3,dtVI4,dtVI5,dtVI6],[HS,HI,HR,VS_1,VS_2,VS_3,VS_4,VS_5,VS_6,VI_1,VI_2,VI_3,VI_4,VI_5,VI_6])

SOLVER

sym.print_latex(SOLVER)



"""$\left[ \left( 0.0, \  0.0, \  0.0, \  - 625.0 T, \  - 1428.57142857143 T, \  - 555.555555555556 T, \  - 588.235294117647 T, \  - 1666.66666666667 T, \  - 666.666666666667 T, \  0.0, \  0.0, \  0.0, \  0.0, \  0.0, \  0.0\right)\right]$"""'''