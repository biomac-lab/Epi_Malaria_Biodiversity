# Modelo de Ecuaciones - Gabriel Arango - 23 oct 2023
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def model(variables_population, t):

    H_s = variables_population[0] 
    H_i = variables_population[1]

    VS = [variables_population[i] for i in range(2, 2+num_species)]
    VI = [variables_population[i] for i in range(2+num_species, 2+num_species+num_species)]

    derivadasVS = list()
    derivadasVI = list()

    for i in range(num_species):
        interaction_V = 0
        for j in range(num_species):
            if i != j:
                result = gamma[i][j] * (VS[i] + VI[i])  * (VS[i] + VI[i])
                interaction_V += result

        dtVS = r[i] * (VS[i] + VI[i]) - b[i] * p[i] * VS[i] * H_i - u_s[i]* VS[i] + a * (VS[i] + VI[i])  * H_i - interaction_V
        dtVI = b[i] * p[i] * VS[i] * H_i - u_i[i] * VI[i]

        derivadasVS.append(dtVS)
        derivadasVI.append(dtVI)

    infection_H = 0

    for i in range(num_species):
        infection_H = b[i] * p[i] * VI[i]

    dtH_s = r_h * (H_s + H_i) - infection_H * H_s - g_s * H_s
    dtH_i = infection_H * H_s - g_i * H_i

    derivadas = list()
    derivadas.append(dtH_s)
    derivadas.append(dtH_i)
    
    for i in range(len(VS)):
        derivadas.append(derivadasVS[i])

    for i in range(len(VI)):
        derivadas.append(derivadasVI[i])

    return tuple(derivadas)


#NO SE SI ESTO VA ACA O VA EN OTRA CARPETA :)

variables_population = (1000,100,                       # Variables iniciales de HOST
                        10000,10000,10000,10000,0,0,    # Variables iniciales de Vectores Suceptible
                        1000,1000,1000,1000,0,0)        # Variables inicialed de Vectores Infectados

r = [0/1000,0/1000,0/1000,0/1000,0/1000,0/1000]         # Variables de crecimiento de Vectores
b = [0.0001,0.0001,0.0001,0.0001,0.0001,0.0001]         # Variables de biting rate de Vectores
p = [0.01,0.01,0.01,0.01,0.01,0.01]                     # Variables de probabilidad de infectar de los Vectores
u_s = [0/1000,0/1000,0/1000,0/1000,0/1000,0/1000]       # Variables de muerte de Vectores Suceptibles
u_i = [0/100,0/100,0/100,0/100,0/100,0/100]             # Variables de muerte de Vectores Infectados
a = 0


gamma =np.zeros((6,6))                                  # Matriz de interaccion NULA

'''
gamma = np.random.random((6,6))/100000000               # Matriz de interaccion 
'''

r_h = 0/1000        # Variable de crecimiento de Host
g_s = 0/1000        # Variable de muerte de Host Suceptible
g_i = 0/1000        # Variable de muerte de HOst Infectado

num_species = 6     # Numero de especies en el modelo

vecTime = np.linspace(0,1000,1000)      # Vector de Tiempo

solve = odeint(model,variables_population, vecTime)

fig, ax = plt.subplots(1,4, figsize = (19,3))

ax[0].plot(solve[:,0], 'b', label = 'Suceptible population')
ax[0].plot(solve[:,0] + solve[:,1], 'k', label = 'max pupolation')
ax[1].plot(solve[:,1], 'r',label = 'Infected population')

ax[2].plot(solve[:,2], 'b', label = 'Vector 1')
ax[2].plot(solve[:,3], 'y', label = 'Vector 2')
ax[2].plot(solve[:,4], 'r', label = 'Vector 3')
ax[2].plot(solve[:,5], 'g', label = 'Vector 4')
ax[2].plot(solve[:,6], 'purple', label = 'Vector 5')
ax[2].plot(solve[:,7], 'cyan', label = 'Vector 6')

ax[3].plot(solve[:,8], 'b', label = 'Vector 1')
ax[3].plot(solve[:,9], 'y', label = 'Vector 2')
ax[3].plot(solve[:,10], 'r', label = 'Vector 3')
ax[3].plot(solve[:,11], 'g', label = 'Vector 4')
ax[3].plot(solve[:,12], 'purple', label = 'Vector 5')
ax[3].plot(solve[:,13], 'cyan', label = 'Vector 6')

ax[0].grid(1)
ax[1].grid(1)
ax[2].grid(1)
ax[3].grid(1)

ax[0].legend(loc = 'best')
ax[1].legend(loc = 'best')
ax[2].legend(loc = 'best')
ax[3].legend(loc = 'best')

ax[0].set_title('HOST')
ax[1].set_title('HOST')
ax[2].set_title('Suceptible Vectors')
ax[3].set_title('Infected Vectors')
plt.show()