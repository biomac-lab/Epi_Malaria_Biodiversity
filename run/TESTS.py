from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sympy as sym
import cv2

def model(variables_population, t):

    HS = variables_population[0]
    HI = variables_population[1]
    HR = variables_population[2]

    VS = [variables_population[i] for i in range(3, 3+num_species)]
    VI = [variables_population[i] for i in range(3+num_species, 3+num_species+num_species)]

    derivadasVS = list()
    derivadasVI = list()

    for i in range(num_species):
        interaction_V = 0
        for j in range(num_species):
            if i != j:
                result = gamma[i][j] * (VS[i] + VI[i])  * (VS[j] + VI[j])
                interaction_V += result

        dtVS = r[i] * (VS[i] + VI[i]) - beta[i] * rho[i] * VS[i] * HI - mu_S[i]* VS[i] + alpha * (VS[i] + VI[i])  * HI - interaction_V
        dtVI = beta[i] * rho[i] * VS[i] * HI - mu_I[i] * VI[i]

        derivadasVS.append(dtVS)
        derivadasVI.append(dtVI)

        infection_H = 0

    for i in range(num_species):
        infection_H = beta[i] * rho[i] * VI[i]

    dtHS = r_H * (HS + HI + HR) + tau * HR - infection_H * HS - g_S * HS
    dtHI = infection_H * HS - g_I * HI
    dtHR = landa * HI - tau * HR - g_R * HR

    derivadas = list()
    derivadas.append(dtHS)
    derivadas.append(dtHI)
    derivadas.append(dtHR)

    for i in range(len(VS)):
        derivadas.append(derivadasVS[i])

    for i in range(len(VI)):
        derivadas.append(derivadasVI[i])

    return tuple(derivadas)
#VARIABLES
variables_population = (1000,100,0,                               # Variables iniciales de HOST (S,I,R)
                        10000,10000,10000,10000,10000,10000,      # Variables iniciales de Vectores Suceptible (V1,V2,...,V5,V6)
                        1000,1000,1000,1000,1000,1000)            # Variables inicialed de Vectores Infectados (V1,V2,...,V5,V6)

r = [1.4/1000,1.3/1000,1.2/1000,1.3/1000,1.4/1000,1.5/1000]        # Variables de crecimiento de Vectores (V1,V2,...,V5,V6)
beta = [0.9/1000,0.3/1000,0.6/1000,0.4/1000,0.2/1000,0.1/1000]                # Variables de biting rate de Vectores (V1,V2,...,V5,V6)
rho = [0.001,0.02,0.02,0.02,0.02,0.02]                            # Variables de probabilidad de infectar de los Vectores (V1,V2,...,V5,V6)
mu_S = [3/1000,2/1000,3/1000,3/1000,2/1000,3/1000]               # Variables de muerte de Vectores Suceptibles (V1,V2,...,V5,V6)
mu_I = [1/1000,2/1000,3/1000,3/1000,2/1000,1/1000]               # Variables de muerte de Vectores Infectados (V1,V2,...,V5,V6)
alpha = 0


#gamma =np.zeros((6,6))                                  # Matriz de interaccion NULA


#Toca definir com se hara la matriz, opcion 1 iun cociente y de ahi la interaccion

gamma = np.random.random((6,6))/50000000           # Matriz de interaccion
'''
gamma = [[0,1,1,1,1,1],
         [1,0,1,1,1,1],
         [1,1,0,1,1,1],
         [1,1,1,0,1,1],
         [1,1,1,1,0,1],
         [1,1,1,1,1,0]]
'''

r_H = 0      # Variable de crecimiento de Host
g_S = 0       # Variable de muerte de Host Suceptible
g_I = 0      # Variable de muerte de Host Infectado
g_R = 0      # Variable de muerte de Host Recuperado
tau = 1/100       # Variable de suceptibilisacion ¿?
landa = 1/1000  # Variable de recuperacion de Host Infectado

num_species = 6      # Numero de especies en el modelo

vecTime = np.linspace(0,5000,5000)      # Vector de Tiempo
#solver
solve = odeint(model,variables_population, vecTime)
#Mapeo del sistema
fig, ax = plt.subplots(1,3, figsize = (18,3))
ax[0].plot()


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
'''
ax[3].plot(vecTime, solve[:,9]+solve[:,3], 'b', label = 'Vector 1')
ax[3].plot(vecTime, solve[:,10]+solve[:,4], 'y', label = 'Vector 2')
ax[3].plot(vecTime, solve[:,11]+solve[:,5], 'r', label = 'Vector 3')
ax[3].plot(vecTime, solve[:,12]+solve[:,6], 'g', label = 'Vector 4')
ax[3].plot(vecTime, solve[:,13]+solve[:,7], 'm', label = 'Vector 5')
ax[3].plot(vecTime, solve[:,14]+solve[:,8], 'c', label = 'Vector 6')

'''
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

plt.savefig('results1_beta_model.jpg')
plt.show()
#print(gamma)
'''
solve_data = pd.DataFrame(solve, columns = ['SUCEPTIBLE HOST','INFECTED HOST','RECOVERED HOST',
                                            'SUCEPTIBLE VECTOR 1','SUCEPTIBLE VECTOR 2',
                                            'SUCEPTIBLE VECTOR 3','SUCEPTIBLE VECTOR 4',
                                            'SUCEPTIBLE VECTOR 5','SUCEPTIBLE VECTOR 6',
                                            'INFECTED VECTOR 1','INFECTED VECTOR 2',
                                            'INFECTED VECTOR 3','INFECTED VECTOR 4',
                                            'INFECTED VECTOR 5','INFECTED VECTOR 6'])
#Poner los resultados en un csv
solve_data.to_csv('results1_beta_model.csv')

#Puntos de equilibrio del sistema con valores 
HS,HI,HR = sym.symbols('H^S H^I H^R')

VS_1,VS_2,VS_3,VS_4,VS_5,VS_6 = sym.symbols('V^S_1 V^S_2 V^S_3 V^S_4 V^S_5 V^S_6')
VI_1,VI_2,VI_3,VI_4,VI_5,VI_6 = sym.symbols('V^I_1 V^I_2 V^I_3 V^I_4 V^I_5 V^I_6')

interaction_V = sym.symbols('T')
infection_H = sym.symbols('Y')

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

print(SOLVER)

sym.print_latex(SOLVER)
'''
