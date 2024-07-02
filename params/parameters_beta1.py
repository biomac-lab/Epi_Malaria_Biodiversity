variables_population = (1000,100,                       # Variables iniciales de HOST
                        10000,10000,10000,10000,0,0,    # Variables iniciales de Vectores Suceptible
                        1000,1000,1000,1000,0,0)        # Variables inicialed de Vectores Infectados

r = [1/1000,1.1/1000,1.2/1000,1.3/1000,1/1000,1/1000]         # Variables de crecimiento de Vectores
b = [0.001,0.0008,0.0006,0.0004,0.0002,0.0001]         # Variables de biting rate de Vectores
p = [0.02,0.02,0.02,0.02,0.02,0.02]                     # Variables de probabilidad de infectar de los Vectores
u_s = [1/1000,1/1000,1/1000,1/1000,1/1000,1/1000]       # Variables de muerte de Vectores Suceptibles
u_i = [5/1000,4/1000,3/1000,2/1000,1/1000,1/1000]             # Variables de muerte de Vectores Infectados
a = 0


gamma =np.zeros((6,6))                                  # Matriz de interaccion NULA

'''
gamma = np.random.random((6,6))/100000000               # Matriz de interaccion 
'''

r_h = 1/10000       # Variable de crecimiento de Host
g_s = 1/10000     # Variable de muerte de Host Suceptible
g_i = 1.2/1000       # Variable de muerte de HOst Infectado

num_species = 6     # Numero de especies en el modelo

vecTime = np.linspace(0,1000,1000)      # Vector de Tiempo
