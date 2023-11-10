#Codigo de solo el modelo
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

        dtVS = r[i] * (VS[i] + VI[i]) - beta[i] * rho[i] * VS[i] * HI - mu_S[i]* VS[i] - interaction_V
        dtVI = beta[i] * rho[i] * VS[i] * HI - mu_I[i] * VI[i]

        derivadasVS.append(dtVS)
        derivadasVI.append(dtVI)

        infection_H = 0

        for i in range(num_species):
            infection_H = beta[i] * rho[i] * VI[i]

        dtHS = r_H * (HS + HI + HR) + tau * HR - infection_H * HS - g_S * HS
        dtHI = infection_H * HS - g_I * HI - lambda_ * HI
        dtHR = lambda_ * HI - tau * HR - g_R * HR

    derivadas = list()
    derivadas.append(dtHS)
    derivadas.append(dtHI)
    derivadas.append(dtHR)

    for i in range(len(VS)):
        derivadas.append(derivadasVS[i])

    for i in range(len(VI)):
        derivadas.append(derivadasVI[i])

    return tuple(derivadas)
