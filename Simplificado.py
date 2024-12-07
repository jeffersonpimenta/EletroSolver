from EletroSolver import Barra, SistemaPotencia
import numpy as np
 
# Uso simplificado

# Define o sistema de barras na sequência: indice, tipo, V, theta, P e Q
barras = [
    Barra(0, 3, 1.02, 0.0, -0.5327, 5),
    Barra(1, 1, 1.0, 0.0, 0.6, 0.25),
    Barra(2, 1, 1.0, 0.0, 0.4, 0.15),
    Barra(3, 2, 1.01, 0.0, -0.21, 0.15),
    Barra(4, 1, 1.0, 0.0, 1.0, 0.4)
]

# Define matriz Y e admitâncias
Y = np.zeros((len(barras), len(barras)), dtype=complex)
Y[0, 1] = Y[1, 0] = -1 / (0.05 + 0.2j)
Y[0, 3] = Y[3, 0] = -1 / (0.05 + 0.25j)
Y[1, 2] = Y[2, 1] = -1 / (0.05 + 0.15j)
Y[1, 4] = Y[4, 1] = -1 / (0.1 + 0.3j)
Y[2, 3] = Y[3, 2] = -1 / (0.05 + 0.2j)
Y[3, 4] = Y[4, 3] = -1 / (0.05 + 0.1j)

for i in range(len(barras)):
    Y[i, i] = -np.sum(Y[i, :])

# Criar o sistema e calcular fluxo
sistema = SistemaPotencia(barras, Y)
sistema.calcular_fluxo()

# Calculando o trânsito de potência
transito_12 = sistema.transito(1, 2)
print(f"Trânsito de potência de Barra 1 para Barra 2: {transito_12["S_ij"]:.4f}")

import cmath
magnitude = abs(sistema.transito(3, 4)["S_ij"])  # Magnitude
angulo = cmath.phase(sistema.transito(3, 4)["S_ij"])  # ângulo
print(f"Trânsito de 3 para 4: {magnitude:.4f} /_ {angulo:.4f} pu")

# Calculando perdas de potência
perdas_12 = sistema.losses(1, 2)
print(f"Perdas de potência entre Barra 1 e Barra 2: {perdas_12["S_loss"]:.4f}")

print(f"Convergência em {sistema.convergencia} iterações com tempo de {sistema.tempo*1000:.4f} ms")
print(f"Tensão na Barra 3: {sistema.v(3)} pu")
print(f"Ângulo na Barra 3: {sistema.theta(3)} rad")
sistema.imprimir_estado()
sistema.imprimir_transito()
sistema.imprimir_perdas()
