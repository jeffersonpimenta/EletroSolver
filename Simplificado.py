from EletroSolver import Barra, Linha, SistemaPotencia
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

# Define os ramos da rede (impedância série; b e tap opcionais)
linhas = [
    Linha(1, 2, z=0.05 + 0.2j),
    Linha(1, 4, z=0.05 + 0.25j),
    Linha(2, 3, z=0.05 + 0.15j),
    Linha(2, 5, z=0.1 + 0.3j),
    Linha(3, 4, z=0.05 + 0.2j),
    Linha(4, 5, z=0.05 + 0.1j),
]

# Monta a matriz Y a partir dos ramos
Y = np.zeros((len(barras), len(barras)), dtype=complex)
for lin in linhas:
    i, j = lin.de - 1, lin.para - 1
    ys = 1 / lin.z
    Y[i, i] += ys
    Y[j, j] += ys
    Y[i, j] -= ys
    Y[j, i] -= ys

# Criar o sistema e calcular fluxo (linhas= habilita trânsito exato por ramo)
sistema = SistemaPotencia(barras, Y, linhas=linhas)
sistema.calcular_fluxo()

# Calculando o trânsito de potência
transito_12 = sistema.transito(1, 2)
print(f"Trânsito de potência de Barra 1 para Barra 2: {transito_12['S_ij']:.4f}")

import cmath
magnitude = abs(sistema.transito(3, 4)["S_ij"])  # Magnitude
angulo = cmath.phase(sistema.transito(3, 4)["S_ij"])  # ângulo
print(f"Trânsito de 3 para 4: {magnitude:.4f} /_ {angulo:.4f} pu")

# Calculando perdas de potência
perdas_12 = sistema.losses(1, 2)
print(f"Perdas ativas de potência entre Barra 1 e Barra 2: {perdas_12['P_loss']:.4f}")

print(f"Perdas ativas totais do sistema: {sistema.totlosses()['P_loss']:.4f} pu")
print(f"Perdas reativas totais do sistema: {sistema.totlosses()['Q_loss']:.4f} pu")

print(f"Convergência em {sistema.convergencia} iterações com tempo de {sistema.tempo*1000:.4f} ms")
print(f"Tensão na Barra 3: {sistema.v(3)} pu")
print(f"Ângulo na Barra 3: {sistema.theta(3)} rad")
sistema.imprimir_estado()
sistema.imprimir_transito()
sistema.imprimir_perdas()

sistema.exportar()