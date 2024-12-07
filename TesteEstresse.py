from EletroSolver import Barra, SistemaPotencia
import numpy as np

n_barras = 200

# Gera as barras com diferentes tipos
barras = [
    Barra(indice=i, 
          tipo=(3 if i == 0 else np.random.choice([1, 2])),  # Slack é a primeira, outras são PQ ou PV
          V=1.0 + np.random.uniform(-0.05, 0.05),  # Tensão inicial com pequenas variações
          theta=0.0,  # Ângulo inicial
          P=np.random.uniform(-1.0, 1.0) if i > 0 else 0.0,  # Potência ativa
          Q=np.random.uniform(-0.5, 0.5) if i > 0 else 0.0)  # Potência reativa
    for i in range(n_barras)
]

# Gera uma matriz de admitância Y com conexões aleatórias
Y = np.zeros((n_barras, n_barras), dtype=complex)

# Conexões aleatórias entre as barras
for i in range(n_barras):
    for j in range(i + 1, n_barras):
        if np.random.rand() < 0.4:  # 40% de chance de haver conexão
            z = complex(np.random.uniform(0.05, 0.2), np.random.uniform(0.1, 0.3))  # Impedância
            Y[i, j] = Y[j, i] = -1 / z  # Admitância inversa da impedância

# Ajusta as admitâncias diagonais
for i in range(n_barras):
    Y[i, i] = -np.sum(Y[i, :])
    
print("Sistema Desafio Gerado:")
for barra in barras:
    print(f"Barra {barra.indice + 1}: Tipo={barra.tipo}, V={barra.V:.4f}, P={barra.P:.4f}, Q={barra.Q:.4f}")
print("\nMatriz Y:")
print(Y)  
    
    
sistema = SistemaPotencia(barras, Y)
sistema.calcular_fluxo()

print(f"Convergência em {sistema.convergencia} iterações com tempo de {sistema.tempo*1000:.4f} ms")

sistema.imprimir_estado()
#sistema.imprimir_transito()
#sistema.imprimir_perdas()