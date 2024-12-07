# EletroSolver: Toolkit de Análise de Sistemas Elétricos

## Visão Geral

EletroSolver é uma biblioteca Python para análise de sistemas de energia, oferecendo ferramentas para cálculo de fluxo de potência, modelagem de redes e avaliação do desempenho de sistemas elétricos. O toolkit implementa o método de Newton-Raphson para resolução de problemas de fluxo de potência e fornece informações sobre o comportamento de sistemas elétricos.

## Funcionalidades

### Modelagem de Sistemas de Potência
- Suporte para tipos de barras PQ, PV e Slack
- Configuração flexível de barras e redes
- Geração de matriz de admitância complexa

### Cálculo de Fluxo de Potência
- Implementação do método iterativo de Newton-Raphson
- Cálculos de potência ativa e reativa
- Rastreamento de convergência e métricas de desempenho
- Cálculos abrangentes da matriz Jacobiana

### Capacidades de Análise de Rede
- Cálculo de trânsito de potência entre barras
- Estimativa de perdas de potência
- Determinação de tensão e ângulo de fase
- Relatórios sobre o estado do sistema

## Componentes Principais

### Classe `Barra`
- Representa barras individuais do sistema de energia
- Atributos incluem:
  - Índice da barra
  - Tipo de barra
  - Magnitude de tensão
  - Ângulo de fase
  - Injeção de potência ativa e reativa

### Classe `SistemaPotencia`
- Motor principal de análise
- Métodos para:
  - Cálculo de fluxo de potência
  - Computação de trânsito e perdas
  - Visualização do estado do sistema

## Exemplos de Uso

### Análise Simples de Sistema de Potência
```python
from EletroSolver import Barra, SistemaPotencia
import numpy as np

# Define configurações de barras
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

# Cria e analisa sistema de potência
sistema = SistemaPotencia(barras, Y)
sistema.calcular_fluxo()

# Acessa resultados
print(f"Convergência em {sistema.convergencia} iterações")
sistema.imprimir_estado()
```

### Teste de Estresse para avaliação de desempenho
A biblioteca suporta geração e análise de sistemas de energia com n barras, possibilitando avaliação abrangente de desempenho do algoritmo.

```python

n_barras = 1000

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

```

## Instalação

```bash
# Clona o repositório
git clone https://github.com/jeffersonpimenta/EletroSolver.git

# Instala dependências
pip install numpy
```

## Dependências
- NumPy
- Python 3.7+

## Licença
GNU GENERAL PUBLIC LICENSE

## Contribuições
Contribuições são bem-vindas!

## Reconhecimentos
- Inspirado em princípios de engenharia de sistemas de potência
- Implementação do método de Newton-Raphson
