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
- Cálculo de trânsito de potência entre barras (convenção: `S_ij` é a potência
  que sai da barra `de` em direção a `para`, medida na extremidade `de`)
- Trânsito exato com taps de transformador e carregamento shunt de linha,
  passando os dados de ramo via `linhas=` (dataclass `Linha`); sem eles, o
  cálculo via Ybus é exato para ramos série puros
- Estimativa de perdas de potência com sinal: `Q_loss` negativo indica linha
  gerando reativo (carregamento capacitivo)
- Após `calcular_fluxo()`, P e Q da slack e Q das barras PV são atualizados
  com os valores resolvidos (refletidos também em `exportar()`)
- Determinação de tensão e ângulo de fase
- Relatórios sobre o estado do sistema

### Análise de Curto-Circuito (Corrente de Faltas)
- Faltas trifásica simétrica, monofásica (fase-terra), bifásica (fase-fase) e
  bifásica-terra, via componentes simétricas (sequências positiva, negativa e zero)
- Modelagem de geradores por reatâncias subtransitórias (X''d) e aterramento de neutro
- Roteamento de sequência zero por tipo de ligação de transformador (YNyn, Dyn, YNd, Dd…)
- Cargas modeladas como impedância constante, com sua influência na impedância de
  Thévenin e contribuição à corrente de falta
- Carregamento capacitivo das linhas (susceptância shunt `b1`/`b0`), incluído por
  padrão, com sua contribuição capacitiva à corrente de falta
- Tensão pré-falta selecionável (flat 1∠0 pu, padrão, ou solução do fluxo)
- Tensões pós-falta em todas as barras e contribuições por gerador, linha e carga
- Conversão para kA e potência de curto-circuito (MVA)

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

### Dataclass `Linha`
- Dados de ramo para trânsito/perdas exatos: `de`, `para`, `z` (impedância
  série), `b` (susceptância shunt total, modelo pi) e `tap` (lado `de`)
- Passada ao construtor via `SistemaPotencia(barras, Y, linhas=[...])`

### Módulo `Faltas` (curto-circuito)
- `Gerador`, `Ramo`, `Carga`: dados de sequência das fontes, linhas/trafos e cargas
- `EstudoCurtoCircuito`: monta as redes de sequência (Z1, Z2, Z0) e calcula faltas
  - `falta_trifasica(barra, Zf=0)`
  - `falta_monofasica(barra, Zf=0, Zg=0)`
  - `falta_bifasica(barra, Zf=0)`
  - `falta_bifasica_terra(barra, Zf=0, Zg=0)`
  - `corrente_kA`, `potencia_curto_MVA`, `imprimir_falta`

## Exemplos de Uso

### Análise Simples de Sistema de Potência
```python
from EletroSolver import Barra, Linha, SistemaPotencia
import numpy as np

# Define configurações de barras
barras = [
    Barra(0, 3, 1.02, 0.0, -0.5327, 5),
    Barra(1, 1, 1.0, 0.0, 0.6, 0.25),
    Barra(2, 1, 1.0, 0.0, 0.4, 0.15),
    Barra(3, 2, 1.01, 0.0, -0.21, 0.15),
    Barra(4, 1, 1.0, 0.0, 1.0, 0.4)
]

# Define os ramos da rede (b = carregamento shunt total; tap p/ trafos)
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

# Cria e analisa sistema de potência
# (linhas= é opcional: habilita trânsito/perdas exatos por ramo,
#  inclusive com tap e carregamento; sem ele, usa-se a Ybus)
sistema = SistemaPotencia(barras, Y, linhas=linhas)
sistema.calcular_fluxo()

# Acessa resultados
print(f"Convergência em {sistema.convergencia} iterações")
sistema.imprimir_estado()
sistema.imprimir_transito()
sistema.imprimir_perdas()
```

### Análise de Curto-Circuito
```python
from EletroSolver import Barra, SistemaPotencia
from Faltas import Carga, EstudoCurtoCircuito, Gerador, Ramo
import numpy as np

# Dados de sequência: geradores, ramos (linhas/trafos) e cargas
geradores = [
    Gerador(1, X1=0.15, X2=0.15, X0=0.05, aterrado=True),
    Gerador(3, X1=0.20, X2=0.20, X0=0.08, aterrado=True),
]
ramos = [
    # b1/b0: carregamento capacitivo shunt da linha (sequência +/- e zero)
    Ramo(1, 2, z1=0.02 + 0.08j, z0=0.06 + 0.24j, b1=0.10, b0=0.06),
    Ramo(2, 3, z1=0.02 + 0.08j, z0=0.06 + 0.24j, b1=0.10, b0=0.06),
]
cargas = [Carga(2, P=1.0, Q=0.4)]

# SistemaPotencia fornece nº de barras e Sbase (pré-falta flat dispensa o fluxo)
barras = [Barra(0, 3, 1.0, 0.0, 0.0, 0.0),
          Barra(1, 1, 1.0, 0.0, -1.0, -0.4),
          Barra(2, 2, 1.0, 0.0, 0.5, 0.0)]
sistema = SistemaPotencia(barras, np.zeros((3, 3), dtype=complex), Sbase=100)

est = EstudoCurtoCircuito(sistema, geradores, ramos, cargas=cargas)

# Falta trifásica na barra 2
res = est.falta_trifasica(2)
est.imprimir_falta(res)
print(f"Corrente: {abs(res['I_fase']['a']):.4f} pu = "
      f"{est.corrente_kA(res['I_fase']['a'], kV_base=138):.4f} kA")

# Contribuições para a corrente de falta
for g in res["contrib_gerador"]:
    print(f"Gerador {g['barra']}: {abs(g['I_fase']['a']):.4f} pu")
for c in res["contrib_carga"]:
    print(f"Carga barra {c['barra']}: {abs(c['I_fase']['a']):.4f} pu")
```

Veja `Curtos.py` para um exemplo completo com as quatro faltas.

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
pip install -r requirements.txt
```

## Dependências
- NumPy
- Python 3.7+

## Testes

```bash
pip install -r requirements-dev.txt
pytest
```

A suíte cobre o Jacobiano (analítico vs. diferenças finitas), invariantes físicos
(mismatch nulo e balanço de potência), validação de entrada, detecção de
(não-)convergência e os casos IEEE padrão (WSCC 9 barras e IEEE 14 barras,
dados MATPOWER) — com a solução do IEEE 14 conferida contra a referência publicada.
A análise de curto-circuito é validada por valores analíticos exatos (correntes 3φ,
FT e FF em sistemas pequenos), propriedades das componentes simétricas, balanço de
KCL das contribuições e tratamento de erros (redes de sequência singulares, dados
de sequência zero ausentes).
## Licença
GNU GENERAL PUBLIC LICENSE

## Contribuições
Contribuições são bem-vindas!

## Reconhecimentos
- Inspirado em princípios de engenharia de sistemas de potência
- Implementação do método de Newton-Raphson
