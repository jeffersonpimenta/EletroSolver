import numpy as np
import time
import json

class Barra:
    def __init__(self, indice, tipo, V, theta, P, Q):
        self.indice = indice  # Índice da barra
        self.tipo = tipo      # Tipo de barra: 1 = PQ, 2 = PV, 3 = Slack
        self.V = V            # Tensão na barra
        self.theta = theta    # Ângulo de fase
        self.P = P            # Potência ativa injetada
        self.Q = Q            # Potência reativa injetada

class SistemaPotencia:
    def __init__(self, barras, Y, tolerancia=1e-6, max_iter=100, Sbase=100):
        self.barras = barras  # Lista de objetos Barra
        self.Y = Y            # Matriz de admitâncias
        self.n_barras = len(barras)
        self.tolerancia = tolerancia
        self.max_iter = max_iter
        self.sbase = Sbase
        self.convergencia = None

    def v(self, indice):
        """Retorna a tensão em pu da barra `indice`."""
        return self.barras[indice - 1].V

    def theta(self, indice):
        """Retorna o ângulo em radianos da barra `indice`."""
        return self.barras[indice - 1].theta

    def inicializar_estado(self):
        """Inicializa os vetores de tensões e ângulos."""
        V = np.array([barra.V for barra in self.barras])
        theta = np.array([barra.theta for barra in self.barras])
        return V, theta
        
    def alterar_barra(self, indice, tipo=None, V=None, theta=None, P=None, Q=None):
        """Altera as propriedades de uma barra."""
        barra = self.barras[indice - 1]
        if tipo is not None:
            barra.tipo = tipo
        if V is not None:
            barra.V = V
        if theta is not None:
            barra.theta = theta
        if P is not None:
            barra.P = P
        if Q is not None:
            barra.Q = Q
        
    def calcular_fluxo(self):
        """Resolve o fluxo de potência usando o método de Newton-Raphson."""
        start_time = time.time()  # Marca o tempo inicial
        V, theta = self.inicializar_estado()

        for k in range(self.max_iter):
            P_calc, Q_calc = self.fluxo_potencia(V, theta)
            dP, dQ = self.calcular_desvios(P_calc, Q_calc)
            dX = np.concatenate([dP, dQ])

            if np.linalg.norm(dX) < self.tolerancia:
                self.convergencia = k + 1
                #print(f"Convergência atingida em {k + 1} iterações.")
                self.ultimajacobiana = J
                break

            J = self.calcular_jacobiano(V, theta, P_calc, Q_calc)
            delta_X = np.linalg.solve(J, dX)

            # Atualiza theta e V
            theta[1:] += delta_X[:self.n_barras - 1]
            indices_pq = [i for i, barra in enumerate(self.barras) if barra.tipo == 1]
            delta_values = delta_X[self.n_barras - 1:]
            for idx, i in enumerate(indices_pq):
                V[i] += delta_values[idx]

        # Atualiza os valores de V e theta nas barras
        for i, barra in enumerate(self.barras):
            barra.V = V[i]
            barra.theta = theta[i]
        
        self.tempo = time.time() - start_time #calcula o tempo de solução do sistema.

    def calcular_desvios(self, P_calc, Q_calc):
        dP = np.array([barra.P - P_calc[i] for i, barra in enumerate(self.barras) if barra.tipo != 3])
        dQ = np.array([barra.Q - Q_calc[i] for i, barra in enumerate(self.barras) if barra.tipo == 1])
        return dP, dQ

    def calcular_jacobiano(self, V, theta, P_calc, Q_calc):
        n_pq = len([b for b in self.barras if b.tipo == 1])
        H = np.zeros((self.n_barras - 1, self.n_barras - 1))
        N = np.zeros((self.n_barras - 1, n_pq))
        M = np.zeros((n_pq, self.n_barras - 1))
        L = np.zeros((n_pq, n_pq))

        # Cálculo das partes da matriz Jacobiana
        pq_indices = [i for i, b in enumerate(self.barras) if b.tipo == 1]
        pv_pq_indices = [i for i, b in enumerate(self.barras) if b.tipo == 1 or b.tipo == 2]
        
        # Cálculo de H (Derivada de P em relação a theta)
        for i in range(1, self.n_barras):  # Ignorando a barra slack (barra 1)
            for j in range(1, self.n_barras):
                if i == j:
                    H[i-1, j-1] = -Q_calc[i] - (V[i] ** 2) * self.Y[i, i].imag#diagonal principal
                else:
                    H[i-1, j-1] = V[i] * V[j] * (self.Y[i, j].real * np.sin(self.theta(i) - self.theta(j)) - self.Y[i, j].imag * np.cos(self.theta(i) - self.theta(j)))#demais elementos
        
        # Cálculo de N (Derivada de P em relação a V)
        for i in range(1, self.n_barras):
            for idx, j in enumerate(pq_indices):
                if i == j:
                    soma = 0
                    for k in range(self.n_barras):
                        if k != i:
                            soma += V[k] * (self.Y[i, k].real * np.cos(self.theta(i) - self.theta(k)) + self.Y[i, k].imag * np.sin(self.theta(i) - self.theta(k)))
                    N[i-1, idx] = 2 * V[i] * self.Y[i, i].real + soma
                else:
                    N[i-1, idx] = V[i] * (self.Y[i, j].real * np.cos(self.theta(i) - self.theta(j)) + self.Y[i, j].imag * np.sin(self.theta(i) - self.theta(j)))

        # Cálculo de M (Derivada de Q em relação a theta)
        for idx_i, i in enumerate(pq_indices):
            for j in range(1, self.n_barras):
                if i == j:
                    M[idx_i, j-1] = P_calc[i] - (V[i] ** 2) * self.Y[i, i].real
                else:
                    M[idx_i, j-1] = -V[i] * V[j] * (self.Y[i, j].real * np.cos(self.theta(i) - self.theta(j)) + self.Y[i, j].imag * np.sin(self.theta(i) - self.theta(j)))

        # Cálculo de L (Derivada de Q em relação a V)
        for idx_i, i in enumerate(pq_indices):
            for idx_j, j in enumerate(pq_indices):
                if i == j:
                    soma = 0
                    for k in range(self.n_barras):
                        if k != i:
                            soma += V[k] * (self.Y[i, k].real * np.sin(self.theta(i) - self.theta(k)) - self.Y[i, k].imag * np.cos(self.theta(i) - self.theta(k)))
                    L[idx_i, idx_j] = -2 * V[i] * self.Y[i, i].imag + soma
                else:
                    L[idx_i, idx_j] = V[i] * (self.Y[i, j].real * np.sin(self.theta(i) - self.theta(j)) - self.Y[i, j].imag * np.cos(self.theta(i) - self.theta(j)))
                    
        return np.block([[H, N], [M, L]])#monta a matriz jacobiana com as submatrizes das derivadas de P e Q em relação a v e theta

    def fluxo_potencia(self, V, theta):
        """Calcula a potência ativa e reativa nas barras."""
        P_calc = np.zeros(self.n_barras)
        Q_calc = np.zeros(self.n_barras)

        for i in range(self.n_barras):
            for j in range(self.n_barras):
                P_calc[i] += V[i] * V[j] * (
                    self.Y[i, j].real * np.cos(theta[i] - theta[j]) +
                    self.Y[i, j].imag * np.sin(theta[i] - theta[j])
                )
                Q_calc[i] += V[i] * V[j] * (
                    self.Y[i, j].real * np.sin(theta[i] - theta[j]) -
                    self.Y[i, j].imag * np.cos(theta[i] - theta[j])
                )
        return P_calc, Q_calc

    def transito(self, de, para):
        """Calcula o trânsito de potência ativa e reativa entre duas barras."""
        i = de - 1  # Índice da barra 'de'
        j = para - 1  # Índice da barra 'para'

        # Tensão complexa nas barras
        Vi = self.barras[i].V * (np.cos(self.barras[i].theta) + 1j * np.sin(self.barras[i].theta))
        Vj = self.barras[j].V * (np.cos(self.barras[j].theta) + 1j * np.sin(self.barras[j].theta))

        # Corrente complexa de i para j
        I_ij = self.Y[i, j] * (Vi - Vj)

        # Potência complexa de i para j
        S_ij = Vi * np.conj(I_ij)
        P_ij = S_ij.real  # Potência ativa
        Q_ij = S_ij.imag  # Potência reativa

        return {"S_ij": S_ij, "P_ij": P_ij, "Q_ij": Q_ij}


    def losses(self, de, para):
        """Calcula as perdas de potência ativa e reativa na linha entre duas barras."""
        transito_ij = self.transito(de, para)
        transito_ji = self.transito(para, de)

        # Perdas de potência ativa e reativa
        S_loss = abs(transito_ij["S_ij"] + transito_ji["S_ij"])  # Potência complexa
        P_loss = abs(transito_ij["P_ij"] + transito_ji["P_ij"])  # Potência ativa
        Q_loss = abs(transito_ij["Q_ij"] + transito_ji["Q_ij"])  # Potência reativa

        return {"S_loss": S_loss, "P_loss": P_loss, "Q_loss": Q_loss}

    def totlosses(self):
        """Calcula as perdas totais de potência ativa e reativa em todas as linhas do sistema."""
        total_P_loss = 0.0  # Inicializa a soma das perdas de potência ativa
        total_Q_loss = 0.0  # Inicializa a soma das perdas de potência reativa

        # Itera sobre todas as combinações de barras conectadas
        for i in range(self.n_barras):
            for j in range(i + 1, self.n_barras):  # Considera apenas pares únicos (i, j)
                if self.Y[i, j] != 0:  # Verifica se há uma conexão entre as barras
                    perdas = self.losses(i + 1, j + 1)  # Calcula as perdas entre as barras
                    total_P_loss += perdas["P_loss"]  # Soma as perdas de potência ativa
                    total_Q_loss += perdas["Q_loss"]  # Soma as perdas de potência reativa

        # retorna as perdas totais
        return {"P_loss":total_P_loss,"Q_loss":total_Q_loss}
    
    def calcular_sensibilidade(self):
        """Calcula as sensibilidades de potência ativa e reativa com base na matriz Jacobiana."""
        if hasattr(self, 'ultimajacobiana'):
            return np.linalg.inv(self.ultimajacobiana)  # Retorna a matriz de sensibilidade
        else:
            raise ValueError("A Jacobiana ainda não foi calculada.")    

    def exportar(self, arquivo="sistema.json"):
        """Exporta todas as informações do sistema para um arquivo JSON."""
        def serializar_valor(valor):
            """Converte valores complexos em um formato serializável."""
            if isinstance(valor, complex):
                return {"real": valor.real, "imag": valor.imag}
            if isinstance(valor, np.ndarray):  # Converte arrays para listas
                return valor.tolist()
            return valor

        dados = {
            "barras": [
                {
                    "indice": barra.indice + 1,
                    "tipo": barra.tipo,
                    "tensao": barra.V,
                    "angulo": barra.theta,
                    "potencia_ativa": barra.P,
                    "potencia_reativa": barra.Q
                }
                for barra in self.barras
            ],
            "fluxos": [
                {
                    "linha": f"{i + 1} -> {j + 1}",
                    **{chave: serializar_valor(valor) for chave, valor in self.transito(i + 1, j + 1).items()}
                }
                for i in range(self.n_barras)
                for j in range(i + 1, self.n_barras)
                if self.Y[i, j] != 0
            ] + [
                {
                    "linha": f"{j + 1} -> {i + 1}",
                    **{chave: serializar_valor(valor) for chave, valor in self.transito(j + 1, i + 1).items()}
                }
                for i in range(self.n_barras)
                for j in range(i + 1, self.n_barras)
                if self.Y[j, i] != 0
            ],
            "perdas": self.totlosses(),
            "parametros": {
                "tolerancia": self.tolerancia,
                "max_iteracoes": self.max_iter,
                "sbase": self.sbase,
                "tempo_solucao": getattr(self, 'tempo', None)
            },
            "matriz_admitancia": [
                [serializar_valor(self.Y[i, j]) for j in range(self.n_barras)] 
                for i in range(self.n_barras)
            ]
        }

        with open(arquivo, 'w') as f:
            json.dump(dados, f, indent=4)
        print(f"Informações exportadas para {arquivo}.")

    
    def imprimir_estado(self,cd=4):
        print("Estado do Sistema:")
        print("-" * 40)   
        for barra in self.barras:
            print(f"Barra {barra.indice + 1}: V = {barra.V:.{cd}f} pu, theta = {barra.theta:.{cd}f} rad")
        print("-" * 40)

    def imprimir_transito(self,cd=4):
        """Imprime o trânsito de potência ativa e reativa em todas as linhas do sistema."""
        print("Trânsito de Potências:")
        print("-" * 40)
        for i in range(self.n_barras):
            for j in range(i + 1, self.n_barras):  # Considera apenas pares únicos (i, j)
                if self.Y[i, j] != 0:  # Verifica se há uma conexão entre as barras
                    transito = self.transito(i + 1, j + 1)  # Calcula trânsito entre barras
                    print(f"Linha {i + 1} -> {j + 1}: "
                          f"P = {transito['P_ij']:.{cd}f} pu, Q = {transito['Q_ij']:.{cd}f} pu")
        print("-" * 40)

    def imprimir_perdas(self,cd=4):
        """Imprime as perdas de potência ativa e reativa em todas as linhas do sistema."""
        print("Perdas de Potência nas Linhas:")
        print("-" * 40)
        for i in range(self.n_barras):
            for j in range(i + 1, self.n_barras):  # Considera apenas pares únicos (i, j)
                if self.Y[i, j] != 0:  # Verifica se há uma conexão entre as barras
                    perdas = self.losses(i + 1, j + 1)  # Calcula as perdas entre as barras
                    print(f"Linha {i + 1} -> {j + 1}: "
                          f"P_loss = {perdas['P_loss']:.{cd}f} pu, "
                          f"Q_loss = {perdas['Q_loss']:.{cd}f} pu")
        print("-" * 40)