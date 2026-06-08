"""Valida calcular_jacobiano contra um Jacobiano numerico (diferencas finitas).

Este teste prova diretamente o bug B1: a versao original montava o Jacobiano com
self.theta(i) (o metodo, que le barras[i-1].theta) em vez do array de trabalho theta[i].
Em um ponto de operacao nao-flat (angulos != 0), o Jacobiano analitico errado nao bate
com o numerico. Corrigido, ambos coincidem.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from EletroSolver import Barra, SistemaPotencia


def _sistema_5_barras():
    barras = [
        Barra(0, 3, 1.02, 0.0, -0.5327, 5),
        Barra(1, 1, 1.0, 0.0, 0.6, 0.25),
        Barra(2, 1, 1.0, 0.0, 0.4, 0.15),
        Barra(3, 2, 1.01, 0.0, -0.21, 0.15),
        Barra(4, 1, 1.0, 0.0, 1.0, 0.4),
    ]
    Y = np.zeros((len(barras), len(barras)), dtype=complex)
    Y[0, 1] = Y[1, 0] = -1 / (0.05 + 0.2j)
    Y[0, 3] = Y[3, 0] = -1 / (0.05 + 0.25j)
    Y[1, 2] = Y[2, 1] = -1 / (0.05 + 0.15j)
    Y[1, 4] = Y[4, 1] = -1 / (0.1 + 0.3j)
    Y[2, 3] = Y[3, 2] = -1 / (0.05 + 0.2j)
    Y[3, 4] = Y[4, 3] = -1 / (0.05 + 0.1j)
    for i in range(len(barras)):
        Y[i, i] = -np.sum(Y[i, :])
    return SistemaPotencia(barras, Y)


def _jacobiano_numerico(sistema, V, theta, eps=1e-6):
    n = sistema.n_barras
    nao_slack = [i for i in range(n) if sistema.barras[i].tipo != 3]
    pq = [i for i in range(n) if sistema.barras[i].tipo == 1]

    def F(Vx, thx):
        P_calc, Q_calc = sistema.fluxo_potencia(Vx, thx)
        return np.concatenate([[P_calc[i] for i in nao_slack],
                               [Q_calc[i] for i in pq]])

    n_vars = len(nao_slack) + len(pq)
    J = np.zeros((n_vars, n_vars))
    col = 0
    # derivadas em relacao a theta (barras nao-slack)
    for i in nao_slack:
        th_p = theta.copy(); th_p[i] += eps
        th_m = theta.copy(); th_m[i] -= eps
        J[:, col] = (F(V, th_p) - F(V, th_m)) / (2 * eps)
        col += 1
    # derivadas em relacao a V (barras PQ)
    for i in pq:
        V_p = V.copy(); V_p[i] += eps
        V_m = V.copy(); V_m[i] -= eps
        J[:, col] = (F(V_p, theta) - F(V_m, theta)) / (2 * eps)
        col += 1
    return J


def test_jacobiano_bate_com_diferencas_finitas():
    sistema = _sistema_5_barras()
    # Ponto de operacao NAO-flat: e aqui que o bug B1 aparece.
    V = np.array([1.02, 1.0, 1.0, 1.01, 1.0])
    theta = np.array([0.0, 0.10, -0.05, 0.08, -0.12])

    P_calc, Q_calc = sistema.fluxo_potencia(V, theta)
    J_analitico = sistema.calcular_jacobiano(V, theta, P_calc, Q_calc)
    J_numerico = _jacobiano_numerico(sistema, V, theta)

    assert np.allclose(J_analitico, J_numerico, atol=1e-4, rtol=1e-4), (
        "Jacobiano analitico difere do numerico — calcular_jacobiano usa "
        "angulos errados (bug B1)."
    )
