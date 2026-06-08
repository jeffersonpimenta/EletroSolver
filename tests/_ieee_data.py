"""Casos de teste IEEE padrao (dados MATPOWER) e construtor de sistema.

Fonte: MATPOWER (https://github.com/MATPOWER/matpower), case9 (WSCC 9 barras) e
case14 (IEEE 14 barras). baseMVA = 100. As barras sao numeradas 1..n de forma
consecutiva, com a slack na barra 1 (indice 0), compativel com o solver.

build_ybus monta a Ybus padrao incluindo carregamento de linha (b/2 em cada
extremidade), taps de transformador (modelo MATPOWER, defasagem nula) e shunts
de barra (Gs, Bs em MW/MVAr na base).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from EletroSolver import Barra, SistemaPotencia

BASE_MVA = 100.0

# --- WSCC 9 barras (MATPOWER case9) ---
# buses: (bus, type, Pd, Qd, Gs, Bs)
CASE9 = {
    "buses": [
        (1, 3, 0, 0, 0, 0),
        (2, 2, 0, 0, 0, 0),
        (3, 2, 0, 0, 0, 0),
        (4, 1, 0, 0, 0, 0),
        (5, 1, 90, 30, 0, 0),
        (6, 1, 0, 0, 0, 0),
        (7, 1, 100, 35, 0, 0),
        (8, 1, 0, 0, 0, 0),
        (9, 1, 125, 50, 0, 0),
    ],
    # gens: (bus, Pg, Vg)  -- Pg da slack (bus 1) e ignorado pelo solver
    "gens": [(1, 72.3, 1.04), (2, 163, 1.025), (3, 85, 1.025)],
    # branches: (fbus, tbus, r, x, b, ratio)  -- ratio 0 => linha (tap 1.0)
    "branches": [
        (1, 4, 0.0, 0.0576, 0.0, 0),
        (4, 5, 0.017, 0.092, 0.158, 0),
        (5, 6, 0.039, 0.17, 0.358, 0),
        (3, 6, 0.0, 0.0586, 0.0, 0),
        (6, 7, 0.0119, 0.1008, 0.209, 0),
        (7, 8, 0.0085, 0.072, 0.149, 0),
        (8, 2, 0.0, 0.0625, 0.0, 0),
        (8, 9, 0.032, 0.161, 0.306, 0),
        (9, 4, 0.01, 0.085, 0.176, 0),
    ],
}

# --- IEEE 14 barras (MATPOWER case14) ---
CASE14 = {
    "buses": [
        (1, 3, 0, 0, 0, 0),
        (2, 2, 21.7, 12.7, 0, 0),
        (3, 2, 94.2, 19, 0, 0),
        (4, 1, 47.8, -3.9, 0, 0),
        (5, 1, 7.6, 1.6, 0, 0),
        (6, 2, 11.2, 7.5, 0, 0),
        (7, 1, 0, 0, 0, 0),
        (8, 2, 0, 0, 0, 0),
        (9, 1, 29.5, 16.6, 0, 19),
        (10, 1, 9, 5.8, 0, 0),
        (11, 1, 3.5, 1.8, 0, 0),
        (12, 1, 6.1, 1.6, 0, 0),
        (13, 1, 13.5, 5.8, 0, 0),
        (14, 1, 14.9, 5, 0, 0),
    ],
    "gens": [
        (1, 232.4, 1.06),
        (2, 40, 1.045),
        (3, 0, 1.01),
        (6, 0, 1.07),
        (8, 0, 1.09),
    ],
    "branches": [
        (1, 2, 0.01938, 0.05917, 0.0528, 0),
        (1, 5, 0.05403, 0.22304, 0.0492, 0),
        (2, 3, 0.04699, 0.19797, 0.0438, 0),
        (2, 4, 0.05811, 0.17632, 0.034, 0),
        (2, 5, 0.05695, 0.17388, 0.0346, 0),
        (3, 4, 0.06701, 0.17103, 0.0128, 0),
        (4, 5, 0.01335, 0.04211, 0.0, 0),
        (4, 7, 0.0, 0.20912, 0.0, 0.978),
        (4, 9, 0.0, 0.55618, 0.0, 0.969),
        (5, 6, 0.0, 0.25202, 0.0, 0.932),
        (6, 11, 0.09498, 0.1989, 0.0, 0),
        (6, 12, 0.12291, 0.25581, 0.0, 0),
        (6, 13, 0.06615, 0.13027, 0.0, 0),
        (7, 8, 0.0, 0.17615, 0.0, 0),
        (7, 9, 0.0, 0.11001, 0.0, 0),
        (9, 10, 0.03181, 0.0845, 0.0, 0),
        (9, 14, 0.12711, 0.27038, 0.0, 0),
        (10, 11, 0.08205, 0.19207, 0.0, 0),
        (12, 13, 0.22092, 0.19988, 0.0, 0),
        (13, 14, 0.17093, 0.34802, 0.0, 0),
    ],
    # solucao MATPOWER (coluna Vm, Va[graus] da matriz de barras): (bus: (Vm, Va))
    "solucao": {
        1: (1.060, 0.00), 2: (1.045, -4.98), 3: (1.010, -12.72),
        4: (1.019, -10.33), 5: (1.020, -8.78), 6: (1.070, -14.22),
        7: (1.062, -13.37), 8: (1.090, -13.36), 9: (1.056, -14.94),
        10: (1.051, -15.10), 11: (1.057, -14.79), 12: (1.055, -15.07),
        13: (1.050, -15.16), 14: (1.036, -16.04),
    },
}


def build_ybus(n, branches, shunts):
    """Monta a Ybus padrao (modelo MATPOWER, defasagem nula)."""
    Y = np.zeros((n, n), dtype=complex)
    for fbus, tbus, r, x, b, ratio in branches:
        i, j = fbus - 1, tbus - 1
        ys = 1.0 / complex(r, x)
        bc = 1j * b / 2.0
        tap = ratio if ratio != 0 else 1.0
        Yff = (ys + bc) / (tap * tap)
        Ytt = ys + bc
        Yft = -ys / tap
        Ytf = -ys / tap
        Y[i, i] += Yff
        Y[j, j] += Ytt
        Y[i, j] += Yft
        Y[j, i] += Ytf
    for bus, gs, bs in shunts:
        Y[bus - 1, bus - 1] += complex(gs, bs) / BASE_MVA
    return Y


def build_system(case, tolerancia=1e-8, max_iter=50):
    """Constroi um SistemaPotencia a partir de um dict de caso IEEE."""
    n = len(case["buses"])
    tipo = [0] * n
    Pd = np.zeros(n)
    Qd = np.zeros(n)
    Pg = np.zeros(n)
    Vset = np.ones(n)
    for bus, tp, pd, qd, gs, bs in case["buses"]:
        i = bus - 1
        tipo[i] = tp
        Pd[i] = pd
        Qd[i] = qd
    for bus, pg, vg in case["gens"]:
        i = bus - 1
        Pg[i] += pg
        Vset[i] = vg

    barras = []
    for i in range(n):
        # PQ parte de flat start (V=1.0); slack/PV mantem o setpoint de tensao.
        V = 1.0 if tipo[i] == 1 else Vset[i]
        P = (Pg[i] - Pd[i]) / BASE_MVA       # injecao liquida ativa
        Q = (0.0 - Qd[i]) / BASE_MVA         # PQ: Qg=0; PV/slack: Q ignorado pelo solver
        barras.append(Barra(i, tipo[i], V, 0.0, P, Q))

    shunts = [(bus, gs, bs) for bus, tp, pd, qd, gs, bs in case["buses"]
              if gs != 0 or bs != 0]
    Y = build_ybus(n, case["branches"], shunts)
    return SistemaPotencia(barras, Y, tolerancia=tolerancia, max_iter=max_iter)
