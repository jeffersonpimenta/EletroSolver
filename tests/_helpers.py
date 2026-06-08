"""Builders compartilhados pelos testes (sistema de 5 barras do Simplificado.py)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from EletroSolver import Barra, SistemaPotencia


def barras_5():
    return [
        Barra(0, 3, 1.02, 0.0, -0.5327, 5),
        Barra(1, 1, 1.0, 0.0, 0.6, 0.25),
        Barra(2, 1, 1.0, 0.0, 0.4, 0.15),
        Barra(3, 2, 1.01, 0.0, -0.21, 0.15),
        Barra(4, 1, 1.0, 0.0, 1.0, 0.4),
    ]


def Y_5():
    Y = np.zeros((5, 5), dtype=complex)
    Y[0, 1] = Y[1, 0] = -1 / (0.05 + 0.2j)
    Y[0, 3] = Y[3, 0] = -1 / (0.05 + 0.25j)
    Y[1, 2] = Y[2, 1] = -1 / (0.05 + 0.15j)
    Y[1, 4] = Y[4, 1] = -1 / (0.1 + 0.3j)
    Y[2, 3] = Y[3, 2] = -1 / (0.05 + 0.2j)
    Y[3, 4] = Y[4, 3] = -1 / (0.05 + 0.1j)
    for i in range(5):
        Y[i, i] = -np.sum(Y[i, :])
    return Y


def sistema_5(**kwargs):
    return SistemaPotencia(barras_5(), Y_5(), **kwargs)
