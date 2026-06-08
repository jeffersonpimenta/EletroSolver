"""Invariantes fisicos do fluxo de potencia no sistema de 5 barras.

Nao dependem de valores 'golden' arbitrarios: validam que a solucao satisfaz as
equacoes de fluxo (mismatch ~ 0) e o balanco de potencia (injecoes = perdas).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from _helpers import sistema_5


def test_converge():
    sistema = sistema_5()
    sistema.calcular_fluxo()
    assert sistema.convergiu is True
    assert sistema.convergencia is not None
    # Newton verdadeiro: convergencia rapida (bem abaixo do max_iter padrao).
    assert sistema.convergencia <= 6


def test_mismatch_nulo_na_solucao():
    sistema = sistema_5()
    sistema.calcular_fluxo()

    V = np.array([b.V for b in sistema.barras])
    theta = np.array([b.theta for b in sistema.barras])
    P_calc, Q_calc = sistema.fluxo_potencia(V, theta)
    dP, dQ = sistema.calcular_desvios(P_calc, Q_calc)
    mismatch = np.concatenate([dP, dQ])

    assert np.linalg.norm(mismatch) < sistema.tolerancia


def test_balanco_de_potencia():
    # Soma das injecoes liquidas de potencia ativa = perdas ativas totais da rede.
    sistema = sistema_5()
    sistema.calcular_fluxo()

    V = np.array([b.V for b in sistema.barras])
    theta = np.array([b.theta for b in sistema.barras])
    P_calc, _ = sistema.fluxo_potencia(V, theta)

    perdas = sistema.totlosses()
    assert np.isclose(np.sum(P_calc), perdas["P_loss"], atol=1e-3)
