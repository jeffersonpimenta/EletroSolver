"""Slack em posicao arbitraria: invariancia fisica por permutacao das barras.

Resolve o sistema de 5 barras com a slack no indice 0 (referencia) e uma versao
permutada que move a slack para outra posicao, permutando barras E linhas/colunas
de Y de forma consistente. Cada barra fisica deve convergir para o mesmo (V, theta)
nas duas solucoes — prova que a slack pode estar em qualquer indice.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from EletroSolver import Barra, SistemaPotencia
from _helpers import barras_5, Y_5


def _permutar(barras, Y, perm):
    """Reconstroi (barras, Y) na nova ordem: nova posicao k <- antiga perm[k]."""
    novas = []
    for k, o in enumerate(perm):
        b = barras[o]
        novas.append(Barra(k, b.tipo, b.V, b.theta, b.P, b.Q))
    Y = np.asarray(Y)
    Y_perm = Y[np.ix_(perm, perm)]
    return novas, Y_perm


def test_slack_em_qualquer_posicao_da_o_mesmo_resultado():
    # Solucao de referencia: slack no indice 0.
    ref = SistemaPotencia(barras_5(), Y_5(), tolerancia=1e-10)
    ref.calcular_fluxo()
    assert ref.convergiu
    V_ref = [b.V for b in ref.barras]
    th_ref = [b.theta for b in ref.barras]

    # Permutacao que leva a slack (antigo indice 0) para o indice 2.
    perm = [1, 2, 0, 3, 4]
    barras_p, Y_p = _permutar(barras_5(), Y_5(), perm)

    sistema = SistemaPotencia(barras_p, Y_p, tolerancia=1e-10)
    assert sistema.slack_idx == 2  # slack agora no indice 2
    sistema.calcular_fluxo()
    assert sistema.convergiu

    # A barra na nova posicao k corresponde a antiga perm[k]: mesmo V e theta.
    for k, o in enumerate(perm):
        assert np.isclose(sistema.barras[k].V, V_ref[o], atol=1e-8), (
            f"V difere na barra fisica {o} (nova posicao {k})"
        )
        assert np.isclose(sistema.barras[k].theta, th_ref[o], atol=1e-8), (
            f"theta difere na barra fisica {o} (nova posicao {k})"
        )
