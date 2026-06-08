"""Validacao contra casos IEEE/MATPOWER padrao (case9 WSCC, case14 IEEE).

- O mismatch ~ 0 prova que o solver resolve corretamente as equacoes de fluxo
  para a Ybus construida.
- A comparacao do case14 com a solucao publicada do MATPOWER prova que a
  construcao da Ybus (carregamento de linha, taps de transformador e shunt de
  barra) esta correta.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from _ieee_data import CASE9, CASE14, build_system


def _mismatch(sistema):
    V = np.array([b.V for b in sistema.barras])
    theta = np.array([b.theta for b in sistema.barras])
    P_calc, Q_calc = sistema.fluxo_potencia(V, theta)
    dP, dQ = sistema.calcular_desvios(P_calc, Q_calc)
    return np.linalg.norm(np.concatenate([dP, dQ]))


def test_case9_resolve_e_consistente():
    sistema = build_system(CASE9)
    sistema.calcular_fluxo()

    assert sistema.convergiu is True
    assert sistema.convergencia <= 6              # Newton verdadeiro: convergencia rapida
    assert _mismatch(sistema) < 1e-6

    # Setpoints de tensao mantidos: slack (bus 1) e PV (buses 2, 3).
    assert abs(sistema.barras[0].V - 1.04) < 1e-9
    assert abs(sistema.barras[1].V - 1.025) < 1e-9
    assert abs(sistema.barras[2].V - 1.025) < 1e-9
    assert abs(sistema.barras[0].theta) < 1e-12   # referencia angular

    # Tensoes plausiveis em todas as barras.
    assert all(0.9 < b.V < 1.1 for b in sistema.barras)


def test_case14_resolve():
    sistema = build_system(CASE14)
    sistema.calcular_fluxo()
    assert sistema.convergiu is True
    assert sistema.convergencia <= 6
    assert _mismatch(sistema) < 1e-6


def test_case14_bate_com_matpower():
    sistema = build_system(CASE14)
    sistema.calcular_fluxo()

    for bus, (vm, va) in CASE14["solucao"].items():
        barra = sistema.barras[bus - 1]
        assert abs(barra.V - vm) < 3e-3, f"|V| da barra {bus} fora de tolerancia"
        assert abs(np.degrees(barra.theta) - va) < 0.1, (
            f"angulo da barra {bus} fora de tolerancia"
        )
