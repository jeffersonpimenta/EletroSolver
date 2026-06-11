"""Transito de potencia e perdas: convencao de sinal, modelo de ramo e balanco.

Cobre a correcao do sinal de transito() (Y[i, j] = -y_serie na convencao Ybus),
as perdas com sinal (Q_loss negativo = linha gera reativo) e o modo exato com
dados de ramo (tap de transformador e carregamento shunt, modelo MATPOWER).
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from EletroSolver import Barra, Linha, SistemaPotencia
from _helpers import sistema_5
from _ieee_data import BASE_MVA, CASE14, build_system


def _sistema_2barras(linhas=None):
    """Slack alimentando carga de 0.5 + j0.2 pu por uma linha 0.01 + j0.1."""
    barras = [
        Barra(0, 3, 1.0, 0.0, 0.0, 0.0),
        Barra(1, 1, 1.0, 0.0, -0.5, -0.2),
    ]
    z = 0.01 + 0.1j
    y = 1 / z
    Y = np.array([[y, -y], [-y, y]], dtype=complex)
    sistema = SistemaPotencia(barras, Y, linhas=linhas)
    sistema.calcular_fluxo()
    return sistema


# --- Convencao de sinal (sem dados de ramo) ---

def test_sinal_transito_slack_para_carga():
    sistema = _sistema_2barras()
    t12 = sistema.transito(1, 2)
    t21 = sistema.transito(2, 1)

    # A potencia sai da slack (positiva) e chega na carga (negativa).
    assert t12["P_ij"] > 0
    assert t21["P_ij"] < 0
    # Na barra de carga o fluxo que sai e exatamente a injecao especificada
    # (nao ha outro ramo): P_21 = -0.5, Q_21 = -0.2.
    assert np.isclose(t21["P_ij"], -0.5, atol=1e-6)
    assert np.isclose(t21["Q_ij"], -0.2, atol=1e-6)
    # O lado da fonte entrega a carga mais as perdas.
    assert t12["P_ij"] > 0.5


def test_perda_eh_soma_com_sinal_e_bate_com_rI2():
    sistema = _sistema_2barras()
    perdas = sistema.losses(1, 2)
    t12 = sistema.transito(1, 2)
    t21 = sistema.transito(2, 1)

    assert np.isclose(perdas["P_loss"], t12["P_ij"] + t21["P_ij"])
    assert perdas["P_loss"] > 0

    # Referencia independente: perda serie = |I|^2 * R, com I = y (V1 - V2).
    V1 = sistema.barras[0].V * np.exp(1j * sistema.barras[0].theta)
    V2 = sistema.barras[1].V * np.exp(1j * sistema.barras[1].theta)
    I = (V1 - V2) / (0.01 + 0.1j)
    assert np.isclose(perdas["P_loss"], 0.01 * abs(I) ** 2)
    assert np.isclose(perdas["Q_loss"], 0.1 * abs(I) ** 2)


def test_balanco_total_sem_abs():
    # Soma das injecoes = perdas totais (com sinal, sem abs mascarando).
    sistema = sistema_5()
    sistema.calcular_fluxo()
    V = np.array([b.V for b in sistema.barras])
    theta = np.array([b.theta for b in sistema.barras])
    P_calc, Q_calc = sistema.fluxo_potencia(V, theta)
    perdas = sistema.totlosses()
    assert np.isclose(np.sum(P_calc), perdas["P_loss"], atol=1e-6)
    assert np.isclose(np.sum(Q_calc), perdas["Q_loss"], atol=1e-6)


# --- Modo exato com dados de ramo (tap e carregamento) ---

def test_transito_com_linhas_equivale_ao_ybus_em_linha_simples():
    # Sem tap nem charging, o modelo de ramo deve reproduzir o calculo via Ybus.
    linhas = [Linha(1, 2, z=0.01 + 0.1j)]
    com = _sistema_2barras(linhas=linhas).transito(1, 2)
    sem = _sistema_2barras().transito(1, 2)
    assert np.isclose(com["S_ij"], sem["S_ij"])


def test_qloss_negativo_com_carregamento_e_carga_leve():
    # Linha com charging alto e quase sem carga: gera mais reativo do que perde.
    barras = [
        Barra(0, 3, 1.0, 0.0, 0.0, 0.0),
        Barra(1, 1, 1.0, 0.0, -0.01, 0.0),
    ]
    linhas = [Linha(1, 2, z=0.01 + 0.1j, b=0.4)]
    z = 0.01 + 0.1j
    y = 1 / z
    bc = 1j * 0.4 / 2
    Y = np.array([[y + bc, -y], [-y, y + bc]], dtype=complex)
    sistema = SistemaPotencia(barras, Y, linhas=linhas)
    sistema.calcular_fluxo()

    perdas = sistema.losses(1, 2)
    assert perdas["Q_loss"] < 0
    assert perdas["P_loss"] > 0


def test_balanco_por_barra_ieee14_com_taps_e_shunts():
    """Invariante forte: em cada barra, injecao = soma dos fluxos + shunt da barra.

    So fecha com o modelo de ramo (tap e charging por linha); o calculo via
    Ybus nao fecharia nas barras com transformador ou carregamento.
    """
    base = build_system(CASE14)
    linhas = [
        Linha(f, t, z=complex(r, x), b=b, tap=(ratio if ratio != 0 else 1.0))
        for f, t, r, x, b, ratio in CASE14["branches"]
    ]
    sistema = SistemaPotencia(base.barras, base.Y, tolerancia=1e-8,
                              max_iter=50, linhas=linhas)
    sistema.calcular_fluxo()
    assert sistema.convergiu

    V = np.array([b.V for b in sistema.barras])
    theta = np.array([b.theta for b in sistema.barras])
    P_calc, Q_calc = sistema.fluxo_potencia(V, theta)

    # Shunts de barra (gs, bs em MW/MVAr a V=1): S_sh = |V|^2 * conj(y_sh).
    shunt = {bus: complex(gs, bs) / BASE_MVA
             for bus, tp, pd, qd, gs, bs in CASE14["buses"]
             if gs != 0 or bs != 0}

    vizinhos = {}
    for lin in linhas:
        vizinhos.setdefault(lin.de, set()).add(lin.para)
        vizinhos.setdefault(lin.para, set()).add(lin.de)

    for barra in range(1, sistema.n_barras + 1):
        S_flux = sum(sistema.transito(barra, nb)["S_ij"]
                     for nb in vizinhos[barra])
        S_sh = (abs(V[barra - 1]) ** 2) * np.conj(shunt.get(barra, 0.0))
        assert np.isclose(P_calc[barra - 1], S_flux.real + S_sh.real,
                          atol=1e-8), f"balanco de P falhou na barra {barra}"
        assert np.isclose(Q_calc[barra - 1], S_flux.imag + S_sh.imag,
                          atol=1e-8), f"balanco de Q falhou na barra {barra}"


def test_transito_sem_ramo_entre_barras_levanta_erro():
    linhas = [Linha(1, 2, z=0.01 + 0.1j)]
    barras = [
        Barra(0, 3, 1.0, 0.0, 0.0, 0.0),
        Barra(1, 1, 1.0, 0.0, -0.1, 0.0),
        Barra(2, 1, 1.0, 0.0, 0.0, 0.0),
    ]
    z = 0.01 + 0.1j
    y = 1 / z
    Y = np.zeros((3, 3), dtype=complex)
    Y[0, 0] += y
    Y[1, 1] += y
    Y[0, 1] -= y
    Y[1, 0] -= y
    Y[2, 2] = 1.0  # barra 3 isolada com shunt so para nao singularizar
    sistema = SistemaPotencia(barras, Y, linhas=linhas)
    with pytest.raises(ValueError):
        sistema.transito(1, 3)


# --- Valores resolvidos gravados nas barras (slack e PV) ---

def test_fluxo_atualiza_pq_da_slack_e_q_do_pv():
    sistema = sistema_5()
    # Entradas sabidamente "lixo" nas incognitas: Q=5 na slack (helper).
    assert sistema.barras[0].Q == 5
    sistema.calcular_fluxo()

    V = np.array([b.V for b in sistema.barras])
    theta = np.array([b.theta for b in sistema.barras])
    P_calc, Q_calc = sistema.fluxo_potencia(V, theta)

    slack = sistema.barras[sistema.slack_idx]
    assert np.isclose(slack.P, P_calc[sistema.slack_idx], atol=1e-9)
    assert np.isclose(slack.Q, Q_calc[sistema.slack_idx], atol=1e-9)
    for i, barra in enumerate(sistema.barras):
        if barra.tipo == 2:
            assert np.isclose(barra.Q, Q_calc[i], atol=1e-9)
