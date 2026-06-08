"""Validacao numerica do curto-circuito contra valores analiticos exatos."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Faltas import Carga, EstudoCurtoCircuito, Gerador, Ramo
from _seq_data import estudo_1barra, estudo_2barras, injecao_kcl, sistema_dummy


# --- Caso 1 barra (gerador atras de X''d), valores analiticos ---

def test_1barra_trifasica():
    est = estudo_1barra()
    res = est.falta_trifasica(1)
    assert np.isclose(abs(res["I_fase"]["a"]), 5.0)


def test_1barra_monofasica():
    est = estudo_1barra()
    res = est.falta_monofasica(1)
    assert np.isclose(abs(res["I_fase"]["a"]), 3 / 0.45)   # 6.6667


def test_1barra_bifasica():
    est = estudo_1barra()
    res = est.falta_bifasica(1)
    assert np.isclose(abs(res["I_fase"]["b"]), np.sqrt(3) / 0.4)  # 4.3301
    assert np.isclose(res["I_fase"]["b"], -res["I_fase"]["c"])


def test_1barra_corrente_kA_e_potencia():
    est = estudo_1barra()
    res = est.falta_trifasica(1)
    # base 100 MVA, 13.8 kV: Ibase = 100e3/(sqrt(3)*13.8) A ~ 4.184 kA
    ika = est.corrente_kA(res["I_fase"]["a"], kV_base=13.8)
    assert np.isclose(ika, 5.0 * 100 / (np.sqrt(3) * 13.8))
    assert np.isclose(est.potencia_curto_MVA(1), 5.0 * 100)  # 500 MVA


# --- Caso 2 barras (gerador + linha), Thevenin Z1=0.2 na barra 2 ---

def test_2barras_trifasica_e_contribuicao_da_linha():
    est = estudo_2barras(incluir_cargas=False)
    res = est.falta_trifasica(2)
    assert np.isclose(abs(res["I_fase"]["a"]), 5.0)        # 1/0.2
    # rede radial: a linha 1->2 conduz toda a corrente de falta
    Ilinha = res["contrib_linha"][0]["I_seq"][1]
    assert np.isclose(Ilinha, res["I_seq"][1])
    # o gerador a fornece integralmente
    Iger = res["contrib_gerador"][0]["I_seq"][1]
    assert np.isclose(Iger, res["I_seq"][1])


def test_2barras_tensao_intermediaria():
    est = estudo_2barras(incluir_cargas=False)
    res = est.falta_trifasica(2)
    V0, V1, V2 = res["V_seq"]
    # barra 1 cai para a fracao da reatancia do gerador: 0.1/0.2 = 0.5
    assert np.isclose(V1[0], 0.5)
    assert np.isclose(V1[1], 0.0, atol=1e-12)


# --- Influencia das cargas ---

def test_carga_altera_thevenin_e_corrente():
    sem = estudo_2barras(incluir_cargas=False)
    com = estudo_2barras(incluir_cargas=True, cargas=[Carga(2, P=1.0, Q=0.5)])
    # carga shunt em paralelo reduz a impedancia de Thevenin -> maior corrente
    assert abs(com.Z1[1, 1]) < abs(sem.Z1[1, 1])
    i_sem = abs(sem.falta_trifasica(2)["I_fase"]["a"])
    i_com = abs(com.falta_trifasica(2)["I_fase"]["a"])
    assert i_com > i_sem


def test_contribuicao_de_carga_presente():
    est = estudo_2barras(incluir_cargas=True, cargas=[Carga(1, P=0.8, Q=0.3)])
    res = est.falta_trifasica(2)
    barras_carga = [c["barra"] for c in res["contrib_carga"]]
    assert 1 in barras_carga
    # carga na barra 1 (que nao colapsa totalmente) contribui corrente nao nula
    contrib = next(c for c in res["contrib_carga"] if c["barra"] == 1)
    assert abs(contrib["I_seq"][1]) > 0


def test_desligar_cargas_recupera_sem_carga():
    com_flag = estudo_2barras(incluir_cargas=False,
                              cargas=[Carga(2, P=1.0, Q=0.5)])
    sem = estudo_2barras(incluir_cargas=False)
    assert np.allclose(com_flag.Z1, sem.Z1)


# --- Carregamento capacitivo de linha (shunt) ---

def test_carregamento_linha_altera_corrente():
    ger = [Gerador(1, X1=0.1, X2=0.1, X0=0.05)]
    sem = EstudoCurtoCircuito(
        sistema_dummy(2), ger, [Ramo(1, 2, z1=0.1j, z2=0.1j, z0=0.3j)])
    com = EstudoCurtoCircuito(
        sistema_dummy(2), ger,
        [Ramo(1, 2, z1=0.1j, z2=0.1j, z0=0.3j, b1=0.5, b0=0.2)])
    i_sem = abs(sem.falta_trifasica(2)["I_fase"]["a"])
    i_com = abs(com.falta_trifasica(2)["I_fase"]["a"])
    assert not np.isclose(i_sem, i_com)


def test_carregamento_linha_aparece_nas_contribuicoes_e_mantem_kcl():
    ger = [Gerador(1, X1=0.1, X2=0.1, X0=0.05)]
    est = EstudoCurtoCircuito(
        sistema_dummy(2), ger,
        [Ramo(1, 2, z1=0.1j, z2=0.1j, z0=0.3j, b1=0.5, b0=0.2)])
    res = est.falta_monofasica(2)
    # o shunt de linha contribui corrente nao nula
    assert res["contrib_shunt_linha"]
    assert any(abs(s["I_seq"][1]) > 0 for s in res["contrib_shunt_linha"])
    # KCL continua exata em sequencia positiva e zero
    for seq in (0, 1, 2):
        inj = injecao_kcl(est, res, seq=seq)
        esperado = np.zeros(2, dtype=complex)
        esperado[1] = res["I_seq"][seq]
        assert np.allclose(inj, esperado, atol=1e-9)


def test_desligar_shunt_linha_ignora_b():
    ger = [Gerador(1, X1=0.1, X2=0.1, X0=0.05)]
    com_b_off = EstudoCurtoCircuito(
        sistema_dummy(2), ger,
        [Ramo(1, 2, z1=0.1j, z2=0.1j, z0=0.3j, b1=0.5, b0=0.2)],
        incluir_shunt_linha=False)
    sem_b = EstudoCurtoCircuito(
        sistema_dummy(2), ger, [Ramo(1, 2, z1=0.1j, z2=0.1j, z0=0.3j)])
    assert np.allclose(com_b_off.Z1, sem_b.Z1)
    assert com_b_off.falta_trifasica(2)["contrib_shunt_linha"] == []


# --- Pre-falta a partir do fluxo ---

def test_prefault_fluxo_usa_tensao_convergida():
    from EletroSolver import Barra, SistemaPotencia
    from Faltas import EstudoCurtoCircuito, Gerador, Ramo
    # sistema 2 barras com fluxo trivial; marca convergiu e tensoes.
    sis = SistemaPotencia(
        [Barra(0, 3, 1.05, 0.0, 0.0, 0.0), Barra(1, 1, 1.0, 0.0, 0.0, 0.0)],
        np.zeros((2, 2), dtype=complex),
    )
    sis.convergiu = True
    sis.barras[1].V = 0.97
    sis.barras[1].theta = -0.1
    est = EstudoCurtoCircuito(
        sis, [Gerador(1, X1=0.1, X2=0.1, X0=0.05)],
        [Ramo(1, 2, z1=0.1j, z2=0.1j, z0=0.3j)], prefault="fluxo",
    )
    assert np.isclose(est.Vpre[1], 0.97 * np.exp(-0.1j))
    res = est.falta_trifasica(2)
    assert np.isclose(res["Vf"], 0.97 * np.exp(-0.1j))
