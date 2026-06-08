"""Testes de propriedades do motor de curto-circuito (sem referencia externa)."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Faltas import A, EstudoCurtoCircuito
from _seq_data import estudo_1barra, estudo_2barras, injecao_kcl


def test_matriz_A_inversivel():
    # A @ inv(A) == I
    assert np.allclose(A @ np.linalg.inv(A), np.eye(3))


def test_trifasica_igual_thevenin_positiva():
    est = estudo_2barras(incluir_cargas=False)
    res = est.falta_trifasica(2)
    I1 = res["I_seq"][1]
    assert np.isclose(I1, est.Vpre[1] / est.Z1[1, 1])
    # so sequencia positiva
    assert np.isclose(res["I_seq"][0], 0.0)
    assert np.isclose(res["I_seq"][2], 0.0)


def test_trifasica_bolted_zera_tensao_na_falta():
    est = estudo_2barras(incluir_cargas=False)
    res = est.falta_trifasica(2)
    V0, V1, V2 = res["V_seq"]
    assert np.isclose(V1[1], 0.0, atol=1e-12)
    # tensoes de fase na barra em falta tambem nulas
    assert np.allclose(res["V_fase"][1], 0.0, atol=1e-12)


def test_monofasica_sequencias_iguais():
    est = estudo_2barras(incluir_cargas=False)
    res = est.falta_monofasica(2)
    I0, I1, I2 = res["I_seq"]
    assert np.isclose(I0, I1) and np.isclose(I1, I2)
    # Ia = 3*I1 ; Ib = Ic = 0
    assert np.isclose(res["I_fase"]["a"], 3 * I1)
    assert np.isclose(res["I_fase"]["b"], 0.0, atol=1e-12)
    assert np.isclose(res["I_fase"]["c"], 0.0, atol=1e-12)


def test_bifasica_fase_a_nula_e_b_oposta_c():
    est = estudo_2barras(incluir_cargas=False)
    res = est.falta_bifasica(2)
    assert np.isclose(res["I_seq"][0], 0.0)            # sem sequencia zero
    assert np.isclose(res["I_fase"]["a"], 0.0, atol=1e-12)
    assert np.isclose(res["I_fase"]["b"], -res["I_fase"]["c"])


def test_bifasica_terra_fase_a_nula():
    est = estudo_2barras(incluir_cargas=False)
    res = est.falta_bifasica_terra(2)
    assert np.isclose(res["I_fase"]["a"], 0.0, atol=1e-12)
    # corrente de terra = Ia+Ib+Ic = 3*I0
    Ig = sum(res["I_fase"][f] for f in ("a", "b", "c"))
    assert np.isclose(Ig, 3 * res["I_seq"][0])


def test_zbus_simetrica():
    est = estudo_2barras(incluir_cargas=False)
    assert np.allclose(est.Z1, est.Z1.T)
    assert np.allclose(est.Z2, est.Z2.T)


@pytest.mark.parametrize("tipo", ["trifasica", "monofasica", "bifasica",
                                  "bifasica_terra"])
def test_kcl_contribuicoes_somam_corrente_de_falta(tipo):
    # Soma das contribuicoes (gerador + carga + ramos) = corrente de falta na
    # barra k e ~0 nas demais (sequencia positiva).
    est = estudo_2barras(incluir_cargas=True,
                         cargas=_carga_simples())
    res = getattr(est, f"falta_{tipo}")(2)
    inj = injecao_kcl(est, res, seq=1)
    esperado = np.zeros(est.n, dtype=complex)
    esperado[1] = res["I_seq"][1]
    assert np.allclose(inj, esperado, atol=1e-9)


def test_prefault_flat_e_um():
    est = estudo_1barra()
    assert np.allclose(est.Vpre, 1.0 + 0j)


def _carga_simples():
    from Faltas import Carga
    return [Carga(2, P=0.5, Q=0.2)]
