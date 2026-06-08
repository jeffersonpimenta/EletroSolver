"""Testes de bordas e tratamento de erros do curto-circuito."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Faltas import EstudoCurtoCircuito, Gerador, Ramo
from _seq_data import estudo_1barra, estudo_2barras, sistema_dummy


def test_exige_gerador():
    with pytest.raises(ValueError, match="gerador"):
        EstudoCurtoCircuito(sistema_dummy(2), [], [Ramo(1, 2, z1=0.1j)])


def test_prefault_invalido():
    with pytest.raises(ValueError, match="prefault"):
        EstudoCurtoCircuito(sistema_dummy(1), [Gerador(1, X1=0.2)], [],
                            prefault="xyz")


def test_indice_de_falta_invalido():
    est = estudo_1barra()
    with pytest.raises(ValueError, match="faixa"):
        est.falta_trifasica(5)


def test_prefault_fluxo_sem_convergir():
    # sistema com convergiu == False -> Vpre via fluxo deve falhar na construcao.
    with pytest.raises(ValueError, match="convergido"):
        EstudoCurtoCircuito(sistema_dummy(1), [Gerador(1, X1=0.2)], [],
                            prefault="fluxo")


def test_ramo_de_igual_para():
    with pytest.raises(ValueError, match="'de' == 'para'"):
        EstudoCurtoCircuito(sistema_dummy(2), [Gerador(1, X1=0.1)],
                            [Ramo(1, 1, z1=0.1j)])


def test_trifasica_e_bifasica_funcionam_sem_dados_de_seq_zero():
    # ramo sem z0; trifasica e bifasica nao tocam a rede de sequencia zero.
    est = EstudoCurtoCircuito(
        sistema_dummy(2), [Gerador(1, X1=0.1, X2=0.1)],
        [Ramo(1, 2, z1=0.1j, z2=0.1j)],  # z0=None
    )
    assert np.isclose(abs(est.falta_trifasica(2)["I_fase"]["a"]), 5.0)
    assert est.falta_bifasica(2) is not None


def test_monofasica_sem_z0_na_linha_levanta_erro():
    est = EstudoCurtoCircuito(
        sistema_dummy(2), [Gerador(1, X1=0.1, X2=0.1, X0=0.05)],
        [Ramo(1, 2, z1=0.1j, z2=0.1j)],  # z0 ausente
    )
    with pytest.raises(ValueError, match="z0 ausente"):
        est.falta_monofasica(2)


def test_sequencia_zero_sem_terra_e_singular():
    # gerador nao aterrado e nenhum outro caminho de terra -> Y0 singular.
    est = estudo_1barra(aterrado=False)
    # trifasica continua funcionando
    assert np.isclose(abs(est.falta_trifasica(1)["I_fase"]["a"]), 5.0)
    with pytest.raises(np.linalg.LinAlgError, match="singular"):
        est.falta_monofasica(1)


def test_gerador_nao_aterrado_nao_contribui_seq_zero():
    # Com um caminho de terra externo (trafo), o gerador nao aterrado nao deve
    # aparecer na sequencia zero. Aqui basta verificar Y0 finita e Ig0 = 0.
    est = EstudoCurtoCircuito(
        sistema_dummy(2),
        # gerador 1 nao aterrado; gerador 2 aterrado da caminho de terra a barra 2
        [Gerador(1, X1=0.1, X2=0.1, X0=0.05, aterrado=False),
         Gerador(2, X1=0.1, X2=0.1, X0=0.05, aterrado=True)],
        # trafo YNd aterra o lado 'de' (barra 1) -> caminho de terra em seq zero
        [Ramo(1, 2, z1=0.1j, z2=0.1j, z0=0.3j, ligacao="YNd")],
    )
    res = est.falta_monofasica(1)
    # gerador 1 (nao aterrado) nao contribui em sequencia zero
    contrib_g1 = next(g for g in res["contrib_gerador"] if g["barra"] == 1)
    assert np.isclose(contrib_g1["I_seq"][0], 0.0)
