"""Deteccao de (nao-)convergencia e ausencia do NameError do bug B2."""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from _helpers import sistema_5


def test_nao_convergencia_emite_warning():
    sistema = sistema_5(max_iter=1)
    with pytest.warns(UserWarning):
        sistema.calcular_fluxo()
    assert sistema.convergiu is False
    assert sistema.convergencia is None


def test_convergencia_primeira_iteracao_sem_nameerror():
    # Bug B2: convergir na 1a iteracao acessava J antes de defini-lo (NameError).
    sistema = sistema_5()
    sistema.calcular_fluxo()  # resolve o sistema
    assert sistema.convergiu

    # Re-rodar partindo da solucao converge na iteracao 1 — sem NameError.
    sistema.calcular_fluxo()
    assert sistema.convergiu is True
    assert sistema.convergencia == 1
    # ultimajacobiana foi atribuida nessa iteracao convergente.
    assert hasattr(sistema, "ultimajacobiana")
