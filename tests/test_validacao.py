"""Validacao de entrada em Barra e SistemaPotencia."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from EletroSolver import Barra, SistemaPotencia
from _helpers import barras_5, Y_5, sistema_5


# --- Barra ---

def test_barra_tipo_invalido():
    with pytest.raises(ValueError):
        Barra(0, 5, 1.0, 0.0, 0.0, 0.0)


def test_barra_tensao_nao_positiva():
    with pytest.raises(ValueError):
        Barra(0, 1, -1.0, 0.0, 0.0, 0.0)


# --- SistemaPotencia.__init__ ---

def test_slack_fora_do_indice_zero():
    barras = barras_5()
    barras[0].tipo = 1  # tira a slack da posicao 0
    barras[1].tipo = 3  # coloca em outra posicao
    with pytest.raises(ValueError):
        SistemaPotencia(barras, Y_5())


def test_sem_slack():
    barras = barras_5()
    barras[0].tipo = 1
    with pytest.raises(ValueError):
        SistemaPotencia(barras, Y_5())


def test_multiplas_slacks():
    barras = barras_5()
    barras[1].tipo = 3
    with pytest.raises(ValueError):
        SistemaPotencia(barras, Y_5())


def test_Y_dimensao_errada():
    with pytest.raises(ValueError):
        SistemaPotencia(barras_5(), np.zeros((4, 4), dtype=complex))


def test_tolerancia_invalida():
    with pytest.raises(ValueError):
        sistema_5(tolerancia=0)


def test_max_iter_invalido():
    with pytest.raises(ValueError):
        sistema_5(max_iter=0)


def test_sbase_invalido():
    with pytest.raises(ValueError):
        sistema_5(Sbase=0)


# --- metodos 1-based ---

def test_transito_indice_fora_de_faixa():
    sistema = sistema_5()
    with pytest.raises(ValueError):
        sistema.transito(1, 99)


def test_transito_de_igual_para():
    sistema = sistema_5()
    with pytest.raises(ValueError):
        sistema.transito(2, 2)
