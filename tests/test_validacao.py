"""Validacao de entrada em Barra e SistemaPotencia."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from EletroSolver import Barra, Linha, SistemaPotencia
from _helpers import barras_5, Y_5, sistema_5


# --- Barra ---

def test_barra_tipo_invalido():
    with pytest.raises(ValueError):
        Barra(0, 5, 1.0, 0.0, 0.0, 0.0)


def test_barra_tensao_nao_positiva():
    with pytest.raises(ValueError):
        Barra(0, 1, -1.0, 0.0, 0.0, 0.0)


# --- SistemaPotencia.__init__ ---

def test_slack_em_indice_arbitrario():
    # A slack pode estar em qualquer posicao: construir nao deve lancar erro
    # (a correcao fisica e verificada em test_slack_arbitraria.py).
    barras = barras_5()
    barras[0].tipo = 1  # tira a slack da posicao 0
    barras[2].tipo = 3  # coloca a slack no indice 2
    sistema = SistemaPotencia(barras, Y_5())
    assert sistema.slack_idx == 2


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


# --- linhas (dados de ramo) ---

def test_linha_indice_fora_de_faixa():
    with pytest.raises(ValueError):
        SistemaPotencia(barras_5(), Y_5(), linhas=[Linha(1, 99, z=0.1j)])


def test_linha_de_igual_para():
    with pytest.raises(ValueError):
        SistemaPotencia(barras_5(), Y_5(), linhas=[Linha(2, 2, z=0.1j)])


def test_linha_z_nula():
    with pytest.raises(ValueError):
        SistemaPotencia(barras_5(), Y_5(), linhas=[Linha(1, 2, z=0)])


def test_linha_tap_nulo():
    with pytest.raises(ValueError):
        SistemaPotencia(barras_5(), Y_5(), linhas=[Linha(1, 2, z=0.1j, tap=0)])


# --- alterar_barra preserva o invariante da slack ---

def test_alterar_barra_nao_remove_unica_slack():
    sistema = sistema_5()
    with pytest.raises(ValueError):
        sistema.alterar_barra(1, tipo=1)  # barra 1 e a slack
    assert sistema.barras[0].tipo == 3    # nada mudou
    assert sistema.slack_idx == 0


def test_alterar_barra_nao_cria_segunda_slack():
    sistema = sistema_5()
    with pytest.raises(ValueError):
        sistema.alterar_barra(2, tipo=3)
    assert sistema.barras[1].tipo == 1
    assert sistema.slack_idx == 0


def test_alterar_barra_tipo_valido_mantem_slack_idx():
    sistema = sistema_5()
    sistema.alterar_barra(4, tipo=1)      # PV -> PQ, permitido
    assert sistema.barras[3].tipo == 1
    assert sistema.slack_idx == 0


# --- metodos 1-based ---

def test_transito_indice_fora_de_faixa():
    sistema = sistema_5()
    with pytest.raises(ValueError):
        sistema.transito(1, 99)


def test_transito_de_igual_para():
    sistema = sistema_5()
    with pytest.raises(ValueError):
        sistema.transito(2, 2)
