"""Builders e casos de sequencia para os testes de curto-circuito.

Sistemas pequenos com resposta analitica (verificavel a mao), usados como
referencia para os testes de falta. Tudo em pu, base 100 MVA.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from EletroSolver import Barra, SistemaPotencia
from Faltas import Carga, EstudoCurtoCircuito, Gerador, Ramo


def sistema_dummy(n):
    """SistemaPotencia minimo (slack na barra 1) so para fornecer n e Sbase.

    A Ybus aqui nao e usada pela analise de falta (que monta as suas proprias
    redes de sequencia); serve apenas para satisfazer o construtor.
    """
    barras = [Barra(0, 3, 1.0, 0.0, 0.0, 0.0)]
    barras += [Barra(i, 1, 1.0, 0.0, 0.0, 0.0) for i in range(1, n)]
    Y = np.zeros((n, n), dtype=complex)
    return SistemaPotencia(barras, Y)


def estudo_1barra(aterrado=True, X0=0.05):
    """Gerador unico atras de X''d: caso analitico classico.

    X1 = X2 = 0.2, X0 = 0.05, neutro solidamente aterrado (Zn=0).
    3f:  |I| = 1/0.2          = 5.0
    FT:  |Ia| = 3/(0.2+0.2+0.05) = 6.6667
    FF:  |Ib| = sqrt(3)/(0.2+0.2) = 4.3301
    """
    sis = sistema_dummy(1)
    ger = [Gerador(1, X1=0.2, X2=0.2, X0=X0, aterrado=aterrado)]
    return EstudoCurtoCircuito(sis, ger, ramos=[])


def estudo_2barras(cargas=None, incluir_cargas=True):
    """Gerador na barra 1, linha 1-2. Thevenin em 2: Z1 = 0.1 + 0.1 = 0.2.

    Gerador: X1=X2=0.1, X0=0.05 (aterrado). Linha: z1=z2=0.1j, z0=0.3j.
    """
    sis = sistema_dummy(2)
    ger = [Gerador(1, X1=0.1, X2=0.1, X0=0.05)]
    ramos = [Ramo(1, 2, z1=0.1j, z2=0.1j, z0=0.3j)]
    return EstudoCurtoCircuito(sis, ger, ramos, cargas=cargas,
                               incluir_cargas=incluir_cargas)


def injecao_kcl(estudo, res, seq=1):
    """Reconstroi a injecao de corrente de falta por barra a partir das
    contribuicoes (gerador + carga + ramos), para a sequencia `seq`.

    Deve resultar em I_seq na barra em falta e ~0 nas demais (KCL).
    """
    n = estudo.n
    inj = np.zeros(n, dtype=complex)
    for g in res["contrib_gerador"]:
        inj[g["barra"] - 1] += g["I_seq"][seq]
    for c in res["contrib_carga"]:
        inj[c["barra"] - 1] += c["I_seq"][seq]
    for s in res.get("contrib_shunt_linha", []):
        inj[s["barra"] - 1] += s["I_seq"][seq]
    for lin, r in zip(res["contrib_linha"], estudo.ramos):
        I = lin["I_seq"][seq]   # corrente de->para (sai de 'de', chega em 'para')
        inj[r.de - 1] -= I
        inj[r.para - 1] += I
    return inj
