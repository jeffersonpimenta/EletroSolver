"""Exemplo de uso do modulo de curto-circuito (Faltas.py).

Sistema didatico de 3 barras: geradores nas barras 1 e 3, linhas 1-2 e 2-3,
carga na barra 2. Calcula as quatro faltas na barra 2 e mostra as
contribuicoes de geradores, linhas e carga.

Base 100 MVA. Tudo em pu. Tensao base de linha assumida 138 kV para a
conversao em kA.
"""
import numpy as np

from EletroSolver import Barra, SistemaPotencia
from Faltas import Carga, EstudoCurtoCircuito, Gerador, Ramo

BASE_MVA = 100.0
KV_BASE = 138.0

# --- Dados de sequencia ------------------------------------------------------
geradores = [
    Gerador(1, X1=0.15, X2=0.15, X0=0.05, aterrado=True),
    Gerador(3, X1=0.20, X2=0.20, X0=0.08, aterrado=True),
]
ramos = [
    # b1/b0 = carregamento capacitivo shunt da linha (sequencia +/- e zero)
    Ramo(1, 2, z1=0.02 + 0.08j, z0=0.06 + 0.24j, b1=0.10, b0=0.06),
    Ramo(2, 3, z1=0.02 + 0.08j, z0=0.06 + 0.24j, b1=0.10, b0=0.06),
]
cargas = [Carga(2, P=1.0, Q=0.4)]

# --- SistemaPotencia (fornece n e Sbase; pre-falta flat dispensa o fluxo) ----
barras = [
    Barra(0, 3, 1.0, 0.0, 0.0, 0.0),   # slack (barra 1)
    Barra(1, 1, 1.0, 0.0, -1.0, -0.4),  # carga
    Barra(2, 2, 1.0, 0.0, 0.5, 0.0),    # gerador
]
Y = np.zeros((3, 3), dtype=complex)
for r in ramos:
    i, j = r.de - 1, r.para - 1
    ys = 1.0 / r.z1
    Y[i, i] += ys
    Y[j, j] += ys
    Y[i, j] -= ys
    Y[j, i] -= ys
sistema = SistemaPotencia(barras, Y, Sbase=BASE_MVA)

# --- Estudo de curto-circuito ------------------------------------------------
est = EstudoCurtoCircuito(sistema, geradores, ramos, cargas=cargas,
                          incluir_cargas=True, prefault="flat")

BARRA_FALTA = 2
faltas = {
    "trifasica": est.falta_trifasica(BARRA_FALTA),
    "monofasica": est.falta_monofasica(BARRA_FALTA),
    "bifasica": est.falta_bifasica(BARRA_FALTA),
    "bifasica_terra": est.falta_bifasica_terra(BARRA_FALTA),
}

for res in faltas.values():
    est.imprimir_falta(res)
    print()

# --- Detalhe da falta trifasica: contribuicoes e correntes em kA -------------
res = faltas["trifasica"]
print(f"Falta trifasica na barra {BARRA_FALTA} - detalhe das contribuicoes")
print("-" * 60)
Ia = res["I_fase"]["a"]
print(f"Corrente de falta: {abs(Ia):.4f} pu = "
      f"{est.corrente_kA(Ia, KV_BASE):.4f} kA")
print(f"Potencia de curto: {est.potencia_curto_MVA(BARRA_FALTA):.2f} MVA")
print("Contribuicao dos geradores (fase a):")
for g in res["contrib_gerador"]:
    Ig = g["I_fase"]["a"]
    print(f"  G{g['barra']}: {abs(Ig):.4f} pu = "
          f"{est.corrente_kA(Ig, KV_BASE):.4f} kA")
print("Contribuicao das linhas (fase a):")
for lin in res["contrib_linha"]:
    Il = lin["I_fase"]["a"]
    print(f"  {lin['linha']}: {abs(Il):.4f} pu")
print("Contribuicao das cargas (fase a):")
for c in res["contrib_carga"]:
    Ic = c["I_fase"]["a"]
    print(f"  carga barra {c['barra']}: {abs(Ic):.4f} pu")
print("Contribuicao do carregamento de linha (fase a):")
for s in res["contrib_shunt_linha"]:
    Is = s["I_fase"]["a"]
    print(f"  shunt barra {s['barra']}: {abs(Is):.4f} pu")
print("-" * 60)
