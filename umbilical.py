"""
Este arquivo é um exemplo de como fazer a simulação de um umbilical tripolar
composto de três fases, três blindagens e uma armadura.
"""
# %% Imports
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

# Interfacear com julia
from juliacall import Main as jl

jl.include("param_cabo.jl")
jl.include("simulador_circuito.jl")

# O primeiro passo é definir o circuito através de uma netlist.
# Veja "simulador_circuito.jl" para detalhes.

# %% Parâmetros do tempo
dt = 2.5e-9  # passo no tempo [s]
# lembre-se que ~20% dos tempos finais serão inutilizados por causa de erro numérico acumulado
nt = 131072
tmax = (nt - 1) * dt
# `nt` é parâmetro ao invés de `tmax` pois nem sempre `tmax` será múltiplo inteiro de `dt`

# %% Fonte de excitação
# criar arquivo que descreve a fonte de corrente no tempo. Os valores que não
# estão presentes nesse arquivo serão interpolados considerando `v(tmax) = 0`.

# pulso quadrado com largura 1 us e amplitude 50 A
with open("ifonte.csv", "w", encoding="utf8") as file:
    i0 = 50  # amplitude
    ts = 5e-9  # tempo de subida e descida
    t0 = 1e-6  # largura
    file.write("t,v" + "\n")  # colunas precisam ter nome "t" e "v"
    file.write("0,0" + "\n")
    file.write(str(ts) + "," + str(i0) + "\n")
    file.write(str(t0) + "," + str(i0) + "\n")
    file.write(str(t0 + ts) + ",0")
    # Cuidado para não colocar uma linha vazia no final desse arquivo, senão
    # causará um erro por dados estarem faltando!


# %% Umbilical tripolar de 10 km com emendas a cada 2.5 km, netlist base, sem falha

# Nomeia os nós de uma seção do umbilical tripolar: fase_1, blindagem_1, ..., armadura
def nos_umbilical(secao):
    nos = str.split("f1e b1e f2e b2e f3e b3e ae f1r b1r f2r b2r f3r b3r ar")
    for n in range(len(nos)):
        nos[n] = "M" + str(secao) + nos[n]
    return nos


netlist = [
    "dt " + str(dt),
    "nt " + str(nt),
    "I1 ifonte.csv 0 f1",
    "Z1 50,0.1e-6,0 0 f1",  # RLC série: 50 Ohms, 0.1 uH, 0 F; separados apenas por vírgula! Espaços causarão erro.
]
Remenda = 0.01
num_secoes = 4
d = 7  # número de terminais em uma ponta
for i in range(num_secoes):
    nos = nos_umbilical(i)
    if i == 0:
        for n in range(d):
            nos[n] = "0"  # aterrar terminais do emissor
        nos[0] = "f1"  # conectar à fonte
    netlist += ["M" + str(i) + " tripolar2500.h5 " + ' '.join(nos)]
    if i > 0:  # emendas conectando o emissor da seção atual ao receptor da seção anterior
        nos0 = nos_umbilical(i - 1)
        nos1 = nos_umbilical(i)
        for k in range(d):
            netlist += ["Remenda " + str(Remenda) + " " + nos0[k + d] + " " + nos1[k]]


netlist_file = "netlist_umbilical.txt"
with open(netlist_file, "w", encoding="utf8") as file:
    for e in netlist:
        file.write(e + "\n")

# inserir uma falha é fácil: basta alterar o R da emenda defeituosa
# ou colocar um resistor (R > 0) conectando dois terminais quaisquer
# ou, ainda, uma fonte de tensão nula.

# %% Calcular YN dos segmentos
## Geometria do umbilical tripolar
rc = 1e-3 * np.array([9.6, 17.054, 18.054, 19.50])  # raios dos cabos SC [m]: condutor central, isolante, blindagem e isolante externo
ra = 1e-3 * np.array([48, 59, 65])  # raios da armadura [m]: interno, externo, capa
nx = 2  # número de condutores "ativos" em cada cabo SC (single core)

# distância do centro do umbilical ao centro das veias
d01 = rc[3] / np.cos(np.deg2rad(30))
d02 = d01
di = [d01, d02, d02]
theta_jk = (2 * np.pi) / 3  # ângulo entre cabos do umbilical

# permissividades elétricas
epsr_c = 3.31  # condutor central
epsr_b = 2.3  # blindagem, jaqueta da veia
epsr_a = 10  # última camada, que envolve a armadura

# resistividades elétricas [Ohm.m]
rho_c = 1.7241e-8  # condutor central
rho_b = 2.2e-7  # blindagem, jaqueta da veia
rho_a = 2.86e-8  # armadura

mur_a = 300  # permeabilidade magnética relativa da armadura
sig_s = 5.0  # condutividade elétrica do mar [S/m]
eps_s = 81.0  # permissividade elétrica relativa do mar

L = 2500  # comprimento do segmento

# obter frequência complexa `sk = c + jω` necessárias para o cálculo
_, sk = jl.laplace(jl.ones(nt), tmax, nt)
nf = np.size(sk)  # número de frequências

# %% Calcular parâmetros por unidade de comprimento do umbilical
zc = np.empty((7, 7, nf), dtype=np.complex128, order='F')
yc = np.empty((7, 7, nf), dtype=np.complex128, order='F')
for f in range(nf):
    w = -1.0j * sk[f]
    zc[:,:,f], yc[:,:,f] = jl.zy_cabo(w, rc, ra, rho_a, rho_b, rho_c, epsr_a, epsr_b, epsr_c, mur_a, theta_jk, di, sig_s, eps_s, nx)


# admitância nodal yn equivalente para o comprimento
yn = np.empty((14, 14, nf), dtype=np.complex128, order='F')
for f in range(nf):
    yn[:,:,f] = jl.ynodal(zc[:,:,f], yc[:,:,f], L)


# salvar yn como arquivo HDF5
arquivo_yn = "tripolar" + str(L) + ".h5"
with h5py.File(arquivo_yn, "w") as f:
    f.create_dataset("s", data=sk)
    f.create_dataset("v", data=np.transpose(yn))  # Python está salvando o array transposto...


# calcular zc e yc é demorado. Pode ser uma boa ideia salvá-los num arquivo
# para depois reutilizar para calcular o yn(L) equivalente...

# %% Simular a netlist
vout = jl.simular_netlist(netlist_file)
nt_trunc = int(np.ceil(nt * 0.8))  # tempos finais serão inutilizados por causa de erro numérico acumulado
# converter o Dict julia para python e então criar um DataFrame
df = pd.DataFrame(jl.pydict(vout))[0:nt_trunc]
tempo = np.array([dt * k for k in range(nt_trunc)]) * 1e6

# %% Plotar resultados
terminais = str.split("M0f1r M0b1r M0f2r M0b2r M0f3r M0b3r M0ar")

plt.figure()
plt.plot(tempo, df[terminais])
plt.xlabel("Tempo [us]")
plt.ylabel("Potencial [V]")
plt.grid()
plt.show()
