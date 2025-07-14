# /2424!/usr/bin/env python
# coding: utf-8

# # Generate trace sampling input from RadicalPy simulator

# In[1]:


import radicalpy as rp
from radicalpy.simulation import State
import numpy as np
import matplotlib.pyplot as plt

import sys
import pathlib

sys.path.append(str(pathlib.Path().resolve()))

from utils import dump_input, parse_output  # ty: ignore


# In[2]:


is_small_case = False

if is_small_case:
    # You can use following block instead
    n_nuc_spins = 1

    flavin = rp.simulation.Molecule.fromisotopes(
        isotopes=["1H", "14N"], hfcs=[0.4, 0.5]
    )
    Z = rp.simulation.Molecule.fromisotopes(
        isotopes=["1H"] * n_nuc_spins,
        hfcs=[0.5] * n_nuc_spins,
        # isotopes=[], hfcs=[]
    )
    sim = rp.simulation.LiouvilleSimulation([Z, flavin])

    # Parameters
    A = {}  # mT
    isotropic = True

    # Isotropic
    for i in range(len(sim.radicals)):
        for j, nuc in enumerate(sim.molecules[i].nuclei):
            if isotropic:
                A[(i, j)] = np.eye(3) * nuc.hfc.isotropic
            else:
                A[(i, j)] = nuc.hfc.anisotropic
    B0 = 0.2  # 2J
    B = np.array((0.0, 0.0, 1.0)) * B0  # mT
    J = 0.1  # Typically 1.0e+03 scale # mT
    D = 0.1  # mT
    kS = 1.0e06  # Exponential model in s-1
    kT = 1.0e06
    if isinstance(D, float):
        D = 2 / 3 * np.diag((-1.0, -1.0, 2.0)) * D * sim.radicals[0].gamma_mT
else:
    flavin = rp.simulation.Molecule.all_nuclei("flavin_anion")
    trp = rp.simulation.Molecule.all_nuclei("tryptophan_cation")
    for mol in [flavin, trp]:
        nucs = []
        for nuc in mol.nuclei:
            # if nuc.name == "1H":
            #     nucs.append(nuc)
            if abs(nuc.hfc.isotropic) > 1e-01:
                nucs.append(nuc)
        mol.nuclei = nucs
    sim = rp.simulation.LiouvilleSimulation([flavin, trp])
    A = {}
    isotropic = True
    for i in range(len(sim.radicals)):
        for j, nuc in enumerate(sim.molecules[i].nuclei):
            if isotropic:
                A[(i, j)] = np.eye(3) * nuc.hfc.isotropic
            else:
                A[(i, j)] = nuc.hfc.anisotropic
    B0 = 0.050  # Fay 2020
    B = np.array((0.0, 0.0, 1.0)) * B0
    J = 0.224  # Fay 2020

    D = (
        np.array(
            [
                [-0.225, 0.156, 0.198],
                [0.156, 0.117, -0.082],
                [0.198, -0.082, 0.107],
            ]
        )
        * sim.radicals[0].gamma_mT
    )
    D = -0.38  # Fay 2020
    kS = 1.0e06  # Exponential model in s-1
    kT = 1.0e06
    if isinstance(D, float):
        D = 2 / 3 * np.diag((-1.0, -1.0, 2.0)) * D * sim.radicals[0].gamma_mT
sim


# In[3]:


if has_rp_result := (is_small_case and len(sim.particles) < 7):
    # when D is Float, a bug appears
    assert isinstance(D, np.ndarray)
    H = sim.total_hamiltonian(B0=B0, D=D, J=J)
    time = np.arange(0, 5.0e-8 + 5e-10, 5e-10)
    kinetics = [
        rp.kinetics.Haberkorn(kS, rp.simulation.State.SINGLET),
        rp.kinetics.Haberkorn(kT, rp.simulation.State.TRIPLET),
    ]
    sim.apply_liouville_hamiltonian_modifiers(H, kinetics)
    rhos = sim.time_evolution(State.SINGLET, time, H)

    time_evol_s = sim.product_probability(State.SINGLET, rhos)
    time_evol_tp = sim.product_probability(State.TRIPLET_PLUS, rhos)
    time_evol_tz = sim.product_probability(State.TRIPLET_ZERO, rhos)
    time_evol_tm = sim.product_probability(State.TRIPLET_MINUS, rhos)
    x = time * 1e6

    plt.plot(x, time_evol_tp, linewidth=2, label="T+")
    plt.plot(x, time_evol_tz, linewidth=2, label="T0")
    plt.plot(x, time_evol_s, linewidth=2, label="S")
    plt.plot(x, time_evol_tm, linewidth=2, label="T-")
    plt.plot(
        x,
        time_evol_tp + time_evol_tz + time_evol_s + time_evol_tm,
        linewidth=2,
        label="trace",
    )
    plt.legend()
    plt.title(f"RadicalPy: Density Matrix Approach {B0=}, {J=}")
    plt.xlabel(r"Time ($\mu s$)")
    plt.ylabel("Probability")
    # plt.ylim([0, 1])
    plt.grid()
    plt.show()


# In[4]:


input_path, output_path = dump_input(
    sim=sim,
    output_folder="out",
    simulation_type="trace_sampling",  # "exact_dynamics", "trace_sampling"
    J=J,
    D=D,
    kS=kS,
    kT=kT,
    B=B0,
    dt=5e-10,
    simulation_time=5e-08 - 5e-10,
    N_krylov=7,
    integrator_tolerance=1e-08,
    N_samples=128,
    M1=1,
    M2=2,
)
print(input_path)
# cat contents of input file
with open(input_path) as f:
    print(f.read())


# In[5]:


# Execute sphinchem
import subprocess

subprocess.run(["../bin/spinchem", input_path])


# In[6]:


df = parse_output(output_path="out", subdir=".000", dt=5e-10)


# In[7]:


df


# In[8]:


plt.plot(df["time"], df["Tp_prob"], label="Tp")
plt.plot(df["time"], df["T0_prob"], label="T0")
plt.plot(df["time"], df["S_prob"], label="S")
plt.plot(df["time"], df["Tm_prob"], label="Tm")
plt.plot(
    df["time"],
    df["S_prob"] + df["T0_prob"] + df["Tp_prob"] + df["Tm_prob"],
    label="trace",
)
if has_rp_result:
    plt.plot(x, time_evol_tp, linewidth=2, label="T+ (RadicalPy)", ls="--")
    plt.plot(x, time_evol_tz, linewidth=2, label="T0 (RadicalPy)", ls="--")
    plt.plot(x, time_evol_s, linewidth=2, label="S (RadicalPy)", ls="--")
    plt.plot(x, time_evol_tm, linewidth=2, label="T- (RadicalPy)", ls="--")
plt.xlabel("Time (microsec)")
plt.ylabel("Probability")
plt.grid()
plt.legend()
plt.show()


# In[ ]:
