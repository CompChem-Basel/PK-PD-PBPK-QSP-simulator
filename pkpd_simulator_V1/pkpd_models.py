"""
pkpd_models.py
==============
Defines and solves PK, PBPK, PD, and QSP models.

Models:
- PK: One-Compartment & Two-Compartment
- PBPK: Simplified 3-compartment (Blood, Liver, Kidney)
- PD: Emax model
- QSP: Toy receptor-ligand activation cascade

All solved with SciPy ODE solvers.
"""

import numpy as np
from scipy.integrate import solve_ivp
from pkpd_simulator_V1.utils import log_step


# =====================================================================
# PK MODELS
# =====================================================================
def pk_one_compartment(t, y, ke):
    """
    One-compartment PK model (IV bolus or lumped oral).
    Args:
        t (float): Time (h)
        y (list): [A] drug amount in central compartment
        ke (float): Elimination rate constant (1/h)
    Returns:
        list: dA/dt
    """
    A = y[0]  # drug amount in central compartment
    dA_dt = -ke * A
    return [dA_dt]


def pk_two_compartment(t, y, k12, k21, ke):
    """
    Two-compartment PK model (central + peripheral).
    Args:
        t (float): Time (h)
        y (list): [Ac, Ap] drug amounts
        k12 (float): rate constant central -> peripheral
        k21 (float): rate constant peripheral -> central
        ke (float): elimination rate constant (1/h)
    Returns:
        list: [dAc/dt, dAp/dt]
    """
    Ac, Ap = y
    dAc_dt = -(ke + k12) * Ac + k21 * Ap
    dAp_dt = k12 * Ac - k21 * Ap
    return [dAc_dt, dAp_dt]


# =====================================================================
# PBPK MODEL
# =====================================================================
def pbpk_model(t, y, Ql, Qk, Cl_hepatic, Cl_renal):
    """
    Simplified 3-compartment PBPK model.
    Compartments:
        y[0] = Blood concentration
        y[1] = Liver concentration
        y[2] = Kidney concentration
    Args:
        Ql (float): Liver blood flow (L/h)
        Qk (float): Kidney blood flow (L/h)
        Cl_hepatic (float): Hepatic clearance rate (1/h)
        Cl_renal (float): Renal clearance rate (1/h)
    Returns:
        list: [dCb/dt, dCl/dt, dCk/dt]
    """
    Cb, Cl, Ck = y

    dCb_dt = -(Ql * (Cb - Cl) + Qk * (Cb - Ck))
    dCl_dt = Ql * (Cb - Cl) - Cl_hepatic * Cl
    dCk_dt = Qk * (Cb - Ck) - Cl_renal * Ck

    return [dCb_dt, dCl_dt, dCk_dt]


# =====================================================================
# PD MODEL
# =====================================================================
def pd_emax(conc, Emax=1.0, EC50=1.0):
    """
    Emax pharmacodynamic (PD) model.
    Args:
        conc (ndarray): drug concentration (mg/L)
        Emax (float): maximum effect
        EC50 (float): concentration at 50% effect
    Returns:
        ndarray: PD effect over time
    """
    return (Emax * conc) / (EC50 + conc + 1e-12)


# =====================================================================
# QSP MODEL
# =====================================================================
def qsp_model(t, y, kon, koff, kact):
    """
    Simple receptor-ligand activation cascade.
    Compartments:
        y[0] = free receptor
        y[1] = activated receptor
    Args:
        kon (float): binding rate constant
        koff (float): unbinding rate constant
        kact (float): activation rate constant
    Returns:
        list: [dR/dt, dA/dt]
    """
    R, A = y
    dR_dt = -kon * R + koff * A
    dA_dt = kon * R - (koff + kact) * A
    return [dR_dt, dA_dt]


# =====================================================================
# MASTER FUNCTION
# =====================================================================
def run_all_simulations(
    smiles: str,
    dose_mg: float,
    route: str,
    population: str,
    measured_pka: float = None,
    measured_type: str = None,
    extra_params: dict = None,
    duration: int = 24,
) -> dict:
    """
    Run PK, PBPK, PD, QSP simulations for a given drug.

    Args:
        smiles (str): Drug SMILES
        dose_mg (float): Dose (mg)
        route (str): Administration route ("oral" or "iv")
        population (str): Population type
        measured_pka (float): Optional experimental pKa
        measured_type (str): Optional pKa type
        extra_params (dict): Extra parameters (e.g., ML predictions)
        duration (int): Simulation duration (hours)

    Returns:
        dict: Simulation results with ODE solutions and PD effect
    """

    log_step("Simulation", f"Starting for {smiles}")

    # -----------------------------
    # Time settings
    # -----------------------------
    t_end = duration  # hours
    t_eval = np.linspace(0, t_end, 200)

    # -----------------------------
    # Parameters
    # -----------------------------
    if extra_params is None:
        extra_params = {}

    # PK parameters
    Vd = float(extra_params.get("Vd_L_per_kg", 0.5))   # L/kg
    Cl = float(extra_params.get("CL_pred", 0.1))       # L/h
    ka = float(extra_params.get("ka", 1.0))            # absorption rate (oral)
    ke = Cl / max(Vd, 1e-6)                            # elimination rate constant

    k12 = float(extra_params.get("k12", 0.3))
    k21 = float(extra_params.get("k21", 0.2))

    # PBPK flows
    Ql = float(extra_params.get("Ql", 1.2))            # liver flow
    Qk = float(extra_params.get("Qk", 0.8))            # kidney flow
    Cl_hepatic = Cl * 0.6
    Cl_renal = Cl * 0.4

    # PD parameters
    EC50 = float(extra_params.get("EC50_pred", 1.0))

    # QSP parameters
    kon = float(extra_params.get("kon", 0.5))
    koff = float(extra_params.get("koff", 0.2))
    kact = float(extra_params.get("kact", 0.1))

    # -----------------------------
    # Initial conditions
    # -----------------------------
    dose = dose_mg
    if route.lower() == "iv":
        A0_1C = dose
        A0_2C = [dose, 0.0]
    else:  # oral (lumped into central compartment)
        A0_1C = dose
        A0_2C = [dose, 0.0]

    y0_pbpk = [dose / max(Vd, 1e-6), 0.0, 0.0]
    y0_qsp = [1.0, 0.0]

    # -----------------------------
    # Solve ODE systems
    # -----------------------------
    one_c = solve_ivp(
        pk_one_compartment,
        [0, t_end],
        [A0_1C],
        t_eval=t_eval,
        args=(ke,),
    )

    two_c = solve_ivp(
        pk_two_compartment,
        [0, t_end],
        A0_2C,
        t_eval=t_eval,
        args=(k12, k21, ke),
    )

    pbpk = solve_ivp(
        pbpk_model,
        [0, t_end],
        y0_pbpk,
        t_eval=t_eval,
        args=(Ql, Qk, Cl_hepatic, Cl_renal),
    )

    pd_effect = pd_emax(pbpk.y[0], Emax=1.0, EC50=EC50)

    qsp = solve_ivp(
        qsp_model,
        [0, t_end],
        y0_qsp,
        t_eval=t_eval,
        args=(kon, koff, kact),
    )

    # -----------------------------
    # Package results
    # -----------------------------
    results = {
        "t": t_eval,
        "one_compartment": one_c,
        "two_compartment": two_c,
        "pbpk": pbpk,
        "pd_effect": pd_effect,
        "qsp": qsp,
        "duration": duration,
    }

    log_step("Simulation", f"Finished for {smiles}")
    return results
