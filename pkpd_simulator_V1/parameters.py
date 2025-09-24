"""
parameters.py
=============
Stores and retrieves physiological and population-specific
parameters for PK/PD/PBPK/QSP simulations.
"""

import logging
from pkpd_simulator_V1.utils import log_step


def get_baseline_physiology() -> dict:
    """Return baseline physiological parameters (adult healthy male)."""
    return {
        "Vd": 0.5,    # Volume of distribution (L/kg)
        "Cl": 0.1,    # Clearance (L/hr)
        "ka": 1.0,    # Absorption rate constant (1/hr)
        "t_abs": 2.0  # Absorption duration (hr, oral)
    }


def choose_population_parameters(population: str,
                                 route: str,
                                 baseline: dict = None) -> dict:
    """
    Adjust PK parameters for:
    - Population: healthy_male/female, diabetic, liver impairment
    - Route: oral vs IV
    """
    if baseline is None:
        baseline = get_baseline_physiology()

    params = baseline.copy()
    population = population.lower().strip()
    route = route.lower().strip()

    log_step("Population", f"{population}, Route={route}")

    # Population adjustments
    if population == "healthy_male":
        pass
    elif population == "healthy_female":
        params["Vd"] *= 0.9
        params["Cl"] *= 0.9
    elif population == "type1_diabetic":
        params["Vd"] *= 1.1
        params["Cl"] *= 1.1
    elif population == "liver_impairment_moderate":
        params["Vd"] *= 1.2
        params["Cl"] *= 0.7
    elif population == "custom":
        logging.info("Custom parameters requested (manual input).")
        try:
            params["Vd"] = float(input(f"Vd (L/kg) [{params['Vd']}]: ") or params["Vd"])
            params["Cl"] = float(input(f"Cl (L/hr) [{params['Cl']}]: ") or params["Cl"])
            params["ka"] = float(input(f"ka (1/hr) [{params['ka']}]: ") or params["ka"])
            params["t_abs"] = float(input(f"t_abs (hr) [{params['t_abs']}]: ") or params["t_abs"])
        except Exception as e:
            logging.error(f"Invalid custom input, using baseline: {e}")

    # Route adjustments
    if route == "oral":
        params["ka"] = params.get("ka", 1.0)
    elif route == "iv":
        params["ka"] = 0.0
    else:
        raise ValueError(f"Unsupported route: {route}")

    return params
