"""
run_simulation.py
=================
Main entry point for PK/PD/PBPK/QSP simulations.
"""

import os
import sys
import logging
from rdkit import Chem
from pkpd_simulator_V1.descriptors import featurize_smiles, determine_absorption_site, _guess_pka_type
from typing import Optional
from pkpd_simulator_V1.descriptors import featurize_smiles, determine_absorption_site
from pkpd_simulator_V1.parameters import choose_population_parameters, get_baseline_physiology
from pkpd_simulator_V1.pkpd_models import run_all_simulations
from pkpd_simulator_V1.ml_model import predict_pka, predict_ec50, predict_clearance, train_all_models
from pkpd_simulator_V1.utils import (
    setup_logging,
    save_all_results,
    log_step,
    smiles_to_2d,
    safe_filename,
    make_comparison_overlay,
    make_comparison_report,
    make_detailed_report,
)


# ==============================================================
# Helpers
# ==============================================================
def ask_user_inputs():
    """Prompt user for SMILES, population, and dose settings (CLI)."""
    print("\nChoose input mode:")
    print("1. Single SMILES")
    print("2. Two SMILES (comparison)")
    print("3. Three SMILES (comparison)")
    print("4. List of SMILES (batch, no comparison)")
    mode = input("Enter choice [1-4]: ").strip()

    smiles_list = []
    if mode == "1":
        smiles_list = [input("Enter SMILES: ").strip()]
    elif mode == "2":
        smiles_list = [
            input("Enter SMILES #1: ").strip(),
            input("Enter SMILES #2: ").strip(),
        ]
    elif mode == "3":
        smiles_list = [
            input("Enter SMILES #1: ").strip(),
            input("Enter SMILES #2: ").strip(),
            input("Enter SMILES #3: ").strip(),
        ]
    elif mode == "4":
        n = int(input("How many SMILES in list? "))
        for i in range(n):
            smiles_list.append(input(f"Enter SMILES #{i+1}: ").strip())
    else:
        print("Invalid choice. Defaulting to single SMILES.")
        smiles_list = [input("Enter SMILES: ").strip()]

    print("\nChoose patient population:")
    print("1. Healthy Male")
    print("2. Healthy Female")
    print("3. Type 1 Diabetic")
    print("4. Liver Impairment (Moderate)")
    print("5. Custom parameters")
    pop_choice = input("Enter choice [1-5]: ").strip()
    pop_map = {
        "1": "healthy_male",
        "2": "healthy_female",
        "3": "type1_diabetic",
        "4": "liver_impairment_moderate",
        "5": "custom",
    }
    population = pop_map.get(pop_choice, "healthy_male")

    route = input("Route (oral/iv) [oral]: ").strip().lower() or "oral"
    dose_mg = float(input("Dose in mg [100]: ").strip() or "100")
    duration = int(input("Duration (hours) [24]: ").strip() or "24")

    return smiles_list, population, route, dose_mg, duration


def resource_path(relative_path: str) -> str:
    """Resolve path for bundled resources (PyInstaller support)."""
    if hasattr(sys, "_MEIPASS"):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)


def get_absorption_image(site: str) -> Optional[str]:
    """
    Map absorption site → static image path in /images.
    Returns None for "Poor absorption" or "Unknown".
    """
    mapping = {
        "Stomach": "images/Stomach.png",
        "Duodenum": "images/Duodenum.png",
        "Ileum": "images/ileum.png",
    }
    return resource_path(mapping.get(site)) if site in mapping else None


# ==============================================================
# Simulation Driver
# ==============================================================
def simulate_drug(smiles_list, population, route, dose_mg, duration=24):
    """
    Main simulation driver (called by CLI or GUI).
    """
    output_dir = "pkpd_output"
    os.makedirs(output_dir, exist_ok=True)

    log_file = setup_logging("sim_run", output_dir)
    log_step("Logging initialized", f"Log file: {log_file}")
    logging.info(
        f"Population: {population}, Route: {route}, "
        f"Dose: {dose_mg} mg, Duration: {duration} h"
    )

    try:
        _ = predict_pka("CCO")
    except Exception:
        logging.warning("ML models not trained. Training now...")
        train_all_models()

    results = {}

    for smiles in smiles_list:
        safe_smiles = safe_filename(smiles)
        drug_dir = os.path.join(output_dir, safe_smiles)
        os.makedirs(drug_dir, exist_ok=True)

        log_step("Featurization", f"Processing {smiles}")
        features = featurize_smiles(smiles)

        # ML predictions
        raw_pka = predict_pka(smiles)
        try:
            pka_val = float(raw_pka) if raw_pka is not None else None
        except Exception:
            pka_val = None 
			
        features["pKa_pred"] = pka_val
        features["EC50_pred"] = predict_ec50(smiles)
        features["CL_pred"] = predict_clearance(smiles)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            features["pKa_type"] = _guess_pka_type(mol)

        # Absorption site (unpack tuple)
        site, site_img = determine_absorption_site(
            pka=pka_val, heuristic_type=features.get("pKa_type", "neutral"), images_dir="images"
        )
        features["absorption_site"] = site
        features["absorption_image"] = site_img
        logging.info(f"Absorption site: {site}, pKa: {pka_val}")

        # Parameters
        pop_params = choose_population_parameters(population, route)
        physiology = get_baseline_physiology()
        logging.info(f"PK parameters: {pop_params}")
        logging.info(f"Physiology: {physiology}")

        # 2D structure
        smiles_img = smiles_to_2d(smiles, drug_dir, safe_smiles)
        features["structure_img"] = smiles_img.get("png") if isinstance(smiles_img, dict) else None

        # Simulation
        sim_results = run_all_simulations(
            smiles=smiles,
            dose_mg=dose_mg,
            route=route,
            population=population,
            measured_pka=pka_val,
            measured_type="predicted",
            extra_params=features,
            duration=duration,
        )

        # Save outputs
        files = save_all_results(sim_results, smiles, drug_dir, features, smiles_img)

        results[smiles] = {
            "features": features,
            "files": files,
            "sim": sim_results,
            "absorption_site": site,
            "absorption_image": site_img,
        }

    return results


# ==============================================================
# CLI Entry
# ==============================================================
def main():
    smiles_list, population, route, dose_mg, duration = ask_user_inputs()
    results = simulate_drug(smiles_list, population, route, dose_mg, duration)

    if len(smiles_list) in [2, 3]:
        overlay_file = make_comparison_overlay(results, "pkpd_output")
        pdf_file = make_comparison_report(results, "pkpd_output", overlay_file)
        print(f"\n📊 Comparison overlay: {overlay_file}")
        print(f"📑 Comparison report: {pdf_file}")

    detailed_pdf = make_detailed_report(results, "pkpd_output")
    print(f"\n📑 Detailed report: {detailed_pdf}")

    print("\n✅ Simulation complete. Results stored in pkpd_output/")


if __name__ == "__main__":
    main()
