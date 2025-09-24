"""
ml_model.py
===========
ML models for predicting pKa, EC50, and Clearance.

Production-Ready:
- Reads `.dat` files (whitespace or comma-separated, ignores comments).
- Extracts descriptors via RDKit.
- Trains RandomForest regressors.
- Saves trained models as joblib files in `models/`.
- Generates scatter plots with R² values.
"""

import os
import re
import logging
import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
import matplotlib.pyplot as plt

from pkpd_simulator_V1.utils import resource_path, log_step, safe_filename

# --------------------------------------------------
# Feature extraction
# --------------------------------------------------
def smiles_to_features(smiles: str) -> np.ndarray:
    """Convert SMILES to numeric feature vector (6 descriptors)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(6)
    return np.array([
        Descriptors.MolWt(mol),
        Crippen.MolLogP(mol),
        rdMolDescriptors.CalcTPSA(mol),
        rdMolDescriptors.CalcNumHBD(mol),
        rdMolDescriptors.CalcNumHBA(mol),
        rdMolDescriptors.CalcNumRotatableBonds(mol),
    ])


def load_dataset(path: str, target_col: str):
    """Load dataset from .dat file, clean comments, extract SMILES + values."""
    path = resource_path(path)
    data = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[,\s]+", line)
            if len(parts) < 2:
                continue
            smi, val = parts[0], parts[1]
            try:
                val = float(val)
                data.append((smi, val))
            except ValueError:
                continue

    smiles = [d[0] for d in data]
    y = np.array([d[1] for d in data])
    X = np.array([smiles_to_features(s) for s in smiles])
    return X, y, smiles


# --------------------------------------------------
# Training and persistence
# --------------------------------------------------
def train_model(dataset_file: str, target_col: str, model_name: str,
                out_dir="models", plot_dir="pkpd_output/ml_plots"):
    """Train model, save it, return (model, R², plot_path)."""
    os.makedirs(resource_path(out_dir), exist_ok=True)
    os.makedirs(resource_path(plot_dir), exist_ok=True)

    X, y, smiles = load_dataset(dataset_file, target_col)
    if len(y) < 5:
        raise ValueError(f"Not enough samples in {dataset_file} to train.")

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42
    )

    model = RandomForestRegressor(n_estimators=300, random_state=42)
    model.fit(X_train, y_train)

    model_path = os.path.join(resource_path(out_dir), f"{model_name}.joblib")
    joblib.dump(model, model_path)

    # Predictions + R²
    y_pred = model.predict(X_test)
    r2 = r2_score(y_test, y_pred)

    # Scatter plot
    plt.figure(figsize=(6, 6))
    plt.scatter(y_test, y_pred, alpha=0.6, edgecolor="k")
    plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], "r--", lw=2)
    plt.xlabel("True Values")
    plt.ylabel("Predicted Values")
    plt.title(f"{target_col} Prediction (R²={r2:.3f})")
    plot_path = os.path.join(resource_path(plot_dir), f"{model_name}_scatter.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()

    logging.info(f"✔ Trained {model_name} (R²={r2:.3f})")
    return model, r2, plot_path


def _load_model(name: str, out_dir="models"):
    """Load trained model, raise if missing."""
    path = os.path.join(resource_path(out_dir), f"{name}.joblib")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Model {name} not trained yet. Run train_all_models().")
    return joblib.load(path)


def train_all_models():
    """Train all models: pKa, EC50, Clearance."""
    results, r2_scores, ml_plots = {}, {}, {}
    models = {
        "pka": ("data/pka_dataset.dat", "pKa", "pka_model"),
        "ec50": ("data/ec50_dataset.dat", "EC50", "ec50_model"),
        "clearance": ("data/clearance_dataset.dat", "CL", "clearance_model"),
    }

    for key, (path, target, name) in models.items():
        model, r2, plot_path = train_model(path, target, name)
        results[key] = model
        r2_scores[key] = r2
        ml_plots[key] = plot_path

    return results, r2_scores, ml_plots


# --------------------------------------------------
# Predictors
# --------------------------------------------------
def predict_pka(smiles: str) -> float:
    try:
        model = _load_model("pka_model")
    except FileNotFoundError:
        logging.warning("pKa model not found. Training...")
        train_all_models()
        model = _load_model("pka_model")
    return float(model.predict([smiles_to_features(smiles)])[0])


def predict_ec50(smiles: str) -> float:
    try:
        model = _load_model("ec50_model")
    except FileNotFoundError:
        logging.warning("EC50 model not found. Training...")
        train_all_models()
        model = _load_model("ec50_model")
    return float(model.predict([smiles_to_features(smiles)])[0])


def predict_clearance(smiles: str) -> float:
    try:
        model = _load_model("clearance_model")
    except FileNotFoundError:
        logging.warning("Clearance model not found. Training...")
        train_all_models()
        model = _load_model("clearance_model")
    return float(model.predict([smiles_to_features(smiles)])[0])
