"""
train_pka_model.py
==================
Train a regression model to predict pKa values from SMILES.

Outputs (in user-defined folder inside pkpd_output/):
- Trained model (.joblib)
- Scatter plot (true vs predicted pKa, with R² score)
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import joblib

from descriptors import featurize_smiles


def train_pka_model(db_csv: str = "pka_dataset.csv",
                    smiles_col: str = "smiles",
                    target_col: str = "pKa") -> tuple[str, float]:
    # --- Folder setup ---
    base_dir = "pkpd_output"
    os.makedirs(base_dir, exist_ok=True)
    folder_name = input("Enter folder name for this pKa training run: ").strip() or "pka_run"
    out_dir = os.path.join(base_dir, folder_name)
    os.makedirs(out_dir, exist_ok=True)

    model_out = os.path.join(out_dir, "pka_model.joblib")
    plot_out = os.path.join(out_dir, "pka_scatter.png")

    # Load dataset
    df = pd.read_csv(db_csv)

    # Featurize SMILES
    X, y = [], []
    for _, row in df.iterrows():
        feats = featurize_smiles(row[smiles_col])
        X.append([feats["MWT"], feats["LogP"], feats["TPSA"], feats["HBD"], feats["HBA"]])
        y.append(row[target_col])

    X, y = np.array(X), np.array(y)

    # Train/test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Train model
    model = RandomForestRegressor(n_estimators=200, random_state=42)
    model.fit(X_train, y_train)

    # Predictions
    y_pred = model.predict(X_test)

    # R² score
    r2 = r2_score(y_test, y_pred)
    print(f"✅ pKa Model R² Score: {r2:.3f}")

    # Scatter plot
    plt.figure(figsize=(6, 6))
    plt.scatter(y_test, y_pred, alpha=0.7, edgecolors="k")
    plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], "r--", lw=2)
    plt.xlabel("True pKa")
    plt.ylabel("Predicted pKa")
    plt.title(f"pKa Prediction (R² = {r2:.3f})")
    plt.savefig(plot_out, dpi=300)
    plt.close()
    print(f"📊 Scatter plot saved: {plot_out}")

    # Save model
    joblib.dump(model, model_out)
    print(f"💾 Model saved: {model_out}")
    return model_out, r2


if __name__ == "__main__":
    train_pka_model()
