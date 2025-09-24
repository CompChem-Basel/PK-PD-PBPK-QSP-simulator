"""
descriptors.py
==============
Transforms SMILES into molecular descriptors and ADME/T predictions.

⚡ Features:
- Full featurization with MWT, LogP, TPSA, HBD, HBA, RotB.
- pKa prediction (heuristic + ML hook).
- Absorption site prediction based on pKa rules:
    * pKa 1–3 → Stomach
    * pKa 5–6 → Duodenum
    * pKa 6–8 → Ileum
    * Else → Poor absorption (unionized/ionized, low permeability)
- Each valid absorption site has a linked organ image from `images/` folder.
- Returns dict of descriptors + PK parameters for simulations.
"""

from typing import Dict, Optional, Tuple
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, AllChem
from pkpd_simulator_V1.utils import resource_path


# ============================================================
# Henderson–Hasselbalch helpers
# ============================================================
def _fraction_unionized_acid(pH: float, pKa: float) -> float:
    return 1.0 / (1.0 + 10.0 ** (pH - pKa))


def _fraction_unionized_base(pH: float, pKa: float) -> float:
    return 1.0 / (1.0 + 10.0 ** (pKa - pH))


def _guess_pka_type(mol: Chem.Mol, predicted_pka: float = None) -> Optional[str]:
    """
    Classify molecule as acidic, basic, or neutral.
    Uses SMARTS for known functional groups, but corrects misclassifications
    using predicted pKa value if available.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule object.
    predicted_pka : float, optional
        Predicted pKa value (from ML model). Used to override obvious mismatches.

    Returns
    -------
    str
        "acidic", "basic", or "neutral".
    """

    # --------------------------
    # SMARTS patterns
    # --------------------------
    smarts_acidic = [
        Chem.MolFromSmarts("C(=O)[OH]"),      # carboxylic acids
        Chem.MolFromSmarts("S(=O)(=O)[OH]"),  # sulfonic acids
        Chem.MolFromSmarts("c1ccc(O)cc1"),    # phenols (weak acids)
        Chem.MolFromSmarts("C(=O)[O;H1]"),   # protonated COOH
        Chem.MolFromSmarts("C(=O)O"),        # generic carboxyl
    ]
    smarts_basic = [
        Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]"),  # primary/secondary amines
        Chem.MolFromSmarts("[NX3;H0;!$(NC=O)]"),     # tertiary amines
        Chem.MolFromSmarts("n")                      # heteroaromatic N
    ]

    # --------------------------
    # Initial assignment
    # --------------------------
    for patt in smarts_acidic:
        if mol.HasSubstructMatch(patt):
            return "acidic"

    for patt in smarts_basic:
        if mol.HasSubstructMatch(patt):
            return "basic"
                
    return "neutral"

    # --------------------------
    # Correction using predicted pKa
    # --------------------------
    if predicted_pka is not None:
        try:
            pka_val = float(predicted_pka)

            # If SMARTS says "basic" but pKa is clearly acidic (<6) → override
            if pka_type == "basic" and pka_val < 6:
                pka_type = "acidic"

            # If SMARTS says "acidic" but pKa is clearly basic (>8) → override
            elif pka_type == "acidic" and pka_val > 8:
                pka_type = "basic"

            # If SMARTS says "neutral" but pKa is in a strong acidic or basic range
            elif pka_type == "neutral":
                if pka_val < 5.5:
                    pka_type = "acidic"
                elif pka_val > 7.5:
                    pka_type = "basic"

        except Exception:
            pass  # keep original SMARTS type

    return pka_type

# ============================================================
# Optional ML pKa model hook
# ============================================================
_pka_model = None


def load_pka_model(path: str) -> None:
    """Load ML model for pKa prediction."""
    global _pka_model
    import joblib
    _pka_model = joblib.load(path)


def _predict_pka_ml(smiles: str, fallback: float, pka_type: Optional[str]) -> Tuple[float, Optional[str]]:
    """Predict pKa with ML model if loaded, otherwise fallback heuristic."""
    if _pka_model is None:
        return fallback, pka_type
    try:
        mol = Chem.MolFromSmiles(smiles)
        feats = [
            Descriptors.MolWt(mol),
            Crippen.MolLogP(mol),
            rdMolDescriptors.CalcTPSA(mol),
            rdMolDescriptors.CalcNumHBD(mol),
            rdMolDescriptors.CalcNumHBA(mol),
            rdMolDescriptors.CalcNumRotatableBonds(mol),
        ]
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=512)
        arr = np.zeros((512,), dtype=int)
        from rdkit import DataStructs
        DataStructs.ConvertToNumpyArray(fp, arr)
        x = np.concatenate([feats, arr]).reshape(1, -1)
        pred = float(_pka_model.predict(x)[0])
        return pred, pka_type
    except Exception:
        return fallback, pka_type


# ============================================================
# Main featurization
# ============================================================
def featurize_smiles(smiles: str) -> Dict:
    """Generate descriptors, predicted PK parameters, and pKa values."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # Base descriptors
    MWT = Descriptors.MolWt(mol)
    LogP = Crippen.MolLogP(mol)
    TPSA = rdMolDescriptors.CalcTPSA(mol)
    HBD = rdMolDescriptors.CalcNumHBD(mol)
    HBA = rdMolDescriptors.CalcNumHBA(mol)
    RotB = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # pKa prediction
    pka_type = _guess_pka_type(mol, predicted_pka=None)  # at first
    heuristic_pka = 7.4 if LogP > 3.5 else (3.5 if TPSA < 20 else 5.0)
    pKa, pka_type = _predict_pka_ml(smiles, heuristic_pka, pka_type)

    # Unionized fractions at GI pH
    ph_st, ph_si, ph_co = 1.5, 6.5, 7.5
    if pka_type == "acidic":
        fu_st = _fraction_unionized_acid(ph_st, pKa)
        fu_si = _fraction_unionized_acid(ph_si, pKa)
        fu_co = _fraction_unionized_acid(ph_co, pKa)
    elif pka_type == "basic":
        fu_st = _fraction_unionized_base(ph_st, pKa)
        fu_si = _fraction_unionized_base(ph_si, pKa)
        fu_co = _fraction_unionized_base(ph_co, pKa)
    else:
        fu_st = fu_si = fu_co = 0.8

    # PK mappings
    fu_plasma = float(max(0.01, min(0.95, 0.6 - 0.1 * max(LogP, 0) - 0.03 * HBD)))
    Vd = float(max(0.1, min(6.0, 0.7 + 0.25 * LogP - 0.002 * TPSA - 0.03 * HBD)))
    CLint = float(max(0.01, 0.3 + 0.02 * (HBD + HBA) - 0.03 * max(LogP, 0)))
    Qh = 1.35
    CLh = (Qh * fu_plasma * CLint) / (Qh + fu_plasma * CLint + 1e-12)
    extraction = CLh / (Qh + 1e-12)
    F_h = 1.0 - extraction
    k_elim = CLh / max(1e-6, Vd)

    if TPSA < 60 and 0 <= LogP <= 5:
        ka = 1.5 * fu_si
    elif TPSA < 90:
        ka = 0.8 * fu_si
    else:
        ka = 0.3 * fu_si

    CLl_pbpk = CLh * (LogP / 3.0 if LogP > 0 else 0.5)
    CLk_pbpk = 0.2 * fu_plasma * (1.0 + 0.1 * HBA)

    return {
        "MWT": MWT, "LogP": LogP, "TPSA": TPSA,
        "HBD": HBD, "HBA": HBA, "RotatableBonds": RotB,
        "pKa": pKa, "pKa_type": pka_type,
        "unionized_stomach": fu_st,
        "unionized_intestine": fu_si,
        "unionized_colon": fu_co,
        "fu_plasma": fu_plasma,
        "Vd_L_per_kg": Vd,
        "CLint_L_h_kg": CLint,
        "CLh_L_h_kg": CLh,
        "hepatic_extraction": extraction,
        "F_hepatic_only": F_h,
        "k_elim": k_elim,
        "ka": ka,
        "CLl_pbpk": CLl_pbpk,
        "CLk_pbpk": CLk_pbpk,
    }


# ============================================================
# Absorption Site Prediction
# ============================================================
def determine_absorption_site(
    pka: Optional[float] = None,
    heuristic_type: str = None,
    images_dir: str = "images"
) -> Tuple[str, Optional[str]]:
    """
    Predict absorption site based on pKa rules + acid/base type.

    Rules (approximate physiology):
    - Acidic drugs:
        pKa < 4      → Stomach
        4 ≤ pKa ≤ 6  → Duodenum
        > 6          → Poor absorption
    - Basic drugs:
        pKa < 6      → Poor absorption
        6 ≤ pKa ≤ 9  → Ileum
        > 9          → Duodenum
    - Neutral / unknown:
        5 ≤ pKa ≤ 7  → Duodenum
        7 < pKa ≤ 8  → Ileum
        Else         → Poor absorption
    """
    site, img_file = "Poor absorption (ionized)", None

    try:
        if pka is not None:
            pka = float(pka)

            # Acidic case
            if heuristic_type == "acidic" and 2 <= pka <= 4:
                site, img_file = "Stomach", "Stomach.png"
				
            elif heuristic_type == "acidic" and 4 <= pka <= 6:
                site, img_file = "Duodenum", "Duodenum.png"

            # Basic case
            elif heuristic_type == "basic" and 6 <= pka <= 9:
                site, img_file = "Ileum", "ileum.png"
            elif heuristic_type == "basic" and pka > 9:
                site, img_file = "Duodenum", "Duodenum.png"

            # Neutral / fallback case
            elif heuristic_type == "neutral" and 5 <= pka <= 7:
                site, img_file = "Duodenum", "Duodenum.png"
            elif heuristic_type == "neutral" and 7 < pka <= 8:
                site, img_file = "Ileum", "ileum.png"
            else:
                site, img_file = "Poor absorption (ionized)", None
        else:
            # No pKa value, use type hint
            if heuristic_type == "acidic":
                site, img_file = "Stomach", "Stomach.png"
            elif heuristic_type == "basic":
                site, img_file = "Ileum", "ileum.png"

    except Exception:
        site, img_file = "Poor absorption (ionized)", None

    # Build full path if image is available
    if img_file:
        img_path = resource_path(os.path.join(images_dir, img_file))
        return site, img_path
    return site, None