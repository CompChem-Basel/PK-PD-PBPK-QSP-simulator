"""
ddi.py
======
Handles drug-drug interaction (DDI) trials between compounds.
Currently a placeholder — you can expand later with mechanistic models.
"""

def generate_ddi_trials(drug_smiles_list: list[str],
                        mode: str = "inhibition") -> list[dict]:
    """
    For each ordered pair of drugs in list, generate a DDI scenario:
    - perpetrator, victim
    - effect on victim’s hepatic clearance (factor)
    - interaction score
    """
    trials = []
    for i, smi1 in enumerate(drug_smiles_list):
        for smi2 in drug_smiles_list[i+1:]:
            trial = {
                "perpetrator": smi1,
                "victim": smi2,
                "effect": 0.5,   # Placeholder (e.g., inhibition factor)
                "score": 0.8     # Placeholder interaction score
            }
            trials.append(trial)
    return trials
