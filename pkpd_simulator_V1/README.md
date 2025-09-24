# PK/PD/PBPK/QSP Simulator

A modular Python-based simulator for pharmacokinetics (PK), physiologically-based PK (PBPK), pharmacodynamics (PD), and quantitative systems pharmacology (QSP).  
Includes ADME/T descriptors, machine learning predictions (pKa, EC50, Clearance), and visualization.

---

## 📦 Installation

```bash
# Create a clean environment (recommended)
conda create -n pkpd python=3.10 -y
conda activate pkpd

# Install dependencies
pip install -r requirements.txt


1️⃣ descriptors.py

Purpose:
This module handles all molecular-level calculations starting from the SMILES string. It extracts chemical descriptors and predicts the absorption site based on pKa or heuristic rules.

Workflow:

Take a SMILES string as input.

Convert it to an RDKit molecule object.

Compute molecular descriptors:

MWT → Molecular weight (g/mol)

LogP → Octanol-water partition coefficient (hydrophobicity)

TPSA → Topological polar surface area (Å², influences permeability)

HBD → Hydrogen bond donors

HBA → Hydrogen bond acceptors

RotatableBonds → Number of rotatable bonds (affects flexibility & absorption)

Assign a pKa (acid dissociation constant) if available; classify molecule as acidic, basic, or neutral.

Determine the primary absorption site in the gastrointestinal (GI) tract based on pKa or heuristic type:

Stomach → acidic drugs, pKa < 3

Small intestine → basic/neutral drugs, pKa 3–7 or >7

Abbreviations:

SMILES → Simplified Molecular Input Line Entry System

pKa → Acid dissociation constant, indicates how easily a molecule donates or accepts a proton

HBD/HBA → Hydrogen bond donor/acceptor

TPSA → Topological polar surface area

Use Case: This module is the first step in your pipeline for predicting PK parameters from chemical structure.

2️⃣ parameters.py

Purpose:
Stores and retrieves physiological and population-specific parameters needed for PK/PD simulations.

Workflow:

Define baseline parameters for a “generic” human:

Vd → Volume of distribution (L/kg)

Cl → Clearance (L/hr)

ka → Absorption rate constant (1/hr)

Adjust parameters for specific population groups:

healthy_male, healthy_female

type1_diabetic

liver_impairment_moderate

Adjust parameters for drug administration route:

Oral → requires ka (absorption)

IV → no absorption (ka = 0)

Abbreviations:

PK → Pharmacokinetics, study of drug absorption, distribution, metabolism, excretion

Vd → Volume of distribution, indicates how widely the drug distributes in body tissues

Cl → Clearance, volume of plasma cleared of drug per unit time

ka → Absorption rate constant

Use Case: Provides input parameters for PK models depending on patient population and administration route.

3️⃣ pkpd_models.py

Purpose:
Contains the core PK/PD/QSP/PBPK mathematical models.

Workflow:

Define PK models:

1-compartment: drug distributes in a single compartment; simple first-order elimination

2-compartment: central + peripheral compartments with inter-compartment transfer

PBPK 3-organ: physiologically-based PK model, placeholder for organ-specific drug kinetics

Define PD model:

Emax model → Effect increases with concentration until maximum effect Emax

Define QSP model (toy version):

Captures simplified systems-level dynamics of drug-target interactions

Provide high-level simulation runner:

Takes SMILES, dose, population, route

Calls descriptor functions, selects parameters, runs PK/PD/QSP models

Abbreviations:

PBPK → Physiologically-Based Pharmacokinetics, models drug distribution in organs explicitly

PD → Pharmacodynamics, study of drug effect on the body

QSP → Quantitative Systems Pharmacology, combines PK/PD with systems biology models

Use Case: The computational engine of your simulator.

4️⃣ ml_model.py

Purpose:
Handles machine learning-based dose prediction using chemical descriptors.

Workflow:

Featurize all molecules in a dataset using descriptors.py.

Train a model (RandomForest, but can be swapped) to predict dose from features.

Save trained model for future prediction.

Provide single-SMILES prediction with uncertainty.

Abbreviations:

ML → Machine Learning

RF → Random Forest (ensemble learning method)

Use Case: Predict starting dose for an unknown compound based on chemical descriptors.

5️⃣ ddi.py

Purpose:
Handles drug-drug interaction (DDI) trials between compounds.

Workflow:

Accept a list of SMILES drugs.

Generate all pairwise combinations (perpetrator/victim).

Assign effect on victim clearance or other PK parameters.

Return list of DDI scenarios.

Abbreviations:

DDI → Drug-Drug Interaction

Perpetrator → Drug causing the interaction

Victim → Drug affected by the interaction

Use Case: Simulate potential interactions in combination therapies.

6️⃣ utils.py

Purpose:
Contains helper functions for logging, exporting, plotting, and visualization.

Workflow:

setup_logging() → create log file & console logger.

export_to_excel() → write simulation results to Excel.

plot_time_course() → generate PK/PD/QSP plots.

make_timecourse_gif() → optional animation of concentration/time.

Abbreviations: None new; mainly plotting and file utilities.

Use Case: Support module to store, visualize, and track simulations.

7️⃣ run_simulation.py

Purpose:
The main entry point. Glues all modules together for an end-to-end simulation.

Workflow:

Accept input: SMILES, dose, population, route, etc.

Compute molecular descriptors (descriptors.py).

Select parameters (parameters.py).

Determine absorption site.

Run PK/PD/PBPK/QSP simulation (pkpd_models.py).

Export results (utils.py) and optionally plot or animate.

Abbreviations: Already covered (SMILES, PK, PD, PBPK, QSP).

Use Case: Acts as the CLI or script runner, orchestrating a full simulation for a compound.

✅ Summary of Workflow Across Files

Input: SMILES, dose, population, administration route → run_simulation.py

Featurization & Absorption Site → descriptors.py

Select PK/PD parameters → parameters.py

Run PK/PD/PBPK/QSP models → pkpd_models.py

Machine Learning Dose Prediction (optional) → ml_model.py

DDI scenarios (optional) → ddi.py

Logging, Excel export, plots → utils.py