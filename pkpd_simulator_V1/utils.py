"""
utils.py
========
Helper functions for:
- Logging
- Safe filenames
- Resource path resolution
- SMILES rendering (2D, PNG & SVG)
- Saving simulation outputs (Excel, plots, PDFs)
- Comparison overlays and reports

This is the verbose, long-form version with expanded docstrings and
step-by-step logging for clarity.
"""

import os
import re
import sys
import logging
import datetime
import pandas as pd
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Draw

from openpyxl import load_workbook
from openpyxl.drawing.image import Image as XLImage

from reportlab.lib.pagesizes import A4
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image as RLImage,
    Table, TableStyle
)
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet


# =====================================================================
# Resource path helper
# =====================================================================
def resource_path(relative_path: str) -> str:
    """
    Resolve absolute path for bundled resources (PyInstaller support).

    Parameters
    ----------
    relative_path : str
        The relative path to the resource (e.g., "images/icon.png").

    Returns
    -------
    str
        The absolute path to the resource.
    """
    try:
        # PyInstaller creates a temp _MEIPASS folder
        base_path = sys._MEIPASS  # type: ignore
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)


# =====================================================================
# Safe filename
# =====================================================================
def safe_filename(name: str) -> str:
    """
    Sanitize a string into a safe filename.

    Parameters
    ----------
    name : str
        Input name (e.g., SMILES string).

    Returns
    -------
    str
        Safe string that can be used as a filename.
    """
    return re.sub(r'[\\/:\*\?"<>\|]', "_", name)


# =====================================================================
# Logging setup
# =====================================================================
def setup_logging(identifier: str = "sim", output_dir: str = "pkpd_output") -> str:
    """
    Configure logging for both console and file output.

    Parameters
    ----------
    identifier : str
        Unique label for log file (e.g., "sim_run").
    output_dir : str
        Directory to store log files.

    Returns
    -------
    str
        Full path to the created log file.
    """
    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"{identifier}_{timestamp}.log")

    # Reset any existing logging handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    # Ensure UTF-8 encoding for stdout/stderr
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode="w", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )

    logging.info("=" * 60)
    logging.info("PK/PD Simulation Log Started")
    logging.info("=" * 60)

    return log_file


def log_step(step_name: str, details: str = ""):
    """
    Log a step in the simulation pipeline.

    Parameters
    ----------
    step_name : str
        Short description of the step (e.g., "Featurization").
    details : str
        Optional details to append.
    """
    arrow = "->"
    if details:
        logging.info(f"[STEP] {step_name} {arrow} {details}")
    else:
        logging.info(f"[STEP] {step_name}")


# =====================================================================
# SMILES rendering
# =====================================================================
def smiles_to_2d(smiles: str, output_dir: str, prefix: str) -> dict:
    """
    Generate 2D images (PNG + SVG) for a SMILES structure.

    Parameters
    ----------
    smiles : str
        Input SMILES string.
    output_dir : str
        Directory where files should be saved.
    prefix : str
        Prefix for filenames.

    Returns
    -------
    dict
        File paths for {"png", "svg", "2d_image"}.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logging.error(f"Invalid SMILES: {smiles}")
            return {}

        safe_prefix = safe_filename(prefix)
        png_file = os.path.join(output_dir, f"{safe_prefix}_2D.png")
        svg_file = os.path.join(output_dir, f"{safe_prefix}_2D.svg")

        # PNG image
        img = Draw.MolToImage(mol, size=(400, 400))
        img.save(png_file)

        # SVG image
        drawer = Draw.MolDraw2DSVG(400, 400)
        Draw.rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        with open(svg_file, "w", encoding="utf-8") as f:
            f.write(drawer.GetDrawingText())

        logging.info(f"✔ 2D structure saved: {png_file}, {svg_file}")
        return {"png": png_file, "svg": svg_file, "2d_image": png_file}

    except Exception as e:
        logging.error(f"Error generating 2D structure: {e}")
        return {}


# =====================================================================
# Excel export
# =====================================================================
def export_to_excel(sim_results: dict, features: dict, filename: str, smiles_image: str = None):
    """
    Save simulation results + parameters to Excel, embed SMILES image.

    Parameters
    ----------
    sim_results : dict
        Dictionary with ODE simulation results.
    features : dict
        Drug features and descriptors.
    filename : str
        Output Excel filename.
    smiles_image : str
        Path to PNG image of molecule (optional).
    """
    log_step("Export", f"Saving results to {filename}")

    # Time-course DataFrame
    df = pd.DataFrame({
        "Time": sim_results["t"],
        "PK_1C": sim_results["one_compartment"].y[0],
        "PK_2C_Central": sim_results["two_compartment"].y[0],
        "PK_2C_Periph": sim_results["two_compartment"].y[1],
        "PBPK_Blood": sim_results["pbpk"].y[0],
        "PBPK_Liver": sim_results["pbpk"].y[1],
        "PBPK_Kidney": sim_results["pbpk"].y[2],
        "PD_Effect": sim_results["pd_effect"],
        "QSP_Receptor": sim_results["qsp"].y[0],
        "QSP_Activated": sim_results["qsp"].y[1],
    })

    with pd.ExcelWriter(filename, engine="openpyxl") as writer:
        # Results
        df.to_excel(writer, sheet_name="Sim_Results", index=False)

        # Parameters
        pd.DataFrame(list(features.items()), columns=["Parameter", "Value"]) \
            .to_excel(writer, sheet_name="Parameters", index=False)

    # Embed molecule image
    if smiles_image and os.path.exists(smiles_image):
        wb = load_workbook(filename)
        ws = wb["Parameters"]
        img = XLImage(smiles_image)
        img.anchor = "E2"
        ws.add_image(img)
        wb.save(filename)
        logging.info(f"✔ 2D structure embedded in Excel: {filename}")


# =====================================================================
# Plotting
# =====================================================================
def plot_time_course(sim_results: dict, output_png: str):
    """
    Save standard PK/PD/PBPK/QSP plots as a 2x2 grid.

    Parameters
    ----------
    sim_results : dict
        Simulation results with ODE solutions.
    output_png : str
        Path to save PNG plot.
    """
    log_step("Plotting", f"Saving plots to {output_png}")
    t = sim_results["t"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("PK/PD/PBPK/QSP Results", fontsize=18)

    # PK
    axes[0, 0].plot(t, sim_results["one_compartment"].y[0], label="1C")
    axes[0, 0].plot(t, sim_results["two_compartment"].y[0], label="2C Central")
    axes[0, 0].plot(t, sim_results["two_compartment"].y[1], label="2C Peripheral")
    axes[0, 0].set_title("PK Models")
    axes[0, 0].set_xlabel("Time (h)")
    axes[0, 0].set_ylabel("Conc (mg/L)")
    axes[0, 0].legend()

    # PBPK
    axes[0, 1].plot(t, sim_results["pbpk"].y[0], label="Blood")
    axes[0, 1].plot(t, sim_results["pbpk"].y[1], label="Liver")
    axes[0, 1].plot(t, sim_results["pbpk"].y[2], label="Kidney")
    axes[0, 1].set_title("PBPK Model")
    axes[0, 1].set_xlabel("Time (h)")
    axes[0, 1].set_ylabel("Conc (mg/L)")
    axes[0, 1].legend()

    # PD
    axes[1, 0].plot(t, sim_results["pd_effect"], color="purple", label="PD Effect")
    axes[1, 0].set_title("PD (Emax)")
    axes[1, 0].set_xlabel("Time (h)")
    axes[1, 0].set_ylabel("Effect")
    axes[1, 0].legend()

    # QSP
    axes[1, 1].plot(t, sim_results["qsp"].y[0], label="Receptor")
    axes[1, 1].plot(t, sim_results["qsp"].y[1], label="Activated")
    axes[1, 1].set_title("QSP Cascade")
    axes[1, 1].set_xlabel("Time (h)")
    axes[1, 1].set_ylabel("Level")
    axes[1, 1].legend()

    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    plt.close()
    logging.info(f"✔ Plots saved: {output_png}")


# =====================================================================
# PDF Reports
# =====================================================================
def make_detailed_report(results: dict, output_dir: str, overlay_file: str = None) -> str:
    """
    Generate detailed PDF report for one or more drugs.
    """
    pdf_path = os.path.join(output_dir, "detailed_report.pdf")
    doc = SimpleDocTemplate(pdf_path, pagesize=A4)
    styles = getSampleStyleSheet()
    elements = []

    absorption_expl = (
        "For weak acids (2.5 ≤ pKa ≤ 7.5, e.g., aspirin), the non-ionized form "
        "predominates in the acidic stomach, favoring absorption there. However, most weak "
        "acids are still primarily absorbed in the small intestine due to its superior surface area. "
        "Weak bases (5 ≤ pKa ≤ 11) are absorbed mainly in the small intestine. "
        "Strong acids (pKa < 2.5) and strong bases (pKa > 11) remain ionized and have poor absorption."
    )

    for smi, data in results.items():
        feats = data["features"]
        absorption = data.get("absorption_site", "Not determined")
        img_path = data.get("absorption_image", None)

        elements.append(Paragraph(f"<b>Drug:</b> {smi}", styles["Title"]))
        elements.append(Spacer(1, 12))

        elements.append(Paragraph(f"<b>Absorption Site:</b> {absorption}", styles["Normal"]))
        if img_path and os.path.exists(img_path):
            elements.append(RLImage(img_path, width=250, height=250))
            elements.append(Spacer(1, 12))

        elements.append(Paragraph(absorption_expl, styles["Normal"]))
        elements.append(Spacer(1, 12))

        # 2D structure
        if "2d_image" in data["files"]:
            elements.append(RLImage(data["files"]["2d_image"], width=200, height=200))
            elements.append(Spacer(1, 12))

        # Graphs
        if "plot" in data["files"]:
            elements.append(RLImage(data["files"]["plot"], width=400, height=250))
            elements.append(Spacer(1, 12))

        elements.append(Spacer(1, 24))

    if overlay_file:
        elements.append(Paragraph("<b>Comparison Overlay</b>", styles["Heading2"]))
        elements.append(RLImage(overlay_file, width=400, height=250))

    doc.build(elements)
    logging.info(f"✔ PDF report saved: {pdf_path}")
    return pdf_path


# =====================================================================
# Comparison Overlays
# =====================================================================
def make_comparison_overlay(results: dict, output_dir: str) -> str:
    """
    Generate PBPK overlay plot for multiple drugs.
    """
    plt.figure(figsize=(10, 6))
    for smiles, data in results.items():
        sim = data.get("sim")
        if sim is None:
            continue
        t = sim["t"]
        blood = sim["pbpk"].y[0]
        plt.plot(t, blood, label=smiles[:15])
    plt.title("PBPK Blood Comparison")
    plt.xlabel("Time (h)")
    plt.ylabel("Conc (mg/L)")
    plt.legend()
    overlay_file = os.path.join(output_dir, "comparison_overlay.png")
    plt.savefig(overlay_file, dpi=300)
    plt.close()
    logging.info(f"✔ Comparison overlay saved: {overlay_file}")
    return overlay_file


def make_comparison_report(results: dict, output_dir: str, overlay_file: str) -> str:
    """
    Generate PDF comparison report across multiple drugs.
    """
    pdf_path = os.path.join(output_dir, "comparison_report.pdf")
    doc = SimpleDocTemplate(pdf_path, pagesize=A4)
    styles = getSampleStyleSheet()
    story = []

    story.append(Paragraph("PK/PD/PBPK/QSP Comparison Report", styles["Title"]))
    story.append(Spacer(1, 12))

    for smiles, data in results.items():
        f = data["features"]
        story.append(Paragraph(f"<b>Drug:</b> {smiles}", styles["Heading2"]))
        rows = [[k, f"{v:.3f}" if isinstance(v, (int, float)) else str(v)]
                for k, v in f.items()]
        table = Table(rows, colWidths=[150, 200])
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.grey),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
            ("ALIGN", (0, 0), (-1, -1), "LEFT"),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.black),
        ]))
        story.append(table)
        story.append(Spacer(1, 12))

    story.append(Paragraph("Overlay:", styles["Heading2"]))
    story.append(RLImage(overlay_file, width=400, height=250))

    doc.build(story)
    logging.info(f"✔ Comparison PDF saved: {pdf_path}")
    return pdf_path


# =====================================================================
# Save-all wrapper
# =====================================================================
def save_all_results(sim_results: dict,
                     smiles: str,
                     output_dir: str,
                     features: dict,
                     smiles_image: dict = None) -> dict:
    """
    Save Excel, plots, and 2D image for a single drug.

    Returns
    -------
    dict
        Paths to generated files.
    """
    safe_smiles = safe_filename(smiles)
    os.makedirs(output_dir, exist_ok=True)

    excel_path = os.path.join(output_dir, f"{safe_smiles}_results.xlsx")
    png_path   = os.path.join(output_dir, f"{safe_smiles}_timecourse.png")

    export_to_excel(sim_results, features, excel_path,
                    smiles_image["png"] if smiles_image else None)
    plot_time_course(sim_results, png_path)

    results = {"excel": excel_path, "plot": png_path}
    if smiles_image:
        results.update(smiles_image)
    results["2d_image"] = smiles_image.get("png")
	
    logging.info("✔ All results saved")
    return results
