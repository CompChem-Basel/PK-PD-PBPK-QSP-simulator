"""
main_window.py
==============
Qt MainWindow for PK/PD Simulator v1.1
Expanded interface:
- Sidebar with drug inputs, population, route, dose, duration, advanced options.
- Results tab with sub-tabs for each drug.
    - Each drug tab contains PK, PBPK, PD, QSP, Parameters, Absorption Site.
- Logs tab for runtime messages.
"""

import os
import logging
from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QLineEdit, QTextEdit, QLabel,
    QPushButton, QComboBox, QTabWidget, QAction, QFileDialog,
    QMessageBox, QTableWidget, QTableWidgetItem, QSplitter, QStatusBar,
    QCheckBox, QSpinBox, QHBoxLayout, QScrollArea
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from pkpd_simulator_V1.run_simulation import simulate_drug


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("PK/PD Simulator v1.1")
        self.setGeometry(200, 100, 1600, 950)

        self._create_menu()
        self._create_status_bar()

        # Splitter for sidebar + workspace
        splitter = QSplitter(Qt.Horizontal)

        # Sidebar (input panel)
        self.sidebar = self._build_sidebar()
        splitter.addWidget(self.sidebar)

        # Workspace (tabs for results + logs)
        self.tabs = QTabWidget()
        splitter.addWidget(self.tabs)
        self.setCentralWidget(splitter)

        # Tabs
        self.results_tab = QTabWidget()
        self.tabs.addTab(self.results_tab, "Results")
        self._add_logs_tab()

        # Internal state
        self.results = {}

    # =========================================================
    # Menu bar
    # =========================================================
    def _create_menu(self):
        menubar = self.menuBar()

        # ------------------- FILE MENU ------------------------
        file_menu = menubar.addMenu("&File")

        open_smiles = QAction("Open SMILES File…", self)
        open_smiles.triggered.connect(self.load_smiles_file)
        file_menu.addAction(open_smiles)

        file_menu.addSeparator()

        save_results = QAction("Save Simulation Results…", self)
        save_results.triggered.connect(self.save_results)
        file_menu.addAction(save_results)

        export_report = QAction("Export Detailed Report (PDF)…", self)
        export_report.triggered.connect(self.export_report)
        file_menu.addAction(export_report)

        file_menu.addSeparator()

        exit_action = QAction("Exit", self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # ------------------- TOOLS MENU -----------------------
        tools_menu = menubar.addMenu("&Tools")

        train_action = QAction("Train ML Models", self)
        train_action.triggered.connect(self.train_models)
        tools_menu.addAction(train_action)

        clear_logs_action = QAction("Clear Logs", self)
        clear_logs_action.triggered.connect(self.clear_logs)
        tools_menu.addAction(clear_logs_action)

        # ------------------- HELP MENU ------------------------
        help_menu = menubar.addMenu("&Help")

        about_action = QAction("About PK/PD Simulator v1.1", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

    def _create_status_bar(self):
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)

    # =========================================================
    # Sidebar (Input Panel)
    # =========================================================
    def _build_sidebar(self):
        panel = QWidget()
        layout = QVBoxLayout()

        # ---------------- Drug Input ----------------
        layout.addWidget(QLabel("SMILES (one per line):"))
        self.smiles_input = QTextEdit()
        layout.addWidget(self.smiles_input)

        load_btn = QPushButton("Load SMILES File")
        load_btn.clicked.connect(self.load_smiles_file)
        layout.addWidget(load_btn)

        # ---------------- Simulation Controls ----------------
        layout.addWidget(QLabel("Population:"))
        self.pop_combo = QComboBox()
        self.pop_combo.addItems([
            "healthy_male", "healthy_female",
            "type1_diabetic", "liver_impairment_moderate"
        ])
        layout.addWidget(self.pop_combo)

        layout.addWidget(QLabel("Route:"))
        self.route_combo = QComboBox()
        self.route_combo.addItems(["oral", "iv"])
        layout.addWidget(self.route_combo)

        layout.addWidget(QLabel("Dose (mg):"))
        self.dose_input = QLineEdit("100")
        layout.addWidget(self.dose_input)

        layout.addWidget(QLabel("Duration (hours):"))
        self.duration_input = QSpinBox()
        self.duration_input.setRange(1, 240)
        self.duration_input.setValue(24)
        layout.addWidget(self.duration_input)

        # ---------------- Advanced Options ----------------
        layout.addWidget(QLabel("Advanced Options:"))
        self.chk_pbpk = QCheckBox("Enable PBPK")
        self.chk_pbpk.setChecked(True)
        layout.addWidget(self.chk_pbpk)

        self.chk_pd = QCheckBox("Enable PD")
        self.chk_pd.setChecked(True)
        layout.addWidget(self.chk_pd)

        self.chk_qsp = QCheckBox("Enable QSP")
        self.chk_qsp.setChecked(True)
        layout.addWidget(self.chk_qsp)

        # ---------------- Run / Reset ----------------
        run_layout = QHBoxLayout()
        self.run_btn = QPushButton("Run Simulation")
        self.run_btn.clicked.connect(self.run_simulation)
        reset_btn = QPushButton("Reset Inputs")
        reset_btn.clicked.connect(self.reset_inputs)
        run_layout.addWidget(self.run_btn)
        run_layout.addWidget(reset_btn)

        layout.addLayout(run_layout)
        panel.setLayout(layout)
        return panel

    # =========================================================
    # Tabs
    # =========================================================
    def _add_logs_tab(self):
        self.log_view = QTextEdit()
        self.log_view.setReadOnly(True)
        self.tabs.addTab(self.log_view, "Logs")

    # =========================================================
    # Actions
    # =========================================================
    def run_simulation(self):
        smiles_text = self.smiles_input.toPlainText().strip()
        if not smiles_text:
            QMessageBox.warning(self, "Error", "Please enter at least one SMILES.")
            return

        smiles_list = [s.strip() for s in smiles_text.splitlines() if s.strip()]
        population = self.pop_combo.currentText()
        route = self.route_combo.currentText()
        try:
            dose = float(self.dose_input.text())
        except ValueError:
            QMessageBox.warning(self, "Error", "Invalid dose value.")
            return
        duration = self.duration_input.value()

        self.status_bar.showMessage("Running simulation…")
        self.log_view.append(">>> Running simulation…")

        try:
            results = simulate_drug(smiles_list, population, route, dose, duration=duration)
            self.results = results
            self.update_output(results)
            self.status_bar.showMessage("Simulation finished.")
        except Exception as e:
            logging.error(f"Simulation failed: {e}")
            QMessageBox.critical(self, "Error", f"Simulation failed: {e}")
            self.status_bar.showMessage("Simulation failed.")

    def update_output(self, results):
        self.results_tab.clear()

        for smi, res in results.items():
            sim = res["sim"]
            feats = res["features"]

            drug_tab = QTabWidget()
            self.results_tab.addTab(drug_tab, smi[:15])

            t = sim["t"]

            # ---------------- PK Tab ----------------
            pk_tab = QWidget()
            pk_layout = QVBoxLayout()
            fig = Figure(figsize=(7, 5))
            canvas = FigureCanvas(fig)
            ax = fig.add_subplot(111)
            ax.plot(t, sim["one_compartment"].y[0], label="1C")
            ax.plot(t, sim["two_compartment"].y[0], label="2C Central")
            ax.plot(t, sim["two_compartment"].y[1], label="2C Peripheral")
            ax.set_title(f"PK Models – {smi}")
            ax.set_xlabel("Time (h)")
            ax.set_ylabel("Conc (mg/L)")
            ax.legend()
            canvas.draw()
            pk_layout.addWidget(canvas)
            pk_tab.setLayout(pk_layout)
            drug_tab.addTab(pk_tab, "PK")

            # ---------------- PBPK Tab ----------------
            pbpk_tab = QWidget()
            pbpk_layout = QVBoxLayout()
            fig = Figure(figsize=(7, 5))
            canvas = FigureCanvas(fig)
            ax = fig.add_subplot(111)
            ax.plot(t, sim["pbpk"].y[0], label="Blood")
            ax.plot(t, sim["pbpk"].y[1], label="Liver")
            ax.plot(t, sim["pbpk"].y[2], label="Kidney")
            ax.set_title(f"PBPK – {smi}")
            ax.set_xlabel("Time (h)")
            ax.set_ylabel("Conc (mg/L)")
            ax.legend()
            canvas.draw()
            pbpk_layout.addWidget(canvas)
            pbpk_tab.setLayout(pbpk_layout)
            drug_tab.addTab(pbpk_tab, "PBPK")

            # ---------------- PD Tab ----------------
            pd_tab = QWidget()
            pd_layout = QVBoxLayout()
            fig = Figure(figsize=(7, 5))
            canvas = FigureCanvas(fig)
            ax = fig.add_subplot(111)
            ax.plot(t, sim["pd_effect"], color="purple", label="PD Effect")
            ax.set_title(f"PD (Emax) – {smi}")
            ax.set_xlabel("Time (h)")
            ax.set_ylabel("Effect")
            ax.legend()
            canvas.draw()
            pd_layout.addWidget(canvas)
            pd_tab.setLayout(pd_layout)
            drug_tab.addTab(pd_tab, "PD")

            # ---------------- QSP Tab ----------------
            qsp_tab = QWidget()
            qsp_layout = QVBoxLayout()
            fig = Figure(figsize=(7, 5))
            canvas = FigureCanvas(fig)
            ax = fig.add_subplot(111)
            ax.plot(t, sim["qsp"].y[0], label="Receptor")
            ax.plot(t, sim["qsp"].y[1], label="Activated")
            ax.set_title(f"QSP Cascade – {smi}")
            ax.set_xlabel("Time (h)")
            ax.set_ylabel("Level")
            ax.legend()
            canvas.draw()
            qsp_layout.addWidget(canvas)
            qsp_tab.setLayout(qsp_layout)
            drug_tab.addTab(qsp_tab, "QSP")

            # ---------------- Parameters Tab ----------------
            param_tab = QWidget()
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            param_inner = QWidget()
            param_layout = QVBoxLayout(param_inner)

            # add 2D structure if available
            struct_img = feats.get("structure_img")
            if struct_img and os.path.exists(struct_img):
                struct_label = QLabel()
                pixmap = QPixmap(struct_img).scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation)
                struct_label.setPixmap(pixmap)
                struct_label.setAlignment(Qt.AlignCenter)
                param_layout.addWidget(struct_label)

            # parameters table
            table = QTableWidget()
            table.setColumnCount(2)
            table.setHorizontalHeaderLabels(["Parameter", "Value"])
            table.setRowCount(len(feats))
            for i, (k, v) in enumerate(feats.items()):
                table.setItem(i, 0, QTableWidgetItem(str(k)))
                table.setItem(i, 1, QTableWidgetItem(f"{v:.3f}" if isinstance(v, (int, float)) else str(v)))
            param_layout.addWidget(table)

            scroll.setWidget(param_inner)
            vbox = QVBoxLayout()
            vbox.addWidget(scroll)
            param_tab.setLayout(vbox)
            drug_tab.addTab(param_tab, "Parameters")

            # ---------------- Absorption Site Tab ----------------
            abs_tab = QWidget()
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            abs_inner = QWidget()
            abs_layout = QVBoxLayout(abs_inner)

            site = feats.get("absorption_site", "Unknown")
            abs_label = QLabel(f"Predicted absorption site: {site}")
            abs_label.setAlignment(Qt.AlignCenter)
            abs_layout.addWidget(abs_label)

            abs_img = res["features"].get("absorption_image")
            if abs_img and os.path.exists(abs_img):
                img_label = QLabel()
                img_pixmap = QPixmap(abs_img).scaled(400, 400, Qt.KeepAspectRatio)
                img_label.setPixmap(img_pixmap)
                img_label.setAlignment(Qt.AlignCenter)
                abs_layout.addWidget(img_label)
            else:
                abs_text = QLabel("(No absorption image available)")
                abs_text.setAlignment(Qt.AlignCenter)
                abs_layout.addWidget(abs_text)

            abs_layout.addStretch(1)
            scroll.setWidget(abs_inner)
            vbox = QVBoxLayout()
            vbox.addWidget(scroll)
            abs_tab.setLayout(vbox)
            drug_tab.addTab(abs_tab, "Absorption Site")

        self.log_view.append(">>> Simulation finished.")

    # =========================================================
    # Utilities
    # =========================================================
    def load_smiles_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open SMILES File", "", "Text Files (*.txt)")
        if path:
            with open(path, "r") as f:
                self.smiles_input.setText(f.read())

    def save_results(self):
        if not self.results:
            QMessageBox.warning(self, "Error", "No results to save.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save Results", "results.xlsx", "Excel Files (*.xlsx)")
        if path:
            smi, res = list(self.results.items())[0]
            excel_file = res["files"]["excel"]
            os.replace(excel_file, path)
            QMessageBox.information(self, "Saved", f"Results saved to {path}")

    def export_report(self):
        if not self.results:
            QMessageBox.warning(self, "Error", "No results to export.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export Report", "report.pdf", "PDF Files (*.pdf)")
        if path:
            smi, res = list(self.results.items())[0]
            report_file = res["files"].get("report")
            if report_file and os.path.exists(report_file):
                os.replace(report_file, path)
                QMessageBox.information(self, "Exported", f"Report exported to {path}")
            else:
                QMessageBox.warning(self, "Error", "No report available.")

    def train_models(self):
        self.log_view.append(">>> Training ML models (this may take time)…")
        from pkpd_simulator_V1.ml_model import train_all_models
        try:
            models, r2, plots = train_all_models()
            self.log_view.append(f"✔ Training finished. R² scores: {r2}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Training failed: {e}")

    def clear_logs(self):
        self.log_view.clear()

    def reset_inputs(self):
        self.smiles_input.clear()
        self.dose_input.setText("100")
        self.duration_input.setValue(24)
        self.pop_combo.setCurrentIndex(0)
        self.route_combo.setCurrentIndex(0)
        self.chk_pbpk.setChecked(True)
        self.chk_pd.setChecked(True)
        self.chk_qsp.setChecked(True)
        self.status_bar.showMessage("Inputs reset.")

    def show_about(self):
        QMessageBox.information(
            self,
            "About",
            "PK/PD Simulator v1.1\n\nMulti-scale PK/PD, PBPK, QSP simulations.\n"
            "Built with Python, RDKit, SciPy, Qt, and Matplotlib."
        )
