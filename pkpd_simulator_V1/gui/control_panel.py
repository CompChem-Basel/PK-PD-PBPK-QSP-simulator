from PySide6.QtWidgets import QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QComboBox, QSpinBox
from pkpd_simulator_V1.run_simulation import run_single_simulation

class ControlPanel(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window

        layout = QVBoxLayout()

        # SMILES
        layout.addWidget(QLabel("Drug SMILES:"))
        self.smiles_input = QLineEdit()
        layout.addWidget(self.smiles_input)

        # Dose
        layout.addWidget(QLabel("Dose (mg):"))
        self.dose_input = QSpinBox()
        self.dose_input.setRange(1, 10000)
        self.dose_input.setValue(100)
        layout.addWidget(self.dose_input)

        # Population
        layout.addWidget(QLabel("Population:"))
        self.pop_input = QComboBox()
        self.pop_input.addItems(["Healthy Male", "Healthy Female", "Type 1 Diabetic", "Liver Impairment"])
        layout.addWidget(self.pop_input)

        # Route
        layout.addWidget(QLabel("Route:"))
        self.route_input = QComboBox()
        self.route_input.addItems(["Oral", "IV"])
        layout.addWidget(self.route_input)

        # Run button
        self.run_button = QPushButton("Run Simulation")
        self.run_button.clicked.connect(self.run_simulation)
        layout.addWidget(self.run_button)

        self.setLayout(layout)

    def run_simulation(self):
        smiles = self.smiles_input.text().strip()
        dose = self.dose_input.value()
        pop = self.pop_input.currentText()
        route = self.route_input.currentText()

        results = run_single_simulation(smiles, dose, pop, route)
        self.main_window.add_results_tab(smiles, results)
