from PySide6.QtWidgets import QWidget, QVBoxLayout, QTextEdit, QSplitter
from PySide6.QtCore import Qt
from pkpd_simulator_V1.gui.plots import MatplotlibWidget

class ResultsTab(QWidget):
    def __init__(self, results: dict):
        super().__init__()
        layout = QVBoxLayout()

        splitter = QSplitter(Qt.Vertical)

        # Text results
        text_box = QTextEdit()
        text_box.setReadOnly(True)
        for k, v in results.items():
            text_box.append(f"{k}: {v}")
        splitter.addWidget(text_box)

        # PK curve plot
        self.plot_widget = MatplotlibWidget()
        self.plot_widget.plot_pk_curve(results.get("time", []), results.get("concentration", []))
        splitter.addWidget(self.plot_widget)

        layout.addWidget(splitter)
        self.setLayout(layout)
