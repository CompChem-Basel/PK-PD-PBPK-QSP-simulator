from PySide6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class MatplotlibWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(6, 4))
        self.canvas = FigureCanvas(self.figure)

        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def plot_pk_curve(self, t, conc):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        if len(t) > 0 and len(conc) > 0:
            ax.plot(t, conc, label="Concentration-Time")
            ax.set_xlabel("Time (h)")
            ax.set_ylabel("Concentration (mg/L)")
            ax.legend()
        self.canvas.draw()
