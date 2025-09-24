# gui/app.py
# ==========
# Application launcher for Qt GUI

import sys
import logging
from PyQt5.QtWidgets import QApplication
from .main_window import MainWindow


def main():
    """Launch the Qt GUI application."""
    app = QApplication(sys.argv)

    # Configure logging for dev + GUI log tab
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    window = MainWindow()
    window.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
