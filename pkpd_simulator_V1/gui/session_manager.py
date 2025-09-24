from PySide6.QtWidgets import QFileDialog, QMessageBox
import json

class SessionManager:
    def __init__(self, main_window):
        self.main_window = main_window

    def save_session(self):
        filename, _ = QFileDialog.getSaveFileName(self.main_window, "Save Session", "", "JSON Files (*.json)")
        if filename:
            data = {"sessions": []}  # Collect session data
            # TODO: pull real session info from tabs
            with open(filename, "w") as f:
                json.dump(data, f, indent=2)
            QMessageBox.information(self.main_window, "Saved", f"Session saved to {filename}")

    def load_session(self):
        filename, _ = QFileDialog.getOpenFileName(self.main_window, "Open Session", "", "JSON Files (*.json)")
        if filename:
            with open(filename, "r") as f:
                data = json.load(f)
            QMessageBox.information(self.main_window, "Loaded", f"Loaded session: {filename}")
            # TODO: recreate tabs from saved data
