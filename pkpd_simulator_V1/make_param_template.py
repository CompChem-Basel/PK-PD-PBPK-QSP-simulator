"""
make_param_template.py
======================
Utility to generate a JSON file with default physiological parameters
for custom PK/PD/PBPK simulations.

Run:
    python make_param_template.py --output custom_params.json

Then edit the JSON and run simulation with:
    python run_simulation.py --smiles "CCO" --population custom --param_file custom_params.json
"""

import json
import argparse
import os


def main():
    parser = argparse.ArgumentParser(description="Generate template JSON for custom physiology parameters")
    parser.add_argument("--output", type=str, default="custom_params.json",
                        help="Output JSON file (default: custom_params.json)")
    args = parser.parse_args()

    default_params = {
        "Vd": 0.5,     # Volume of distribution (L/kg)
        "Cl": 0.1,     # Clearance (L/h)
        "ka": 1.0,     # Absorption rate constant (1/h)
        "Ql": 1.5,     # Liver blood flow (L/h)
        "Qk": 1.0,     # Kidney blood flow (L/h)
        "Vb": 5.0,     # Blood volume (L)
        "Vl": 1.5,     # Liver volume (L)
        "Vk": 0.3,     # Kidney volume (L)
        "CLl": 0.5,    # Intrinsic clearance in liver (L/h)
        "CLk": 0.3     # Intrinsic clearance in kidney (L/h)
    }

    with open(args.output, "w") as f:
        json.dump(default_params, f, indent=4)

    print(f"✅ Template written to {os.path.abspath(args.output)}")
    print("   Edit this JSON and use it with --population custom --param_file <file>")


if __name__ == "__main__":
    main()
