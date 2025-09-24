# pkpd.spec
# =========
# PyInstaller spec for PK/PD Simulator
# Forces dist folder → app/dist

import os
from PyInstaller.utils.hooks import collect_submodules
from PyInstaller.building.build_main import Analysis, PYZ, EXE, COLLECT

here = os.path.abspath(os.getcwd())

hiddenimports = [
    'pkpd_simulator_V1.run_simulation',
    'pkpd_simulator_V1.pkpd_models',
    'pkpd_simulator_V1.utils',
    'pkpd_simulator_V1.descriptors',
    'pkpd_simulator_V1.parameters',
    'pkpd_simulator_V1.ml_model',
    'pkpd_simulator_V1.gui.app',
    'pkpd_simulator_V1.gui.main_window',
]

hiddenimports += collect_submodules('rdkit')
hiddenimports += collect_submodules('matplotlib')
hiddenimports += collect_submodules('scipy')
hiddenimports += collect_submodules('PyQt5')

datas = [
    (os.path.join(here, 'data'), 'data'),
    (os.path.join(here, 'models'), 'models'),
    (os.path.join(here, 'images'), 'images'),
    (os.path.join(here, 'pkpd_simulator_V1', 'gui'), 'gui'),
]

a = Analysis(
    ['main.py'],
    pathex=[here],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=[],
    noarchive=False
)

pyz = PYZ(a.pure, a.zipped_data, cipher=None)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='PKPD_Simulator',
    debug=False,
    strip=False,
    upx=True,
    console=True,
    distpath=os.path.join(here, 'dist'),
    workpath=os.path.join(here, 'build'),
    icon=os.path.join(here, 'icon.ico') if os.path.exists(os.path.join(here, 'icon.ico')) else None
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    name='PKPD_Simulator'
)
