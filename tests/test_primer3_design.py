"""
pytest configuration.

This lets the tests import modules from the project's src/ directory.
Expected project layout:

DNASeqA/
├── src/
│   ├── main.py
│   ├── region_scanner.py
│   └── primer3_design.py
└── tests/
"""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIRECTORY = PROJECT_ROOT / "src"

if str(SRC_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(SRC_DIRECTORY))
