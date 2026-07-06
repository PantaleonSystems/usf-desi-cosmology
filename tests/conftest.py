"""
Pytest configuration for the USF pipeline validation suite.

Ensures the repository root is importable so that `theory_module` and
`usf_constants` resolve exactly as they do when Cobaya runs from the root.
"""

import os
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
