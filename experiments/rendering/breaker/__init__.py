"""Grouped breaker instruments; retain their existing sibling import names."""
from pathlib import Path
import sys
_directory = str(Path(__file__).resolve().parent)
if _directory not in sys.path:
    sys.path.insert(0, _directory)
