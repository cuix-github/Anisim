"""Compatibility entry point; implementation lives in rendering/breaker/."""
from pathlib import Path as _Path
import sys as _sys
__file__ = str(_Path(__file__).resolve().parent / 'rendering' / 'breaker' / _Path(__file__).name)
_sys.path.insert(0, str(_Path(__file__).parent))
exec(compile(_Path(__file__).read_bytes(), __file__, 'exec'), globals())
