from pathlib import Path
import sys

# Makes relative imports work consistently
here = Path(__file__).resolve().parent
sys.path.append(str(here))

# Ground source of truth for version information
version_file = here / ".variants-version"

try:
    version = version_file.read_text().splitlines()[0].strip()
except FileNotFoundError:
    # When namespace is __main__ (script executed from subdir)
    version = (here.parent / ".variants-version").read_text().splitlines()[0].strip()
