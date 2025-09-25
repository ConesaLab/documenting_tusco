"""
Shared optional dependency handling for tusco_selector.

This module provides a centralized way to handle optional dependencies like 
tqdm and colorama across the package, reducing code duplication.
"""

# Try importing tqdm for progress bars
try:
    from tqdm import tqdm

    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

    def tqdm(iterable, **kwargs):
        """Fallback tqdm that doesn't show progress when tqdm is not available."""
        return iterable


# Try importing colorama for colored output
try:
    import colorama
    from colorama import Fore, Style

    colorama.init()
    HAS_COLORS = True
except ImportError:
    HAS_COLORS = False

    # Create dummy color objects if colorama is not available
    class DummyColor:
        def __getattr__(self, name):
            return ""

    Fore = Style = DummyColor()

__all__ = ["tqdm", "HAS_TQDM", "Fore", "Style", "HAS_COLORS"]
