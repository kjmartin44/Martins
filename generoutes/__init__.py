"""Ways to move one gene's expression, ranked by their evidence.

Several independent routes, from the promoter to the protein. See
docs/how-it-works.md.
"""
__version__ = "0.1.0"
from .pipeline import run, markdown_report, write_outputs  # noqa: F401
