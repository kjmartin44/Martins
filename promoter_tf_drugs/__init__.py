"""Ways to move one gene's expression, ranked by their evidence.

Several independent routes, of which the promoter motif scan the
package is named after is one. See docs/how-it-works.md.
"""
__version__ = "0.1.0"
from .pipeline import run, markdown_report, write_outputs  # noqa: F401
