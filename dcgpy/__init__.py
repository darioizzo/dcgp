__all__ = ['test']
from .core import *

# Add to expression a method to visualize the d-CGP graph
from ._graphviz_visualize import _graphviz_visualize
from ._visualize import _dcgpann_visualize
expression_double.visualize = _graphviz_visualize
expression_gdual_double.visualize = _graphviz_visualize
expression_gdual_vdouble.visualize = _graphviz_visualize
expression_ann_double.visualize = _dcgpann_visualize


# Add to expression a method to simplify the d-CGP expression
from ._sympy_simplify import _sympy_simplify, _ann_simplify
expression_double.simplify = _sympy_simplify
expression_gdual_double.simplify = _sympy_simplify
expression_gdual_vdouble.simplify = _sympy_simplify
expression_ann_double.simplify = _ann_simplify

# Version setup.
from ._version import __version__