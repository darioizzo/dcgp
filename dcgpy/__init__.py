__all__ = ['test']
from ._core import *

# Add to expression a method to visualize the d-CGP graph
from ._graphviz_visualize import _graphviz_visualize
expression_double.visualize = _graphviz_visualize
expression_gdual_double.visualize = _graphviz_visualize
expression_gdual_vdouble.visualize = _graphviz_visualize

# Add to expression a method to simplify the d-CGP expression
from ._sympy_simplify import _sympy_simplify
expression_double.simplify = _sympy_simplify
expression_gdual_double.simplify = _sympy_simplify
expression_gdual_vdouble.simplify = _sympy_simplify

__version__ = {'major': 1, 'minor': 0, 'bugfix': 1}
