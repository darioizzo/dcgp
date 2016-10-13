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

# Add to expression_weighted a method to substitute the weight values in the d-CGP expression string
from ._subs_weights import _subs_weights
expression_weighted_double.subs_weights = _subs_weights
expression_weighted_gdual_double.subs_weights = _subs_weights
expression_weighted_gdual_vdouble.subs_weights = _subs_weights
