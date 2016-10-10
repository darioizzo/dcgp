__all__ = ['test']
from ._core import *

# Add to expression a method to visualize the d-CGP graph
from ._graphviz_visualize import _graphviz_visualize
expression_double.visualize = _graphviz_visualize
expression_gdual_double.visualize = _graphviz_visualize
expression_gdual_vdouble.visualize = _graphviz_visualize
