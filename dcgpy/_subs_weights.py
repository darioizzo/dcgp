def _subs_weights(self, expr_str):
    """
    ex.subs_weights(self, expr_str)

    Substitutes the weights in the d-CGP expression string with their value

    Returns the d-CGP expression with substituted weights

    Note:
        This method requires the sympy module installed in your Python system

    Args:
        expr_str (a ``str``): the expression

    Returns:
        A string containing the expression with substituted weights

    Raises:
        ImportError: if the sympy module is not installed in your Python system

    Examples:
        >>> ex = dcgpy.expression_weighted_double(3,1,3,3,2,2,dcgpy.kernel_set_double(["sum","diff"])(),0)
        >>> expr_str = ex.simplify(['x','c0','c1'],[1,2])
        >>> print(ex.subs_weights(expr_str))
        1.0*x + 6.0
    """

    n = self.get_n()
    r = self.get_rows()
    c = self.get_cols()
    a = self.get_arity()
    
    keys = ['w' + str(n + i) + '_' + str(j) for i in range(r * c) for j in range(a)]
    subs_dict = dict(zip(keys, self.get_weights()))
    expr_str = expr_str.subs(subs_dict)

    return expr_str
