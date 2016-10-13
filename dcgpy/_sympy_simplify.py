def _sympy_simplify(self, in_sym, subs_weights = False, erc = []):
    """
    ex.simplify(self, in_sym, subs_weights, erc)

    Simplifies the d-CGP expression

    Returns the simplified d-CGP expression

    Note:
        This method requires the sympy module installed in your Python system

    Args:
        in_sym (a ``List[str]``): input symbols (its length must match the number of inputs)
        erc (a ``List[float]``): values of the ephemeral random constants (if empty their symbolic representation is used instead)
        subs_weights (a ``bool``): indicates whether to substitute the weights symbols with their values

    Returns:
        A string containing the simplified expression

    Raises:
        ValueError: if the length of in_sym does not match the number of inputs
        ValueError: if the length of erc is larger than the number of inputs
        ImportError: if the sympy module is not installed in your Python system

    Examples:
        >>> ex = dcgpy.expression_weighted_gdual_double(3,1,3,3,2,2,dcgpy.kernel_set_gdual_double(["sum","diff"])(),0)
        >>> print(ex.simplify(['x','c0','c1'],True,[1,2]))
        x + 6
    """

    n = self.get_n()

    if len(in_sym) != n:
        raise ValueError("The length of in_sym must match the number of inputs")

    if len(erc) > n:
        raise ValueError("The length of erc can't be larger than the number of inputs")

    try:
        import sympy
        from sympy.parsing.sympy_parser import parse_expr
    except ImportError:
        print("Failed to import the required module sympy")
        raise

    pe = parse_expr(self(in_sym)[0])

    # substitute the weights symbols with their values
    if subs_weights:
        r = self.get_rows()
        c = self.get_cols()
        a = self.get_arity()

        keys = ['w' + str(n + i) + '_' + str(j) for i in range(r * c) for j in range(a)]
        subs_dict = dict(zip(keys, self.get_weights()))
        pe = pe.subs(subs_dict)

    # substitute the ephemeral random constants symbols with their values
    if len(erc) > 0:
        keys = in_sym[n - len(erc):]
        subs_dict = dict(zip(keys, erc))
        pe = pe.subs(subs_dict)

    simplex = sympy.expand(pe)

    return simplex
