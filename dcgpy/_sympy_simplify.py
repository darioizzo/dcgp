def _ann_simplify(self, *args, **kwargs):
    print("No point in simplifying a dCGPANN expression, sorry :)")


def _sympy_simplify(self, in_sym, subs_weights = False, erc = []):
    """
    simplify(in_sym, subs_weights = False, erc = [])

    Simplifies the symbolic output of dCGP expressions

    Returns the simplified dCGP expression for each output

    Note:
        This method requires ``sympy`` and ``pyaudi`` modules installed in your Python system

    Args:
        in_sym (a ``List[str]``): input symbols (its length must match the number of inputs)
        subs_weights (a ``bool``): indicates whether to substitute the weights symbols with their values
        erc (a ``List[float]``): values of the ephemeral random constants (if empty their symbolic representation is used instead)

    Returns:
        A ``List[str]`` containing the simplified expressions

    Raises:
        ValueError: if the length of in_sym does not match the number of inputs
        ValueError: if the length of erc is larger than the number of inputs
        ImportError: if modules sympy or pyaudi are not installed in your Python system

    Examples:
        >>> ex = dcgpy.expression_weighted_gdual_double(3,2,3,3,2,2,dcgpy.kernel_set_gdual_double(["sum","diff"])(),0)
        >>> print(ex.simplify(['x','c0','c1'],True,[1,2])[0])
        x + 6
    """

    n = self.get_n()
    m = self.get_m()

    if len(in_sym) != n:
        raise ValueError("The length of in_sym must match the number of inputs")

    if len(erc) > n:
        raise ValueError("The length of erc can't be larger than the number of inputs")

    try:
        import sympy
    except ImportError:
        print("Failed to import the required module sympy")
        raise

    try:
        import pyaudi
    except ImportError:
        print("Failed to import the required module pyaudi")
        raise

    # define symbols
    pe = []
    exv = self(in_sym)
    ns = {}
    for i in range(n):
        ns[in_sym[i]] = sympy.Symbol(in_sym[i], real = True)
    if subs_weights:
        r = self.get_rows()
        c = self.get_cols()
        a = self.get_arity()
        for i in range(r * c):
            for j in range(a):
                ws = 'w' + str(n + i) + '_' + str(j)
                ns[ws] = sympy.Symbol(ws, real = True)

    # create Sympy expressions from strings
    for i in range(m):
        pe.append(sympy.sympify(exv[i], locals = ns))

    # substitute the weights symbols with their values
    if subs_weights:
        wl = self.get_weights()
        if type(wl[0]) == pyaudi.core.gdual_vdouble:
            wl = [w.constant_cf[0] for w in wl]
        for i in range(m):
            for j in range(r * c):
                for k in range(a):
                    pe[i] = pe[i].subs(ns['w' + str(n + j) + '_' + str(k)], wl[a * j + k])

    # substitute the ephemeral random constants symbols with their values
    if len(erc) > 0:
        for i in range(m):
            for j in range(len(erc)):
                pe[i] = pe[i].subs(ns[in_sym[n - len(erc) + j]],erc[j])

    # simplifications
    simplex = []
    for i in range(m):
        currex = sympy.expand(pe[i])
        if subs_weights:
            currex = sympy.collect(sympy.expand(pe[i]),list(ns.values()))
        addends = sympy.Add.make_args(currex)
        currsum = 0
        for j in range(len(addends)):
            currsum = sympy.Add(currsum, sympy.simplify(addends[j]))
        simplex.append(currsum)

    return simplex
