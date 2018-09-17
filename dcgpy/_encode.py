def encode_ffnn(inputs, outputs, layers_size, kernels, levels_back):
    """
    encode_ffnn(inputs, outputs, layers, kernels, levels_back)

    Encodes a feed forward neural network as a dCGPANN expression. While there are infinite ways to perform such an encoding
    this function generates one of the possible ones.

    Args:
        inputs (``int``): number of inputs
        outputs (``int``): number of outputs
        layers_size (``List[int]``): size of hidden neural layers 
        kernels (``List[string])``: list containing the non linearity name for each hidden neural layer plus the one for the output layer
        levels_back (``int``): number of levels-back in the cartesian program

    Returns:
        A ``expression_ann_double`` encoding the requested network.

    Raises:
        ValueError: if the kernel list contains unknown kernel names or if layers_size and kernels are malformed
    """

    from dcgpy import expression_ann_double, kernel_set_double

    if (len(kernels) != len(layers_size) + 1):
        raise ValueError("The size of layers_size must be one less than the size of kernels (as this also includes the output nonlinearity) ")

    # We create the list of possible kernels
    kernel_list = kernel_set_double(list(set(kernels)))()
    # We compute the cartesian substrate able to host the network
    arity = [inputs] + layers_size 
    extended_layer_size = layers_size + [outputs]
    cols = len(extended_layer_size)
    rows = max(extended_layer_size)
    # We create a dCGPANN expression that is able to encode the FFNN
    retval = expression_ann_double(inputs, outputs, rows, cols, levels_back, arity, kernel_list)
    # We hand write the chromosome.
    start_prev_col = [0] + [inputs + c * rows for c in range(cols)]
    x = retval.get()
    for c in range(cols):
        kernel_id = list(set(kernels)).index(kernels[c])
        for r in range(extended_layer_size[c]):
            node_id = inputs + c * rows +r
            g_idx = retval.get_gene_idx()[node_id]
            x[g_idx] = kernel_id
            for conn in range(arity[c]):
                x[g_idx+conn+1] = start_prev_col[c] + conn
    # And the output
    for i in range(outputs):
        x[-outputs+i] = start_prev_col[c+1] + i
    retval.set(x)


    return retval
