#include <string>

#include "docstrings.hpp"

namespace dcgpy
{

std::string kernel_init_doc(const std::string &type)
{
    return R"(__init__(callable_f, callable_s, name)

Constructs a kernel function from callables.

Args:
    callable_f (``callable - List[)"
           + type + R"(] -> )" + type + R"(``): a callable taking a list of )" + type + R"( as inputs and returning a )"
           + type + R"( (the value of the kernel function evaluated on the inputs)
    callable_s (``callable - List[string] -> string``): a callable taking a list of string as inputs and returning a string (the symbolic representation of the kernel function evaluated on the input symbols)
    name (``string``): name of the kernel

Examples:

>>> from dcgpy import *
>>> def my_sum(x):
...     return sum(x)
>>> def print_my_sum(x)
...     s = "+"
...     return "(" + s.join(x) + ") "
>>> my_kernel = kernel_)"
           + type + R"((my_sum, print_my_sum, "my_sum")
    )";
}

std::string kernel_set_init_doc(const std::string &type)
{
    return R"(__init__(kernels)

Constructs a set of common kernel functions from their common name. The kernel
functions can be then retrieved via the call operator.

Args:
    kernels (``list`` of ``string``): a list of strings indicating names of kernels to use. The following are available: "sum", "diff", "mul", "div", "sig", "sin", "log", "exp"

Examples:

>>> from dcgpy import *
>>> kernels = kernel_set_)"
           + type + R"((["sum", "diff", "mul", "div"])
>>> kernels()[0](["x", "y"])
    )";
}

std::string expression_init_doc(const std::string &type)
{
    return R"(__init__(inputs, outputs, rows, columns, levels_back, arity, kernels, n_eph, seed = randint)

Constructs a CGP expression operating on )"
           + type + R"(

Args:
    inputs (``int``): number of inputs
    outputs (``int``): number of outputs
    rows (``int``): number of rows in the cartesian program
    columns (``int``): number of columns in the cartesian program
    levels_back (``int``): number of levels-back in the cartesian program
    arity (``int`` on ``list``): arity of the kernels. Assumed equal for all columns unless its specified by a list. The list must contain a number of entries equal to the number of columns.
    kernels (``List[dcgpy.kernel_)"
           + type + R"(]``): kernel functions
    n_eph (int): Number of ephemeral constants. Their values and their symbols can be set via dedicated methods.
    seed (``int``): random seed to generate mutations and chromosomes

Attributes:
    eph_val (list(```))" + type + R"(```): Values of the ephemeral constants.
    eph_symb (list(``str``)): Symbols used for the ephemeral constants.

Examples:

>>> from dcgpy import *
>>> dcgp = expression_)"
           + type + R"((1,1,1,10,11,2,kernel_set(["sum","diff","mul","div"])(), 0u, 32u)
>>> print(dcgp)
...
>>> num_out = dcgp([in])
>>> sym_out = dcgp(["x"])
    )";
}

std::string kernel_set_push_back_str_doc()
{
    return R"(**push_back(kernel_name)**

Adds one more kernel to the set by common name.

Args:
    kernel_name (``string``): a string containing the kernel name
    )";
}

std::string kernel_set_push_back_ker_doc(const std::string &type)
{
    return R"(**push_back(kernel)**

Adds one more kernel to the set.

Args:
    kernel (``dcgpy.kernel_)"
           + type + R"(``): the kernel to add
    )";
}

std::string expression_loss_doc()
{
    return R"(loss(points, labels, loss_type)

Computes the loss of the model on the data

Args:
    points (2D NumPy float array or ``list of lists`` of ``float``): the input data
    labels (2D NumPy float array or ``list of lists`` of ``float``): the output labels (supervised signal)
    loss_type (``str``): the loss, one of "MSE" for Mean Square Error and "CE" for Cross-Entropy.

Raises:
    ValueError: if *points* or *labels* are malformed or if *loss_type* is not one of the available types.
    )";
}

std::string expression_set_doc()
{
    return R"(set(chromosome)

Sets the chromosome.

Args:
    chromosome (``List[int]``): the new chromosome

Raises:
    ValueError: if the the chromosome is incompatible with the expression (n.inputs, n.outputs, levels-back, etc.)
    )";
}

std::string expression_set_f_gene_doc()
{
    return R"(set_f_gene(node_id, f_id)

Sets for a valid node (i.e. not an input node) a new kernel.

Args:
    node_id (``List[int]``): the node id
    f_id (``List[int]``): the kernel id


Raises:
    ValueError: if the node_id or f_id are  incompatible with the expression.
    )";
}

std::string expression_mutate_doc()
{
    return R"(mutate(idxs)

Mutates multiple genes within their allowed bounds.

Args:
    idxs (``List[int]``): indexes of the genes to me mutated

Raises:
    ValueError: if the index of a gene is out of bounds
    )";
}

std::string expression_weighted_set_weight_doc()
{
    return R"(set_weight(node_id, input_id, weight)

Sets a weight.

Note:
    Convention adopted for node numbering: http://ppsn2014.ijs.si/files/slides/ppsn2014-tutorial3-miller.pdf

Args:
    node_id (``int``): the id of the node whose weight is being set
    input_id (``int``): the id of the node input (0 for the first one up to arity-1)
    weight (``float``): the new value of the weight

Raises:
    ValueError: if *node_id* or *input_id* are not valid
    )";
}

std::string expression_weighted_set_weights_doc()
{
    return R"(set_weights(weights)

Sets all weights.

Args:
    weights (``List[float]``): the new values of the weights

Raises:
    ValueError: if the input vector dimension is not valid (r*c*arity)
    )";
}

std::string expression_weighted_get_weight_doc()
{
    return R"(get_weight(node_id, input_id)

Gets a weight.

Note:
    Convention adopted for node numbering: http://ppsn2014.ijs.si/files/slides/ppsn2014-tutorial3-miller.pdf

Args:
    node_id (``int``): the id of the node
    input_id (``int``): the id of the node input (0 for the first one up to arity-1)

Returns:
    The value of the weight (a ``float``)

Raises:
    ValueError: if *node_id* or *input_id* are not valid
    )";
}

std::string expression_ann_set_weight_doc()
{
    return R"(set_weight(node_id, input_id, weight)
set_weight(idx, weight)

Sets a weight. Two overloads are available. You can set the weight specifying the node and the input id (that needs
to be less than the arity), or directly specifying its position in the weight vector.

Note:
    Convention adopted for node numbering: http://ppsn2014.ijs.si/files/slides/ppsn2014-tutorial3-miller.pdf

Args:
    node_id (``int``): the id of the node whose weight is being set
    input_id (``int``): the id of the node input (0 for the first one up to arity-1)
    weight (``float``): the new value of the weight
    idx (``int``): the idx of weight to be set

Raises:
    ValueError: if *node_id* or *input_id* or *idx* are not valid
    )";
}

std::string expression_ann_get_weight_doc()
{
    return R"(get_weight(node_id, input_id)
get_weight(idx)

Gets a weight. Two overloads are available. You can get the weight specifying the node and the input id (that needs
to be less than the arity), or directly specifying its position in the weight vector.

Note:
    Convention adopted for node numbering: http://ppsn2014.ijs.si/files/slides/ppsn2014-tutorial3-miller.pdf

Args:
    node_id (``int``): the id of the node
    input_id (``int``): the id of the node input (0 for the first one up to arity-1)

Returns:
    The value of the weight (``float``)

Raises:
    ValueError: if *node_id* or *input_id* or *idx* are not valid
    )";
}

std::string expression_ann_randomise_weights_doc()
{
    return R"(randomise_weights(mean = 0, std = 0.1, seed = randomint)

Randomises all the values for the weights using a normal distribution.

Args:
    mean (``float``): the mean of the normal distribution.
    std (``float``): the standard deviation of the normal distribution.
    seed (``int``): the random seed to use.
)";
}

std::string expression_ann_set_bias_doc()
{
    return R"(set_bias(node_id, bias)

Sets a bias.

Note:
    Convention adopted for node numbering: http://ppsn2014.ijs.si/files/slides/ppsn2014-tutorial3-miller.pdf

Args:
    node_id (``int``): the id of the node whose weight is being set
    weight (``float``): the new value of the weight

Raises:
    ValueError: if *node_id* is not valid
    )";
}

std::string expression_ann_set_biases_doc()
{
    return R"(set_biases(biases)

Sets all biases.

Args:
    biases (``List[float]``): the new values of the biases

Raises:
    ValueError: if the input vector dimension is not valid (r*c)
    )";
}

std::string expression_ann_get_bias_doc()
{
    return R"(get_bias(node_id)

Gets a bias.

Note:
    Convention adopted for node numbering: http://ppsn2014.ijs.si/files/slides/ppsn2014-tutorial3-miller.pdf

Args:
    node_id (``int``): the id of the node

Returns:
    The value of the bias (a ``float``)

Raises:
    ValueError: if *node_id* is not valid
    )";
}

std::string expression_ann_randomise_biases_doc()
{
    return R"(randomise_biases(mean = 0, std = 0.1, seed = randomint)

Randomises all the values for the biases using a normal distribution.

Args:
    mean (``float``): the mean of the normal distribution.
    std (``float``): the standard deviation of the normal distribution.
    seed (``int``): the random seed to use.
)";
}

std::string expression_ann_sgd_doc()
{
    return R"(sgd(points, labels, lr, batch_size, loss_type, parallel = 0, shuffle = True)

Performs one epoch of mini-batch (stochastic) gradient descent updating the weights and biases using the 
*points* and *labels* to decrease the loss.

Args:
    points (2D NumPy float array or ``list of lists`` of ``float``): the input data
    labels (2D NumPy float array or ``list of lists`` of ``float``): the output labels (supervised signal)
    lr (``float``): the learning generate
    batch_size (``int``): the batch size
    loss_type (``str``): the loss, one of "MSE" for Mean Square Error and "CE" for Cross-Entropy.
    parallel (``int``): sets the grain for parallelism. 0 -> no parallelism n -> divides the data into n parts and processes them in parallel threads 
    shuffle (``bool``): when True it shuffles the points and labels before performing one epoch of training.


Returns:
    The average error across the batches (``float``). Note: this is only a proxy for the real loss on the whole data set.

Raises:
    ValueError: if *points* or *labels* are malformed or if *loss_type* is not one of the available types.
    )";
}

std::string expression_ann_set_output_f_doc()
{
    return R"(set_output_f(name)

Sets the nonlinearities of all nodes connected to the output nodes.
This is useful when, for example, the dCGPANN is used for a regression task where output values are expected in [-1 1]
and hence the output layer should have some sigmoid or tanh nonlinearity, or in a classification task when one wants to have a softmax
layer by having a sum in all output neurons.

Args:
    name (``string``): the kernel name

Raises:
    ValueError: if *name* is not one of the kernels in the expression.
    )";
}

std::string expression_ann_n_active_weights_doc()
{
    return R"(n_active_weights(unique = False)

Computes the number of weights influencing the result. This will also be the number
of weights that are updated when calling sgd. The number of active weights, as well as
the number of active nodes, define the complexity of the expression expressed by the chromosome.

Args:
    unique (``bool``): when True weights are counted only once if connecting the same two nodes.
    )";
}

std::string generate_koza_quintic_doc()
{
    return R"(

Generates the data for the classic Koza quintic regression problem.

.. math::
   y = x^5 - 2 x^3 + x

x is sampled in ten equally spaced points in [-3,3].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> X, Y = generate_koza_quintic()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_P1_doc()
{
    return R"(

Generates the problem P1 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
   y = x^5 - \pi x^3 + x

x is sampled in ten equally spaced points in [1,3].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> X, Y = generate_P1()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_P2_doc()
{
    return R"(

Generates the problem P2 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
   y = x^5 - \pi x^3 + \frac{\pi}{x}

x is sampled in ten equally spaced points in [0.1,5].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> X, Y = generate_P2()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_P3_doc()
{
    return R"(

Generates the problem P3 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
   y = \frac{e x^5 + x^3}{x+1}

x is sampled in ten equally spaced points in [-0.9,1].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> X, Y = generate_P3()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_P4_doc()
{
    return R"(

Generates the problem P4 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
   y = \sin(\pi x) + \frac 1x

x is sampled in ten equally spaced points in [-1,1].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> X, Y = generate_P4()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_P5_doc()
{
    return R"(

Generates the problem P5 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
   y = e x^5 - \pi x^3 + x

x is sampled in ten equally spaced points in [1,3].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> X, Y = generate_P5()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_P6_doc()
{
    return R"(

Generates the problem P6 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
   y = \frac{e x^2 - 1}{\pi (x + 2)}

x is sampled in ten equally spaced points in [-2.1,1].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> X, Y = generate_P6()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_P7_doc()
{
    return R"(

Generates the problem P7 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
   y = \cos(\pi x) + \sin(e x)

x is sampled in ten equally spaced points in [-1,1].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> X, Y = generate_P7()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_kotanchek_doc()
{
    return R"(

Generates one of the problems  from the paper:

Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:

.. math::
   y = \frac{e^{-(x_1-1)^2}}{1.2+(x_2-2.5)^2}

:math:`x_1` and :math:`x_2` are sampled in one hundred randomly selected points in [0.3,4]x[0.3,4].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> X, Y = generate_kotanchek()
>>> fig = plt.figure() 
>>> ax = fig.gca(projection='3d')
>>> ax.scatter(X[:,0], X[:,1], Y)
    )";
}
} // namespace dcgpy