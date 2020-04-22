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
>>> def print_my_sum(x):
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
    eph_val (list(```))"
           + type + R"(```): Values of the ephemeral constants.
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
    idxs (``int``): indexes of the single gene to me mutated

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
>>> import matplotlib.pyplot as plt 
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
>>> import matplotlib.pyplot as plt 
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
>>> import matplotlib.pyplot as plt 
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
>>> import matplotlib.pyplot as plt 
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
>>> import matplotlib.pyplot as plt 
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
>>> import matplotlib.pyplot as plt 
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
>>> import matplotlib.pyplot as plt 
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
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_P7()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_kotanchek_doc()
{
    return R"(

Generates the problem Kotanchek from the paper:

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
>>> from mpl_toolkits.mplot3d import Axes3D  
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_kotanchek()
>>> fig = plt.figure() 
>>> ax = fig.gca(projection='3d')
>>> ax.scatter(X[:,0], X[:,1], Y)
    )";
}

std::string generate_salutowicz_doc()
{
    return R"(

Generates the problem Salutowicz from the paper:

Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:

.. math::
   y = e^{-x} x^3 \cos x\sin x (\cos x \sin^2 x - 1)

x is sampled in one hundred points uniformly sampled in [0.5,10].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_salutowicz()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_salutowicz2d_doc()
{
    return R"(

Generates the problem Salutowicz2D from the paper:

Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:

.. math::
   y = e^{-x} x^3 \cos x_1\sin x_1 (\cos x_1 \sin^2 x_1 - 1) * (x_2 - 5)

:math:`x_1` and :math:`x_2` are sampled in 601 randomly selected points in [0.05,10]x[0.05,10].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> from mpl_toolkits.mplot3d import Axes3D  
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_kotanchek()
>>> fig = plt.figure() 
>>> ax = fig.gca(projection='3d')
>>> ax.scatter(X[:,0], X[:,1], Y)
    )";
}

std::string generate_uball5d_doc()
{
    return R"(

Generates the problem UBall5D from the paper:

Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:

.. math::
   y = \frac{10}{5 + \sum_{i=1}^5 (x_i-3)^2}

:math:`x_i` are sampled in 1024 randomly selected points in :math:`[0.05,6.05]^5`.

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_uball5d()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_ratpol3d_doc()
{
    return R"(

Generates the problem RatPol3D from the paper:

Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:

.. math::
   y = 30 \frac{(x_1 - 3)(x_3 - 1)}{x_2^2(x_1-10)}

:math:`x_1`, :math:`x_2`, :math:`x_3` are sampled in 300 randomly selected points in [0.05,2] x [1,2].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_ratpol3d()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_sinecosine_doc()
{
    return R"(

Generates the problem SineCosine from the paper:

Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:

.. math::
   y = 6 \sin(x_1)\cos(x_2)

:math:`x_1`, :math:`x_2` are sampled in 30 randomly selected points in [0.1,5.9] x [0.1,5.9].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> from mpl_toolkits.mplot3d import Axes3D  
>>> import matplotlib.pyplot as plt 
>>> fig = plt.figure() 
>>> ax = fig.gca(projection='3d')
>>> ax.scatter(X[:,0], X[:,1], Y)
    )";
}

std::string generate_ripple_doc()
{
    return R"(

Generates the problem Ripple from the paper:

Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:

.. math::
   y = (x_1-3)(x_2-3) + 2\sin((x_1-4)(x_2-4))

:math:`x_1`, :math:`x_2` are sampled in 300 randomly selected points in [0.05,6.05] x [0.05,6.05].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> from mpl_toolkits.mplot3d import Axes3D  
>>> import matplotlib.pyplot as plt 
>>> fig = plt.figure() 
>>> ax = fig.gca(projection='3d')
>>> ax.scatter(X[:,0], X[:,1], Y)
    )";
}

std::string generate_ratpol2d_doc()
{
    return R"(

Generates the problem RatPol2D from the paper:

Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:

.. math::
   y = \frac{(x_1-3)^4+(x_2-3)^3-(x_2-3)}{(x_2-2)^4+10}

:math:`x_1`, :math:`x_2` are sampled in 50 randomly selected points in [0.05,6.05] x [0.05,6.05].

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> from mpl_toolkits.mplot3d import Axes3D  
>>> import matplotlib.pyplot as plt 
>>> fig = plt.figure() 
>>> ax = fig.gca(projection='3d')
>>> ax.scatter(X[:,0], X[:,1], Y)
    )";
}

std::string generate_chwirut1_doc()
{
    return R"(

These data are the result of a NIST study involving ultrasonic calibration. The response variable is ultrasonic response, 
and the predictor variable is metal distance. (see https://www.itl.nist.gov/div898/strd/nls/data/chwirut1.shtml)

A proposed good model for such a problem is:

.. math::
   y = \frac{e^{-\beta_1 x}}{\beta_2 + \beta_3 x} + \epsilon

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_chwirut1()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_chwirut2_doc()
{
    return R"(

These data are the result of a NIST study involving ultrasonic calibration. The response variable is ultrasonic response, 
and the predictor variable is metal distance. (see https://www.itl.nist.gov/div898/strd/nls/data/chwirut2.shtml)

A proposed good model for such a problem is:

.. math::
   y = \frac{e^{-\beta_1 x}}{\beta_2 + \beta_3 x} + \epsilon

with respect to the problem chwirut1, less points are included here.

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_chwirut2()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_daniel_wood_doc()
{
    return R"(

These data and model are described in Daniel and Wood (1980), and originally published in E.S.Keeping, 
"Introduction to Statistical Inference," Van Nostrand Company, Princeton, NJ, 1962, p. 354. The response variable is energy
radieted from a carbon filament lamp per cm**2 per second, and the predictor variable is the absolute temperature 
of the filament in 1000 degrees Kelvin. (see https://www.itl.nist.gov/div898/strd/nls/data/daniel_wood.shtml)

A proposed good model for such a problem is:

.. math::
   y = \beta_1 x^{\beta_2} + \epsilon

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_daniel_wood()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_gauss1_doc()
{
    return R"(

The data are two well-separated Gaussians on a decaying exponential baseline plus normally 
distributed zero-mean noise with variance = 6.25. (see https://www.itl.nist.gov/div898/strd/nls/data/gauss1.shtml)

A proposed good model for such a problem is:

.. math::
   y = \beta_1 e^{-\beta_2 x} + \beta_3 e^{-\frac{(x-\beta_4)^2}{\beta_5^2}} + \beta_6 e^{-\frac{(x-\beta_7)^2}{\beta_8^2}} + \epsilon

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_gauss1()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_kirby2_doc()
{
    return R"(

These data are the result of a NIST study involving scanning electron microscope line with standards. 151 
observations are included. (see https://www.itl.nist.gov/div898/strd/nls/data/kirby2.shtml)

A proposed good model for such a problem is:

.. math::
   y = \frac{\beta_1 + \beta_2 x + \beta_3 x^2}{1 + \beta_4 x + \beta_5 x^2} + \epsilon

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_kirby2()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_lanczos2_doc()
{
    return R"(

These data are taken from an example discussed in Lanczos (1956). The data were generated to 6-digits
of accuracy using the formula below. (see https://www.itl.nist.gov/div898/strd/nls/data/lanczos2.shtml)

A good model for such a problem is, trivially:

.. math::
   y = \beta_1 e^{-\beta_2 x} + \beta_3 e^{-\beta_4 x} + \beta_5 e^{-\beta_6 x} + \epsilon

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_lanczos2()
>>> plt.plot(X,Y, '.')
    )";
}

std::string generate_misra1b_doc()
{
    return R"(

These data are the result of a NIST study involving dental research in monomolecular adsorption. 
The response variable is volume, and the predictor variable is pressure. 14 observations are 
available. (see https://www.itl.nist.gov/div898/strd/nls/data/misra1b.shtml)

A good model for such a problem is:

.. math::
   y = \beta_1 \left( 1 - \frac 1{\left(1 + \beta_1 \frac x2\right)^2}\right) + \epsilon

Returns:
    A tuple containing the points (``2D NumPy float array``) and labels (``2D NumPy float array``).

Examples:

>>> from dcgpy import *
>>> import matplotlib.pyplot as plt 
>>> X, Y = generate_misra1b()
>>> plt.plot(X,Y, '.')
    )";
}

std::string symbolic_regression_doc()
{
    return R"(

Symbolic regression is a type of regression analysis that searches the space of mathematical expressions to 
find the model that best fits a given dataset, both in terms of accuracy and simplicity 
(ref: https://en.wikipedia.org/wiki/Symbolic_regression). It also is one of the applications
for Differentiable Cartesian Genetic Programming.

This class provides an easy way to instantiate symbolic regression problems as optimization problems having
a continuous part (i.e. the value of the parameters in the model) and an integer part (i.e. the representation of
the model computational graph). The instantiated object can be used as UDP (User Defined Problem) in the pygmo optimization suite.

The symbolic regression problem can be instantiated both as a single and as a two-objectives problem. In the second
case, aside the chosen loss on the data, the model complexity will be considered as an objective.

    )";
}

std::string symbolic_regression_init_doc()
{
    return R"(__init__(points, labels, rows = 1, columns=16, levels_back=17, arity=2, kernels, n_eph=0, multi_objective=False, parallel_batches=0, loss="MSE")

Constructs a symbolic_regression optimization problem compatible with the pagmo UDP interface.

Args:
    points (2D NumPy float array or ``list of lists`` of ``float``): the input data
    labels (2D NumPy float array or ``list of lists`` of ``float``): the output data (to be predicted)
    rows (``int``): number of rows in the cartesian program
    columns (``int``): number of columns in the cartesian program
    levels_back (``int``): number of levels-back in the cartesian program
    arity (``int`` on ``list``): arity of the kernels. Assumed equal for all columns.
    kernels (``List[dcgpy.kernel_]``): kernel functions
    n_eph (``int``): Number of ephemeral constants. 
    multi_objective (``bool``): when True the problem will be considered as multiobjective (loss and model complexity).
    parallel_batches (``int``): allows to split the data into batches for parallel evaluation.
    loss (``str``): loss type used, one of "MSE" (for mean squared error) or "CE" (for cross entropy).

Raises:
    unspecified: any exception thrown by failures at the intersection between C++ and Python (e.g.,
      type conversion errors, mismatched function signatures, etc.)

Examples:
    >>> import dcgpy
    >>> import pygmo as pg
    >>> X, Y = dcgpy.generate_koza_quintic()
    >>> udp = dcgpy.symbolic_regression(
    ...        points = X,
    ...        labels = Y,
    ...        rows = 1,
    ...        cols = 20,
    ...        levels_back = 21,
    ...        arity = 2, 
    ...        kernels = dcgpy.kernel_set_double(["sum", "diff"])(), 
    ...        n_eph = 1, 
    ...        multi_objective = True, 
    ...        parallel_batches = 0)   
    >>> prob = pg.problem(udp)
    >>> print(prob) # doctest: +NORMALIZE_WHITESPACE
    Problem name: a CGP symbolic regression problem
        Global dimension:			62
        Integer dimension:			61
        Fitness dimension:			2
        Number of objectives:			2
        Equality constraints dimension:		0
        Inequality constraints dimension:	0
        Lower bounds: [-10, 0, 0, 0, 0, ... ]
        Upper bounds: [10, 1, 1, 1, 1, ... ]
        Has batch fitness evaluation: false
        <BLANKLINE>
        Has gradient: true
        User implemented gradient sparsity: true
        Expected gradients: 1
        Has hessians: true
        User implemented hessians sparsity: true
        Expected hessian components: [1, 1]
        <BLANKLINE>
        Fitness evaluations: 0
        Gradient evaluations: 0
        Hessians evaluations: 0
        <BLANKLINE>
        Thread safety: basic
        <BLANKLINE>
    Extra info:
        Data dimension (in): 1
        Data dimension (out): 1
        Data size: 10
        Kernels: [sum, diff]
)";
}

std::string symbolic_regression_predict_doc()
{
    return R"(predict(points, chromosome)

Predicts the labels of *points* using the model encoded in *x*.

Args:
    points (2D NumPy float array or ``list of lists`` of ``float`` or the 1D equivalents): the input point / points
    chromosome (1D NumPy array or ``list`` of ``int``):  the encoded model.

Raises:
    unspecified: any exception thrown by failures at the intersection between C++ and Python (e.g.,
      type conversion errors, mismatched function signatures, etc.)
)";
}

std::string generic_uda_get_seed_doc()
{
    return R"(get_seed()
This method will return the random seed used internally by this uda.

Returns:
    ``int``: the random seed of the population
)";
}

std::string es4cgp_doc()
{
    return R"(__init__(gen = 1, max_mut = 4, ftol = 1e-4, learn_constants = False, seed = random)

Evolutionary strategies are popular global optimization meta-heuristics essentially based
on the following simple pseudo-algorithm:

* Start from a population (pop) of dimension N

*  while i < gen

*  > > Mutation: create a new population pop2 mutating N times the best individual

*  > > Evaluate all new chromosomes in pop2

*  > > Reinsertion: set pop to contain the best N individuals taken from pop and pop2

The key to the success of such a search strategy is in the quality of its mutation operator. In the
case of chrosomoses that encode a Cartesian Genetic Program (CGP), it makes sense to have mutation act
on active genes only (that is on that part of the chromosome that is actually expressed in the
final CGP / formula / model). This introduces a coupling between the optimization problem (say a symbolic
regression problem) and its solution strategy which, although not preventing, makes the use of general purpose
optimization algorithms inefficient (e.g. a generic evolutionary strategy would have a mutation operator which
is agnostic of the existence of active genes).

.. note::
    ES4CGP is tailored to solve :class:`dcgpy.symbolic_regression` problems and will not work on different types.

In this class we provide an evolutionary strategy tailored to solve :class:`dcgpy.symbolic_regression` problems
leveraging the kowledge on the genetic structure of Cartesian Genetic Programs (i.e. able to mutate only active
genes).

Args:
    gen (``int``): number of generations.
    max_mut (``int``): number of active genes to be mutated.
    ftol (``float``): the algorithm will exit when the loss is below this tolerance.
    learn_constants (``bool``): when true a gaussian mutation is applied to the ephemeral constants (std = 0.1).
    use_bfe (``bool``):  when true the fitness evaluation will happen in parallel batches using the pygmo::bfe_mp.
    seed (``int``): seed used by the internal random number generator (default is random).

Raises:
    unspecified: any exception thrown by failures at the intersection between C++ and Python (e.g.,
      type conversion errors, mismatched function signatures, etc.)
    ValueError: if  *max_mut* is 0 or *ftol* is negative.

.. note::
    When the use_bfe argument is True the algorithm cannot be used in a pygmo.archipelago as nested parallelism
    would lead to AssertionError: daemonic processes are not allowed to have children
    )";
}

std::string es4cgp_get_log_doc()
{
    return R"(get_log()
Returns a log containing relevant parameters recorded during the last call to ``evolve()``. The log frequency depends
on the verbosity parameter (by default nothing is logged) which can be set calling the
method :func:`~pygmo.algorithm.set_verbosity()` on an :class:`~pygmo.algorithm`
constructed with a :class:`~dcgpy.es4cgp`. A verbosity of ``N`` implies a log
line each ``N`` generations.

Returns:
    ``list`` of ``tuples``: at each logged epoch, the values ``Gen``, ``Fevals``, ``Current best``, ``Best``, where:

    * ``Gen`` (``int``), generation number.
    * ``Fevals`` (``int``), number of functions evaluation made.
    * ``Best`` (``float``), the best fitness found.
    * ``Constants`` (``list``), the current values for the ephemeral constants.
    * ``Model`` (``string``), the string representation of the current best model
Examples:
    >>> import dcgpy
    >>> from pygmo import *
    >>> 
    >>> algo = algorithm(dcgpy.es4cgp(gen = 2000, max_mut = 4, ftol = 1e-4, learn_constants=True))       
    >>> X, Y = dcgpy.generate_koza_quintic()    
    >>> udp = dcgpy.symbolic_regression(X, Y ,1,20,21,2, dcgpy.kernel_set_double(["sum", "diff", "mul"])(), 1, False, 0)
    >>> pop = population(udp, 4)
    >>> algo.set_verbosity(200)
    >>> pop = algo.evolve(pop) # doctest: +SKIP
    Gen:        Fevals:          Best:   Constants:   Model:
       0              0        7398.14   [-1.22497]   [x0*c1**2 + c1**2] ...
     200            800        233.979   [-1.34118]   [x0*(x0 - x0**2 + x0**3 + (x0 - x0**2)** ...
     400           1600        4.26131   [-1.15376]   [x0*(x0 + x0*(c1 + x0**2) - x0**2 + (x0  ...
     600           2400        4.26126   [-1.15198]   [x0*(x0 + x0*(c1 + x0**2) - x0**2 + (x0  ...
     800           3200        4.26126   [-1.15198]   [x0*(x0 + x0*(c1 + x0**2) - x0**2 + (x0  ...
    1000           4000        4.26126   [-1.15198]   [x0*(x0 + x0*(c1 + x0**2) - x0**2 + (x0  ...
    1200           4800        4.26126   [-1.15198]   [x0*(x0 + x0*(c1 + x0**2) - x0**2 + (x0  ...
    1400           5600        4.26126   [-1.15198]   [x0*(x0 + x0*(c1 + x0**2) - x0**2 + (x0  ...
    1600           6400       0.664691   [-1.12614]   [x0*(x0 + x0*(c1 + x0**2) - (c1 + x0**2) ...
    1800           7200       0.664691   [-1.12614]   [x0*(x0 + x0*(c1 + x0**2) - (c1 + x0**2) ...
    2000           8000       0.664689   [-1.12548]   [x0*(x0 + x0*(c1 + x0**2) - (c1 + x0**2) ...
    Exit condition -- generations = 2000
    >>> uda = algo.extract(dcgpy.es4cgp)
    >>> uda.get_log() # doctest: +SKIP
    [(0, 0, 7398.139620548432, array([-1.22496858]), '[x0*c1**2 + c1**2]'), ...

See also the docs of the relevant C++ method :cpp:func:`dcgp::es4cgp::get_log()`.
)";
}

std::string moes4cgp_doc()
{
    return R"(__init__(gen = 1, max_mut = 4, learn_constants = False, seed = random)

Multi-Objective Evolutionary strategies are able to trade off, in symbolic regression tasks, the 
model complexity with its loss. This particular implementation follows the simple pseudo-code:

* Start from a population (pop) of dimension N

*  while i < gen

*  > > Mutation: create a new population pop2 mutating each individual in pop.

*  > > Evaluate all new chromosomes in pop2.

*  > > Reinsertion: set pop to contain the best N individuals taken from pop and pop2 according to non dominated sorting.

.. note::
    MOES4CGP is tailored to solve :class:`dcgpy.symbolic_regression` problems and will not work on different problem types.

Args:
    gen (``int``): number of generations.
    max_mut (``int``): maximum number of active genes to be mutated.
    learn_constants (``bool``): when true a gaussian mutation is applied to the ephemeral constants (std = 0.1).
    seed (``int``): seed used by the internal random number generator (default is random).

Raises:
    unspecified: any exception thrown by failures at the intersection between C++ and Python (e.g.,
      type conversion errors, mismatched function signatures, etc.)
    ValueError: if  *max_mut* is 0.

    )";
}

std::string moes4cgp_get_log_doc()
{
    return R"(get_log()
Returns a log containing relevant parameters recorded during the last call to ``evolve()``. The log frequency depends
on the verbosity parameter (by default nothing is logged) which can be set calling the
method :func:`~pygmo.algorithm.set_verbosity()` on an :class:`~pygmo.algorithm`
constructed with a :class:`~dcgpy.moes4cgp`. A verbosity of ``N`` implies a log
line each ``N`` generations.

Returns:
    ``list`` of ``tuples``: at each logged epoch, the values ``Gen``, ``Fevals``, ``Current best``, ``Best``, where:

    * ``Gen`` (``int``), generation number.
    * ``Fevals`` (``int``), number of functions evaluation made.
    * ``Best loss`` (``float``), the best fitness found.
    * ``Ndf size`` (``int``), number of models in the non dominated front.
    * ``Compl.`` (``in``), minimum complexity across the models in the non dominated front.
Examples:
    >>> import dcgpy
    >>> from pygmo import *
    >>> 
    >>> algo = algorithm(dcgpy.moes4cgp(gen = 90, max_mut = 2))       
    >>> X, Y = dcgpy.generate_koza_quintic()    
    >>> udp = dcgpy.symbolic_regression(X, Y ,1,20,21,2, dcgpy.kernel_set_double(["sum", "diff", "mul"])(), 1, True, 0)
    >>> pop = population(udp, 100)
    >>> algo.set_verbosity(10)
    >>> pop = algo.evolve(pop) # doctest: +SKIP
    Gen:        Fevals:     Best loss: Ndf size:   Compl.:
       0              0        6.07319         3        92
      10           1000        2.15419         5        10
      20           2000        1.92403         8        33
      30           3000       0.373663        12        72
      40           4000        0.36954        13        72
      50           5000       0.235749        16        73
      60           6000       0.235749        12        73
      70           7000       0.235749        13        73
      80           8000       0.217968        12        75
      90           9000       0.217968        12        75
     100          10000       0.217968        12        75
     110          11000       0.217968        14        75
     120          12000       0.217968        14        75
     130          13000       0.217968        13        75
     140          14000       0.162293        12        52
    Exit condition -- generations = 140
    >>> uda = algo.extract(dcgpy.moes4cgp)
    >>> uda.get_log() # doctest: +SKIP
    [(0, 0, 6.0731942123423, 3, 92), ...

See also the docs of the relevant C++ method :cpp:func:`dcgp::moes4cgp::get_log()`.
)";
}

std::string gd4cgp_doc()
{
    return R"(__init__(max_iter = 1, lr = 1., lr_min = 1e-3)

In a symbolic regression problem, models parameters are typically present in the form of
ephemeral constants (i.e. extra input terminals). The actual values of the constants have a profound
effect on the resulting loss and its an open question how to balance the learning of the model parameters
(continuous optimization) with learning the model itself (integer optimization).

.. note::
    GD4CGP is tailored to solve :class:`dcgpy.symbolic_regression` problems and will not work on different types.

In this class we provide a simple gradient descent algorithm able to tackle :class:`dcgpy.symbolic_regression` problems.
The gradient descent will only modify the continuous part of the chromosome, leaving the integer part (i.e. the 
actual model) unchanged.

Args:
    gen (``int``): maximum number of iterations.
    lr (``float``): initial learning rate (or step size).
    lr_min (``float``): stopping criteria on the minimum value for the learning rate (or step size).

Raises:
    unspecified: any exception thrown by failures at the intersection between C++ and Python (e.g.,
      type conversion errors, mismatched function signatures, etc.)
    ValueError: if  *lr_min* is smaller than 0, or larger than *lr*
    )";
}

std::string gd4cgp_get_log_doc()
{
    return R"(get_log()
Returns a log containing relevant parameters recorded during the last call to ``evolve()``. The log frequency depends
on the verbosity parameter (by default nothing is logged) which can be set calling the
method :func:`~pygmo.algorithm.set_verbosity()` on an :class:`~pygmo.algorithm`
constructed with a :class:`~dcgpy.gd4cgp`. A verbosity of ``N`` implies a log
line each ``N`` iterations.

Returns:
    ``list`` of ``tuples``: at each logged epoch, the values ``Gen``, ``Fevals``, ``Current best``, ``Best``, where:

    * ``Gen`` (``int``), generation number.
    * ``Fevals`` (``int``), number of functions evaluation made.
    * ``Gevals`` (``int``), number of gradient evaluation made.
    * ``grad norm`` (``float``), norm of the loss gradient.
    * ``lr`` (``float``), the current learning rate.
    * ``Best`` (``float``), current fitness value.

Examples:
    >>> import dcgpy
    >>> from pygmo import *
    >>> 
    >>> algo = algorithm(dcgpy.gd4cgp(4, 0.1, 1e-4))       
    >>> X, Y = dcgpy.generate_koza_quintic()    
    >>> udp = dcgpy.symbolic_regression(X, Y ,1,20,21,2, dcgpy.kernel_set_double(["sum", "diff", "mul"])(), 1, False, 0)
    >>> pop = population(udp, 10)
    >>> algo.set_verbosity(1)
    >>> pop = algo.evolve(pop) # doctest: +SKIP
    Iter:        Fevals:        Gevals:     grad norm:            lr:         Best:
        0              0              0              0            0.1         4588.6
        1              1              1        687.738           0.15        4520.41
        2              2              2        676.004          0.225        4420.33
        3              3              3        658.404         0.3375        4275.16
        4              4              4        632.004        0.50625        4068.54
    Exit condition -- max iterations = 4
    >>> uda = algo.extract(gd4cgp)
    >>> uda.get_log() # doctest: +SKIP
    [(0, 0, 0, 0.0, 0.1, 4588.5979303850145), ...

See also the docs of the relevant C++ method :cpp:func:`dcgp::gd4cgp::get_log()`.
)";
}

std::string mes4cgp_doc()
{
    return R"(__init__(gen = 1, max_mut = 1, ftol = 1e-4, learn_constants = False, seed = random)

The term Memetic is widely used, in the context of meta-heuristic search, to indicate a synergy between any
population-based approach with local improvement procedures. The resulting algorithms are also referred to, in the
literature, as Baldwinian evolutionary algorithms (EAs), Lamarckian EAs, cultural algorithms, or genetic local
searches. The very same approach, is seen by many just as an hybridization of a global search technique with a
local search technique. Regardless of the terminology and point of view, a memetic approach is applicable to symbolic
regression tasks and able to improve considerably on the long standing issue of finding constants in
Genetic Programming.

Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
In European Conference on Genetic Programming, pp. 35-51. Springer, 2017.

In this class we offer an UDA (User Defined Algorithm for the pagmo optimization suite) hybridizing the classic
Evolutionary Strategy that is traditionally used in Cartesian Genetic Programming research with a second order Newton
search step able to help finding the best values for the ephemeral constants. The resulting algorithm is
outlined by the following pseudo-algorithm:

* Start from a population (pop) of dimension N

*  while i < gen

*  > > Mutation: create a new population pop2 mutating N times the best individual

*  > > Life long learning: apply a one step of a second order Newton method to each individual (only the continuous part is affected) 

*  > > Reinsertion: set pop to contain the best N individuals taken from pop and pop2

.. note::
    MES4CGP is tailored to solve :class:`dcgpy.symbolic_regression` problems and will not work on different types.

Args:
    gen (``int``): number of generations.
    max_mut (``int``): number of active genes to be mutated.
    ftol (``float``): the algorithm will exit when the loss is below this tolerance.
    seed (``int``): seed used by the internal random number generator (default is random).

Raises:
    unspecified: any exception thrown by failures at the intersection between C++ and Python (e.g.,
      type conversion errors, mismatched function signatures, etc.)
    ValueError: if  *max_mut* is 0 or *ftol* is negative.

    )";
}

std::string mes4cgp_get_log_doc()
{
    return R"(get_log()
Returns a log containing relevant parameters recorded during the last call to ``evolve()``. The log frequency depends
on the verbosity parameter (by default nothing is logged) which can be set calling the
method :func:`~pygmo.algorithm.set_verbosity()` on an :class:`~pygmo.algorithm`
constructed with a :class:`~dcgpy.mes4cgp`. A verbosity of ``N`` implies a log
line each ``N`` generations.

Returns:
    ``list`` of ``tuples``: at each logged epoch, the values ``Gen``, ``Fevals``, ``Current best``, ``Best``, where:

    * ``Gen`` (``int``), generation number.
    * ``Fevals`` (``int``), number of functions evaluation made.
    * ``Best`` (``float``), the best fitness found.
    * ``Constants`` (``list``), the current values for the ephemeral constants.
    * ``Model`` (``string``), the string representation of the current best model
Examples:
    >>> import dcgpy
    >>> from pygmo import *
    >>> 
    >>> algo = algorithm(dcgpy.mes4cgp(gen = 90, max_mut = 4, ftol = 1e-4))       
    >>> X, Y = dcgpy.generate_koza_quintic()    
    >>> udp = dcgpy.symbolic_regression(X, Y ,1,20,21,2, dcgpy.kernel_set_double(["sum", "diff", "mul"])(), 1, False, 0)
    >>> pop = population(udp, 4)
    >>> algo.set_verbosity(10)
    >>> pop = algo.evolve(pop) # doctest: +SKIP
    Gen:        Fevals:          Best:   Constants:   Model:
       0              0        2802.82    [5.35943]    [c1**2] ...
      10             40        948.839    [10.9722]    [x0**2*c1] ...
      20             80        823.816    [8.38173]    [(c1 + x0)*x0**2] ...
      30            120        473.274    [4.48466]    [x0**3*c1] ...
      40            160        338.735    [24.2287]    [-x0 + x0**2*c1 - (c1 + x0*c1) + x0**2] ...
      50            200        107.126    [24.2287]    [x0**2*(-x0 - x0**2 + x0**3)] ...
      60            240        10.2064    [0.844799]   [x0**2*(-(c1 + x0**2) + x0**3)] ...
      70            280        10.2064    [0.844799]   [x0**2*(-(c1 + x0**2) + x0**3)] ...
      80            320         6.3605    [1.03424]    [x0**2*(x0**3*c1 - (c1 + x0**2*c1))] ...
      90            360         6.3605    [1.03424]    [x0**2*(x0**3*c1 - (c1 + x0**2*c1))] ...
    Exit condition -- generations = 90
    >>> uda = algo.extract(dcgpy.mes4cgp)
    >>> uda.get_log() # doctest: +SKIP
    [(0, 0, 2802.8212344354, array([5.35943212]), '[c1**2]'), ...

See also the docs of the relevant C++ method :cpp:func:`dcgp::mes4cgp::get_log()`.
)";
}

std::string momes4cgp_doc()
{
    return R"(__init__(gen = 1, max_mut = 1, seed = random)

Symbolic regression tasks seek for good mathematical models to represent input data. By increasing
the model complexity it is always (theoretically) possible to find almost perfect fits of any input data.
As a consequence, the model complexity must be traded off with its accuracy so that symbolic regression
is, ultimately, a two-objectives optimization problem. 

In this class we offer an UDA (User Defined Algorithm for the pygmo optimization suite) which extends
:class:`dcgpy.mes4cgp` for a multiobjective problem. The resulting algorithm, is
outlined by the following pseudo-algorithm:

* Start from a population (pop) of dimension N

*  while i < gen

*  > > Mutation: create a new population pop2 mutating N times the best individual

*  > > Life long learning: apply a one step of a second order Newton method to each individual (only the continuous part is affected) 

*  > > Reinsertion: set pop to contain the best N individuals taken from pop and pop2 according to non dominated sorting.

.. note::
    MOMES4CGP is tailored to solve :class:`dcgpy.symbolic_regression` problems and will not work on different types.

Args:
    gen (``int``): number of generations.
    max_mut (``int``): maximum number of active genes to be mutated. The minimum is zero (this allow the memetic learning to act in consecutive gens)
    seed (``int``): seed used by the internal random number generator (default is random).

Raises:
    unspecified: any exception thrown by failures at the intersection between C++ and Python (e.g.,
      type conversion errors, mismatched function signatures, etc.)
    ValueError: if  *max_mut* is 0.

    )";
}

std::string momes4cgp_get_log_doc()
{
    return R"(get_log()
Returns a log containing relevant parameters recorded during the last call to ``evolve()``. The log frequency depends
on the verbosity parameter (by default nothing is logged) which can be set calling the
method :func:`~pygmo.algorithm.set_verbosity()` on an :class:`~pygmo.algorithm`
constructed with a :class:`~dcgpy.mes4cgp`. A verbosity of ``N`` implies a log
line each ``N`` generations.

Returns:
    ``list`` of ``tuples``: at each logged epoch, the values ``Gen``, ``Fevals``, ``Current best``, ``Best``, where:

    * ``Gen`` (``int``), generation number.
    * ``Fevals`` (``int``), number of functions evaluation made.
    * ``Best loss`` (``float``), the best fitness found.
    * ``Ndf size`` (``int``), number of models in the non dominated front.
    * ``Compl.`` (``in``), minimum complexity across the models in the non dominated front.
Examples:
    >>> import dcgpy
    >>> from pygmo import *
    >>> 
    >>> algo = algorithm(dcgpy.momes4cgp(gen = 90, max_mut = 2))       
    >>> X, Y = dcgpy.generate_koza_quintic()    
    >>> udp = dcgpy.symbolic_regression(X, Y ,1,20,21,2, dcgpy.kernel_set_double(["sum", "diff", "mul"])(), 1, False, 0)
    >>> pop = population(udp, 100)
    >>> algo.set_verbosity(10)
    >>> pop = algo.evolve(pop) # doctest: +SKIP
    Gen:        Fevals:     Best loss: Ndf size:   Compl.:
       0              0        6.07319         3        92
      10           1000        2.15419         5        10
      20           2000        1.92403         8        33
      30           3000       0.373663        12        72
      40           4000        0.36954        13        72
      50           5000       0.235749        16        73
      60           6000       0.235749        12        73
      70           7000       0.235749        13        73
      80           8000       0.217968        12        75
      90           9000       0.217968        12        75
     100          10000       0.217968        12        75
     110          11000       0.217968        14        75
     120          12000       0.217968        14        75
     130          13000       0.217968        13        75
     140          14000       0.162293        12        52
    Exit condition -- generations = 140
    >>> uda = algo.extract(dcgpy.momes4cgp)
    >>> uda.get_log() # doctest: +SKIP
    [(0, 0, 6.0731942123423, 3, 92), ...

See also the docs of the relevant C++ method :cpp:func:`dcgp::momes4cgp::get_log()`.
)";
}

} // namespace dcgpy
