.. contents::

The algebra of truncated polynomials 
====================================
Here we introduce, formally, a basic algebraic structure over the set of truncated polynomials and we show how such a structure allows to compute the partial derivatives of multivariate functions up to arbitrary order.

Formal definition 
-----------------
Consider the set :math:`\mathcal P_{n,m}` of all polynomials of order :math:`\le m` in :math:`n` variables and having
coefficients in :math:`\mathbf K`. We indicate with the symbols :math:`T_f, T_g, T_h`, etc. the generic members of such a set. Such a set is an algebra over the field :math:`\mathbf K` if we introduce
the truncated multiplication as the standard polynomial multiplication truncated at order :math:`m`.
When needed, we will indicate such a multiplication with the symbol :math:`T_f \cdot T_g`.

This algebra is commonly referred to as the algebra of truncated polynomials. A first important
property of this algebra is that, under the multiplication, all polynomials having a zero constant
coefficient are nil-potent of order :math:`m+1`, as easily verified. We will indicate the generic
truncated polynomial :math:`\in \mathcal P_{n,m}` as :math:`T_f` and often we will consider its constant part
separated from the rest, writing :math:`T_f = f_0 + \hat f`.
It is worth noting at this point how such an algebra is unitary and associative.
The first property, in particular, deserves a few more words as it is a property that the
algebra of (non-truncated) polynomials does not possess. Formally
:math:`\forall T_f \in \mathcal P_{n,m}, \: \exists !\: T_g\in \mathcal P_{n,m}  \Rightarrow T_g\cdot T_f = 1`.
In practice:

.. math::
   T_g = \frac 1{T_f} =  \frac 1{f_0} \left(1 +\sum_{k=1}^m (-1)^k (\hat f / f_0)^k\right)
   :label: reciprocal

as its easily verified by performing the truncated multiplication :math:`T_g \cdot T_f` and accounting for the nilpotency of :math:`\hat f`. 

The link to Taylor expansions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We make use of the multi-index notation according to which

:math:`\alpha = (\alpha_1, ..,\alpha_n)` and :math:`\mathbf x = (x_1, .., x_n)`
are n-tuples and the Taylor expansion around the point :math:`\mathbf a` to order
:math:`m` of a multivariate function :math:`f` of :math:`\mathbf x` is written as:

.. math::
   T_f(\mathbf x) = \sum_{|\alpha| = 0}^m  \frac{(\mathbf x-\mathbf a)^\alpha}{\alpha!}(\partial^\alpha f)(\mathbf a)

where:

.. math::
   \partial^\alpha f = \frac{\partial^{|\alpha|} f}{\partial^{\alpha_1} x_1\partial^{\alpha_2} x_2\dots\partial^{\alpha_n} x_n}

.. math::
   \alpha ! = \prod_{i=j}^n \alpha_j

.. math::
   |\alpha| = \sum_{j=0}^n \alpha_j

The summation :math:`\sum_{|\alpha| = 0}^n` must then be taken over all possible
combinations of :math:`\alpha_j \in N` such that :math:`\sum_{j=1}^n \alpha_j = |\alpha|`. 
The expression above, i.e. the Taylor expansion truncated
at order :math:`m` of a generic function :math:`f`, is a polynomial 
:math:`\in \mathcal P_{n,m}` in the variables :math:`\mathbf{dx} = \mathbf x-\mathbf a`.
We now show that if :math:`T_f, T_g \in \mathcal P_{n,m}` are Taylor expansions
of two functions :math:`f, g` then the Taylor expansion of :math:`f\pm g, fg, f/g`
can be found operating on the algebra :math:`\mathcal P`, thus computing 
:math:`T_f\pm T_g, T_f\cdot T_g, T_f/T_g`. We may thus compute high order
derivatives of multivariate functions computing their Taylor expansions
and then extracting the desired coefficient.

Multiplication
^^^^^^^^^^^^^^

We here prove that the product of two truncated Taylor expansions is the
truncated Taylor expansion of the product. We perform the proof for
:math:`n=1` as the notation is there clearer. The multivariate case is
formally identical, requiring a more complex notation. The truncated
Taylor expansion of the product between two functions :math:`f` and :math:`g` is:

.. math::
   T_{(fg)} = \sum_{k=0}^m \frac{(x-a)^k}{k!}(fg)^{(k)}

where we indicate with the superscript :math:`(i)` the :math:`i`-th derivative with
respect to the independent variable.
We show how the same expression is derived by multiplying the Taylor
expansions of :math:`f` and :math:`g`:

.. math::
   T_f T_g = \sum_{k=0}^m \frac{(x-a)^k}{k!}(f)^{(k)}\sum_{k=0}^m \frac{(x-a)^k}{k!}(g)^{(k)} = \sum_{k=0}^m c_k (x-a)^k

The coefficients :math:`c_k` in the last power series are determined as the
Cauchy product of the two Taylor series (or discrete convolution) and are:

.. math::
   c_k = \sum_{n=0}^k \frac{f^{(n)}}{n!}  \frac{g^{(k-n)}}{(k-n)!}  = \frac{1}{k!}\sum_{n=0}^k {{k}\choose{n}}(f)^{(n)}(g)^{(k-n)}

applying now the general Leibniz rule to the last expression we get:

.. math::
   c_k =  \frac{1}{k!} (fg)^{(k)}

which allows us to conclude:

.. math::
   T_{(fg)} = T_f \cdot T_g.

Reciprocal
^^^^^^^^^^

We here prove that the reciprocal of a truncated Taylor expansion,
as defined in the algebra :math:`\mathcal P_{n,m}` is the Taylor expansion of
the reciprocal. Consider the generic
function :math:`f` and its truncated Taylor expansion :math:`T_f`.
We denote with :math:`T_{(1/f)}` the truncated Taylor expansion of the
reciprocal and apply the multiplication rule to derive that, necessarily,
:math:`T_f  T_{(1/f)} = 1`. We separate the constant part of :math:`T_f` from the
rest writing :math:`T_f = f_0 +\hat f` and we compute the product between :math:`T_f`
and the definition of reciprocal:

.. math::
   \left(f_0 + \hat f\right)\frac 1f_0 \left(1 +\sum_{j=1}^m (-1)^j (\hat f / f_0)^j\right)= \frac 1f_0 \left(f_0 + \hat f\right)\left(1 - \frac{\hat f}{f_0} + \frac{\hat f^2}{f_0^2} - ... \right) = 1

which allows us to conclude:

.. math::
   T_{(1/f)} = \frac 1f_0 \left(1 +\sum_{j=1}^m (-1)^j (\hat f / f_0)^j\right)

------------------------------------------------------------------------------

Elementary functions
====================

Consider the MacLaurin expansion of a generic function :math:`g(x) = \sum g_n x^n`. 
Consider now a multivariate function :math:`\hat f(\mathbf x) = \sum_{|\alpha|=1} f_\alpha \mathbf x^\alpha` whose MacLaurin Taylor expansion does not have a constant term. The composition between these two functions will then be, trivially, 
:math:`(g \circ \hat f) (x) = \sum g_n (\sum_{|\alpha|=1} f_\alpha \mathbf x^\alpha)^n`. 
If we now truncate such an expansion to order :math:`m`, we get 
:math:`T_{g\circ f}= \sum_{n=0}^m g_n (\sum_{|\alpha|=1}^m f_\alpha \mathbf x^\alpha)^n`, 
which can be written as:

.. math::
   T_{g\circ \hat f} = T_g\circ T_{\hat f}

The above equation is called the **composition rule** and is only valid for functions whose Taylor expansion 
does not have a constant term and, is thus nil-potent of order 
:math:`m+1` in :math:`\mathcal P_{n,m}`. In  general, we cannot compute the truncated 
Taylor expansion of a composition function directly composing the truncated 
Taylor expansions. For most elementary functions, though, we can consider 
:math:`T_f = f_0 + \hat f` and use some addition formula to be able to 
''extract`` :math:`\hat f` and thus exploit its nil-potency. The details 
on how this is done differ for each particular :math:`f` considered and are 
reported in the following subsections for some commonly used functions.

 Other functions suc as tan, cosh etc. can also be treated similarly and are not reported for convenience.

Exponential
-----------
Let us consider the case of the exponential:

.. math::
   g(x) = \exp(x) = \sum_{i=0} \frac{x^i}{i!} = 1 + x + \frac {x^2}{2} + ...

We want to compute the truncated Taylor expansion of :math:`\exp(f(\mathbf x))` 
starting from the truncated Taylor expansion :math:`T_f = f_0 + \hat f`. 
We thus write:

.. math::
   (g \circ f) (\mathbf x) = \exp(f(\mathbf x)) =  \exp f_0 \exp (f(\mathbf x) - f_0)

note that, now, we can apply the **composition rule** to :math:`\exp (f(\mathbf x) - f_0)` 
since the MacLaurin Taylor expansion of :math:`f(\mathbf x) - f_0` does not have 
a constant term. Hence:

.. math::
   T_{g \circ f} = \exp f_0 T_g \circ T_{\hat f}

and, finally:

.. math::
   T_{(\exp f)} = \exp f_0 \sum_{i=0}^m \frac{\hat f^i}{i!} = \exp f_0 \left( 1 + \hat f + \frac {\hat f^2}{2!} + ... \right)
   :label: exp

Logarithm
---------
Let us consider the case of the natural logarithm:

.. math::
   g(x) = \log(x)

We want to compute the truncated Taylor expansion of 
:math:`\log(f(\mathbf x))` starting from the truncated Taylor expansion 
:math:`T_f = f_0 + \hat f`. We thus write:

.. math::
   (g \circ f) (\mathbf x) = \log(f(\mathbf x)) =  \log (f_0 + (f(\mathbf x) - f_0)) = \log f_0 + \log(1 + \frac{f(\mathbf x) - f_0}{f_0})

We can now apply the **composition rule** to get:

.. math::
   T_{g \circ f} = \log f_0 + T_{\log(1+x)} \circ \frac{\hat f}{f_0}

and, using the known expression for MacLaurin expansion of :math:`\log(1+x)`, we get:

.. math::
   T_{(\log f)} = \log f_0 + \sum_{i=1}^m (-1)^{i+1} \frac 1i \left(\frac{\hat f}{f_0}\right)^i = \log f_0 + \frac{\hat f}{f_0} - \frac 12 \left(\frac{\hat f}{f_0}\right)^2 + ...
   :label: log

Note that the above expression is only defined if :math:`f_0 \ge 0`.

Sine and cosine
---------------
Let us consider the case of the sine and cosine functions:

.. math::
   g_1(x) = \sin(x) = \sum_{i=0} (-1)^{i} \frac{x^{2i+1}}{(2i+1)!} = x - \frac{x^3}{3!} + \frac{x^5}{5!} - ... 

.. math::
   g_2(x) = \cos(x) = \sum_{i=0} (-1)^{i} \frac{x^{2i}}{(2i)!} = 1 - \frac{x^2}{2!} + \frac{x^4}{4!} - ... 


We want to compute the truncated Taylor expansion of :math:`\sin(f(\mathbf x))`, 
:math:`\cos(f(\mathbf x))` starting from the truncated Taylor expansion 
:math:`T_f = f_0 + \hat f`. We thus write:

.. math::
   (g_1 \circ f) (\mathbf x) = \sin(f(\mathbf x)) =  \sin f_0 \cos(f(\mathbf x) - f_0) + \cos f_0 \sin(f(\mathbf x) - f_0) 

.. math::
   (g_2 \circ f) (\mathbf x) = \cos(f(\mathbf x)) =  \cos f_0 \cos(f(\mathbf x) - f_0) - \sin f_0 \sin(f(\mathbf x) - f_0) 

and, applying the **composition rule** to :math:`\cos(f(\mathbf x) - f_0)` and
:math:`\sin(f(\mathbf x) - f_0)`, we get:

.. math::
   T_{(\sin f)} = \sin f_0 \left(\sum_{i=0}^{2i\le m} (-1)^{i} \frac{\hat f^{2i}}{(2i)!}\right) + \cos f_0 \left(\sum_{i=0}^{(2i+1)\le m} (-1)^{i} \frac{\hat f^{2i+1}}{(2i+1)!}\right) \\
   T_{(\cos f)} = \cos f_0 \left(\sum_{i=0}^{2i\le m} (-1)^{i} \frac{\hat f^{2i}}{(2i)!}\right) - \sin f_0 \left(\sum_{i=0}^{(2i+1)\le m} (-1)^{i} \frac{\hat f^{2i+1}}{(2i+1)!}\right)
   :label: sinandcos


Exponentiation
--------------
Let us consider the case of the power function. 

.. math::
   g(x) = x^\alpha

We want to compute the truncated Taylor expansion of :math:`f(\mathbf x)^\alpha` 
assuming to have access to the truncated Taylor expansion of :math:`f`, 
:math:`T_f = f_0 + \hat f`. We thus write:

.. math::
   (g \circ f) (\mathbf x) = f(\mathbf x) ^ \alpha =  (f_0 + (f(\mathbf x) - f_0))^\alpha = f_0^\alpha \left( 1+ \frac{f(x) - f_0}{f_0}\right)^\alpha

We can now apply the **composition rule** to get:

.. math::
   T_{f(\mathbf x)^\alpha} = f_0^\alpha \left(T_{(1+x)^\alpha}\circ \frac{\hat f}{f_0}\right) = 

.. math::
   = f_0^\alpha \sum_{k=0}^m {\alpha \choose k} \left(\frac{\hat f}{f_0}\right)^k = f_0^\alpha\left(1 + \alpha \frac{\hat f}{f_0} + \frac{\alpha (\alpha - 1)}{2}\left(\frac{\hat f}{f_0}\right)^2 + ... \right)
   :label: pow

-----------------------------------------------------------------------

Inverse functions
=================

The case for the inverse trigonometric functions and hyperbolic functions is more complex and deserves some extended discussion. 

Inverse hyperbolic tangent
--------------------------

We consider the inverse hyperbolic tangent function:

.. math::
   g(x) = \mbox{atanh} (x)

As before, we want to compute the truncated Taylor expansion of :math:`\mbox{atanh} (f(\mathbf x))` assuming to have access to the truncated Taylor expansion of :math:`f`, :math:`T_f = f_0 + \hat f`. Since there is not a convenient addition formula we must seek a different expression. We start from the identity:

.. math::
   \mbox{atanh} (x) = \frac 12 \left(\log (1 + x) - \log (1 - x)\right)

which allows us to write:

.. math::
   (g \circ f) (\mathbf x) =  \mbox{atanh} (f(\mathbf x)) = \frac 12 \left(\log (1 + f(\mathbf x)) - \log (1 - f(\mathbf x))\right)

and, using the addition formula for the logarithms, we get:

.. math::
   \mbox{atanh} (f(\mathbf x)) = \mbox{atanh} f_0 + \frac 12\left(\log\left(1+\frac{f-f_0}{1+f_0} \right) - \log\left(1-\frac{f-f_0}{1-f_0} \right)  \right)

we may now apply the composition theorem and find:

.. math::
   T_{(\mbox{atanh} f)} =  \mbox{atanh} f_0 + \frac 12 \left(T_{\log (1+x)} \circ \frac{\hat f}{1+f_0} - T_{\log (1-x)} \circ \frac{\hat f}{1-f_0}\right)

Using the Taylor expansions for the logaritmic functions we get the final expression used in AuDi:

.. math::
   T_{(\mbox{atanh} f)} =  \mbox{atanh} f_0 +\frac 12 \sum_{n=1}^m \left(\frac{1}{(1-f_0)^n} + \frac{(-1)^{n+1}}{(1+f_0)^n}\right) \frac {\hat f^n}{n}


Inverse tangent
---------------

We consider the inverse tangent function:

.. math::
   g(x) = \mbox{atan} (x)

As before, we want to compute the truncated Taylor expansion of :math:`\mbox{atanh} (f(\mathbf x))` assuming to have access to the truncated Taylor expansion of :math:`f`, :math:`T_f = f_0 + \hat f`. We can apply the following identity involving the imaginary unit :math:`i`:

.. math::
   \mbox{atan}(z) = i \mbox{atanh}(-iz)

and therefore write:

.. math::
   T_{(\mbox{atan} f)} =  \mbox{atan} f_0 +\frac i2 \sum_{k=1}^m \left(\frac{1}{(1+if_0)^k} + \frac{(-1)^{k+1}}{(1-if_0)^k}\right) \frac {\hat f^k}{k}(-i)^k

We may avoid to compute the above expression using complex arithmetic by expanding explicitly its terms. Separating the even from the odd powers of $\hat f$ we get:

.. math::
   \begin{multline*}
   T_{(\mbox{atan} f)} =  \mbox{atan} f_0 +\frac 12 \sum_{k=1}^{2k-1\le m} \left(\frac{1}{(1+if_0)^{2k-1}} + \frac{1}{(1-if_0)^{2k-1}}\right) \frac {\hat f^{2k-1}}{2k-1}(-1)^{k+1} + \\ + \frac i2 \sum_{k=1}^{2k\le m} \left(\frac{1}{(1+if_0)^{2k}} - \frac{1}{(1-if_0)^{2k}}\right) \frac {\hat f^{2k}}{2k}(-1)^k
   \end{multline*}

and, expanding further:

.. math::
   \begin{multline*}
   T_{(\mbox{atan} f)} =  \mbox{atan} f_0 +\frac 12 \sum_{k=1}^{2k-1\le m} \left(\frac{(1-if_0)^{2k-1} + (1+if_0)^{2k-1}}{(1+f_0^2)^{2k-1}}\right) \frac {\hat f^{2k-1}}{2k-1}(-1)^{k+1} + \\ + \frac i2 \sum_{k=1}^{2k\le m} \left(\frac{(1-if_0)^{2k} - (1+if_0)^{2k}}{(1+f_0^2)^{2k}}\right) \frac {\hat f^{2k}}{2k}(-1)^k
   \end{multline*}

Using the binomial theorem it is possible to show:

.. math::
   (1+if_0)^n + (1-if_0)^n = 2 + 2\sum_{j=1}^{2j\le n} {n \choose 2j} f_0^{2j}(-1)^j

and 

.. math::
   (1-if_0)^n - (1+if_0)^n = - 2i\sum_{j=1}^{2j-1\le n} {n \choose 2j-1} f_0^{2j-1}(-1)^{j+1}

which, substituted in the expression above, return a formula not involving any more imaginary units:

.. math::
   \begin{multline*}
   T_{(\mbox{atan} f)} =  \mbox{atan} f_0 + \sum_{k=1}^{2k-1\le m} \left(\frac{1 + \sum_{j=1}^{2j\le 2k-1} {2k-1 \choose 2j} f_0^{2j}(-1)^j}{(1+f_0^2)^{2k-1}}\right) \frac {\hat f^{2k-1}}{2k-1}(-1)^{k+1} + \\ + \sum_{k=1}^{2k\le m} \left(\frac{\sum_{j=1}^{2j-1\le 2k} {2k \choose 2j-1} f_0^{2j-1}(-1)^{j+1}}{(1+f_0^2)^{2k}}\right) \frac {\hat f^{2k}}{2k}(-1)^k
   \end{multline*}

This formula is rather friendly to a computer implementation and is used in AuDi.

The other inverse functions
---------------------------
To compute the other inverse functions we may now use the identities:

.. math::
   \mbox{arccosh} x = \log(d + \sqrt{d^2 - 1})


.. math::
   \mbox{arcsinh} x = \log(d + \sqrt{1 + d^2})


.. math::
   \arcsin x = \arctan(d / \sqrt{1 - d^2})


.. math::
   \arccos x = \arctan(d / \sqrt{1 - d^2})

which apply also for the Taylor expansions.

-----------------------------------------------------------------------

Practical Examples (to be done by hand)
=======================================
In the above sections we derived a number of results that allow operating 
on simple Taylor expansions to compute Taylor expansions of increasingly 
complex expressions. We summarize here those results (keep in mind that 
:math:`T_f = f_0 + \hat f`) :

.. math::
   T_{f\pm g} = T_f \pm T_g\\
   T_{fg} = T_f \cdot T_g\\
   T_{(1/f)} = \frac 1f_0 \left(1 +\sum_{k=1}^m (-1)^k (\hat f / f_0)^k\right)\\
   T_{(\exp f)} = \exp f_0 \sum_{k=0}^m \frac{\hat f^k}{k!} \\
   T_{(\log f)} = \log f_0 - \sum_{k=1}^m \frac{(-1)^k}k \left(\hat f / f_0\right)^k \\
   T_{(\sin f)} = \sin f_0 \left(\sum_{k=0}^{2k\le m} (-1)^{k} \frac{\hat f^{2k}}{(2k)!}\right) + \cos f_0 \left(\sum_{k=0}^{(2k+1)\le m} (-1)^k \frac{\hat f^{2k+1}}{(2k+1)!}\right) \\
   T_{(\cos f)} = \cos f_0 \left(\sum_{k=0}^{2k\le m} (-1)^{k} \frac{\hat f^{2k}}{(2k)!}\right) - \sin f_0 \left(\sum_{k=0}^{(2k+1)\le m} (-1)^k \frac{\hat f^{2k+1}}{(2k+1)!}\right) \\
   T_{(f^\alpha)} = f_0^\alpha \sum_{k=0}^m {\alpha \choose k} \left(\hat f / f_0\right)^k\\
   T_{(\mbox{atanh} f)} =  \mbox{atanh} f_0 +\frac 12 \sum_{n=1}^m \left(\frac{1}{(1-f_0)^n} + \frac{(-1)^{n+1}}{(1+f_0)^n}\right) \frac {\hat f^n}{n}\\
   T_{(\mbox{atan} f)} =  \mbox{atan} f_0 + \sum_{k=1}^{2k-1\le m} \left(\frac{1 + \sum_{j=1}^{2j\le 2k-1} {2k-1 \choose 2j} f_0^{2j}(-1)^j}{(1+f_0^2)^{2k-1}}\right) \frac {\hat f^{2k-1}}{2k-1}(-1)^{k+1} + \sum_{k=1}^{2k\le m} \left(\frac{\sum_{j=1}^{2j-1\le 2k} {2k \choose 2j-1} f_0^{2j-1}(-1)^{j+1}}{(1+f_0^2)^{2k}}\right) \frac {\hat f^{2k}}{2k}(-1)^k
   :label: all

It is worth mentioning here that other functions such as the inverse functions, 
the hyperbolic functions etc. can also be treated in this way. 
The above equations can be used to find Taylor expansions of increasingly 
complex functions by simply operating on the algebra :math:`\mathcal P_{n,m}`. 
Once a Taylor expansion is computed, its coefficients can be extracted to 
obtain the value of any desired derivative. We have thus built an automated 
differentiation system. While the formalism presented can, at first, appear 
complex, the system is rather simple as we hope will appear from the following 
examples. 

Example 1 - A multiplication
----------------------------

Consider the simple function of two variables:

.. math::
   f(x,y) = x + 3xy + y^2
 
Its Taylor expansion :math:`T_f \in \mathcal P_{2,2}` can be computed as:

.. math::
   T_f = T_x + 3T_x \cdot T_y + T_y\cdot T_y
 
Let us explicitly compute such an expression at the point :math:`x=3`, :math:`y=7`. The exact sequence of computations to be performed is:

.. math::
   T_x = 3 + 1 dx + 0 dy  + 0 dxdy + 0 dx^2 + 0 dy^2 
 
.. math::
   T_y = 7 + 0 dx + 1 dy  + 0 dxdy + 0 dx^2 + 0 dy^2 
 
.. math::
   T_x \cdot T_y = 21 + 7 d x + 3 d y  + 1 dxdy + 0 dx^2 + 0 dy^2 
 
and

.. math::
   T_y \cdot T_y = 49 + 0 dx + 14 dy  + 0 dxdy + 0 dx^2 + 1 dy^2 
 
We can then derive the final expression:

.. math::
   T_f = 115 + 22 dx + 23 dy +3 dxdy + 0 dx^2 + 1 dy^2 
 
and we may easily extract the derivatives comparing this expression to the generic form of a Taylor expansion:

.. math::
   f = 115, 
   \partial_x f = 22,
   \partial_y f = 23,
   \partial_{xy} f = 3,
   \partial_{xx} f = 0,
   \partial_{yy} f = 2,
 

Example 2 - A division
----------------------

Consider the simple function of two variables:

.. math::
   f = 1 / (x + 2xy + y^2) = 1 / p
 
Its Taylor expansion :math:`T_f \in \mathcal P_{2,2}` in (say) :math:`x=0`, :math:`y=1` 
can be computed as follows:

.. math::
   T_x = 0 + 1 dx + 0 dy  + 0 dxdy + 0 dx^2 + 0 dy^2 
 
.. math::
   T_y = 1 + 0 dx + 1 dy  + 0 dxdy + 0 dx^2 + 0 dy^2 
 
.. math::
   T_p =  1 + 3 dx + 2 dy +2 dxdy + 0 dx^2 + 1 dy^2 
 
and, applying the reciprocal rule, we conclude

.. math::
   T_f = ( 1 - \hat p + \hat p ^ 2 )
 
where :math:`\hat p = 3 dx + 2 dy +2 dxdy + 0 dx^2 + 1 dy^2`, hence:

.. math::
   T_f = 1 -3 dx -2 dy + 10dxdy + 9dx^2 + 3dy^2
 
which allows, as in the previous example, to compute all derivatives up to order two:

.. math::
   f = 1, 
   \partial_x f = -3,
   \partial_y f = -2,
   \partial_{xy} f = 10,
   \partial_{xx} f = 18,
   \partial_{yy} f = 6,
 
Example 3 - An exponential 
--------------------------

Consider the elementary function of two variables:

.. math::
   f = \exp(xy)
 
Its Taylor expansion :math:`T_f \in \mathcal P_{2,2}` in (say) :math:`x=1`, :math:`y=0` 
can be computed as follows:

.. math::
   T_x = 1 + 1 dx + 0 dy  + 0 dxdy + 0 dx^2 + 0 dy^2 
 
.. math::
  T_y = 0 + 0 d x + 1 dy  + 0 dxdy + 0 dx^2 + 0 dy^2 
 
.. math::
   T_x \cdot T_y = 0 + 0 dx + 1 dy  + 1 dxdy + 0 dx^2 + 0dy^2 
 
and, applying the rule for the exponential of Taylor series, we conclude:

.. math::
   T_f = 1 + dy  + dxdy + \frac 12 dy^2
 
and,

.. math::
   f = 1, 
   \partial_x f = 0,
   \partial_y f = 1,
   \partial_{xy} f = 1,
   \partial_{xx} f = 0,
   \partial_{yy} f = 1,
 

