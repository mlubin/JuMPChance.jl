.. _quick-start:

-----------------
Quick Start Guide
-----------------

This quick start guide will introduce the syntax of CCJuMP, again assuming
familiarity with JuMP.


Creating a Model
^^^^^^^^^^^^^^^^

CCJuMP models should be created by using the following constructor::

    m = CCModel()

All variables and constraints are associated with this model object.
As in JuMP, solvers can be specified by using the ``solver=`` argument to the constructor.
For example::

    using CPLEX
    m = CCModel(solver=CplexSolver())

will set the solver to CPLEX, assuming that both CPLEX and the corresponding
Julia package are properly installed.

Defining Variables
^^^^^^^^^^^^^^^^^^

In CCJuMP, you can mix decision variables and random variables in expressions.
Decision variables are declared by using JuMP's ``@defVar`` syntax.
Random variables are declared by using a similar syntax::

    @defIndepNormal(m, x, mean=0, var=1)

creates a single independent normal random variable with the specified
mean and variance. The ``mean`` and ``var`` arguments are always
required. Variables indexed over a given set are supported,
and the means and variances may depend on the given indices. For example::

    @defIndepNormal(m, x[i=1:N,j=1:M], mean=i*j, var = 2j^2)

creates an ``N`` by ``M`` matrix of independent normally distributed
random variables where the variable in index ``(i,j)`` has mean ``i*j``
and variance ``2j^2``.

Index sets do not need to be ranges; they may be arbitrary Julia lists::

    S = [:cat, :dog]
    @defIndepNormal(m, x[S], mean=0, var=1)

defines two variables ``x[:cat]`` and ``x[:dog]``.

Chance Constraints
^^^^^^^^^^^^^^^^^^

A CCJuMP model may contain a combination of standard JuMP constraints
(linear and quadratic) and chance constraints.

Chance constraints are constraints which contain a mix of decision
variables and random variables. Products with decision variables
and random variables are allowed, but products between two decision
variables or two random variables are not currently supported. This
restriction ensures that the resulting chance constraint is a
convex constraint on the decision variables.

Mathematically, the types of constraints supported are

.. math::

    P\left(\sum_{i=1}^k \left(c_i^Tx +d_i\right)z_i \geq b\right) \leq \epsilon

and

.. math::

    P\left(\sum_{i=1}^k \left(c_i^Tx +d_i\right)z_i \leq b\right) \leq \epsilon

where :math:`x` are the decision variables, :math:`c_i` are coefficient vectors, :math:`d_i` and :math:`b` are scalars, :math:`z_i` are independent jointly normal random variables with provided means and variances, and :math:`\epsilon \in (0,1)`.
