.. _quick-start:

-----------------
Quick Start Guide
-----------------

This quick start guide will introduce the syntax of JuMPChance, again assuming
familiarity with JuMP.


Creating a model
^^^^^^^^^^^^^^^^

JuMPChance models should be created by using the following constructor::

    m = ChanceModel()

All variables and constraints are associated with this model object.
As in JuMP, solvers can be specified by using the ``solver=`` argument to the constructor.
For example::

    using CPLEX
    m = ChanceModel(solver=CplexSolver())

will set the solver to CPLEX, assuming that both CPLEX and the corresponding
Julia package are properly installed.

By default, JuMPChance will use `ECOS <https://github.com/JuliaOpt/ECOS.jl>`_,
a lightweight open-source solver which supports the conic constraints needed for the
reformulation method for solving chance-constrained problems.

Defining variables
^^^^^^^^^^^^^^^^^^

In JuMPChance, you can mix decision variables and random variables in expressions.
Decision variables are declared by using JuMP's ``@variable`` syntax.
Random variables are declared by using a similar syntax::

    @indepnormal(m, x, mean=0, var=1)

creates a single independent normal random variable with the specified
mean and variance. The ``mean`` and ``var`` arguments are always
required. Variables indexed over a given set are supported,
and the means and variances may depend on the given indices. For example::

    @indepnormal(m, x[i=1:N,j=1:M], mean=i*j, var = 2j^2)

creates an ``N`` by ``M`` matrix of independent normally distributed
random variables where the variable in index ``(i,j)`` has mean ``i*j``
and variance ``2j^2``.

Index sets do not need to be ranges; they may be arbitrary Julia lists::

    S = [:cat, :dog]
    @indepnormal(m, x[S], mean=0, var=1)

defines two variables ``x[:cat]`` and ``x[:dog]``.

Chance constraints
^^^^^^^^^^^^^^^^^^

A JuMPChance model may contain a combination of standard JuMP constraints
(linear and quadratic) and chance constraints.

Chance constraints are constraints which contain a mix of decision
variables and random variables. Products with decision variables
and random variables are allowed, but products between two decision
variables or two random variables are not currently supported. This
restriction ensures that the resulting chance constraint is a
convex constraint on the decision variables.

Mathematically, the types of constraints supported are

.. math::

    P\left(\sum_{i=1}^k \left(c_i^Tx +d_i\right)z_i \geq b\right) \geq 1- \epsilon

and

.. math::

    P\left(\sum_{i=1}^k \left(c_i^Tx +d_i\right)z_i \leq b\right) \geq 1-\epsilon

where :math:`x` are the decision variables, :math:`c_i` are coefficient vectors, :math:`d_i` and :math:`b` are scalars, :math:`z_i` are independent jointly normal random variables with provided means and variances for :math:`i=1,\ldots,k`, and :math:`\epsilon \in (0,1)`.

Chance constraints of the above form are added by using the ``@constraint`` macro. For example::

    @indepnormal(m, x, mean=0,var=1)
    @variable(m, z)

    @constraint(m, z*x >= -1, with_probability=0.95)

Adds the constraint :math:`P(z*x \geq -1) \geq 0.95`. Note that the ``with_probability`` argument specifies the *minimum* probability :math:`\epsilon` with which the constraint may be satisfied, and so should be a number close to 1.


Distributionally robust chance constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One may also specify normally distributed random variables whose parameters
(mean and variance) are uncertain, that is, known to fall within a certain interval.
These random variables with uncertain distribution are declared as follows::

    @indepnormal(m, x, mean=(-1,1), var=(20,30))

Any combination of the mean, variance, or both may be uncertain.
When these variables appear in constraints, the constraint is
interpreted to be robust, and implies that the the chance constraint
must hold for *all possible* distributions, where the set of possible
distributions will be defined more precisely below. Mathematically,
this is

.. math::

    P\left(\sum_{i=1}^k \left(c_i^Tx +d_i\right)z_i \leq b\right) \geq 1-\epsilon, \forall \text{ distributions of } z_1,\ldots,z_n

Using the above notation, let the uncertainty interval on the mean of :math:`z_i` be :math:`[\hat\mu_i - \alpha_i,\hat\mu_i + \alpha_i]` and on the variance :math:`[\hat\sigma_i^2 - \beta_i, \hat\sigma_i^2 + \beta_i]` where :math:`\alpha_i \geq 0` and :math:`\beta_i \geq 0`.

Currently JuMPChance supports only the following uncertainty sets on the means and variances:

.. math::
    M = \left\{ (\mu_1,\ldots,\mu_k) : \exists (s_1,\ldots,s_k) \text{ such that }\mu_i = \hat\mu_i + s_i, |s_i| \leq \alpha_i, \sum_{i=1}^k \frac{|s_i|}{\alpha_i} \leq \Gamma_\mu \right\}

    V = \left\{ (v_1,\ldots,v_k) : \exists (s_1,\ldots,s_k) \text{ such that }v_i = \hat\sigma_i + s_i, |s_i| \leq \beta_i, \sum_{i=1}^k \frac{|s_i|}{\beta_i} \leq \Gamma_\sigma \right\}

where :math:`\Gamma_\mu` and :math:`\Gamma_\sigma` and given (integer) constants, known as the uncertainty budgets. The interpretation of these sets is that at most :math:`\Gamma` out of :math:`k` uncertain parameters are allows to vary from their nominal values :math:`\hat\mu_i` and :math:`\hat\sigma_i^2`. This is the uncertainty set proposed by Bertsimas and Sim (2004). Note that the means and variances are allowed to vary independently.

The uncertainty budgets :math:`\Gamma_\mu` and :math:`\Gamma_\sigma` are specified as parameters to ``@constraint`` as follows::

    @constraint(m, z*x >= -1, with_probability=0.95,
        uncertainty_budget_mean=1, uncertainty_budget_variance=1)

Solving the model
^^^^^^^^^^^^^^^^^

After the model `m` has been created and all constraints added, calling::

    solve(m,method=:Cuts)

or::

    solve(m,method=:Reformulate)

will tell JuMPChance to solve the model. The available solution methods are described
in the following section.

The ``solve`` function also returns a solution status. This should be checked
to confirm that the model was successfully solved to optimality, for example::

    status = solve(m)
    if status == :Optimal
        println("Solved to optimality")
    else
        println("Not optimal, termination status $status")
    end

Optimal values of the decision variables are available by using
``getvalue``, as with JuMP.

