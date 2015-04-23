
-------------------------------
Solution methods and parameters
-------------------------------

Standard reformulation
^^^^^^^^^^^^^^^^^^^^^^

Consider the chance constraint

.. math::

    P\left(\sum_{i=1}^k \left(c_i^Tx +d_i\right)z_i \leq b\right) \geq 1-\epsilon

where :math:`z \sim \mathcal{N}(\mu,\Sigma)` is a vector in :math:`\mathbb{R}^n` of jointly normal random
variables with mean :math:`\mu` and covariance matrix :math:`\Sigma`. JuMPChance currently only supports a diagonal covariance matrix :math:`\Sigma`, i.e., all variables are independent, but we present the more general case here. For simplicity, we can introduce a new set of variables :math:`y_i = c_i^Tx + d_i` and reduce the constraint to:

.. math::

    P\left(y^Tz \leq b\right) \geq 1-\epsilon

Recall that :math:`y^Tz` is normally distributed with mean :math:`y^T\mu` and variance :math:`y^T\Sigma y`. Then

.. math::

    P\left(y^Tz \leq b\right) = P\left(y^Tz - y^T\mu \leq b - y^T\mu\right) = P\left( \frac{y^Tz - \mu^Tz}{\sqrt{y^T\Sigma y}} \leq \frac{b - y^T\mu}{\sqrt{y^T\Sigma y}}\right)
    
    
    = \Phi\left(\frac{b - y^T\mu}{\sqrt{y^T\Sigma y}}\right)

where :math:`\Phi` is the standard normal cumulative distribution function.

Therefore the chance constraint is satisfied if and only if

.. math::

    \Phi\left(\frac{b - y^T\mu}{\sqrt{y^T\Sigma y}}\right) \geq 1- \epsilon

or, since :math:`\Phi^{-1}` is monotonic increasing,

.. math::

    \frac{b - y^T\mu}{\sqrt{y^T\Sigma y}} \geq \Phi^{-1}(1-\epsilon)

which is

.. math::

    y^T\mu + \Phi^{-1}(1-\epsilon)\sqrt{y^T\Sigma y} \leq b.

For :math:`\epsilon \leq 0`, :math:`\Phi^{-1}(1-\epsilon) > 0`, so the above constraint is convex and equivalent to

.. math::

    ||\Sigma^{\frac{1}{2}}y|| \leq (b-\mu^Ty)/\Phi^{-1}(1-\epsilon)

which is a `second-order conic <http://en.wikipedia.org/wiki/Second-order_cone_programming>`_ constraint, where :math:`\Sigma^{\frac{1}{2}}` is the `square root <http://en.wikipedia.org/wiki/Square_root_of_a_matrix>`_ of :math:`\Sigma`.

Methods for distributionally robust constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Following the notation in the quick start guide, a distributionally robust
chance constraint can be formulated as

.. math::

    ||\Sigma^{\frac{1}{2}}y|| \leq (b-\mu^Ty)/\Phi^{-1}(1-\epsilon)\quad \forall\, \mu \in M, \Sigma \in V

This is a convex constraint because it is the intersection of a large (possibly infinite) set of convex constraints. These are challenging to reformulate into an explicit conic form. Instead, we approximate the constraint by a sequence of linear tangents, i.e., given a point :math:`y`, we detect if the constraint is violated for any choice of :math:`\mu` or :math:`\Sigma`, and if so we add a separating hyperplane which is simple to compute.

solve parameters
^^^^^^^^^^^^^^^^

The ``solve`` method has the following optional keyword parameters when invoked on a JuMPChance model:

    - ``method::Symbol``, either ``:Reformulate`` to use the second-order conic formulation or ``:Cuts`` to approximate the constraints by a sequence of linear outer-approximations. Defaults to ``:Reformulate``.
    - ``linearize_objective::Bool``, either ``true`` or ``false`` indicating whether to provide a convex quadratic objective directly to the solver or to use linear outer approximations. Defaults to ``false``.
    - ``probability_tolerance::Float64``, chance constraints are considered satisfied if they actually hold with probability :math:`1-\epsilon` minus the given tolerance. Defaults to ``0.001``.
    - ``debug::Bool``, enables debugging output for the outer approximation algorithm. Defaults to ``false``.
    - ``iteration_limit::Int``, limits the number of iterations performed by the outer approximation algorithm. (In each iteration, a single linearization is added for each violated constraint.) Defaults to ``60``.
    - ``objective_linearization_tolerance::Float64``, absolute term-wise tolerance used when linearizing quadratic objectives. Defaults to ``1e-6``.
    - ``reformulate_quadobj_to_conic::Bool``, if ``true``, automatically reformulates a quadratic objective into second-order conic form. This is necessary for some solvers like Mosek or ECOS which don't support mixing quadratic and conic constraints. Defaults to ``false`` except if the solver is ECOS.
    - ``lazy_constraints::Bool``, if ``true`` and ``method==:Cuts``, use `lazy constraints <http://jump.readthedocs.org/en/latest/callbacks.html#lazy-constraints>`_ instead of re-solving the mixed-integer relaxation in a loop. This option is experimental and not yet implemented for distributionally robust problems. Defaults to ``false``.
