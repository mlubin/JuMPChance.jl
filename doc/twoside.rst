----------------------------
Two-sided chance constraints
----------------------------

JuMPChance supports two-sided chance constraints of the form

.. math::

    P(l \leq x^Tz \leq u) \geq 1- \epsilon

where :math:`l`, :math:`x`, and :math:`u` are decision variables or affine functions of decision variables, :math:`z` is a vector of independent jointly normal random variables with provided means and variances, and :math:`0 < \epsilon < \frac{1}{2}`. This constraint holds iff there exists :math:`t` such that :math:`t \geq \sqrt{\sum_i (\sigma_ix_i)^2}` and :math:`(l-\mu^Tx,u-\mu^Tx,t) \in \bar S_\epsilon` where

.. math::

    \bar S_\epsilon = \operatorname{closure} \{ (l,u,t) : \Phi(u/t) - \Phi(l/t) \geq 1-\epsilon, t > 0 \}

Note that :math:`\bar S_\epsilon` is the `conic hull <http://en.wikipedia.org/wiki/Conical_combination>`_ of :math:`S_\epsilon` where

.. math::

    S_\epsilon = \{ (l,u) : \Phi(u) - \Phi(l) \geq 1-\epsilon \}

A report outlining these results is available on `arXiv <http://arxiv.org/abs/1507.01995>`_.

Using two-sided constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The syntax for two-sided constraints is as follows::

    @constraint(m, l <= x*z <= u, with_probability = 0.95, approx="1.25")
    @constraint(m, l <= x*z <= u, with_probability = 0.95, approx="2.0")

Any affine expression of the decision variables can appear as lower or upper bounds. Random variables may only appear in the expression in the middle. You can use ``sum()`` as in standard JuMP constraints, e.g.::

    @constraint(m, l <= sum( x[i]*z[i] for i=1:n ) + sum( c[j] for j=1:m ) <= u)

Given a chance constraint with probability :math:`1-\epsilon`, the current implementation provides two different formulations, indicated by the ``approx`` parameter. The ``approx`` parameter may be set to ``"1.25"`` or ``"2.0"``. The formulation guarantees that that the constraint will be satisfied with probability :math:`1-approx*\epsilon`. This is *not* a conservative approximation. After a model is solved, you can check the probability level at which a constraint holds as follows::

    constraint_ref = @constraint(m, l <= x*z <= u, with_probability = 1-0.05, approx="1.25")
    solve(m, method=:Reformulate)
    satisfied_with = JuMPChance.satisfied_with_probability(constraint_ref)
    println("The chance constraint is satisfied with probability $satisfied_with.")
