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

A report outlining these results in more detail is currently in preparation.
