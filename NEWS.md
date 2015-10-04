JuMPChance release notes
========================

Version 0.2.2 (October 4, 2015)
-------------------------------

  * Fix deprecation warnings on Julia 0.4

Version 0.2.1 (September 3, 2015)
---------------------------------

  * Update to JuMP 0.10, new minimum required version.
  * ``getVar`` renamed to ``getVariance``
  * New ``satisfied_with_probability`` to check the probability violation at the current solution
  * Performance improvements
  * Special preprocessing for cases that can be reduced to linear constraints (thanks Kaarthik Sundar)
  * For two-sided chance constraints, choice between 2-approximation and 1.25 approximation

Version 0.2.0 (May 5, 2015)
---------------------------

  * Preliminary support for two-sided chance constraints.
  * ``solvechance`` method renamed to ``solve``
  * The meaning of ``with_probability`` **has changed**. Constraints must now hold with the given probability or greater, which is tractable for 1/2 or greater. A deprecation warning is in place when a small value is provided.
  * The ``@addConstraint`` macro now accepts chance constraints.
  * Following JuMP, the use of operators <= and >= *outside of macros* has been deprecated. Use ``@addConstraint`` instead of ``addConstraint``.

Version 0.1.1 (January 10, 2015)
--------------------------------

  * Clarify and tidy up support for integer variables.

Version 0.1.0 (January 4, 2015)
-------------------------------

  * Initial release
