=====================================================
JuMPChance --- JuMP extensions for chance constraints
=====================================================

.. module:: JuMPChance
   :synopsis: JuMP extensions for robust chance constraints

`JuMPChance <https://github.com/mlubin/JuMPChance.jl>`_ is an extension for
`JuMP <https://github.com/JuliaOpt/JuMP.jl>`_. We assume that the reader
is familiar with the syntax of JuMP, which is described
in its `documentation <http://jump.readthedocs.org/en/latest/>`_.

JuMPChance supports a particular class of chance constraints, namely those involving affine combinations of jointly normal random variables. Distributionally robust chance constraints are supported in the form of intervals on the mean and variance of the normally distributed random variables.

JuMPChance is research code and not officially supported as part of JuMP. JuMPChance is released under the terms of the `LGPL <http://www.gnu.org/licenses/lgpl-3.0.txt>`_ license version 3:

.. code-block:: none

    JuMPChance is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License,
    Version 3, as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

Contents
--------

.. toctree::
    :maxdepth: 2

    installation.rst
    quickstart.rst
    solution.rst
    twoside.rst
    references.rst

