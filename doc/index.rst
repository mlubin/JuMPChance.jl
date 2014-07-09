========================================================
CCJuMP --- JuMP extensions for robust chance constraints
========================================================

.. module:: CCJuMP
   :synopsis: JuMP extensions for robust chance constraints

`CCJuMP <https://github.com/mlubin/CCJuMP.jl>`_ is an extension for 
`JuMP <https://github.com/JuliaOpt/JuMP.jl>`_. We assume that the reader
is familiar with the syntax of JuMP, which is described
in its `documentation <http://jump.readthedocs.org/en/latest/>`_.

CCJuMP supports a particular class of chance constraints, namely those involving affine combinations of jointly normal random variables. Distributionally robust chance constraints are supported in the form of intervals on the mean and variance of the normally distributed random variables.

CCJuMP is currently experimental code and is not publicly released (despite the code being available). Contact the authors if you are interested in using CCJuMP.

Contents
--------

.. toctree::
    :maxdepth: 2

    installation.rst
    quickstart.rst
    solution.rst

