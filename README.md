# JuMPChance

[![Build Status](https://travis-ci.org/mlubin/JuMPChance.jl.svg?branch=master)](https://travis-ci.org/mlubin/JuMPChance.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mlubin/JuMPChance.jl?branch=master&svg=true)](https://ci.appveyor.com/project/mlubin/jumpchance-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/mlubin/JuMPChance.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/mlubin/JuMPChance.jl?branch=master)

JuMPChance (formerly CCJuMP) is an extension to [JuMP](https://github.com/JuliaOpt/JuMP.jl) for formulating and solving optimization problems with chance constraints (also known as probabilistic constraints). JuMPChance currently supports only a particular class of chance constraints involving affine combinations of jointly normal random variables, a classical formulation that's known to be efficiently solvable by using second-order conic programming (SOCP) (although JuMPChance also provides an outer-approximation algorithm which solves a sequence of linear problems).

JuMPChance supports an extension of the classical model to distributionally robust (or ambiguous) chance constraints where the parameters of the normal distributions are known to fall in a symmetric interval or more general uncertainty set.

See the [documentation](http://jumpchance.readthedocs.org/) for installation installation instructions, a quick start guide, and a more detailed discussion of the methods implemented.

Please cite JuMPChance using [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.13740.svg)](http://dx.doi.org/10.5281/zenodo.13740).
