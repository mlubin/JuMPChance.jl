.. _jumpchance-installation:

------------------
Installation Guide
------------------

JuMPChance requires Julia, JuMP, and a solver already installed; for instructions, see the corresponding `JuMP installation guide <http://jump.readthedocs.org/en/latest/installation.html>`_. After these are set up, follow the instructions below to install JuMPChance.

Getting JuMPChance
^^^^^^^^^^^^^^^^^^

JuMPChance is available in the Julia package manager. To install it, run::

    julia> Pkg.add("JuMPChance")

This command will install JuMPChance with `ECOS <https://github.com/JuliaOpt/ECOS.jl>`_ as the default solver.

To start using JuMPChance, it should be imported together with JuMP into the local scope::
    
    julia> using JuMP, JuMPChance
