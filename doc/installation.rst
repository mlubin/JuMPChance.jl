.. _jumpchance-installation:

------------------
Installation Guide
------------------

JuMPChance requires Julia, JuMP, and a solver already installed; for instructions, see the corresponding `JuMP installation guide <http://jump.readthedocs.org/en/latest/installation.html>`_. After these are set up, follow the instructions below to install JuMPChance.

Getting JuMPChance
^^^^^^^^^^^^^^^^^^

JuMPChance is not not yet available in the Julia package manager. To install it, run::

    julia> Pkg.clone("https://github.com/mlubin/JuMPChance.jl.git")

This command will pull the latest version of the code from GitHub. Note that subsequent calls to ``Pkg.update()`` will also update JuMPChance, potentially with breaking changes. To prevent this, you may call ``Pkg.pin("JuMPChance")``, and later ``Pkg.free("JuMPChance")`` to allow JuMPChance to be updated when desired.

To start using JuMPChance, it should be imported together with JuMP into the local scope::
    
    julia> using JuMP, JuMPChance
