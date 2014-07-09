.. _ccjump-installation:

------------------
Installation Guide
------------------

CCJuMP requires Julia, JuMP, and a solver already installed; for instructions, see the corresponding `JuMP installation guide <http://jump.readthedocs.org/en/latest/installation.html>`_. After these are set up, follow the instructions below to install CCJuMP.

Getting CCJuMP
^^^^^^^^^^^^^^

CCJuMP is not not yet available in the Julia package manager. To install it, run::

    julia> Pkg.clone("https://github.com/mlubin/CCJuMP.jl.git")

This command will pull the latest version of the code from GitHub. Note that subsequent calls to ``Pkg.update()`` will also update CCJuMP, potentially with breaking changes. To prevent this, you may call ``Pkg.pin("CCJuMP")``, and later ``Pkg.free("CCJuMP")`` to allow CCJuMP to be updated when desired.

To start using CCJuMP, it should be imported together with JuMP into the local scope::
    
    julia> using JuMP, CCJuMP
