using CCJuMP, JuMP
using Base.Test
using Distributions


let
    m = CCModel()
    @defIndepNormal(m, x, mean=1, var=1)
    @test affToStr(1+x) == "(1)*x + 1"

    @defVar(m, v)

    @test affToStr((3v+1)*x+10) == "(3 v + 1)*x + 10"

    @defIndepNormal(m,z[i=1:10],mean=i,var=1)
    @test affToStr(v*z[1]+3.5) == "(v)*z[1] + 3.5"

    @test getMean(z[5]) == 5
    @test getVar(z[5]) == 1

    c = (3v+1)*x + 10 <= 20
    @test conToStr(c) == "(3 v + 1)*x + -10 <= 0"
    addConstraint(m, c, with_probability=0.05)
    @test conToStr(CCJuMP.getCCData(m).chanceconstr[1]) == "(3 v + 1)*x + -10 <= 0, with probability 0.05"

    @test_throws @defIndepNormal(m, q, mean=1, var=-1)
    @test_throws @defIndepNormal(m, q, mean=1, var=(-1,1))


end

let
    for method in [:Reformulate,:Cuts]
        m = CCModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        addConstraint(m, z*x <= -1, with_probability=0.05)

        status = solvecc(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
    end
end

# robust but no uncertainty budget
let
    m = CCModel()
    @defIndepNormal(m, x, mean=(-1,1),var=1)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, z*x <= -1, with_probability=0.05, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
    status = solvecc(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
end

# flipped signs

let
    m = CCModel()
    @defIndepNormal(m, x, mean=(-1,1),var=1)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, -z*x >= 1, with_probability=0.05, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
    status = solvecc(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
end

# quadratic objective
let
    m = CCModel()
    @defIndepNormal(m, x, mean=0,var=1)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z+2z^2)

    addConstraint(m, z*x <= -1, with_probability=0.05)
    status = solvecc(m, method=:Cuts, debug=true)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/4 1e-4
end
let
    m = CCModel()
    @defIndepNormal(m, x, mean=(-1,1),var=1)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z+2z^2)

    addConstraint(m, z*x <= -1, with_probability=0.05, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
    status = solvecc(m, method=:Cuts, debug=true)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/4 1e-4
end

# uncertainty budget for mean
let
    m = CCModel()
    @defIndepNormal(m, x, mean=(-1,1),var=1)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, z*x <= -1, with_probability=0.05, uncertainty_budget_mean=1, uncertainty_budget_variance=0)
    status = solvecc(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/(1+quantile(Normal(0,1),0.95)) 1e-6
end

# uncertain variance, but no budget
let
    m = CCModel()
    @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, z*x <= -1, with_probability=0.05, uncertainty_budget_mean=1, uncertainty_budget_variance=0)
    status = solvecc(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/(1+quantile(Normal(0,1),0.95)) 1e-6
end

# uncertain variance, with budget
let
    m = CCModel()
    @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, z*x <= -1, with_probability=0.05, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solvecc(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)) 1e-6
end

# more than one R.V.
let
    m = CCModel()
    @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))
    @defIndepNormal(m, y, mean=(-0.01,0.01),var=0.01)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, z*x + y <= -1, with_probability=0.05, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solvecc(m, method=:Cuts)
    @test status == :Optimal
    # In mathematica: Minimize[{z, z - \[Nu]*Sqrt[1.05*z^2 + 0.01] >= -1}, z]
    @test_approx_eq_eps getValue(z) -0.36431227017642165 1e-5
end
