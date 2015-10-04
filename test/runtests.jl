using JuMPChance, JuMP
using FactCheck
using Distributions
using GLPKMathProgInterface
import ECOS
using Compat

facts("Operator overloads and printing") do
    m = ChanceModel()
    @defIndepNormal(m, x, mean=1, var=1)
    @defIndepNormal(m, y, mean=1, var=1)
    @defVar(m, v)
    @defVar(m, q)
    @defIndepNormal(m,z[i=1:10],mean=i,var=1)

    # Number
    # Number--IndepNormal
    @fact affToStr(1+x) --> "(1.0)*x + 1.0"
    @fact affToStr(1-x) --> "(-1.0)*x + 1.0"
    @fact affToStr(2x) --> "(2.0)*x + 0.0"
    # Number--RandomAffExpr
    r = 2x+1
    @fact affToStr(1+r) --> "(2.0)*x + 2.0"
    @fact affToStr(1-r) --> "(-2.0)*x + 0.0"
    @fact affToStr(2r) --> "(4.0)*x + 2.0"
    # Number--CCAffExpr
    r2 = (2v+1)*x+1
    @fact affToStr(1+r2) --> "(2 v + 1)*x + 2"
    @fact affToStr(1-r2) --> "(-2 v - 1)*x + 0"
    @fact affToStr(2r2) --> "(4 v + 2)*x + 2"

    # Variable
    # Variable--IndepNormal
    @fact affToStr(v+x) --> "(1)*x + v"
    @fact affToStr(v-x) --> "(-1)*x + v"
    @fact affToStr(v*x) --> "(v)*x + 0"
    # Variable--RandomAffExpr
    @fact affToStr(v+r) --> "(2)*x + v + 1"
    @fact affToStr(v-r) --> "(-2)*x + v - 1"
    @fact affToStr(v*r) --> "(2 v)*x + v"
    # Variable--CCAffExpr
    @fact affToStr(v+r2) --> "(2 v + 1)*x + v + 1"
    @fact affToStr(v-r2) --> "(-2 v - 1)*x + v - 1"
    # Variable*CCAffExpr not valid

    # AffExpr
    # AffExpr--IndepNormal
    a = 2q+3
    @fact affToStr(a+x) --> "(1)*x + 2 q + 3"
    @fact affToStr(a-x) --> "(-1)*x + 2 q + 3"
    @fact affToStr(a*x) --> "(2 q + 3)*x + 0"
    # AffExpr--RandomAffExpr
    @fact affToStr(a+r) --> "(2)*x + 2 q + 4"
    @fact affToStr(a-r) --> "(-2)*x + 2 q + 2"
    @fact affToStr(a*r) --> "(4 q + 6)*x + 2 q + 3"
    # AffExpr--CCAffExpr
    @fact affToStr(a+r2) --> "(2 v + 1)*x + 2 q + 4"
    @fact affToStr(a-r2) --> "(-2 v - 1)*x + 2 q + 2"
    # AffExpr*CCAffExpr not valid

    # IndepNormal
    @fact affToStr(-x) --> "(-1.0)*x + 0.0"
    # IndepNormal--Number
    @fact affToStr(x+1) --> "(1.0)*x + 1.0"
    @fact affToStr(x-1) --> "(1.0)*x + -1.0"
    @fact affToStr(x*3) --> "(3.0)*x + 0.0"
    @fact affToStr(x/2) --> "(0.5)*x + 0.0"
    # IndepNormal--Variable
    @fact affToStr(x+v) --> "(1)*x + v"
    @fact affToStr(x-v) --> "(1)*x + -v"
    @fact affToStr(x*v) --> "(v)*x + 0"
    # IndepNormal--AffExpr
    @fact affToStr(x+a) --> "(1)*x + 2 q + 3"
    @fact affToStr(x-a) --> "(1)*x + -2 q - 3"
    @fact affToStr(x*a) --> "(2 q + 3)*x + 0"
    # IndepNormal--IndepNormal
    @fact affToStr(y+x) --> "(1.0)*y + (1.0)*x + 0.0"
    @fact affToStr(y-x) --> "(1.0)*y + (-1.0)*x + 0.0"
    # IndepNormal*IndepNormal not valid
    # IndepNormal--RandomAffExpr
    @fact affToStr(x+r) --> "(3.0)*x + 1.0"
    @fact affToStr(x-r) --> "(-1.0)*x + -1.0"
    # IndepNormal*RandomAffExpr not valid
    # IndepNormal--CCAffExpr
    @fact affToStr(x+r2) --> "(2 v + 2)*x + 1"
    @fact affToStr(x-r2) --> "(-2 v)*x + -1"
    # IndepNormal*CCAffExpr not valid

    # RandomAffExpr
    # RandomAffExpr--Number
    @fact affToStr(r+1) --> "(2.0)*x + 2.0"
    @fact affToStr(r-1) --> "(2.0)*x + 0.0"
    @fact affToStr(r*2) --> "(4.0)*x + 2.0"
    @fact affToStr(r/2) --> "(1.0)*x + 0.5"
    # RandomAffExpr--Variable
    @fact affToStr(r+v) --> "(2)*x + v + 1"
    @fact affToStr(r-v) --> "(2)*x + -v + 1"
    @fact affToStr(r*v) --> "(2 v)*x + v"
    # RandomAffExpr--AffExpr
    @fact affToStr(r+a) --> "(2)*x + 2 q + 4"
    @fact affToStr(r-a) --> "(2)*x + -2 q - 2"
    @fact affToStr(r*a) --> "(4 q + 6)*x + 2 q + 3"
    # RandomAffExpr--IndepNormal
    @fact affToStr(r+y) --> "(2.0)*x + (1.0)*y + 1.0"
    @fact affToStr(r-y) --> "(2.0)*x + (-1.0)*y + 1.0"
    # RandomAffExpr*IndepNormal not valid
    # RandomAffExpr--RandomAffExpr
    @fact affToStr(r+(3y-1)) --> "(2.0)*x + (3.0)*y + 0.0"
    @fact affToStr(r-(3y-1)) --> "(2.0)*x + (-3.0)*y + 2.0"
    # RandomAffExpr*RandomAffExpr not valid
    # RandomAffExpr--CCAffExpr
    @fact affToStr(r+r2) --> "(2 v + 3)*x + 2"
    @fact affToStr(r-r2) --> "(-2 v + 1)*x + 0"
    # RandomAffExpr*CCAffExpr not valid

    # CCAffExpr
    # CCAffExpr--Number
    @fact affToStr(r2+1) --> "(2 v + 1)*x + 2"
    @fact affToStr(r2-1) --> "(2 v + 1)*x + 0"
    @fact affToStr(r2*2) --> "(4 v + 2)*x + 2"
    @fact affToStr(r2/2) --> "(v + 0.5)*x + 0.5"
    # CCAffExpr--Variable
    @fact affToStr(r2+v) --> "(2 v + 1)*x + v + 1"
    @fact affToStr(r2-v) --> "(2 v + 1)*x + -v + 1"
    # CCAffExpr*Variable not valid
    # CCAffExpr--AffExpr
    @fact affToStr(r2+a) --> "(2 v + 1)*x + 2 q + 4"
    @fact affToStr(r2-a) --> "(2 v + 1)*x + -2 q - 2"
    # CCAffExpr*AffExpr not valid
    # CCAffExpr--IndepNormal
    @fact affToStr(r2+y) --> "(2 v + 1)*x + (1)*y + 1"
    @fact affToStr(r2-y) --> "(2 v + 1)*x + (-1)*y + 1"
    # CCAffExpr*IndepNormal not valid
    # CCAffExpr--RandomAffExpr
    @fact affToStr(r2+r) --> "(2 v + 3)*x + 2"
    @fact affToStr(r2-r) --> "(2 v - 1)*x + 0"
    # CCAffExpr*RandomAffExpr not valid
    # CCAffExpr--CCAffExpr
    @fact affToStr(r2+r2) --> "(4 v + 2)*x + 2"
    @fact affToStr(r2-r2) --> "(0)*x + 0"
    # CCAffExpr*CCAffExpr not valid


    @fact getMean(z[5]) --> 5
    @fact getVariance(z[5]) --> 1

    @addConstraint(m, (3v+1)*x + 10 <= 20, with_probability=0.05)
    @fact conToStr(JuMPChance.getCCData(m).chanceconstr[1]) --> "(3 v + 1)*x + -10 >= 0, with probability 0.95"
    @addConstraint(m, (3v+1)*x + 10 <= 20, with_probability=0.95)
    @fact conToStr(JuMPChance.getCCData(m).chanceconstr[end]) --> "(3 v + 1)*x + -10 <= 0, with probability 0.95"
    @fact_throws ErrorException @addConstraint(m, (3v+1)*x + 10 >= 20, with_probability=10.0)
    @addConstraint(m, (3v+1)*x + 10 >= 20, with_probability=0.96)
    @fact conToStr(JuMPChance.getCCData(m).chanceconstr[end]) --> "(3 v + 1)*x + -10 >= 0, with probability 0.96"

    ccaff = JuMPChance.CCAffExpr()
    ccaff.constant = 2v
    @fact affToStr(ccaff) --> "2 v"

    raff = JuMPChance.RandomAffExpr()
    raff.constant = 10
    @fact affToStr(raff) --> "10.0"


    @fact_throws ErrorException @defIndepNormal(m, q, mean=1, var=-1)
    @fact_throws ErrorException @defIndepNormal(m, q, mean=1, var=(-1,1))
    @fact startswith(macroexpand(:(@defIndepNormal(m, f(x), mean=1, var=1))).args[1].msg,"Syntax error: Expected") --> true

    @fact affToStr(z[1]+z[2]-2z[3]+10) --> "(1.0)*z[1] + (1.0)*z[2] + (-2.0)*z[3] + 10.0"

    @fact affToStr(v*(x-1)) --> "(v)*x + -v"
    @fact affToStr(x*(v-1)) --> "(v - 1)*x + 0"

    jm = Model()
    @fact_throws ErrorException JuMPChance.getCCData(jm)

end

facts("Printing two-sided constraints") do
    m = ChanceModel()
    @defVar(m, x)
    @defIndepNormal(m, ξ, mean = 0, var = 1)
    cc = JuMPChance.getCCData(m)

    @addConstraint(m, -1 ≤ x*ξ ≤ 1, with_probability = 0.95, approx="1.25")
    @fact conToStr(cc.twosidechanceconstr[end]) --> "-1 <= (x)*ξ + 0 <= 1, with probability 0.95"
    @addConstraint(m, -1 ≤ x*ξ ≤ x, with_probability = 0.95, approx="1.25")
    @fact conToStr(cc.twosidechanceconstr[end]) --> "-1 <= (x)*ξ + 0 <= x, with probability 0.95"
    if VERSION > v"0.4-"
        @defVar(m,y[1:3])
        @addConstraint(m, -1 ≤ x*ξ ≤ sum{y[i],i=1:3}, with_probability = 0.95, approx="1.25")
        @fact conToStr(cc.twosidechanceconstr[end]) --> "-1 <= (x)*ξ + 0 <= y[1] + y[2] + y[3], with probability 0.95"
    end
end

facts("Basic chance constraint model") do
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        cref = @addConstraint(m, z*x >= -1, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/quantile(Normal(0,1),0.95),1e-6)
        @fact JuMPChance.satisfied_with_probability(cref) > 0.95 - 1e-4 --> true
    end
end

facts("Deprecated comparison operator overloads") do
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        addConstraint(m, z*x >= -1, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/quantile(Normal(0,1),0.95),1e-6)
    end
end

facts("Flipped constraint sense") do
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        cref = @addConstraint(m, -z*x <= 1, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/quantile(Normal(0,1),0.95),1e-6)
        @fact JuMPChance.satisfied_with_probability(cref) > 0.95 - 1e-4 --> true
    end
end

facts("Invariance to transformations") do
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=1, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        @addConstraint(m, z*(x-1) >= -1, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/quantile(Normal(0,1),0.95),1e-6)
    end
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=2, var=4)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        @addConstraint(m, z*(x/2-1) >= -1, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/quantile(Normal(0,1),0.95),1e-6)
    end
end

facts("Duplicate terms") do
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        @addConstraint(m, (1/2)z*x + (1/2)z*x >= -1, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/quantile(Normal(0,1),0.95),1e-6)
    end
end

# addConstraints
#=
let
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        @addConstraints m begin
            z >= -100
            z*x >= -1, with_probability->0.95
        end

        status = solve(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
    end
end
=#

facts("Robust but no uncertainty budget") do
    let
        m = ChanceModel()
        @defIndepNormal(m, x, mean=(-1,1),var=1)

        @defVar(m, z >= -100)
        @setObjective(m, Min, z)

        @addConstraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
        status = solve(m, method=:Cuts, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/quantile(Normal(0,1),0.95),1e-6)
    end
    let
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0,var=(0.95,1.05))

        @defVar(m, z >= -100)
        @setObjective(m, Min, z)

        @addConstraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
        status = solve(m, method=:Cuts, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/quantile(Normal(0,1),0.95),1e-6)
    end
    let # flipped signs
        m = ChanceModel()
        @defIndepNormal(m, x, mean=(-1,1),var=1)

        @defVar(m, z >= -100)
        @setObjective(m, Min, z)

        @addConstraint(m, -z*x <= 1, with_probability=0.95, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
        status = solve(m, method=:Cuts,silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/quantile(Normal(0,1),0.95),1e-6)
    end
end



# quadratic objective
facts("Quadratic objective") do
    for method in [:Reformulate, :Cuts], linearize in [true, false]
        (method == :Reformulate && linearize) && continue
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0,var=1)

        @defVar(m, z >= -100)
        @setObjective(m, Min, z+2z^2)

        @addConstraint(m, z*x >= -1, with_probability=0.95)
        status = solve(m, method=method, linearize_objective=linearize, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/4,1e-4)
    end
    for linearize in [true, false]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=(-1,1),var=1)

        @defVar(m, z >= -100)
        @setObjective(m, Min, z+2z^2)

        @addConstraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
        status = solve(m, method=:Cuts, linearize_objective=linearize, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/4,1e-4)
    end
end

facts("Uncertainty budget for mean") do
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(-1,1),var=1)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    cref = @addConstraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=0)
    status = solve(m, method=:Cuts, silent=true)
    @fact status --> :Optimal
    @fact getValue(z) --> roughly(-1/(1+quantile(Normal(0,1),0.95)),1e-6)
    @fact_throws ErrorException JuMPChance.satisfied_with_probability(cref)
end

facts("Uncertain variance, but no budget") do
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    @addConstraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=0)
    status = solve(m, method=:Cuts, silent=true)
    @fact status --> :Optimal
    @fact getValue(z) --> roughly(-1/(1+quantile(Normal(0,1),0.95)),1e-6)
end

facts("Uncertain variance, with budget") do
    let
        m = ChanceModel()
        @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))

        @defVar(m, z >= -100)
        @setObjective(m, Min, z)

        @addConstraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
        status = solve(m, method=:Cuts,silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)),1e-6)
    end
    let # shifted mean
        m = ChanceModel()
        @defIndepNormal(m, x, mean=(0,2),var=(0.95,1.05))

        @defVar(m, z >= -100)
        @setObjective(m, Min, z)

        @addConstraint(m, z*(x-1) >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
        status = solve(m, method=:Cuts,silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-1/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)),1e-6)
    end
    let # rescaled variable
        m = ChanceModel()
        @defIndepNormal(m, x, mean=(0,2),var=(0.95,1.05))

        @defVar(m, z >= -100)
        @setObjective(m, Min, z)

        @addConstraint(m, (z/2)*(x-1) >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
        status = solve(m, method=:Cuts,silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(-2/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)),1e-6)
    end
    let # more than one R.V.
        m = ChanceModel()
        @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))
        @defIndepNormal(m, y, mean=(-0.01,0.01),var=0.01)

        @defVar(m, z >= -100)
        @setObjective(m, Min, z)

        @addConstraint(m, z*x + y >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
        status = solve(m, method=:Cuts,silent=true)
        @fact status --> :Optimal
        # In mathematica: Minimize[{z, z - \[Nu]*Sqrt[1.05*z^2 + 0.01] >= -1}, z]
        @fact getValue(z) --> roughly(-0.36431227017642165,1e-5)
    end
end

facts("Integer variables") do
    m = ChanceModel(solver=GLPKSolverMIP())
    @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))
    @defIndepNormal(m, y, mean=(-0.01,0.01),var=0.01)

    @defVar(m, z >= -100, Int)
    @setObjective(m, Min, z)

    @addConstraint(m, z*x + sum{y,i=1:1} >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solve(m, method=:Cuts,silent=true)
    @fact status --> :Optimal
    @fact getValue(z) --> roughly(0.0,1e-5)
end

# variance == 0 corner case
facts("Variance == 0 corner case") do
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z <= 100) # so original problem is bounded
        @defVar(m, y == 0)

        @setObjective(m, Max, z)
        @addConstraint(m, z + x*y <= 20, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(20,1e-6)
    end
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z <= 100)
        @defVar(m, y == 0)

        @setObjective(m, Max, z)
        @addConstraint(m, -z - x*y >= -20, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(20,1e-6)
    end
    let
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z <= 100)
        @defVar(m, y == 0)

        @setObjective(m, Max, z)
        @addConstraint(m, -z - x*y >= -20, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)

        status = solve(m, method=:Cuts, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(20,1e-6)
    end
    let
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z <= 100)
        @defVar(m, y == 0)

        @setObjective(m, Max, z)
        @addConstraint(m, z + x*y <= 20, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)

        status = solve(m, method=:Cuts, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(20,1e-6)
    end
end


facts("Special cases where chance constraint becomes linear") do
    for method in [:Reformulate, :Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x[1:2], mean=1, var=1)
        @defVar(m, z <= 100)
        @defVar(m, y[1:2])

        @setObjective(m, Max, z)
        @addConstraint(m, y[1] == y[2])
        @addConstraint(m, z + sum{x[i]*y[i], i=1:2} <= 20, with_probability=0.95)


        status = solve(m, method=method, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(20,1e-5)
    end
    for method in [:Reformulate, :Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x[1:2], mean=1, var=1)
        @defVar(m, z >= 100)
        @defVar(m, y[1:2])

        @setObjective(m, Min, z)
        @addConstraint(m, y[1] == y[2])
        @addConstraint(m, z + sum{x[i]*y[i], i=1:2} >= 200, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @fact status --> :Optimal
        @fact getValue(z) --> roughly(200,1e-5)
    end
end

facts("Compare with manual reformulation") do
    for method in [:Reformulate, :Cuts], ϵ in [0.001, 0.01, 0.005, 0.05], σ² in [0, 1]
        m = ChanceModel()
        @defIndepNormal(m, ω₁, mean=0, var=σ²)
        @defIndepNormal(m, ω₂, mean=1, var=σ²)
        @defVar(m, z)
        @defVar(m, x[1:2])
        @addConstraint(m, x[1] >= 3)
        @addConstraint(m, x[2] >= 2)
        @addConstraint(m, x[1]*ω₁ + x[2]*ω₂ <= z, with_probability=1-ϵ)
        @setObjective(m, Min, z)
        
        status = solve(m, probability_tolerance=1e-12, method=method, silent=true)

        objval = getObjectiveValue(m)
        @fact status --> :Optimal

        m = Model(solver=ECOS.ECOSSolver(verbose=0))
        nu = quantile(Normal(0,1),1-ϵ)
        @defVar(m, z)
        @defVar(m, x[1:2])
        @defVar(m, varterm[1:2] >= 0)
        @defVar(m, t >= 0)
        @addConstraint(m, x[1] >= 3)
        @addConstraint(m, x[2] >= 2)
        @addConstraint(m, defvar[i=1:2], varterm[i] == sqrt(σ²)*x[i])
        @addConstraint(m, x[2] + nu*t <= z)
        @addConstraint(m, sum{ varterm[i]^2, i=1:2} <= t^2)
        @setObjective(m, Min, z)
        
        status = solve(m)
        @fact status --> :Optimal
        @fact getObjectiveValue(m) --> roughly(objval,1e-5)

    end
end


# test for treating a two-sided constraints as two one sided constraints
facts("Two-sided constraints vs. two one-sided constraints") do 
    # (a) one random variable model, with mean = 0 and var = 1
    for method in [:Reformulate, :Cuts]
        ϵ = 0.01
        m = ChanceModel()
        @defIndepNormal(m, ω, mean=0, var=1)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x == 1)
        @addConstraint(m, l <= x*ω <= u, with_probability=1-ϵ, approx="2.0")
        @setObjective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        objval = getObjectiveValue(m)
        
        m = ChanceModel()
        @defIndepNormal(m, ω, mean=0, var=1)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x == 1)
        @addConstraint(m, x*ω >= l, with_probability=1-ϵ)
        @addConstraint(m, x*ω <= u, with_probability=1-ϵ)
        @setObjective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        @fact getObjectiveValue(m) --> roughly(objval,1e-5)
    end
    
    # (b) one random variable model, with mean = 1 and var = 1
    for method in [:Reformulate, :Cuts]
        ϵ = 0.005
        m = ChanceModel()
        @defIndepNormal(m, ω, mean=1, var=1)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x == 1)
        @addConstraint(m, l <= x*ω <= u, with_probability=1-ϵ, approx="2.0")
        @setObjective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        objval = getObjectiveValue(m)
        
        m = ChanceModel()
        @defIndepNormal(m, ω, mean=1, var=1)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x == 1)
        @addConstraint(m, x*ω >= l, with_probability=1-ϵ)
        @addConstraint(m, x*ω <= u, with_probability=1-ϵ)
        @setObjective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        @fact getObjectiveValue(m) --> roughly(objval,1e-5)
    end

    # (c) multiple random variable model (μ=1, σ²=1), with same coefficients - linear case
    for method in [:Reformulate, :Cuts]
        ϵ = 0.005
        m = ChanceModel()
        @defIndepNormal(m, ω[1:4], mean=1, var=1)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x[1:4])
        @addConstraint(m, xcons[i=1:4], x[i] == 1)
        @addConstraint(m, l <= sum{x[i]*ω[i], i=1:4} <= u, with_probability=1-ϵ, approx="2.0")
        @setObjective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        objval = getObjectiveValue(m)
        
        m = ChanceModel()
        @defIndepNormal(m, ω[1:4], mean=1, var=1)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x[1:4])
        @addConstraint(m, xcons[i=1:4], x[i] == 1)
        @addConstraint(m, sum{x[i]*ω[i], i=1:4}>= l, with_probability=1-ϵ)
        @addConstraint(m, sum{x[i]*ω[i], i=1:4} <= u, with_probability=1-ϵ)
        @setObjective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        @fact getObjectiveValue(m) --> roughly(objval,1e-5)
    end

    # (d) multiple random variable model (μ=1, σ²=1), with different coefficients
    for method in [:Reformulate, :Cuts]
        ϵ = 0.005
        m = ChanceModel()
        @defIndepNormal(m, ω[1:4], mean=1, var=1)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x[1:4])
        @addConstraint(m, xcons[i=1:4], x[i] == i)
        @addConstraint(m, l <= sum{i*x[i]*ω[i], i=1:4} <= u, with_probability=1-ϵ, approx="2.0")
        @setObjective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        objval = getObjectiveValue(m)
        
        m = ChanceModel()
        @defIndepNormal(m, ω[1:4], mean=1, var=1)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x[1:4])
        @addConstraint(m, xcons[i=1:4], x[i] == i)
        @addConstraint(m, sum{i*x[i]*ω[i], i=1:4}>= l, with_probability=1-ϵ)
        @addConstraint(m, sum{i*x[i]*ω[i], i=1:4} <= u, with_probability=1-ϵ)
        @setObjective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        @fact getObjectiveValue(m) --> roughly(objval,1e-5)
    end

    # (e) multiple random variable model (μ=1, σ²=0), with different coefficients
    for method in [:Reformulate, :Cuts]
        ϵ = 0.005
        m = ChanceModel()
        @defIndepNormal(m, ω[1:4], mean=1, var=0)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x[1:4])
        @addConstraint(m, xcons[i=1:4], x[i] == i)
        @addConstraint(m, l <= sum{i*x[i]*ω[i], i=1:4} <= u, with_probability=1-ϵ, approx="2.0")
        @setObjective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        objval = getObjectiveValue(m)
        
        m = ChanceModel()
        @defIndepNormal(m, ω[1:4], mean=1, var=0)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x[1:4])
        @addConstraint(m, xcons[i=1:4], x[i] == i)
        @addConstraint(m, sum{i*x[i]*ω[i], i=1:4}>= l, with_probability=1-ϵ)
        @addConstraint(m, sum{i*x[i]*ω[i], i=1:4} <= u, with_probability=1-ϵ)
        @setObjective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        @fact getObjectiveValue(m) --> roughly(objval,1e-5)
    end
end


const prob_guarantee = 1.25
facts("Basic two-sided constraints") do
    for ϵ in (0.1, 0.05, 0.005, 0.0005), method in (:Cuts, :Reformulate)
        m = ChanceModel()
        @defIndepNormal(m, ξ, mean=0, var=1)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x == 1)
        cref = @addConstraint(m, l ≤ x*ξ ≤ u, with_probability=1-ϵ, approx="1.25")
        @setObjective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        violation = 1- JuMPChance.satisfied_with_probability(cref)
        @fact violation ≤ prob_guarantee*ϵ + 1e-5 --> true
        @fact violation ≥ ϵ - 1e-5  --> true # should be tight

        lval = getValue(l)
        uval = getValue(u)
        m = ChanceModel()
        @defIndepNormal(m, ξ, mean=1, var=1) # translate with mean
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x == 1)
        cref = @addConstraint(m, l ≤ x*ξ ≤ u, with_probability=1-ϵ, approx="1.25")
        @setObjective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        violation = 1- JuMPChance.satisfied_with_probability(cref)
        @fact violation ≤ prob_guarantee*ϵ + 1e-5 --> true
        @fact violation ≥ ϵ - 1e-5 --> true

        @fact getObjectiveValue(m) --> roughly((uval+1) - 2(lval+1))
    end

    # duplicates
    for ϵ in (0.1, 0.05, 0.005, 0.0005), method in (:Cuts, :Reformulate)
        m = ChanceModel()
        @defIndepNormal(m, ξ, mean=0, var=1)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x == 1)
        cref = @addConstraint(m, l ≤ x*ξ + x*ξ ≤ u, with_probability=1-ϵ, approx="1.25")
        @setObjective(m, Min, u-l)

        solve(m, method=method, silent=true)

        violation = 1- JuMPChance.satisfied_with_probability(cref)
        @fact violation ≤ prob_guarantee*ϵ + 1e-5 --> true
        @fact violation ≥ ϵ - 1e-5 --> true
    end

    # constant bounds
    let
        ϵ = 0.05
        m = ChanceModel()
        @defIndepNormal(m, ξ, mean=0, var=1)
        @defVar(m, x >= 0)
        cref = @addConstraint(m, -1 ≤ x*ξ ≤ 1, with_probability=1-ϵ, approx="1.25")
        @setObjective(m, Max, x)
        solve(m, method=:Reformulate)
        violation = 1- JuMPChance.satisfied_with_probability(cref)
        @fact violation ≤ prob_guarantee*ϵ + 1e-5 --> true
        @fact violation ≥ ϵ - 1e-5 --> true
    end

    # absolute value, random objective
    let
        ϵ = 0.05
        srand(998)
        for q in 1:100
            m = ChanceModel()
            @defIndepNormal(m, ξ, mean=0, var=1)
            @defVar(m, -10 ≤ x ≤ 10)
            @defVar(m, 0 ≤ c ≤ 10)
            @defVar(m, t ≤ 100)
            cref = @addConstraint(m, -t ≤ x*ξ + c ≤ t, with_probability=1-ϵ, approx="1.25")
            @setObjective(m, Min, (rand()-0.5)*x + (rand()-0.5)*c + 0.01*t)
            solve(m, method=:Reformulate)
            violation = 1- JuMPChance.satisfied_with_probability(cref)
            @fact violation ≤ prob_guarantee*ϵ + 1e-5 --> true
        end
    end
end

FactCheck.exitstatus()
