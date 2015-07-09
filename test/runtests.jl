using JuMPChance, JuMP
using Base.Test
using Distributions
using GLPKMathProgInterface


let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=1, var=1)
    @defIndepNormal(m, y, mean=1, var=1)
    @defVar(m, v)
    @defVar(m, q)
    @defIndepNormal(m,z[i=1:10],mean=i,var=1)

    # Number
    # Number--IndepNormal
    @test affToStr(1+x) == "(1.0)*x + 1.0"
    @test affToStr(1-x) == "(-1.0)*x + 1.0"
    @test affToStr(2x) == "(2.0)*x + 0.0"
    # Number--RandomAffExpr
    r = 2x+1
    @test affToStr(1+r) == "(2.0)*x + 2.0"
    @test affToStr(1-r) == "(-2.0)*x + 0.0"
    @test affToStr(2r) == "(4.0)*x + 2.0"
    # Number--CCAffExpr
    r2 = (2v+1)*x+1
    @test affToStr(1+r2) == "(2 v + 1)*x + 2"
    @test affToStr(1-r2) == "(-2 v - 1)*x + 0"
    @test affToStr(2r2) == "(4 v + 2)*x + 2"

    # Variable
    # Variable--IndepNormal
    @test affToStr(v+x) == "(1)*x + v"
    @test affToStr(v-x) == "(-1)*x + v"
    @test affToStr(v*x) == "(v)*x + 0"
    # Variable--RandomAffExpr
    @test affToStr(v+r) == "(2)*x + v + 1"
    @test affToStr(v-r) == "(-2)*x + v - 1"
    @test affToStr(v*r) == "(2 v)*x + v"
    # Variable--CCAffExpr
    @test affToStr(v+r2) == "(2 v + 1)*x + v + 1"
    @test affToStr(v-r2) == "(-2 v - 1)*x + v - 1"
    # Variable*CCAffExpr not valid

    # AffExpr
    # AffExpr--IndepNormal
    a = 2q+3
    @test affToStr(a+x) == "(1)*x + 2 q + 3"
    @test affToStr(a-x) == "(-1)*x + 2 q + 3"
    @test affToStr(a*x) == "(2 q + 3)*x + 0"
    # AffExpr--RandomAffExpr
    @test affToStr(a+r) == "(2)*x + 2 q + 4"
    @test affToStr(a-r) == "(-2)*x + 2 q + 2"
    @test affToStr(a*r) == "(4 q + 6)*x + 2 q + 3"
    # AffExpr--CCAffExpr
    @test affToStr(a+r2) == "(2 v + 1)*x + 2 q + 4"
    @test affToStr(a-r2) == "(-2 v - 1)*x + 2 q + 2"
    # AffExpr*CCAffExpr not valid

    # IndepNormal
    @test affToStr(-x) == "(-1.0)*x + 0.0"
    # IndepNormal--Number
    @test affToStr(x+1) == "(1.0)*x + 1.0"
    @test affToStr(x-1) == "(1.0)*x + -1.0"
    @test affToStr(x*3) == "(3.0)*x + 0.0"
    @test affToStr(x/2) == "(0.5)*x + 0.0"
    # IndepNormal--Variable
    @test affToStr(x+v) == "(1)*x + v"
    @test affToStr(x-v) == "(1)*x + -v"
    @test affToStr(x*v) == "(v)*x + 0"
    # IndepNormal--AffExpr
    @test affToStr(x+a) == "(1)*x + 2 q + 3"
    @test affToStr(x-a) == "(1)*x + -2 q - 3"
    @test affToStr(x*a) == "(2 q + 3)*x + 0"
    # IndepNormal--IndepNormal
    @test affToStr(y+x) == "(1.0)*y + (1.0)*x + 0.0"
    @test affToStr(y-x) == "(1.0)*y + (-1.0)*x + 0.0"
    # IndepNormal*IndepNormal not valid
    # IndepNormal--RandomAffExpr
    @test affToStr(x+r) == "(3.0)*x + 1.0"
    @test affToStr(x-r) == "(-1.0)*x + -1.0"
    # IndepNormal*RandomAffExpr not valid
    # IndepNormal--CCAffExpr
    @test affToStr(x+r2) == "(2 v + 2)*x + 1"
    @test affToStr(x-r2) == "(-2 v)*x + -1"
    # IndepNormal*CCAffExpr not valid

    # RandomAffExpr
    # RandomAffExpr--Number
    @test affToStr(r+1) == "(2.0)*x + 2.0"
    @test affToStr(r-1) == "(2.0)*x + 0.0"
    @test affToStr(r*2) == "(4.0)*x + 2.0"
    @test affToStr(r/2) == "(1.0)*x + 0.5"
    # RandomAffExpr--Variable
    @test affToStr(r+v) == "(2)*x + v + 1"
    @test affToStr(r-v) == "(2)*x + -v + 1"
    @test affToStr(r*v) == "(2 v)*x + v"
    # RandomAffExpr--AffExpr
    @test affToStr(r+a) == "(2)*x + 2 q + 4"
    @test affToStr(r-a) == "(2)*x + -2 q - 2"
    @test affToStr(r*a) == "(4 q + 6)*x + 2 q + 3"
    # RandomAffExpr--IndepNormal
    @test affToStr(r+y) == "(2.0)*x + (1.0)*y + 1.0"
    @test affToStr(r-y) == "(2.0)*x + (-1.0)*y + 1.0"
    # RandomAffExpr*IndepNormal not valid
    # RandomAffExpr--RandomAffExpr
    @test affToStr(r+(3y-1)) == "(2.0)*x + (3.0)*y + 0.0"
    @test affToStr(r-(3y-1)) == "(2.0)*x + (-3.0)*y + 2.0"
    # RandomAffExpr*RandomAffExpr not valid
    # RandomAffExpr--CCAffExpr
    @test affToStr(r+r2) == "(2 v + 3)*x + 2"
    @test affToStr(r-r2) == "(-2 v + 1)*x + 0"
    # RandomAffExpr*CCAffExpr not valid

    # CCAffExpr
    # CCAffExpr--Number
    @test affToStr(r2+1) == "(2 v + 1)*x + 2"
    @test affToStr(r2-1) == "(2 v + 1)*x + 0"
    @test affToStr(r2*2) == "(4 v + 2)*x + 2"
    @test affToStr(r2/2) == "(v + 0.5)*x + 0.5"
    # CCAffExpr--Variable
    @test affToStr(r2+v) == "(2 v + 1)*x + v + 1"
    @test affToStr(r2-v) == "(2 v + 1)*x + -v + 1"
    # CCAffExpr*Variable not valid
    # CCAffExpr--AffExpr
    @test affToStr(r2+a) == "(2 v + 1)*x + 2 q + 4"
    @test affToStr(r2-a) == "(2 v + 1)*x + -2 q - 2"
    # CCAffExpr*AffExpr not valid
    # CCAffExpr--IndepNormal
    @test affToStr(r2+y) == "(2 v + 1)*x + (1)*y + 1"
    @test affToStr(r2-y) == "(2 v + 1)*x + (-1)*y + 1"
    # CCAffExpr*IndepNormal not valid
    # CCAffExpr--RandomAffExpr
    @test affToStr(r2+r) == "(2 v + 3)*x + 2"
    @test affToStr(r2-r) == "(2 v - 1)*x + 0"
    # CCAffExpr*RandomAffExpr not valid
    # CCAffExpr--CCAffExpr
    @test affToStr(r2+r2) == "(4 v + 2)*x + 2"
    @test affToStr(r2-r2) == "(0)*x + 0"
    # CCAffExpr*CCAffExpr not valid


    @test getMean(z[5]) == 5
    @test getVariance(z[5]) == 1

    @addConstraint(m, (3v+1)*x + 10 <= 20, with_probability=0.05)
    @test conToStr(JuMPChance.getCCData(m).chanceconstr[1]) == "(3 v + 1)*x + -10 >= 0, with probability 0.95"
    @addConstraint(m, (3v+1)*x + 10 <= 20, with_probability=0.95)
    @test conToStr(JuMPChance.getCCData(m).chanceconstr[end]) == "(3 v + 1)*x + -10 <= 0, with probability 0.95"
    @test_throws ErrorException @addConstraint(m, (3v+1)*x + 10 >= 20, with_probability=10.0)
    @addConstraint(m, (3v+1)*x + 10 >= 20, with_probability=0.96)
    @test conToStr(JuMPChance.getCCData(m).chanceconstr[end]) == "(3 v + 1)*x + -10 >= 0, with probability 0.96"

    ccaff = JuMPChance.CCAffExpr()
    ccaff.constant = 2v
    @test affToStr(ccaff) == "2 v"

    raff = JuMPChance.RandomAffExpr()
    raff.constant = 10
    @test affToStr(raff) == "10.0"


    @test_throws ErrorException @defIndepNormal(m, q, mean=1, var=-1)
    @test_throws ErrorException @defIndepNormal(m, q, mean=1, var=(-1,1))
    @test beginswith(macroexpand(:(@defIndepNormal(m, f(x), mean=1, var=1))).args[1].msg,"Syntax error: Expected")

    @test affToStr(z[1]+z[2]-2z[3]+10) == "(1.0)*z[1] + (1.0)*z[2] + (-2.0)*z[3] + 10.0"

    @test affToStr(v*(x-1)) == "(v)*x + -v"
    @test affToStr(x*(v-1)) == "(v - 1)*x + 0"

    jm = Model()
    @test_throws ErrorException JuMPChance.getCCData(jm)

end

let
    m = ChanceModel()
    @defVar(m, x)
    @defIndepNormal(m, ξ, mean = 0, var = 1)
    cc = JuMPChance.getCCData(m)

    @addConstraint(m, -1 ≤ x*ξ ≤ 1, with_probability = 0.95)
    @test conToStr(cc.twosidechanceconstr[end]) == "-1 <= (x)*ξ + 0 <= 1, with probability 0.95"
    @addConstraint(m, -1 ≤ x*ξ ≤ x, with_probability = 0.95)
    @test conToStr(cc.twosidechanceconstr[end]) == "-1 <= (x)*ξ + 0 <= x, with probability 0.95"
    if VERSION > v"0.4-"
        @defVar(m,y[1:3])
        @addConstraint(m, -1 ≤ x*ξ ≤ sum{y[i],i=1:3}, with_probability = 0.95)
        @test conToStr(cc.twosidechanceconstr[end]) == "-1 <= (x)*ξ + 0 <= y[1] + y[2] + y[3], with probability 0.95"
    end
end

let
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        cref = @addConstraint(m, z*x >= -1, with_probability=0.95)

        status = solve(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
        @test JuMPChance.satisfied_with_probability(cref) > 0.95 - 1e-4
    end
end

# test deprecated operator overloads
let
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        addConstraint(m, z*x >= -1, with_probability=0.95)

        status = solve(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
    end
end

# flipped constraint sense
let
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        cref = @addConstraint(m, -z*x <= 1, with_probability=0.95)

        status = solve(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
        @test JuMPChance.satisfied_with_probability(cref) > 0.95 - 1e-4
    end
end

# invariance to transformations
let
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=1, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        @addConstraint(m, z*(x-1) >= -1, with_probability=0.95)

        status = solve(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
    end
end
let
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=2, var=4)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        @addConstraint(m, z*(x/2-1) >= -1, with_probability=0.95)

        status = solve(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
    end
end

# duplicate terms
let
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        @addConstraint(m, (1/2)z*x + (1/2)z*x >= -1, with_probability=0.95)

        status = solve(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
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

# robust but no uncertainty budget
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(-1,1),var=1)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    @addConstraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
    status = solve(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
end
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=0,var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    @addConstraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
    status = solve(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
end

# flipped signs
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(-1,1),var=1)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    @addConstraint(m, -z*x <= 1, with_probability=0.95, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
    status = solve(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
end

# quadratic objective
let
    for method in [:Reformulate, :Cuts], linearize in [true, false]
        (method == :Reformulate && linearize) && continue
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0,var=1)

        @defVar(m, z >= -100)
        @setObjective(m, Min, z+2z^2)

        @addConstraint(m, z*x >= -1, with_probability=0.95)
        status = solve(m, method=method, linearize_objective=linearize)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/4 1e-4
    end
end
let
    for linearize in [true, false]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=(-1,1),var=1)

        @defVar(m, z >= -100)
        @setObjective(m, Min, z+2z^2)

        @addConstraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
        status = solve(m, method=:Cuts, linearize_objective=linearize)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/4 1e-4
    end
end

# uncertainty budget for mean
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(-1,1),var=1)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    cref = @addConstraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=0)
    status = solve(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/(1+quantile(Normal(0,1),0.95)) 1e-6
    @test_throws ErrorException JuMPChance.satisfied_with_probability(cref)
end

# uncertain variance, but no budget
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    @addConstraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=0)
    status = solve(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/(1+quantile(Normal(0,1),0.95)) 1e-6
end

# uncertain variance, with budget
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    @addConstraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solve(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)) 1e-6
end

# shifted
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(0,2),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    @addConstraint(m, z*(x-1) >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solve(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)) 1e-6
end

# rescaled variable
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(0,2),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    @addConstraint(m, (z/2)*(x-1) >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solve(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -2/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)) 1e-6
end

# more than one R.V.
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))
    @defIndepNormal(m, y, mean=(-0.01,0.01),var=0.01)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    @addConstraint(m, z*x + y >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solve(m, method=:Cuts)
    @test status == :Optimal
    # In mathematica: Minimize[{z, z - \[Nu]*Sqrt[1.05*z^2 + 0.01] >= -1}, z]
    @test_approx_eq_eps getValue(z) -0.36431227017642165 1e-5
end

# integer variables
let
    m = ChanceModel(solver=GLPKSolverMIP())
    @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))
    @defIndepNormal(m, y, mean=(-0.01,0.01),var=0.01)

    @defVar(m, z >= -100, Int)
    @setObjective(m, Min, z)

    @addConstraint(m, z*x + sum{y,i=1:1} >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solve(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) 0.0 1e-5
end

# variance == 0 corner case
let
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z <= 100) # so original problem is bounded
        @defVar(m, y == 0)

        @setObjective(m, Max, z)
        @addConstraint(m, z + x*y <= 20, with_probability=0.95)

        status = solve(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) 20 1e-6
    end
end

let
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z <= 100)
        @defVar(m, y == 0)

        @setObjective(m, Max, z)
        @addConstraint(m, -z - x*y >= -20, with_probability=0.95)

        status = solve(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) 20 1e-6
    end
end

let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=0, var=1)
    @defVar(m, z <= 100)
    @defVar(m, y == 0)

    @setObjective(m, Max, z)
    @addConstraint(m, -z - x*y >= -20, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)

    status = solve(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) 20 1e-6
end

let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=0, var=1)
    @defVar(m, z <= 100)
    @defVar(m, y == 0)

    @setObjective(m, Max, z)
    @addConstraint(m, z + x*y <= 20, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)

    status = solve(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) 20 1e-6
end

# special cases when cc becomes linear
let 
    m = ChanceModel()
    @defIndepNormal(m, x[1:2], mean=1, var=1)
    @defVar(m, z <= 100)
    @defVar(m, y[1:2])

    @setObjective(m, Max, z)
    @addConstraint(m, y[1] == y[2])
    @addConstraint(m, z + sum{x[i]*y[i], i=1:2} <= 20, with_probability=0.95)


    status = solve(m, method=:Reformulate)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) 20 1e-6
end


let 
    m = ChanceModel()
    @defIndepNormal(m, x[1:2], mean=1, var=1)
    @defVar(m, z >= 100)
    @defVar(m, y[1:2])

    @setObjective(m, Min, z)
    @addConstraint(m, y[1] == y[2])
    @addConstraint(m, z + sum{x[i]*y[i], i=1:2} >= 200, with_probability=0.95)


    status = solve(m, method=:Reformulate)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) 200 1e-6
end

# two-sided constraints
const prob_guarantee = 1.25
let
    for ϵ in (0.1, 0.05, 0.005, 0.0005)
        m = ChanceModel()
        @defIndepNormal(m, ξ, mean=0, var=1)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x == 1)
        cref = @addConstraint(m, l ≤ x*ξ ≤ u, with_probability=1-ϵ)
        @setObjective(m, Min, u-2l)

        solve(m, method=:Reformulate)

        violation = 1- JuMPChance.satisfied_with_probability(cref)
        @test violation ≤ prob_guarantee*ϵ + 1e-5
        @test violation ≥ ϵ - 1e-5 # should be tight

        lval = getValue(l)
        uval = getValue(u)
        m = ChanceModel()
        @defIndepNormal(m, ξ, mean=1, var=1) # translate with mean
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x == 1)
        cref = @addConstraint(m, l ≤ x*ξ ≤ u, with_probability=1-ϵ)
        @setObjective(m, Min, u-2l)

        solve(m, method=:Reformulate)

        violation = 1- JuMPChance.satisfied_with_probability(cref)
        @test violation ≤ prob_guarantee*ϵ + 1e-5
        @test violation ≥ ϵ - 1e-5

        @test_approx_eq getObjectiveValue(m) ((uval+1) - 2(lval+1))
    end
end

# duplicates
let
    for ϵ in (0.1, 0.05, 0.005, 0.0005)
        m = ChanceModel()
        @defIndepNormal(m, ξ, mean=0, var=1)
        @defVar(m, l)
        @defVar(m, u)
        @defVar(m, x == 1)
        cref = @addConstraint(m, l ≤ x*ξ + x*ξ ≤ u, with_probability=1-ϵ)
        @setObjective(m, Min, u-l)

        solve(m, method=:Reformulate)

        violation = 1- JuMPChance.satisfied_with_probability(cref)
        @test violation ≤ prob_guarantee*ϵ + 1e-5
        @test violation ≥ ϵ - 1e-5
    end
end

# constant bounds
let
    ϵ = 0.05
    m = ChanceModel()
    @defIndepNormal(m, ξ, mean=0, var=1)
    @defVar(m, x >= 0)
    cref = @addConstraint(m, -1 ≤ x*ξ ≤ 1, with_probability=1-ϵ)
    @setObjective(m, Max, x)
    solve(m, method=:Reformulate)
    violation = 1- JuMPChance.satisfied_with_probability(cref)
    @test violation ≤ prob_guarantee*ϵ + 1e-5
    @test violation ≥ ϵ - 1e-5
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
        cref = @addConstraint(m, -t ≤ x*ξ + c ≤ t, with_probability=1-ϵ)
        @setObjective(m, Min, (rand()-0.5)*x + (rand()-0.5)*c + 0.01*t)
        solve(m, method=:Reformulate)
        violation = 1- JuMPChance.satisfied_with_probability(cref)
        @test violation ≤ prob_guarantee*ϵ + 1e-5
    end
end
