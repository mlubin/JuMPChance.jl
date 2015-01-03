using CCJuMP, JuMP
using Base.Test
using Distributions


let
    m = CCModel()
    @defIndepNormal(m, x, mean=1, var=1)
    @defIndepNormal(m, y, mean=1, var=1)
    @test affToStr(1+x) == "(1.0)*x + 1.0"
    @test affToStr(y+x) == "(1.0)*y + (1.0)*x + 0.0"
    @test affToStr(y-x) == "(1.0)*y + (-1.0)*x + 0.0"
    @test affToStr(x/2) == "(0.5)*x + 0.0"
    @test affToStr(-x) == "(-1.0)*x + 0.0"

    @defVar(m, v)
    @defVar(m, q)

    @test affToStr((3v+1)*x+10) == "(3 v + 1)*x + 10"

    @defIndepNormal(m,z[i=1:10],mean=i,var=1)
    @test affToStr(v*z[1]+3.5) == "(v)*z[1] + 3.5"
    
    @test affToStr(v*z[1]+2z[2]+3.5) == "(v)*z[1] + (2)*z[2] + 3.5"
    @test affToStr(v+(v*z[1]+2z[2]+3.5)) == "(v)*z[1] + (2)*z[2] + v + 3.5"
    @test affToStr((v*z[1]+2z[2]+3.5)-q) == "(v)*z[1] + (2)*z[2] + -q + 3.5"
    @test affToStr((v*z[1]+2z[2]+3.5)+q) == "(v)*z[1] + (2)*z[2] + q + 3.5"
    @test affToStr(x*v) == "(v)*x + 0"
    @test affToStr(x*10) == "(10.0)*x + 0.0"
    @test affToStr(3v+x) == "(1)*x + 3 v"
    @test affToStr(3v-x) == "(-1)*x + 3 v"
    @test affToStr(x+3v) == "(1)*x + 3 v"
    @test affToStr(x-3v) == "(1)*x + -3 v"
    @test affToStr((3x+1)+v) == "(3)*x + v + 1"
    @test affToStr((3x+1)-v) == "(3)*x + -v + 1"
    @test affToStr((3x+1)*v) == "(3 v)*x + v"

    @test getMean(z[5]) == 5
    @test getVar(z[5]) == 1

    c = (3v+1)*x + 10 <= 20
    @test conToStr(c) == "(3 v + 1)*x + -10 <= 0"
    addConstraint(m, c, with_probability=0.05)
    @test conToStr(CCJuMP.getCCData(m).chanceconstr[1]) == "(3 v + 1)*x + -10 <= 0, with probability 0.05"
    @test_throws ErrorException addConstraint(m, c, with_probability=10.0)
    c = (3v+1)*x + 10 >= 20
    @test conToStr(c) == "(3 v + 1)*x + -10 >= 0"

    ccaff = CCJuMP.CCAffExpr()
    ccaff.constant = 2v
    @test affToStr(ccaff) == "2 v"

    raff = 3x
    @test isa(raff, CCJuMP.RandomAffExpr)
    @test affToStr(raff) == "(3.0)*x + 0.0"
    raff = CCJuMP.RandomAffExpr()
    raff.constant = 10
    @test affToStr(raff) == "10.0"


    @test_throws ErrorException @defIndepNormal(m, q, mean=1, var=-1)
    @test_throws ErrorException @defIndepNormal(m, q, mean=1, var=(-1,1))
    @test beginswith(macroexpand(:(@defIndepNormal(m, f(x), mean=1, var=1))).args[1].msg,"Syntax error: Expected")

    @test affToStr(z[1]+z[2]-2z[3]+10) == "(1.0)*z[1] + (1.0)*z[2] + (-2.0)*z[3] + 10.0"

    @test affToStr(v*(x-1)) == "(v)*x + -v"
    @test affToStr(x*(v-1)) == "(v - 1)*x + 0"

    jm = Model()
    @test_throws ErrorException CCJuMP.getCCData(jm)

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
# flipped constraint sense
let
    for method in [:Reformulate,:Cuts]
        m = CCModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        addConstraint(m, -z*x >= 1, with_probability=0.05)

        status = solvecc(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
    end
end

# invariance to transformations
let
    for method in [:Reformulate,:Cuts]
        m = CCModel()
        @defIndepNormal(m, x, mean=1, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        addConstraint(m, z*(x-1) <= -1, with_probability=0.05)

        status = solvecc(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
    end
end
let
    for method in [:Reformulate,:Cuts]
        m = CCModel()
        @defIndepNormal(m, x, mean=2, var=4)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        addConstraint(m, z*(x/2-1) <= -1, with_probability=0.05)

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
let
    m = CCModel()
    @defIndepNormal(m, x, mean=0,var=(0.95,1.05))

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
    for method in [:Reformulate, :Cuts], linearize in [true, false]
        (method == :Reformulate && linearize) && continue
        m = CCModel()
        @defIndepNormal(m, x, mean=0,var=1)

        @defVar(m, z >= -100)
        @setObjective(m, Min, z+2z^2)

        addConstraint(m, z*x <= -1, with_probability=0.05)
        status = solvecc(m, method=method, linearize_objective=linearize)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/4 1e-4
    end
end
let
    for linearize in [true, false]
        m = CCModel()
        @defIndepNormal(m, x, mean=(-1,1),var=1)

        @defVar(m, z >= -100)
        @setObjective(m, Min, z+2z^2)

        addConstraint(m, z*x <= -1, with_probability=0.05, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
        status = solvecc(m, method=:Cuts, linearize_objective=linearize)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/4 1e-4
    end
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

# shifted
let
    m = CCModel()
    @defIndepNormal(m, x, mean=(0,2),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, z*(x-1) <= -1, with_probability=0.05, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solvecc(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)) 1e-6
end

# rescaled variable
let
    m = CCModel()
    @defIndepNormal(m, x, mean=(0,2),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, (z/2)*(x-1) <= -1, with_probability=0.05, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solvecc(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -2/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)) 1e-6
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
