using JuMPChance, JuMP
using Base.Test
using Distributions


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
    @test getVar(z[5]) == 1

    c = (3v+1)*x + 10 <= 20
    @test conToStr(c) == "(3 v + 1)*x + -10 <= 0"
    addConstraint(m, c, with_probability=0.05)
    @test conToStr(JuMPChance.getCCData(m).chanceconstr[1]) == "(3 v + 1)*x + -10 <= 0, with probability 0.05"
    @test_throws ErrorException addConstraint(m, c, with_probability=10.0)
    c = (3v+1)*x + 10 >= 20
    @test conToStr(c) == "(3 v + 1)*x + -10 >= 0"

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
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        addConstraint(m, z*x <= -1, with_probability=0.05)

        status = solvechance(m, method=method)
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
        addConstraint(m, -z*x >= 1, with_probability=0.05)

        status = solvechance(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
    end
end

# invariance to transformations
let
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @defIndepNormal(m, x, mean=1, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        addConstraint(m, z*(x-1) <= -1, with_probability=0.05)

        status = solvechance(m, method=method)
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
        addConstraint(m, z*(x/2-1) <= -1, with_probability=0.05)

        status = solvechance(m, method=method)
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
        addConstraint(m, (1/2)z*x + (1/2)z*x <= -1, with_probability=0.05)

        status = solvechance(m, method=method)
        @test status == :Optimal
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
    end
end

# robust but no uncertainty budget
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(-1,1),var=1)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, z*x <= -1, with_probability=0.05, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
    status = solvechance(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
end
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=0,var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, z*x <= -1, with_probability=0.05, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
    status = solvechance(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
end

# flipped signs
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(-1,1),var=1)

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, -z*x >= 1, with_probability=0.05, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
    status = solvechance(m, method=:Cuts)
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

        addConstraint(m, z*x <= -1, with_probability=0.05)
        status = solvechance(m, method=method, linearize_objective=linearize)
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

        addConstraint(m, z*x <= -1, with_probability=0.05, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
        status = solvechance(m, method=:Cuts, linearize_objective=linearize)
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

    addConstraint(m, z*x <= -1, with_probability=0.05, uncertainty_budget_mean=1, uncertainty_budget_variance=0)
    status = solvechance(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/(1+quantile(Normal(0,1),0.95)) 1e-6
end

# uncertain variance, but no budget
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, z*x <= -1, with_probability=0.05, uncertainty_budget_mean=1, uncertainty_budget_variance=0)
    status = solvechance(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/(1+quantile(Normal(0,1),0.95)) 1e-6
end

# uncertain variance, with budget
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(-1,1),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, z*x <= -1, with_probability=0.05, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solvechance(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)) 1e-6
end

# shifted
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(0,2),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, z*(x-1) <= -1, with_probability=0.05, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solvechance(m, method=:Cuts)
    @test status == :Optimal
    @test_approx_eq_eps getValue(z) -1/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)) 1e-6
end

# rescaled variable
let
    m = ChanceModel()
    @defIndepNormal(m, x, mean=(0,2),var=(0.95,1.05))

    @defVar(m, z >= -100)
    @setObjective(m, Min, z)

    addConstraint(m, (z/2)*(x-1) <= -1, with_probability=0.05, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solvechance(m, method=:Cuts)
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

    addConstraint(m, z*x + y <= -1, with_probability=0.05, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solvechance(m, method=:Cuts)
    @test status == :Optimal
    # In mathematica: Minimize[{z, z - \[Nu]*Sqrt[1.05*z^2 + 0.01] >= -1}, z]
    @test_approx_eq_eps getValue(z) -0.36431227017642165 1e-5
end
