using Test, Random
using JuMP, JuMPChance
using Distributions
using GLPKMathProgInterface
import ECOS

@testset "Operator overloads and printing" begin
    m = ChanceModel()
    @indepnormal(m, x, mean=1, var=1)
    @indepnormal(m, y, mean=1, var=1)
    @variable(m, v)
    @variable(m, q)
    @indepnormal(m,z[i=1:10],mean=i,var=1)

    # Number
    # Number--IndepNormal
    @test string(1+x) == "(1.0)*x + 1.0"
    @test string(1-x) == "(-1.0)*x + 1.0"
    @test string(2x) == "(2.0)*x + 0.0"
    # Number--RandomAffExpr
    r = 2x+1
    @test string(1+r) == "(2.0)*x + 2.0"
    @test string(1-r) == "(-2.0)*x + 0.0"
    @test string(2r) == "(4.0)*x + 2.0"
    # Number--CCAffExpr
    r2 = (2v+1)*x+1
    @test string(1+r2) == "(2 v + 1)*x + 2"
    @test string(1-r2) == "(-2 v - 1)*x + 0"
    @test string(2r2) == "(4 v + 2)*x + 2"

    # Variable
    # Variable--IndepNormal
    @test string(v+x) == "(1)*x + v"
    @test string(v-x) == "(-1)*x + v"
    @test string(v*x) == "(v)*x + 0"
    # Variable--RandomAffExpr
    @test string(v+r) == "(2)*x + v + 1"
    @test string(v-r) == "(-2)*x + v - 1"
    @test string(v*r) == "(2 v)*x + v"
    # Variable--CCAffExpr
    @test string(v+r2) == "(2 v + 1)*x + v + 1"
    @test string(v-r2) == "(-2 v - 1)*x + v - 1"
    # Variable*CCAffExpr not valid

    # AffExpr
    # AffExpr--IndepNormal
    a = 2q+3
    @test string(a+x) == "(1)*x + 2 q + 3"
    @test string(a-x) == "(-1)*x + 2 q + 3"
    @test string(a*x) == "(2 q + 3)*x + 0"
    # AffExpr--RandomAffExpr
    @test string(a+r) == "(2)*x + 2 q + 4"
    @test string(a-r) == "(-2)*x + 2 q + 2"
    @test string(a*r) == "(4 q + 6)*x + 2 q + 3"
    # AffExpr--CCAffExpr
    @test string(a+r2) == "(2 v + 1)*x + 2 q + 4"
    @test string(a-r2) == "(-2 v - 1)*x + 2 q + 2"
    # AffExpr*CCAffExpr not valid

    # IndepNormal
    @test string(-x) == "(-1.0)*x + 0.0"
    # IndepNormal--Number
    @test string(x+1) == "(1.0)*x + 1.0"
    @test string(x-1) == "(1.0)*x + -1.0"
    @test string(x*3) == "(3.0)*x + 0.0"
    @test string(x/2) == "(0.5)*x + 0.0"
    # IndepNormal--Variable
    @test string(x+v) == "(1)*x + v"
    @test string(x-v) == "(1)*x + -v"
    @test string(x*v) == "(v)*x + 0"
    # IndepNormal--AffExpr
    @test string(x+a) == "(1)*x + 2 q + 3"
    @test string(x-a) == "(1)*x + -2 q - 3"
    @test string(x*a) == "(2 q + 3)*x + 0"
    # IndepNormal--IndepNormal
    @test string(y+x) == "(1.0)*y + (1.0)*x + 0.0"
    @test string(y-x) == "(1.0)*y + (-1.0)*x + 0.0"
    # IndepNormal*IndepNormal not valid
    # IndepNormal--RandomAffExpr
    @test string(x+r) == "(3.0)*x + 1.0"
    @test string(x-r) == "(-1.0)*x + -1.0"
    # IndepNormal*RandomAffExpr not valid
    # IndepNormal--CCAffExpr
    @test string(x+r2) == "(2 v + 2)*x + 1"
    @test string(x-r2) == "(-2 v)*x + -1"
    # IndepNormal*CCAffExpr not valid

    # RandomAffExpr
    # RandomAffExpr--Number
    @test string(r+1) == "(2.0)*x + 2.0"
    @test string(r-1) == "(2.0)*x + 0.0"
    @test string(r*2) == "(4.0)*x + 2.0"
    @test string(r/2) == "(1.0)*x + 0.5"
    # RandomAffExpr--Variable
    @test string(r+v) == "(2)*x + v + 1"
    @test string(r-v) == "(2)*x + -v + 1"
    @test string(r*v) == "(2 v)*x + v"
    # RandomAffExpr--AffExpr
    @test string(r+a) == "(2)*x + 2 q + 4"
    @test string(r-a) == "(2)*x + -2 q - 2"
    @test string(r*a) == "(4 q + 6)*x + 2 q + 3"
    # RandomAffExpr--IndepNormal
    @test string(r+y) == "(2.0)*x + (1.0)*y + 1.0"
    @test string(r-y) == "(2.0)*x + (-1.0)*y + 1.0"
    # RandomAffExpr*IndepNormal not valid
    # RandomAffExpr--RandomAffExpr
    @test string(r+(3y-1)) == "(2.0)*x + (3.0)*y + 0.0"
    @test string(r-(3y-1)) == "(2.0)*x + (-3.0)*y + 2.0"
    # RandomAffExpr*RandomAffExpr not valid
    # RandomAffExpr--CCAffExpr
    @test string(r+r2) == "(2 v + 3)*x + 2"
    @test string(r-r2) == "(-2 v + 1)*x + 0"
    # RandomAffExpr*CCAffExpr not valid

    # CCAffExpr
    # CCAffExpr--Number
    @test string(r2+1) == "(2 v + 1)*x + 2"
    @test string(r2-1) == "(2 v + 1)*x + 0"
    @test string(r2*2) == "(4 v + 2)*x + 2"
    @test string(r2/2) == "(v + 0.5)*x + 0.5"
    # CCAffExpr--Variable
    @test string(r2+v) == "(2 v + 1)*x + v + 1"
    @test string(r2-v) == "(2 v + 1)*x + -v + 1"
    # CCAffExpr*Variable not valid
    # CCAffExpr--AffExpr
    @test string(r2+a) == "(2 v + 1)*x + 2 q + 4"
    @test string(r2-a) == "(2 v + 1)*x + -2 q - 2"
    # CCAffExpr*AffExpr not valid
    # CCAffExpr--IndepNormal
    @test string(r2+y) == "(2 v + 1)*x + (1)*y + 1"
    @test string(r2-y) == "(2 v + 1)*x + (-1)*y + 1"
    # CCAffExpr*IndepNormal not valid
    # CCAffExpr--RandomAffExpr
    @test string(r2+r) == "(2 v + 3)*x + 2"
    @test string(r2-r) == "(2 v - 1)*x + 0"
    # CCAffExpr*RandomAffExpr not valid
    # CCAffExpr--CCAffExpr
    @test string(r2+r2) == "(4 v + 2)*x + 2"
    @test string(r2-r2) == "(0)*x + 0"
    # CCAffExpr*CCAffExpr not valid


    @test getmean(z[5]) == 5
    @test getvariance(z[5]) == 1

    @constraint(m, (3v+1)*x + 10 <= 20, with_probability=0.95)
    @test string(JuMPChance.getCCData(m).chanceconstr[end]) == "(3 v + 1)*x + -10 <= 0, with probability 0.95"
    @test_throws ErrorException @constraint(m, (3v+1)*x + 10 >= 20, with_probability=10.0)
    @constraint(m, (3v+1)*x + 10 >= 20, with_probability=0.96)
    @test string(JuMPChance.getCCData(m).chanceconstr[end]) == "(3 v + 1)*x + -10 >= 0, with probability 0.96"

    ccaff = JuMPChance.CCAffExpr()
    ccaff.constant = 2v
    @test string(ccaff) == "2 v"

    raff = JuMPChance.RandomAffExpr()
    raff.constant = 10
    @test string(raff) == "10.0"


    @test_throws ErrorException @indepnormal(m, q, mean=1, var=-1)
    @test_throws ErrorException @indepnormal(m, q, mean=1, var=(-1,1))

    try @eval @indepnormal(m, f(x), mean=1, var=1)
    catch err
        @test err isa LoadError
        @test err.error isa ErrorException
        @test startswith(err.error.msg, "Syntax error: Expected")
    end

    @test string(z[1]+z[2]-2z[3]+10) == "(1.0)*z[1] + (1.0)*z[2] + (-2.0)*z[3] + 10.0"

    @test string(v*(x-1)) == "(v)*x + -v"
    @test string(x*(v-1)) == "(v - 1)*x + 0"

    jm = Model()
    @test_throws ErrorException JuMPChance.getCCData(jm)

end

@testset "Printing two-sided constraints" begin
    m = ChanceModel()
    @variable(m, x)
    @indepnormal(m, ξ, mean = 0, var = 1)
    cc = JuMPChance.getCCData(m)

    @constraint(m, -1 ≤ x*ξ ≤ 1, with_probability = 0.95, approx="1.25")
    @test string(cc.twosidechanceconstr[end]) == "-1 <= (x)*ξ + 0 <= 1, with probability 0.95"
    @constraint(m, -1 ≤ x*ξ ≤ x, with_probability = 0.95, approx="1.25")
    @test string(cc.twosidechanceconstr[end]) == "-1 <= (x)*ξ + 0 <= x, with probability 0.95"
    @variable(m,y[1:3])
    @constraint(m, -1 ≤ x*ξ ≤ sum(y[i] for i=1:3), with_probability = 0.95, approx="1.25")
    @test string(cc.twosidechanceconstr[end]) == "-1 <= (x)*ξ + 0 <= y[1] + y[2] + y[3], with probability 0.95"
end

@testset "Non-chance constraint model" begin
    m = ChanceModel()
    @indepnormal(m, ξ, mean=0, var=1)
    @variable(m, x)
    @objective(m, Min, x)
    @constraint(m, JuMPChance.CCAffExpr() + x >= 0, with_probability = 0.5)
    status = solve(m, silent=true)
    @show status == :Optimal
    @show getvalue(x) ≈ 0.0 rtol=1e-6 atol=1e-6
end

@testset "Basic chance constraint model" begin
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @indepnormal(m, x, mean=0, var=1)
        @variable(m, z >= -100) # so original problem is bounded

        @objective(m, Min, z)
        cref = @constraint(m, z*x >= -1, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/quantile(Normal(0,1),0.95) rtol=1e-6 atol=1e-6
        @test JuMPChance.satisfied_with_probability(cref) > 0.95 - 1e-4
    end
end

@testset "Flipped constraint sense" begin
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @indepnormal(m, x, mean=0, var=1)
        @variable(m, z >= -100) # so original problem is bounded

        @objective(m, Min, z)
        cref = @constraint(m, -z*x <= 1, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/quantile(Normal(0,1),0.95) rtol=1e-6 atol=1e-6
        @test JuMPChance.satisfied_with_probability(cref) > 0.95 - 1e-4
    end
end

@testset "Invariance to transformations" begin
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @indepnormal(m, x, mean=1, var=1)
        @variable(m, z >= -100) # so original problem is bounded

        @objective(m, Min, z)
        @constraint(m, z*(x-1) >= -1, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/quantile(Normal(0,1),0.95) rtol=1e-6 atol=1e-6
    end
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @indepnormal(m, x, mean=2, var=4)
        @variable(m, z >= -100) # so original problem is bounded

        @objective(m, Min, z)
        @constraint(m, z*(x/2-1) >= -1, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/quantile(Normal(0,1),0.95) rtol=1e-6 atol=1e-6
    end
end

@testset "Duplicate terms" begin
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @indepnormal(m, x, mean=0, var=1)
        @variable(m, z >= -100) # so original problem is bounded

        @objective(m, Min, z)
        @constraint(m, (1/2)z*x + (1/2)z*x >= -1, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/quantile(Normal(0,1),0.95) rtol=1e-6 atol=1e-6
    end
end

@testset "@constraints block" begin
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @indepnormal(m, x, mean=0, var=1)
        @variable(m, z >= -100) # so original problem is bounded

        @objective(m, Min, z)
        @constraints m begin
            z >= -100
            z*x >= -1, (with_probability=0.95)
        end

        status = solve(m, method=method, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/quantile(Normal(0,1),0.95) rtol=1e-6 atol=1e-6
    end
end

@testset "Robust but no uncertainty budget" begin
    let
        m = ChanceModel()
        @indepnormal(m, x, mean=(-1,1),var=1)

        @variable(m, z >= -100)
        @objective(m, Min, z)

        @constraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
        status = solve(m, method=:Cuts, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/quantile(Normal(0,1),0.95) rtol=1e-6 atol=1e-6
    end
    let
        m = ChanceModel()
        @indepnormal(m, x, mean=0,var=(0.95,1.05))

        @variable(m, z >= -100)
        @objective(m, Min, z)

        @constraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
        status = solve(m, method=:Cuts, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/quantile(Normal(0,1),0.95) rtol=1e-6 atol=1e-6
    end
    let # flipped signs
        m = ChanceModel()
        @indepnormal(m, x, mean=(-1,1),var=1)

        @variable(m, z >= -100)
        @objective(m, Min, z)

        @constraint(m, -z*x <= 1, with_probability=0.95, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
        status = solve(m, method=:Cuts,silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/quantile(Normal(0,1),0.95) rtol=1e-6 atol=1e-6
    end
end



# quadratic objective
@testset "Quadratic objective" begin
    for method in [:Reformulate, :Cuts], linearize in [true, false]
        (method == :Reformulate && linearize) && continue
        m = ChanceModel()
        @indepnormal(m, x, mean=0,var=1)

        @variable(m, z >= -100)
        @objective(m, Min, z+2z^2)

        @constraint(m, z*x >= -1, with_probability=0.95)
        status = solve(m, method=method, linearize_objective=linearize, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/4 rtol=1e-4 atol=1e-4
    end
    for linearize in [true, false]
        m = ChanceModel()
        @indepnormal(m, x, mean=(-1,1),var=1)

        @variable(m, z >= -100)
        @objective(m, Min, z+2z^2)

        @constraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=0, uncertainty_budget_variance=0)
        status = solve(m, method=:Cuts, linearize_objective=linearize, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/4 rtol=1e-4 atol=1e-4
    end
end

@testset "Uncertainty budget for mean" begin
    m = ChanceModel()
    @indepnormal(m, x, mean=(-1,1),var=1)

    @variable(m, z >= -100)
    @objective(m, Min, z)

    cref = @constraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=0)
    status = solve(m, method=:Cuts, silent=true)
    @test status == :Optimal
    @test getvalue(z) ≈ -1/(1+quantile(Normal(0,1),0.95)) rtol=1e-6 atol=1e-6
    @test_throws ErrorException JuMPChance.satisfied_with_probability(cref)
end

@testset "Uncertain variance, but no budget" begin
    m = ChanceModel()
    @indepnormal(m, x, mean=(-1,1),var=(0.95,1.05))

    @variable(m, z >= -100)
    @objective(m, Min, z)

    @constraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=0)
    status = solve(m, method=:Cuts, silent=true)
    @test status == :Optimal
    @test getvalue(z) ≈ -1/(1+quantile(Normal(0,1),0.95)) rtol=1e-6 atol=1e-6
end

@testset "Uncertain variance, with budget" begin
    let
        m = ChanceModel()
        @indepnormal(m, x, mean=(-1,1),var=(0.95,1.05))

        @variable(m, z >= -100)
        @objective(m, Min, z)

        @constraint(m, z*x >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
        status = solve(m, method=:Cuts,silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)) rtol=1e-6 atol=1e-6
    end
    let # shifted mean
        m = ChanceModel()
        @indepnormal(m, x, mean=(0,2),var=(0.95,1.05))

        @variable(m, z >= -100)
        @objective(m, Min, z)

        @constraint(m, z*(x-1) >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
        status = solve(m, method=:Cuts,silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -1/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)) rtol=1e-6 atol=1e-6
    end
    let # rescaled variable
        m = ChanceModel()
        @indepnormal(m, x, mean=(0,2),var=(0.95,1.05))

        @variable(m, z >= -100)
        @objective(m, Min, z)

        @constraint(m, (z/2)*(x-1) >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
        status = solve(m, method=:Cuts,silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ -2/(1+sqrt(1.05)*quantile(Normal(0,1),0.95)) rtol=1e-6 atol=1e-6
    end
    let # more than one R.V.
        m = ChanceModel()
        @indepnormal(m, x, mean=(-1,1),var=(0.95,1.05))
        @indepnormal(m, y, mean=(-0.01,0.01),var=0.01)

        @variable(m, z >= -100)
        @objective(m, Min, z)

        @constraint(m, z*x + y >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
        status = solve(m, method=:Cuts,silent=true)
        @test status == :Optimal
        # In mathematica: Minimize[{z, z - \[Nu]*Sqrt[1.05*z^2 + 0.01] >= -1}, z]
        @test getvalue(z) ≈ -0.36431227017642165 rtol=1e-5 atol=1e-5
    end
end

@testset "Integer variables" begin
    m = ChanceModel(solver=GLPKSolverMIP())
    @indepnormal(m, x, mean=(-1,1),var=(0.95,1.05))
    @indepnormal(m, y, mean=(-0.01,0.01),var=0.01)

    @variable(m, z >= -100, Int)
    @objective(m, Min, z)

    @constraint(m, z*x + sum(y for i=1:1) >= -1, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)
    status = solve(m, method=:Cuts,silent=true)
    @test status == :Optimal
    @test getvalue(z) ≈ 0.0 rtol=1e-5 atol=1e-5
end

# variance == 0 corner case
@testset "Variance == 0 corner case" begin
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @indepnormal(m, x, mean=0, var=1)
        @variable(m, z <= 100) # so original problem is bounded
        @variable(m, y == 0)

        @objective(m, Max, z)
        @constraint(m, z + x*y <= 20, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ 20 rtol=1e-6 atol=1e-6
    end
    for method in [:Reformulate,:Cuts]
        m = ChanceModel()
        @indepnormal(m, x, mean=0, var=1)
        @variable(m, z <= 100)
        @variable(m, y == 0)

        @objective(m, Max, z)
        @constraint(m, -z - x*y >= -20, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ 20 rtol=1e-6 atol=1e-6
    end
    let
        m = ChanceModel()
        @indepnormal(m, x, mean=0, var=1)
        @variable(m, z <= 100)
        @variable(m, y == 0)

        @objective(m, Max, z)
        @constraint(m, -z - x*y >= -20, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)

        status = solve(m, method=:Cuts, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ 20 rtol=1e-6 atol=1e-6
    end
    let
        m = ChanceModel()
        @indepnormal(m, x, mean=0, var=1)
        @variable(m, z <= 100)
        @variable(m, y == 0)

        @objective(m, Max, z)
        @constraint(m, z + x*y <= 20, with_probability=0.95, uncertainty_budget_mean=1, uncertainty_budget_variance=1)

        status = solve(m, method=:Cuts, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ 20 rtol=1e-6 atol=1e-6
    end
end


@testset "Special cases where chance constraint becomes linear" begin
    for method in [:Reformulate, :Cuts]
        m = ChanceModel()
        @indepnormal(m, x[1:2], mean=1, var=1)
        @variable(m, z <= 100)
        @variable(m, y[1:2])

        @objective(m, Max, z)
        @constraint(m, y[1] == y[2])
        @constraint(m, z + sum(x[i]*y[i] for i=1:2) <= 20, with_probability=0.95)


        status = solve(m, method=method, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ 20 rtol=1e-5 atol=1e-5
    end
    for method in [:Reformulate, :Cuts]
        m = ChanceModel()
        @indepnormal(m, x[1:2], mean=1, var=1)
        @variable(m, z >= 100)
        @variable(m, y[1:2])

        @objective(m, Min, z)
        @constraint(m, y[1] == y[2])
        @constraint(m, z + sum(x[i]*y[i] for i=1:2) >= 200, with_probability=0.95)

        status = solve(m, method=method, silent=true)
        @test status == :Optimal
        @test getvalue(z) ≈ 200 rtol=1e-5 atol=1e-5
    end
end

@testset "Compare with manual reformulation" begin
    for method in [:Reformulate, :Cuts], ϵ in [0.001, 0.01, 0.005, 0.05], σ² in [0, 1]
        m = ChanceModel()
        @indepnormal(m, ω₁, mean=0, var=σ²)
        @indepnormal(m, ω₂, mean=1, var=σ²)
        @variable(m, z)
        @variable(m, x[1:2])
        @constraint(m, x[1] >= 3)
        @constraint(m, x[2] >= 2)
        @constraint(m, x[1]*ω₁ + x[2]*ω₂ <= z, with_probability=1-ϵ)
        @objective(m, Min, z)

        status = solve(m, probability_tolerance=1e-12, method=method, silent=true)

        objval = getobjectivevalue(m)
        @test status == :Optimal

        m = Model(solver=ECOS.ECOSSolver(verbose=0))
        nu = quantile(Normal(0,1),1-ϵ)
        @variable(m, z)
        @variable(m, x[1:2])
        @variable(m, varterm[1:2] >= 0)
        @variable(m, t >= 0)
        @constraint(m, x[1] >= 3)
        @constraint(m, x[2] >= 2)
        @constraint(m, defvar[i=1:2], varterm[i] == sqrt(σ²)*x[i])
        @constraint(m, x[2] + nu*t <= z)
        @constraint(m, sum( varterm[i]^2 for i=1:2) <= t^2)
        @objective(m, Min, z)

        status = solve(m)
        @test status == :Optimal
        @test getobjectivevalue(m) ≈ objval rtol=1e-5 atol=1e-5

    end
end


# test for treating a two-sided constraints as two one sided constraints
@testset "Two-sided constraints vs. two one-sided constraints" begin
    # (a) one random variable model, with mean = 0 and var = 1
    for method in [:Reformulate, :Cuts]
        ϵ = 0.01
        m = ChanceModel()
        @indepnormal(m, ω, mean=0, var=1)
        @variable(m, l)
        @variable(m, u)
        @variable(m, x == 1)
        @constraint(m, l <= x*ω <= u, with_probability=1-ϵ, approx="2.0")
        @objective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        objval = getobjectivevalue(m)

        m = ChanceModel()
        @indepnormal(m, ω, mean=0, var=1)
        @variable(m, l)
        @variable(m, u)
        @variable(m, x == 1)
        @constraint(m, x*ω >= l, with_probability=1-ϵ)
        @constraint(m, x*ω <= u, with_probability=1-ϵ)
        @objective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        @test getobjectivevalue(m) ≈ objval rtol=1e-5 atol=1e-5
    end

    # (b) one random variable model, with mean = 1 and var = 1
    for method in [:Reformulate, :Cuts]
        ϵ = 0.005
        m = ChanceModel()
        @indepnormal(m, ω, mean=1, var=1)
        @variable(m, l)
        @variable(m, u)
        @variable(m, x == 1)
        @constraint(m, l <= x*ω <= u, with_probability=1-ϵ, approx="2.0")
        @objective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        objval = getobjectivevalue(m)

        m = ChanceModel()
        @indepnormal(m, ω, mean=1, var=1)
        @variable(m, l)
        @variable(m, u)
        @variable(m, x == 1)
        @constraint(m, x*ω >= l, with_probability=1-ϵ)
        @constraint(m, x*ω <= u, with_probability=1-ϵ)
        @objective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        @test getobjectivevalue(m) ≈ objval rtol=1e-5 atol=1e-5
    end

    # (c) multiple random variable model (μ=1, σ²=1), with same coefficients - linear case
    for method in [:Reformulate, :Cuts]
        ϵ = 0.005
        m = ChanceModel()
        @indepnormal(m, ω[1:4], mean=1, var=1)
        @variable(m, l)
        @variable(m, u)
        @variable(m, x[1:4])
        @constraint(m, xcons[i=1:4], x[i] == 1)
        @constraint(m, l <= sum(x[i]*ω[i] for i=1:4) <= u, with_probability=1-ϵ, approx="2.0")
        @objective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        objval = getobjectivevalue(m)

        m = ChanceModel()
        @indepnormal(m, ω[1:4], mean=1, var=1)
        @variable(m, l)
        @variable(m, u)
        @variable(m, x[1:4])
        @constraint(m, xcons[i=1:4], x[i] == 1)
        @constraint(m, sum(x[i]*ω[i] for i=1:4)>= l, with_probability=1-ϵ)
        @constraint(m, sum(x[i]*ω[i] for i=1:4) <= u, with_probability=1-ϵ)
        @objective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        @test getobjectivevalue(m) ≈ objval rtol=1e-5 atol=1e-5
    end

    # (d) multiple random variable model (μ=1, σ²=1), with different coefficients
    for method in [:Reformulate, :Cuts]
        ϵ = 0.005
        m = ChanceModel()
        @indepnormal(m, ω[1:4], mean=1, var=1)
        @variable(m, l)
        @variable(m, u)
        @variable(m, x[1:4])
        @constraint(m, xcons[i=1:4], x[i] == i)
        @constraint(m, l <= sum(i*x[i]*ω[i] for i=1:4) <= u, with_probability=1-ϵ, approx="2.0")
        @objective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        objval = getobjectivevalue(m)

        m = ChanceModel()
        @indepnormal(m, ω[1:4], mean=1, var=1)
        @variable(m, l)
        @variable(m, u)
        @variable(m, x[1:4])
        @constraint(m, xcons[i=1:4], x[i] == i)
        @constraint(m, sum(i*x[i]*ω[i] for i=1:4) >= l, with_probability=1-ϵ)
        @constraint(m, sum(i*x[i]*ω[i] for i=1:4) <= u, with_probability=1-ϵ)
        @objective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        @test getobjectivevalue(m) ≈ objval rtol=1e-5 atol=1e-5
    end

    # (e) multiple random variable model (μ=1, σ²=0), with different coefficients
    for method in [:Reformulate, :Cuts]
        ϵ = 0.005
        m = ChanceModel()
        @indepnormal(m, ω[1:4], mean=1, var=0)
        @variable(m, l)
        @variable(m, u)
        @variable(m, x[1:4])
        @constraint(m, xcons[i=1:4], x[i] == i)
        @constraint(m, l <= sum(i*x[i]*ω[i] for i=1:4) <= u, with_probability=1-ϵ, approx="2.0")
        @objective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        objval = getobjectivevalue(m)

        m = ChanceModel()
        @indepnormal(m, ω[1:4], mean=1, var=0)
        @variable(m, l)
        @variable(m, u)
        @variable(m, x[1:4])
        @constraint(m, xcons[i=1:4], x[i] == i)
        @constraint(m, sum(i*x[i]*ω[i] for i=1:4) >= l, with_probability=1-ϵ)
        @constraint(m, sum(i*x[i]*ω[i] for i=1:4) <= u, with_probability=1-ϵ)
        @objective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        @test getobjectivevalue(m) ≈ objval rtol=1e-5 atol=1e-5
    end
end


const prob_guarantee = 1.25
@testset "Basic two-sided constraints" begin
    for ϵ in (0.1, 0.05, 0.005, 0.0005), method in (:Cuts, :Reformulate)
        m = ChanceModel()
        @indepnormal(m, ξ, mean=0, var=1)
        @variable(m, l)
        @variable(m, u)
        @variable(m, x == 1)
        cref = @constraint(m, l ≤ x*ξ ≤ u, with_probability=1-ϵ, approx="1.25")
        @objective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        violation = 1- JuMPChance.satisfied_with_probability(cref)
        @test violation ≤ prob_guarantee*ϵ + 1e-5
        @test violation ≥ ϵ - 1e-5 # should be tight

        lval = getvalue(l)
        uval = getvalue(u)
        m = ChanceModel()
        @indepnormal(m, ξ, mean=1, var=1) # translate with mean
        @variable(m, l)
        @variable(m, u)
        @variable(m, x == 1)
        cref = @constraint(m, l ≤ x*ξ ≤ u, with_probability=1-ϵ, approx="1.25")
        @objective(m, Min, u-2l)

        solve(m, method=method, silent=true)

        violation = 1- JuMPChance.satisfied_with_probability(cref)
        @test violation ≤ prob_guarantee*ϵ + 1e-5
        @test violation ≥ ϵ - 1e-5

        @test getobjectivevalue(m) ≈ (uval+1) - 2(lval+1)
    end

    # duplicates
    for ϵ in (0.1, 0.05, 0.005, 0.0005), method in (:Cuts, :Reformulate)
        m = ChanceModel()
        @indepnormal(m, ξ, mean=0, var=1)
        @variable(m, l)
        @variable(m, u)
        @variable(m, x == 1)
        cref = @constraint(m, l ≤ x*ξ + x*ξ ≤ u, with_probability=1-ϵ, approx="1.25")
        @objective(m, Min, u-l)

        solve(m, method=method, silent=true)

        violation = 1- JuMPChance.satisfied_with_probability(cref)
        @test violation ≤ prob_guarantee*ϵ + 1e-5
        @test violation ≥ ϵ - 1e-5
    end

    # constant bounds
    let
        ϵ = 0.05
        m = ChanceModel()
        @indepnormal(m, ξ, mean=0, var=1)
        @variable(m, x >= 0)
        cref = @constraint(m, -1 ≤ x*ξ ≤ 1, with_probability=1-ϵ, approx="1.25")
        @objective(m, Max, x)
        solve(m, method=:Reformulate)
        violation = 1- JuMPChance.satisfied_with_probability(cref)
        @test violation ≤ prob_guarantee*ϵ + 1e-5
        @test violation ≥ ϵ - 1e-5
    end

    # absolute value, random objective
    let
        ϵ = 0.05
        Random.seed!(998)
        for q in 1:100
            m = ChanceModel()
            @indepnormal(m, ξ, mean=0, var=1)
            @variable(m, -10 ≤ x ≤ 10)
            @variable(m, 0 ≤ c ≤ 10)
            @variable(m, t ≤ 100)
            cref = @constraint(m, -t ≤ x*ξ + c ≤ t, with_probability=1-ϵ, approx="1.25")
            @objective(m, Min, (rand()-0.5)*x + (rand()-0.5)*c + 0.01*t)
            solve(m, method=:Reformulate)
            violation = 1- JuMPChance.satisfied_with_probability(cref)
            @test violation ≤ prob_guarantee*ϵ + 1e-5
        end
    end
end

@testset "Empty constraints" begin

    m = ChanceModel()
    @indepnormal(m, ξ, mean=0, var=1)
    @variable(m, x)
    @variable(m, u)

    @objective(m, Max, x)
    ex = zero(JuMPChance.CCAffExpr) + x
    @constraint(m, -1 ≤ ex ≤ 1, with_probability=0.95, approx="1.25")
    solve(m, method=:Reformulate)
    @test getvalue(x) ≈ 1.0 rtol=1e-5 atol=1e-5

    m = ChanceModel()
    @indepnormal(m, ξ, mean=0, var=1)
    @variable(m, x)
    @variable(m, u)

    @objective(m, Max, x)
    ex = zero(JuMPChance.CCAffExpr) + x
    @constraint(m, ex ≤ 0.5, with_probability=0.95)
    solve(m, method=:Reformulate)
    @test getvalue(x) ≈ 0.5 rtol=1e-5 atol=1e-5

end

@testset "JuMP macro corner case" begin
    m = ChanceModel()

    @variable(m, x)
    @indepnormal(m, ξ, mean = 0, var = 1)

    @expression(m, ex, (x+1)*(2ξ))
    @test string(ex) == "(2 x + 2)*ξ + 0"

end
