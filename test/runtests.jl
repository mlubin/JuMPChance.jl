using CCJuMP, JuMP
using Base.Test
using Distributions


let
    m = CCModel()
    @defIndepNormal(m, x, mean=1, var=1)
    #x = IndepNormal(m, 1.0,1.0, "x")
    @test affToStr(1+x) == "(1)*x + 1"

    @defVar(m, v)

    @test affToStr((3v+1)*x+10) == "(3 v + 1)*x + 10"

    @defIndepNormal(m,z[i=1:10],mean=i,var=1)
    @test affToStr(v*z[1]+3.5) == "(v)*z[1] + 3.5"

    @test getMean(z[5]) == 5
    @test getVar(z[5]) == 1

    c = (3v+1)*x + 10 <= 20
    @test conToStr(c) == "(3 v + 1)*x + -10 <= 0"
    addConstraint(m, c, with_probability=0.95)
    @test conToStr(CCJuMP.getCCData(m).chanceconstr[1]) == "(3 v + 1)*x + -10 <= 0, with probability 0.95"
end

let
    for method in [:Reformulate,:Cuts]
        m = CCModel()
        @defIndepNormal(m, x, mean=0, var=1)
        @defVar(m, z >= -100) # so original problem is bounded

        @setObjective(m, Min, z)
        addConstraint(m, z*x <= -1, with_probability=0.05)

        solvecc(m, method=method)
        @test_approx_eq_eps getValue(z) -1/quantile(Normal(0,1),0.95) 1e-6
    end
end

