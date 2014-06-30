using CCJuMP, JuMP
using Base.Test



m = CCModel()
x = IndepNormal(m, 1.0,1.0, "x")
@test affToStr(1+x) == "(1)*x + 1"

@defVar(m, v)

@test affToStr((3v+1)*x+10) == "(3 v + 1)*x + 10"
