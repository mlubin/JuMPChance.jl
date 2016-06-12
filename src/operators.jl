
import Base: +, -, /, *, <=, >=


# IndepNormal
(+)(lhs::IndepNormal,rhs::IndepNormal) = RandomAffExpr([lhs,rhs],[1.0,1.0],0.0)
(-)(lhs::IndepNormal,rhs::IndepNormal) = RandomAffExpr([lhs,rhs],[1.0,-1.0],0.0)
(/)(lhs::IndepNormal,rhs::Number) = RandomAffExpr([lhs],[1/rhs],0.0)
(-)(rhs::IndepNormal) = RandomAffExpr([rhs],[-1.0],0.0)


for op in (:+, :-, :*)
    @eval begin
        ($op)(lhs::Number, rhs::IndepNormal) = ($op)(lhs, RandomAffExpr([rhs],[1.0],0.0))
        ($op)(lhs::Variable, rhs::IndepNormal) = ($op)(convert(AffExpr,lhs),rhs)
        ($op)(lhs::Variable, rhs::RandomAffExpr) = ($op)(convert(AffExpr,lhs),rhs)
        ($op)(lhs::Variable, rhs::CCAffExpr) = ($op)(convert(AffExpr,lhs),rhs)
    end
    if op == :-
        @eval ($op)(lhs::CCAffExpr, rhs::Variable) = (+)(lhs,-rhs)
        @eval ($op)(lhs::IndepNormal, rhs::Union{Variable,Number}) = (+)(lhs,-rhs)
    else
        @eval ($op)(lhs::IndepNormal, rhs::Union{Variable,Number}) = ($op)(rhs,lhs)
        @eval ($op)(lhs::CCAffExpr, rhs::Variable) = ($op)(rhs,lhs)
    end
end


# AffExpr
# AffExpr--IndepNormal
(+)(lhs::AffExpr, rhs::IndepNormal) = CCAffExpr([rhs],[convert(AffExpr,1.0)],lhs)
(-)(lhs::AffExpr, rhs::IndepNormal) = CCAffExpr([rhs],[convert(AffExpr,-1.0)],lhs)
(*)(lhs::AffExpr, rhs::IndepNormal) = CCAffExpr([rhs],[lhs],AffExpr())
(+)(lhs::IndepNormal, rhs::AffExpr) = rhs+lhs
(-)(lhs::IndepNormal, rhs::AffExpr) = (+)(lhs,-rhs)
(*)(lhs::IndepNormal, rhs::AffExpr) = rhs*lhs

# AffExpr--CCAffExpr
(+)(lhs::AffExpr,rhs::CCAffExpr) = convert(CCAffExpr,lhs)+rhs
(+)(lhs::CCAffExpr,rhs::AffExpr) = rhs+lhs
(-)(lhs::AffExpr,rhs::CCAffExpr) = convert(CCAffExpr,lhs)-rhs
(-)(lhs::CCAffExpr,rhs::AffExpr) = lhs-convert(CCAffExpr,rhs)
Base.promote_rule(::Type{AffExpr},::Type{CCAffExpr}) = CCAffExpr
Base.promote_rule(::Type{Float64},::Type{CCAffExpr}) = CCAffExpr
Base.convert(::Type{CCAffExpr},a::AffExpr) = CCAffExpr(IndepNormal[],AffExpr[],a)

# AffExpr--RandomAffExpr
(+)(lhs::AffExpr,rhs::RandomAffExpr) = lhs+convert(CCAffExpr,rhs)
(-)(lhs::AffExpr,rhs::RandomAffExpr) = lhs-convert(CCAffExpr,rhs)
(*)(lhs::AffExpr,rhs::RandomAffExpr) = CCAffExpr(rhs.vars, [lhs*c for c in rhs.coeffs], lhs*rhs.constant)

# CCAffExpr--RandomAffExpr
(+)(lhs::CCAffExpr,rhs::RandomAffExpr) = lhs+convert(CCAffExpr,rhs)
(-)(lhs::CCAffExpr,rhs::RandomAffExpr) = lhs-convert(CCAffExpr,rhs)

# CCAffExpr--CCAffExpr
# handled by GenericAffExpr fallback

# RandomAffExpr
# RandomAffExpr--Variable
(+)(lhs::RandomAffExpr, rhs::Variable) = CCAffExpr(lhs.vars,[convert(AffExpr,c) for c in lhs.coeffs],rhs+lhs.constant)
(-)(lhs::RandomAffExpr, rhs::Variable) = CCAffExpr(lhs.vars,[convert(AffExpr,c) for c in lhs.coeffs],lhs.constant-rhs)
(*)(lhs::RandomAffExpr, rhs::Variable) = CCAffExpr(lhs.vars,[c*rhs for c in lhs.coeffs],rhs*lhs.constant)

# RandomAffExpr--AffExpr
(+)(lhs::RandomAffExpr,rhs::AffExpr) = rhs+lhs
(-)(lhs::RandomAffExpr,rhs::AffExpr) = (+)(-rhs,lhs)
(*)(lhs::RandomAffExpr,rhs::AffExpr) = rhs*lhs

# RandomAffExpr--CCAffExpr
(+)(lhs::RandomAffExpr,rhs::CCAffExpr) = convert(CCAffExpr,lhs)+rhs
(-)(lhs::RandomAffExpr,rhs::CCAffExpr) = convert(CCAffExpr,lhs)-rhs
