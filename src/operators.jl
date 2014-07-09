

for op in (:+, :-, :*)
    @eval begin
        ($op)(lhs::Union(Variable,Number), rhs::IndepNormal) = ($op)(convert(AffExpr,lhs),rhs)
        ($op)(lhs::Variable, rhs::CCAffExpr) = ($op)(convert(AffExpr,lhs),rhs)
        ($op)(lhs::IndepNormal, rhs::Union(Variable,Number)) = ($op)(rhs,lhs)
        ($op)(lhs::CCAffExpr, rhs::Variable) = ($op)(rhs,lhs)
    end
end


# AffExpr
# AffExpr--IndepNormal
(+)(lhs::AffExpr, rhs::IndepNormal) = CCAffExpr([rhs],[convert(AffExpr,1.0)],lhs)
(-)(lhs::AffExpr, rhs::IndepNormal) = CCAffExpr([rhs],[convert(AffExpr,1.0)],-lhs)
(*)(lhs::AffExpr, rhs::IndepNormal) = CCAffExpr([rhs],[lhs],AffExpr())

# AffExpr--CCAffExpr
Base.promote_rule(::Type{AffExpr},::Type{CCAffExpr}) = CCAffExpr
Base.convert(::Type{CCAffExpr},a::AffExpr) = CCAffExpr(IndepNormal[],AffExpr[],a)

# CCAffExpr--CCAffExpr
# handled by GenericAffExpr fallback

# comparison operators
(<=)(lhs::CCAffExpr, rhs::Number) = ChanceConstr(lhs-rhs, :(<=))
(>=)(lhs::CCAffExpr, rhs::Number) = ChanceConstr(lhs-rhs, :(>=))
# == not valid
