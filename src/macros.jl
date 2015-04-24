using Base.Meta

import JuMP.@gendict

macro defIndepNormal(m, x, mean, var)
    m = esc(m)
    @assert isexpr(mean,:kw) && mean.args[1] == :mean
    @assert isexpr(var,:kw) && var.args[1] == :var
    mean = esc(mean.args[2])
    var = esc(var.args[2])

    if isa(x,Symbol)
        # easy case
        return quote
            $(esc(x)) = IndepNormal($m,$mean,$var,$(string(x)))
            nothing
        end
    else
        if !isexpr(x,:ref)
            error("Syntax error: Expected $var to be of form var[...]")
        end

        condition = :()
        refcall, idxvars, idxsets, idxpairs = JuMP.buildrefsets(x)
        varname = JuMP.getname(x)

        varstr = :(string($(string(varname)),"["))
        for idxvar in idxvars
            push!(varstr.args,:(string($(esc(idxvar)))))
            push!(varstr.args,",")
        end
        deleteat!(varstr.args,length(varstr.args))
        push!(varstr.args,"]")

        code = :( $(refcall) = IndepNormal($m, $mean, $var, $varstr ) )
        looped = JuMP.getloopedcode(x, code, condition, idxvars, idxsets, idxpairs, :IndepNormal)
        return quote
            $looped
            $(esc(varname))
        end
    end
end

# Extensions to make JuMP macros work with IndepNormals
function JuMP.addToExpression(aff::JuMP.GenericAffExpr, c, x)
    return aff + c*x
end
function JuMP.addToExpression(val::Real, c, x::Union(CCAffExpr,RandomAffExpr))
    return val + c*x
end

function JuMP._construct_constraint!(faff::CCAffExpr, sense::Symbol)
    sense in JuMP.valid_senses || error("Unrecognized sense $sense")
    if sense == :(<=) || sense == :≤
        return ChanceConstr(faff, :(<=))
    elseif sense == :(>=) || sense == :≥
        return ChanceConstr(faff, :(>=))
    elseif sense == :(==)
        error("Equality chance constraints not supported")
    end
    error("Unrecognized constraint type $sense")
end

JuMP._construct_constraint!(faff::CCAffExpr, lb::AffExpr, ub::AffExpr) = TwoSideChanceConstr(faff, lb, ub)
JuMP._construct_constraint!(faff::CCAffExpr, lb, ub) = TwoSideChanceConstr(faff, convert(AffExpr,lb), convert(AffExpr,ub))
