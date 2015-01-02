# extensions to JuMP for (robust) chance constraints


module CCJuMP
importall JuMP
import ECOS # For now, set ECOS as default since it's the only open-source conic solver availible

export CCModel,
    IndepNormal,
    affToStr,
    getMean,
    getVar,
    getStdev,
    solvecc,
    @defIndepNormal


# stores extension data inside JuMP Model
type CCData
    # pure chance constraints
    chanceconstr
    # robust chance constraints
    #robustchanceconstr
    
    numRVs::Int # number of independent normal r.v.'s
    RVmeans
    RVvars
    RVnames
end

function CCModel(;solver=ECOS.ECOSSolver())
    m = Model(solver=solver)
    m.ext[:ChanceConstr] = CCData(Any[],0,Any[],Any[],String[])
    return m
end

function getCCData(m::Model)
    if haskey(m.ext, :ChanceConstr)
        return m.ext[:ChanceConstr]::CCData
    else
        error("This functionality is only available for CCJuMP models")
    end
end


# pointer back to JuMP model and r.v. index
type IndepNormal
    m::Model
    idx::Int
end

function IndepNormal(m::Model, mean, var, name::String)
    ccdata = getCCData(m)
    ccdata.numRVs += 1
    push!(ccdata.RVmeans, mean)
    for v in var
        if v < 0
            error("Invalid value $v for variance of $name")
        end
    end
    push!(ccdata.RVvars, var)
    push!(ccdata.RVnames, name)
    return IndepNormal(m, ccdata.numRVs)
end

function getMean(v::IndepNormal)
    ccdata = getCCData(v.m)
    return ccdata.RVmeans[v.idx]
end

function getVar(v::IndepNormal)
    ccdata = getCCData(v.m)
    return ccdata.RVvars[v.idx]
end

getStdev(v::IndepNormal) = sqrt(getVar(v))

typealias CCAffExpr JuMP.GenericAffExpr{AffExpr,IndepNormal}

CCAffExpr() = CCAffExpr(IndepNormal[],AffExpr[],AffExpr())

Base.print(io::IO, a::CCAffExpr) = print(io, affToStr(a))
Base.show( io::IO, a::CCAffExpr) = print(io, affToStr(a))

# affine expression only involving r.v.'s
typealias RandomAffExpr JuMP.GenericAffExpr{Float64,IndepNormal}

RandomAffExpr() = RandomAffExpr(IndepNormal[],Float64[],0.0)

Base.print(io::IO, a::RandomAffExpr) = print(io, affToStr(a))
Base.show( io::IO, a::RandomAffExpr) = print(io, affToStr(a))

Base.promote_rule(::Type{CCAffExpr},::Type{RandomAffExpr}) = CCAffExpr
Base.promote_rule(::Type{RandomAffExpr},::Type{CCAffExpr}) = CCAffExpr

Base.convert(::Type{CCAffExpr},a::RandomAffExpr) = CCAffExpr(a.vars,[convert(AffExpr,c) for c in a.coeffs],convert(AffExpr,a.constant))



function affToStr(a::CCAffExpr)

    if length(a.vars) == 0
        return JuMP.aff_str(JuMP.REPLMode,a.constant)
    end
    
    m = a.vars[1].m
    ccdata = getCCData(m)

    # don't merge duplicates (yet)
    strs = ["($(JuMP.aff_str(JuMP.REPLMode,a.coeffs[i], show_constant=true)))*$(ccdata.RVnames[a.vars[i].idx])" for i in 1:length(a.vars)]
    return string(join(strs," + "), " + ", JuMP.aff_str(JuMP.REPLMode, a.constant, show_constant=true))
end


function affToStr(a::RandomAffExpr)

    if length(a.vars) == 0
        return string(a.constant)
    end
    
    m = a.vars[1].m
    ccdata = getCCData(m)

    # don't merge duplicates (yet)
    strs = ["($(a.coeffs[i]))*$(ccdata.RVnames[a.vars[i].idx])" for i in 1:length(a.vars)]
    return string(join(strs," + "), " + ", string(a.constant))
end

type ChanceConstr
    ccexpr::CCAffExpr
    sense::Symbol # :(<=) or :(>=), right-hand side assumed to be zero
    with_probability::Float64 # with this probability *or less*
    uncertainty_budget_mean::Int # for now, with Bertsimas-Sim uncertainty
    uncertainty_budget_variance::Int # for now, with Bertsimas-Sim uncertainty
end

ChanceConstr(ccexpr::CCAffExpr,sense::Symbol) = ChanceConstr(ccexpr, sense, NaN, 0, 0)

function addConstraint(m::Model, constr::ChanceConstr; with_probability::Float64=NaN, uncertainty_budget_mean::Int=0, uncertainty_budget_variance::Int=0)
    if !(0 < with_probability < 1)
        error("Must specify with_probability between 0 and 1")
    end
    constr.with_probability = with_probability
    constr.uncertainty_budget_mean = uncertainty_budget_mean
    constr.uncertainty_budget_variance = uncertainty_budget_variance
    
    ccdata = getCCData(m)
    push!(ccdata.chanceconstr, constr)

end



function conToStr(c::ChanceConstr)
    s = "$(affToStr(c.ccexpr)) $(c.sense) 0"
    if isnan(c.with_probability)
        return s
    else
        return s*", with probability $(c.with_probability)"
    end
end

Base.print(io::IO, a::ChanceConstr) = print(io, conToStr(a))
Base.show( io::IO, a::ChanceConstr) = print(io, conToStr(a))

include("operators.jl")
include("macros.jl")
include("solve.jl")

end
