# extensions to JuMP for (robust) chance constraints

__precompile__()

module JuMPChance
using JuMP
using Distributions
import ECOS # For now, set ECOS as default since it's the only open-source conic solver availible

export ChanceModel,
    IndepNormal,
    getmean,
    getvariance,
    getstdev,
    @indepnormal

# stores extension data inside JuMP Model
type CCData
    # pure chance constraints
    chanceconstr
    # robust chance constraints
    #robustchanceconstr
    # two-sided chance constraints
    twosidechanceconstr
    numRVs::Int # number of independent normal r.v.'s
    RVmeans
    RVvars
    RVnames
end

function ChanceModel(;solver=ECOS.ECOSSolver(verbose=false))
    m = Model(solver=solver)
    m.solvehook = solvehook
    m.ext[:ChanceConstr] = CCData(ChanceConstr[],TwoSideChanceConstr[],0,Any[],Any[],AbstractString[])
    return m
end

function getCCData(m::Model)
    if haskey(m.ext, :ChanceConstr)
        return m.ext[:ChanceConstr]::CCData
    else
        error("This functionality is only available for JuMPChance models")
    end
end


# pointer back to JuMP model and r.v. index
type IndepNormal <: JuMP.AbstractJuMPScalar
    m::Model
    idx::Int
end

JuMP.linearindex(v::IndepNormal) = v.idx

function IndepNormal(m::Model, mean, var, name::AbstractString)
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

function getmean(v::IndepNormal)
    ccdata = getCCData(v.m)
    return ccdata.RVmeans[v.idx]
end

@Base.deprecate getMean getmean

function getvariance(v::IndepNormal)
    ccdata = getCCData(v.m)
    return ccdata.RVvars[v.idx]
end

@Base.deprecate getVariance getvariance

getstdev(v::IndepNormal) = sqrt(getvariance(v))

@Base.deprecate getStdev getstdev

const CCAffExpr = JuMP.GenericAffExpr{AffExpr,IndepNormal}

CCAffExpr() = CCAffExpr(IndepNormal[],AffExpr[],AffExpr())

# affine expression only involving r.v.'s
const RandomAffExpr = JuMP.GenericAffExpr{Float64,IndepNormal}

RandomAffExpr() = RandomAffExpr(IndepNormal[],Float64[],0.0)

Base.promote_rule(::Type{CCAffExpr},::Type{RandomAffExpr}) = CCAffExpr
Base.promote_rule(::Type{RandomAffExpr},::Type{CCAffExpr}) = CCAffExpr

Base.convert(::Type{CCAffExpr},a::RandomAffExpr) = CCAffExpr(a.vars,[convert(AffExpr,c) for c in a.coeffs],convert(AffExpr,a.constant))

function Base.show(io::IO,a::CCAffExpr)

    if length(a.vars) == 0
        return print(io,a.constant)
    end
    
    m = a.vars[1].m
    ccdata = getCCData(m)

    v = IndexedVector(JuMP.AffExpr, 0)
    counts = IndexedVector(Int, 0)
    a = merge_duplicates(a, v, counts, m)

    strs = ["($(JuMP.aff_str(JuMP.REPLMode,a.coeffs[i], true)))*$(ccdata.RVnames[a.vars[i].idx])" for i in 1:length(a.vars)]
    print(io,string(join(strs," + "), " + ", JuMP.aff_str(JuMP.REPLMode, a.constant, true)))
end


function Base.show(io::IO,a::RandomAffExpr)

    if length(a.vars) == 0
        return print(io,a.constant)
    end
    
    m = a.vars[1].m
    ccdata = getCCData(m)

    v = IndexedVector(Float64, 0)
    a = merge_duplicates(a, v, m)

    strs = ["($(a.coeffs[i]))*$(ccdata.RVnames[a.vars[i].idx])" for i in 1:length(a.vars)]
    print(io, string(join(strs," + "), " + ", string(a.constant)))
end

type ChanceConstr <: JuMP.AbstractConstraint
    ccexpr::CCAffExpr
    sense::Symbol # :(<=) or :(>=), right-hand side assumed to be zero
    with_probability::Float64 # with this probability *or greater*
    uncertainty_budget_mean::Int # for now, with Bertsimas-Sim uncertainty
    uncertainty_budget_variance::Int # for now, with Bertsimas-Sim uncertainty
end

ChanceConstr(ccexpr::CCAffExpr,sense::Symbol) = ChanceConstr(ccexpr, sense, NaN, 0, 0)

function JuMP.addconstraint(m::Model, constr::ChanceConstr; with_probability::Float64=NaN, uncertainty_budget_mean::Int=0, uncertainty_budget_variance::Int=0)
    if !(0 < with_probability < 1)
        error("Must specify with_probability between 0 and 1")
    end
    if with_probability < 0.5
        error("with_proability < 0.5 is not supported")
    end
    constr.with_probability = with_probability
    constr.uncertainty_budget_mean = uncertainty_budget_mean
    constr.uncertainty_budget_variance = uncertainty_budget_variance
    
    ccdata = getCCData(m)
    push!(ccdata.chanceconstr, constr)
    return ConstraintRef{JuMP.Model,ChanceConstr}(m, length(ccdata.chanceconstr))

end

function satisfied_with_probability(r::ConstraintRef{JuMP.Model,ChanceConstr})
    m = r.m
    ccdata = getCCData(m)
    cc = ccdata.chanceconstr[r.idx]

    no_uncertains = all(x->isa(x,Real), ccdata.RVmeans) && all(x->isa(x,Real), ccdata.RVvars)
    if !no_uncertains
        error("satisfied_with_probability not implemented for distributionally robust constraints")
    end

    nterms = length(cc.ccexpr.vars)
    μ = getvalue(cc.ccexpr.constant)
    σ² = 0.0
    for i in 1:nterms
        μ += getvalue(cc.ccexpr.coeffs[i])*getmean(cc.ccexpr.vars[i])
        σ² += (getstdev(cc.ccexpr.vars[i])*getvalue(cc.ccexpr.coeffs[i]))^2
    end
    σ = sqrt(σ²)
    d = Normal(μ,σ)
    if cc.sense == :(>=)
        return 1 - cdf(d, 0)
    else
        return cdf(d, 0)
    end
end

function JuMP.show(io::IO,c::ChanceConstr)
    s = "$(string(c.ccexpr)) $(c.sense) 0"
    if isnan(c.with_probability)
        return print(io,s)
    else
        return print(io,s*", with probability $(c.with_probability)")
    end
end

type TwoSideChanceConstr <: JuMP.AbstractConstraint
    ccexpr::CCAffExpr
    lb::AffExpr
    ub::AffExpr
    with_probability::Float64 # with this probability *or greater*
    approx::AbstractString # outer approximation model used for this constraint.
end

TwoSideChanceConstr(ccexpr::CCAffExpr,lb::AffExpr,ub::AffExpr) = TwoSideChanceConstr(ccexpr, lb, ub, NaN, "none")

function JuMP.addconstraint(m::Model, constr::TwoSideChanceConstr; with_probability::Float64=NaN, approx::AbstractString="")
    if !(0.5 < with_probability < 1)
        error("Must specify with_probability between 1/2 and 1")
    end
    constr.with_probability = with_probability
    if approx != "1.25" && approx != "2.0"
        error("approx paramater must be provided as \"1.25\" or \"2.0\". See documentation for further details")
    end
    constr.approx = approx

    ccdata = getCCData(m)
    push!(ccdata.twosidechanceconstr, constr)
    return ConstraintRef{JuMP.Model,TwoSideChanceConstr}(m, length(ccdata.twosidechanceconstr))

end

function Base.show(io::IO,c::TwoSideChanceConstr)
    s = "$(string(c.lb)) <= $(string(c.ccexpr)) <= $(string(c.ub))"
    if isnan(c.with_probability)
        return print(io,s)
    else
        return print(io,s*", with probability $(c.with_probability)")
    end
end

function satisfied_with_probability(r::ConstraintRef{JuMP.Model,TwoSideChanceConstr})
    m = r.m
    ccdata = getCCData(m)
    cc = ccdata.twosidechanceconstr[r.idx]

    nterms = length(cc.ccexpr.vars)
    μ = getvalue(cc.ccexpr.constant)
    σ² = 0.0
    for i in 1:nterms
        μ += getvalue(cc.ccexpr.coeffs[i])*getmean(cc.ccexpr.vars[i])
        σ² += (getstdev(cc.ccexpr.vars[i])*getvalue(cc.ccexpr.coeffs[i]))^2
    end
    σ = sqrt(σ²)
    return cdf(Normal(μ,σ), getvalue(cc.ub)) - cdf(Normal(μ,σ), getvalue(cc.lb))

end

include("operators.jl")
include("macros.jl")
include("solve.jl")

end
