# extensions to JuMP for (robust) chance constraints


module CCJuMP
importall JuMP

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

function CCModel(;solver=nothing)
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

function affToStr(a::CCAffExpr)

    if length(a.vars) == 0
        return affToStr(a.constant)
    end
    
    m = a.vars[1].m
    ccdata = getCCData(m)

    # don't merge duplicates (yet)
    strs = ["($(affToStr(a.coeffs[i],true)))*$(ccdata.RVnames[a.vars[i].idx])" for i in 1:length(a.vars)]
    return string(join(strs," "), " + ", affToStr(a.constant,true))
end

type ChanceConstr
    ccexpr::CCAffExpr
    sense::Symbol # :(<=) or :(>=), right-hand side assumed to be zero
    with_probability::Float64 # with this probability *or less*
end

function addConstraint(m::Model, constr::ChanceConstr; with_probability::Float64=NaN)
    if !(0 < with_probability < 1)
        error("Must specify with_probability between 0 and 1")
    end
    constr.with_probability = with_probability
    
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
