include("distributions.jl")
import MathProgBase

function solvechance(m::Model;method=:Refomulate,linearize_objective::Bool=false,probability_tolerance=0.001,debug::Bool = false, iteration_limit::Int=60, objective_linearization_tolerance::Float64=1e-6, reformulate_quadobj_to_conic::Bool=false, lazy_constraints::Bool=false)
    @assert method == :Reformulate || method == :Cuts

    ccdata = getCCData(m)
    no_uncertains = all([isa(x,Real) for x in ccdata.RVmeans]) && all([isa(x,Real) for x in ccdata.RVvars])
    probability_tolerance > 0 || error("Invalid probability tolerance $probability_tolerance")

    # merge possible duplicate terms in constraints
    v = IndexedVector(AffExpr, ccdata.numRVs)
    chanceconstr::Vector{ChanceConstr} = ccdata.chanceconstr
    for i in 1:length(chanceconstr)
        chanceconstr[i].ccexpr = merge_duplicates(chanceconstr[i].ccexpr, v, m)
    end

    if isa(m.solver,ECOS.ECOSSolver) && !linearize_objective
        reformulate_quadobj_to_conic = true
    end
    if reformulate_quadobj_to_conic
        # Use SOC representation of quadratic objective terms
        # for pure conic solvers (Mosek, ECOS).
        # Currently only diagonal terms supported
        qterms = length(m.obj.qvars1)
        quadobj::QuadExpr = m.obj
        if qterms != 0
            # we have quadratic terms
            # assume no duplicates for now
            for i in 1:qterms
                (quadobj.qvars1[i].col == quadobj.qvars2[i].col) || error("Only diagonal quadratic objective terms currently supported")
            end
            m.objSense == :Min || error("Only minimization is currently supported (this is easy to fix)")
            # x^2 <= t iff
            # exists y,z s.t. x^2 + y^2 <= z^2
            # z >= 0, y = z - 1, t = 2z-1, t >= 0
            @defVar(m, quadobj_y[1:qterms])
            @defVar(m, quadobj_z[1:qterms] >= 0)
            @defVar(m, quadobj_t[1:qterms] >= 0)
            @addConstraint(m, quadobj_soc[i=1:qterms], (quadobj.qvars1[i])^2 + quadobj_y[i]^2 <= quadobj_z[i]^2)
            @addConstraint(m, quadobj_eq1[i=1:qterms], quadobj_y[i] == quadobj_z[i] - 1)
            @addConstraint(m, quadobj_eq2[i=1:qterms], quadobj_t[i] == 2quadobj_z[i] - 1)
            @setObjective(m, m.objSense, quadobj.aff + sum{quadobj.qcoeffs[i]*quadobj_t[i], i=1:qterms})
        end
    end

    if method == :Reformulate
        @assert !linearize_objective # not supported
        # check that we have pure chance constraints
        no_uncertains || error("Cannot solve using reformulation, uncertain data are present")
        #display(ccdata.chanceconstr)

        for cc in chanceconstr
            ccexpr = cc.ccexpr
            nu = quantile(Normal(0,1),1-cc.with_probability)
            nterms = length(ccexpr.vars)
            # add auxiliary variables for variance of each term
            @defVar(m, varterm[1:nterms])
            @addConstraint(m, defvar[i=1:nterms], varterm[i] == getStdev(ccexpr.vars[i])*ccexpr.coeffs[i])
            @defVar(m, slackvar >= 0)
            # conic constraint
            @addConstraint(m, sum{ varterm[i]^2, i in 1:nterms } <= slackvar^2)
            if cc.sense == :(>=)
                @addConstraint(m, sum{getMean(ccexpr.vars[i])*ccexpr.coeffs[i], i=1:nterms} + nu*slackvar + ccexpr.constant <= 0)
            else
                @assert cc.sense == :(<=)
                @addConstraint(m, sum{getMean(ccexpr.vars[i])*ccexpr.coeffs[i], i=1:nterms} - nu*slackvar + ccexpr.constant >= 0)
            end
        end
        #println(m)

        return solve(m)

    else
        # check that we have pure chance constraints
        if no_uncertains
            solvecc_cuts(m, probability_tolerance, linearize_objective, debug, iteration_limit, objective_linearization_tolerance, lazy_constraints)
        else
            solverobustcc_cuts(m, probability_tolerance, linearize_objective, debug, iteration_limit, objective_linearization_tolerance, lazy_constraints)
        end
    end



end

function solvecc(m::Model;method=:Refomulate,linearize_objective::Bool=false,probability_tolerance=0.001,debug::Bool = false, iteration_limit::Int=60, objective_linearization_tolerance::Float64=1e-6, reformulate_quadobj_to_conic::Bool=false)

    Base.warn_once("solvecc is deprecated. Use solvechance instead!")
    return solvechance(m, method=method, linearize_objective=linearize_objective, probability_tolerance=probability_tolerance, iteration_limit=iteration_limit, objective_linearization_tolerance=objective_linearization_tolerance, reformulate_quadobj_to_conic=reformulate_quadobj_to_conic)
end


function solvecc_cuts(m::Model, probability_tolerance::Float64, linearize_objective::Bool, debug::Bool, iteration_limit::Int, objective_linearization_tolerance::Float64, lazy_constraints::Bool)

    ccdata = getCCData(m)

    has_integers = any(c-> c != :Cont, m.colCat)

    # set up slack variables and linear constraints
    nconstr = length(ccdata.chanceconstr)
    @defVar(m, slackvar[1:nconstr] >= 0)
    varterm = Dict()
    for i in 1:nconstr
        cc = ccdata.chanceconstr[i]
        ccexpr = cc.ccexpr
        nterms = length(ccexpr.vars)
        nu = quantile(Normal(0,1),1-cc.with_probability)
        if cc.sense == :(>=)
            @addConstraint(m, sum{getMean(ccexpr.vars[k])*ccexpr.coeffs[k], k=1:nterms} + nu*slackvar[i] + ccexpr.constant <= 0)
        else
            @addConstraint(m, sum{getMean(ccexpr.vars[k])*ccexpr.coeffs[k], k=1:nterms} - nu*slackvar[i] + ccexpr.constant >= 0)
        end
        # auxiliary variables
        @defVar(m, varterm[i][1:nterms])
        @addConstraint(m, defvar[k=1:nterms], varterm[i][k] == getStdev(ccexpr.vars[k])*ccexpr.coeffs[k])
    end

    # Optionally linearize quadratic objectives
    # Currently only diagonal terms supported
    qterms = length(m.obj.qvars1)
    quadobj::QuadExpr = m.obj
    if qterms != 0 && linearize_objective
        # we have quadratic terms
        # assume no duplicates for now
        for i in 1:qterms
            (quadobj.qvars1[i].col == quadobj.qvars2[i].col) || error("Only diagonal quadratic objective terms currently supported")
        end
        m.objSense == :Min || error("Only minimization is currently supported (this is easy to fix)")
        # qlinterm[i] >= qcoeffs[i]*qvars1[i]^2
        @defVar(m, qlinterm[1:qterms] >= 0)
        # it would help to add some initial linearizations
        @setObjective(m, m.objSense, quadobj.aff + sum{qlinterm[i], i=1:qterms})
    end

    function addcuts(cb)
        in_callback = isa(cb, MathProgBase.MathProgCallbackData)

        nviol = 0
        nviol_obj = 0
        
        # check violated chance constraints
        for i in 1:nconstr
            cc::ChanceConstr = ccdata.chanceconstr[i]
            mean = 0.0
            var = 0.0
            ccexpr = cc.ccexpr
            nterms = length(ccexpr.vars)
            for k in 1:nterms
                exprval = getValue(ccexpr.coeffs[k])
                mean += getMean(ccexpr.vars[k])*exprval
                var += getVar(ccexpr.vars[k])*exprval^2
            end
            mean += getValue(ccexpr.constant)
            if var < 1e-10 # corner case, need to handle carefully
                if cc.sense == :(<=) # actually this means strict inequality
                    satisfied_prob = (mean >= -1e-7) ? 0.0 : 1.0
                else
                    satisfied_prob = (mean <= 1e-7) ? 0.0 : 1.0
                end
            else
                if cc.sense == :(<=)
                    satisfied_prob = cdf(Normal(mean,sqrt(var)),0.0)
                else
                    satisfied_prob = 1-cdf(Normal(mean,sqrt(var)),0.0)
                end
            end
            if satisfied_prob <= cc.with_probability + probability_tolerance # feasibility tolerance
                # constraint is okay!
                continue
            else
                # check violation of quadratic constraint
                violation = var - getValue(slackvar[i])^2
                debug && println("Violated: $cc")
                debug && println("$satisfied_prob $mean $var")
                debug && println("VIOL $violation")
                nviol += 1
                # add a linearization
                if in_callback
                    @addLazyConstraint(cb, sum{ getValue(varterm[i][k])*varterm[i][k], k in 1:nterms} <= sqrt(var)*slackvar[i])
                else
                    @addConstraint(m, sum{ getValue(varterm[i][k])*varterm[i][k], k in 1:nterms} <= sqrt(var)*slackvar[i])
                end
            end
        end

        if linearize_objective
            # check violated objective linearizations
            for i in 1:qterms
                qval = quadobj.qcoeffs[i]*getValue(quadobj.qvars1[i])^2
                if getValue(qlinterm[i]) <= qval - objective_linearization_tolerance
                    debug && println("Linerization violation: ", qval - getValue(qlinterm[i]))
                    # add another linearization
                    nviol_obj += 1
                    if in_callback
                        @addLazyConstraint(cb, qlinterm[i] >= -qval + 2*quadobj.qcoeffs[i]*getValue(quadobj.qvars1[i])*quadobj.qvars1[i])
                    else
                        @addConstraint(m, qlinterm[i] >= -qval + 2*quadobj.qcoeffs[i]*getValue(quadobj.qvars1[i])*quadobj.qvars1[i])
                    end
                end
            end
        end
        return nviol, nviol_obj
    end
    
    do_lazy = has_integers && lazy_constraints
    (debug && !do_lazy) && println("Solving deterministic model")

    if do_lazy
        setLazyCallback(m, addcuts, fractional=true)
    end
    #tic()
    status = solve(m)

    if do_lazy
        #toc()
        return status
    end
    
    if status != :Optimal
        return status
    end


    niter = 0
    while niter < iteration_limit

        nviol, nviol_obj = addcuts(nothing)
 
        if nviol == 0 && nviol_obj == 0
            println("Done after $niter iterations")
            #toc()
            return :Optimal
        else
            println("Iteration $niter: $nviol constraint violations, $nviol_obj objective linearization violations")
        end
        status = solve(m)
        if status != :Optimal
            return status
        end
        niter += 1
    end

    return :UserLimit # hit iteration limit


end

function solverobustcc_cuts(m::Model, probability_tolerance::Float64, linearize_objective::Bool, debug::Bool, iteration_limit::Int, objective_linearization_tolerance::Float64, lazy_constraints::Bool)

    ccdata = getCCData(m)

    nconstr = length(ccdata.chanceconstr)

    # Optionally linearize quadratic objectives
    # Currently only diagonal terms supported
    # TODO: deduplicate code with solvecc_cuts
    qterms = length(m.obj.qvars1)
    quadobj::QuadExpr = m.obj
    if qterms != 0 && linearize_objective
        # we have quadratic terms
        # assume no duplicates for now
        for i in 1:qterms
            (quadobj.qvars1[i].col == quadobj.qvars2[i].col) || error("Only diagonal quadratic objective terms currently supported")
        end
        m.objSense == :Min || error("Only minimization is currently supported (this is easy to fix)")
        # qlinterm[i] >= qcoeffs[i]*qvars1[i]^2
        @defVar(m, qlinterm[1:qterms] >= 0)
        # it would help to add some initial linearizations
        @setObjective(m, m.objSense, quadobj.aff + sum{qlinterm[i], i=1:qterms})
    end


    # prepare uncertainty set data
    # nominal value is taken to be the center of the given interval
    # and allowed deviation is taken as half the interval length
    means_nominal = zeros(ccdata.numRVs)
    means_deviation = zeros(ccdata.numRVs)
    vars_nominal = zeros(ccdata.numRVs)
    vars_deviation = zeros(ccdata.numRVs)

    for i in 1:ccdata.numRVs
        rvmean = ccdata.RVmeans[i]
        rvvar = ccdata.RVvars[i]
        if isa(rvmean,Real)
            means_nominal[i] = rvmean
            # deviation is zero
        else
            lb,ub = rvmean
            if ub < lb
                error("Variable $(ccdata.RVnames[i]) has invalid mean interval [$lb, $ub]")
            end
            means_nominal[i] = (lb+ub)/2
            means_deviation[i] = (ub-lb)/2
        end
        if isa(rvvar,Real)
            vars_nominal[i] = rvvar
        else
            lb,ub = rvvar
            if ub < lb
                error("Variable $(ccdata.RVnames[i]) has invalid variance interval [$lb, $ub]")
            end
            vars_nominal[i] = (lb+ub)/2
            vars_deviation[i] = (ub-lb)/2
        end

    end

    debug && println("Nominal means: ", means_nominal)
    debug && println("Mean deviations: ", means_deviation)
    debug && println("Nominal variance: ", vars_nominal)
    debug && println("Variance deviations: ", vars_deviation)

    # TODO: special handling for quadratic objectives
    
    debug && println("Solving deterministic model")

    status = solve(m)
    
    if status != :Optimal
        return status
    end

    niter = 0
    while niter < iteration_limit

        nviol = 0
        # check violated chance constraints
        for i in 1:nconstr
            cc::ChanceConstr = ccdata.chanceconstr[i]
            ccexpr = cc.ccexpr
            nterms = length(ccexpr.vars)
            nu = quantile(Normal(0,1),1-cc.with_probability)
            # sort to determine worst case
            meanvals = zeros(nterms)
            varvals = zeros(nterms)
            deviation_sign = zeros(nterms)
            nominal_mean = getValue(ccexpr.constant)
            nominal_var = 0.0
            for k in 1:nterms
                exprval = getValue(ccexpr.coeffs[k])
                #debug && println(ccexpr.coeffs[k], " => ", exprval)
                idx = ccexpr.vars[k].idx
                if cc.sense == :(<=) && exprval < 0.0
                    deviation_sign[k] = 1.0
                elseif cc.sense == :(<=) && exprval >= 0.0
                    deviation_sign[k] = -1.0
                elseif cc.sense == :(>=) && exprval < 0.0
                    deviation_sign[k] = -1.0
                else
                    deviation_sign[k] = 1.0
                end
                meanvals[k] = means_deviation[idx]*abs(exprval)
                varvals[k] = vars_deviation[idx]*exprval^2
                nominal_mean += means_nominal[idx]*exprval
                nominal_var += vars_nominal[idx]*exprval^2
            end
            sorted_mean_idx = sortperm(meanvals,rev=true)
            sorted_var_idx = sortperm(varvals,rev=true)
            @assert cc.uncertainty_budget_mean >= 0
            @assert cc.uncertainty_budget_variance >= 0
            if cc.sense == :(<=)
                worst_mean = nominal_mean - sum(meanvals[sorted_mean_idx[1:cc.uncertainty_budget_mean]])
            else
                worst_mean = nominal_mean + sum(meanvals[sorted_mean_idx[1:cc.uncertainty_budget_mean]])
            end
            worst_var = nominal_var + sum(varvals[sorted_var_idx[1:cc.uncertainty_budget_variance]])
            if worst_var < 1e-10 # corner case, need to handle carefully
                if cc.sense == :(<=) # actually this means strict inequality
                    satisfied_prob = (worst_mean >= -1e-7) ? 0.0 : 1.0
                else
                    satisfied_prob = (worst_mean <= 1e-7) ? 0.0 : 1.0
                end
            else
                if cc.sense == :(<=)
                    satisfied_prob = cdf(Normal(worst_mean,sqrt(worst_var)),0.0)
                else
                    satisfied_prob = 1-cdf(Normal(worst_mean,sqrt(worst_var)),0.0)
                end
            end
            if satisfied_prob <= cc.with_probability
                # constraint is okay!
                continue
            else
                if satisfied_prob >= cc.with_probability + probability_tolerance
                    debug && println(cc)
                    debug && println("nominal_mean = $nominal_mean")
                    debug && println("nominal_var = $nominal_var")
                    debug && println("worst_mean = $worst_mean")
                    debug && println("meanvals = $meanvals")
                    debug && println("worst_var = $worst_var")
                    debug && println("satisfied w.p. $satisfied_prob")
                    debug && println("nu = $nu")
                    debug && println("Extreme indices (mean) ", sorted_mean_idx[1:cc.uncertainty_budget_mean])
                    debug && println("Extreme indices (var) ", sorted_var_idx[1:cc.uncertainty_budget_variance])
                    debug && println(ccexpr.constant, " ===> ", getValue(ccexpr.constant))
                    debug && println(ccexpr.coeffs[1], " ===> ", getValue(ccexpr.coeffs[1]))
                    debug && println("VIOL ", 100*(satisfied_prob - cc.with_probability), "%")
                    nviol += 1
                    var_coeffs = [vars_nominal[ccexpr.vars[k].idx] for k in 1:nterms]
                    for k in 1:cc.uncertainty_budget_variance
                        var_coeffs[sorted_var_idx[k]] += vars_deviation[ccexpr.vars[sorted_var_idx[k]].idx]
                    end
                    debug && println("var_coeffs = $var_coeffs")


                    # add a linearization
                    # f(x') + f'(x')(x-x') <= 0
                    if cc.sense == :(>=)
                        @addConstraint(m, ccexpr.constant + nu*sqrt(worst_var) + 
                            sum{means_nominal[ccexpr.vars[k].idx]*(ccexpr.coeffs[k]), k=1:nterms} +
                            sum{sign(getValue(ccexpr.coeffs[sorted_mean_idx[r]]))*means_deviation[ccexpr.vars[sorted_mean_idx[r]].idx]*(ccexpr.coeffs[sorted_mean_idx[r]]), r=1:cc.uncertainty_budget_mean} +
                            (nu/sqrt(worst_var))*sum{var_coeffs[k]*getValue(ccexpr.coeffs[k])*(ccexpr.coeffs[k]-getValue(ccexpr.coeffs[k])),k=1:nterms}  <= 0)
                        debug && println("ADDED: ", m.linconstr[end])
                    else
                        @addConstraint(m, ccexpr.constant - nu*sqrt(worst_var) + 
                            sum{means_nominal[ccexpr.vars[k].idx]*(ccexpr.coeffs[k]), k=1:nterms} -
                            sum{sign(getValue(ccexpr.coeffs[sorted_mean_idx[r]]))*means_deviation[ccexpr.vars[sorted_mean_idx[r]].idx]*(ccexpr.coeffs[sorted_mean_idx[r]]), r=1:cc.uncertainty_budget_mean} -
                            (nu/sqrt(worst_var))*sum{var_coeffs[k]*getValue(ccexpr.coeffs[k])*(ccexpr.coeffs[k]-getValue(ccexpr.coeffs[k])), k=1:nterms}  >= 0)
                        debug && println("ADDED: ", m.linconstr[end])

                    end

                end
            end
        end

        # check violated objective linearizations
        nviol_obj = 0
        if linearize_objective
            for i in 1:qterms
                qval = quadobj.qcoeffs[i]*getValue(quadobj.qvars1[i])^2
                if getValue(qlinterm[i]) <= qval - objective_linearization_tolerance
                    # add another linearization
                    nviol_obj += 1
                    @addConstraint(m, qlinterm[i] >= -qval + 2*quadobj.qcoeffs[i]*getValue(quadobj.qvars1[i])*quadobj.qvars1[i])
                end
            end
        end

        if nviol == 0 && nviol_obj == 0
            println("Done after $niter iterations")
            return :Optimal
        else
            println("Iteration $niter: $nviol constraint violations, $nviol_obj objective linearization violations")
        end
        status = solve(m)
        if status != :Optimal
            return status
        end
        niter += 1
    end

    return :UserLimit # hit iteration limit


end

# We can remove the code for IndexedVector once JuMP issue #340 is resolved.
type IndexedVector{T}
    elts::Vector{T}
    nzidx::Vector{Int}
    nnz::Int
    empty::BitArray{1}
end

IndexedVector{T}(::Type{T},n::Integer) = IndexedVector(zeros(T,n),zeros(Int,n),0,trues(n))

function addelt!{T}(v::IndexedVector{T},i::Integer,val::T)
    if v.empty[i]  # new index
        v.elts[i] = val
        v.nzidx[v.nnz += 1] = i
        v.empty[i] = false
    else
        v.elts[i] += val
    end
    return nothing
end

function Base.empty!{T}(v::IndexedVector{T})
    elts = v.elts
    nzidx = v.nzidx
    empty = v.empty
    for i in 1:v.nnz
        elts[nzidx[i]] = zero(T)
        empty[nzidx[i]] = true
    end
    v.nnz = 0
end

Base.length(v::IndexedVector) = length(v.elts)
function Base.resize!(v::IndexedVector, n::Integer)
    if n > length(v)
        @assert v.nnz == 0 # only resize empty vector
        resize!(v.elts, n)
        fill!(v.elts,0)
        resize!(v.nzidx, n)
        resize!(v.empty, n)
        fill!(v.empty,true)
    end
end

# Adapted from JuMP:
# returns a GenericAffExpr with terms merged
# assume that v is zero'd
function merge_duplicates{CoefType}(aff::JuMP.GenericAffExpr{CoefType,IndepNormal}, v::IndexedVector{CoefType}, m::Model)
    ccdata = getCCData(m)

    resize!(v, ccdata.numRVs)
    for ind in 1:length(aff.coeffs)
        var = aff.vars[ind]
        is(var.m, m) || error("Variable does not belong to this model")
        addelt!(v, aff.vars[ind].idx, aff.coeffs[ind])
    end
    vars = Array(IndepNormal,v.nnz)
    coeffs = Array(CoefType,v.nnz)
    for i in 1:v.nnz
        idx = v.nzidx[i]
        vars[i] = IndepNormal(m,idx)
        coeffs[i] = v.elts[idx]
    end
    empty!(v)

    return JuMP.GenericAffExpr(vars, coeffs, aff.constant) # do we need to eliminate duplicates in the constant part also?

end
