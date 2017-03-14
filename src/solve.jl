import MathProgBase

Φinv(x) = quantile(Normal(0,1),x)

function solvehook(m::Model; suppress_warnings=false, method=:Refomulate,linearize_objective::Bool=false,probability_tolerance=0.001,debug::Bool = false, iteration_limit::Int=60, objective_linearization_tolerance::Float64=1e-6, reformulate_quadobj_to_conic::Bool=false, lazy_constraints::Bool=false, silent::Bool=false)
    @assert method == :Reformulate || method == :Cuts

    ccdata = getCCData(m)
    no_uncertains = all([isa(x,Real) for x in ccdata.RVmeans]) && all([isa(x,Real) for x in ccdata.RVvars])
    probability_tolerance > 0 || error("Invalid probability tolerance $probability_tolerance")

    # merge possible duplicate terms in constraints
    v = IndexedVector(AffExpr, ccdata.numRVs)
    counts = IndexedVector(Int, ccdata.numRVs)
    chanceconstr::Vector{ChanceConstr} = ccdata.chanceconstr
    for i in 1:length(chanceconstr)
        chanceconstr[i].ccexpr = merge_duplicates(chanceconstr[i].ccexpr, v, counts, m)
    end

    twosidechance::Vector{TwoSideChanceConstr} = ccdata.twosidechanceconstr
    for i in 1:length(twosidechance)
        twosidechance[i].ccexpr = merge_duplicates(twosidechance[i].ccexpr, v, counts, m)
    end
    has_twoside = length(twosidechance) > 0

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
            quadobj_y = @variable(m, [1:qterms])
            quadobj_z = @variable(m, [1:qterms], lowerbound = 0)
            quadobj_t = @variable(m, [1:qterms], lowerbound = 0)
            quadobj_soc = @constraint(m, [i=1:qterms], (quadobj.qvars1[i])^2 + quadobj_y[i]^2 <= quadobj_z[i]^2)
            quadobj_eq1 = @constraint(m, [i=1:qterms], quadobj_y[i] == quadobj_z[i] - 1)
            quadobj_eq2 = @constraint(m, [i=1:qterms], quadobj_t[i] == 2quadobj_z[i] - 1)
            @objective(m, m.objSense, quadobj.aff + sum(quadobj.qcoeffs[i]*quadobj_t[i] for i=1:qterms))
        end
    end

    if method == :Reformulate
        @assert !linearize_objective # not supported
        # check that we have pure chance constraints
        no_uncertains || error("Cannot solve using reformulation, uncertain data are present")
        #display(ccdata.chanceconstr)

        for cc in chanceconstr
            ccexpr = cc.ccexpr
            nu = quantile(Normal(0,1),cc.with_probability)
            nterms = length(ccexpr.vars)
            # case when reformulation results in a linear constraint
            coeffs = ccexpr.coeffs
            if nterms == 0
                # actually not a chance constraint
                if cc.sense == :(<=)
                    @constraint(m, ccexpr.constant <= 0)
                else
                    @constraint(m, ccexpr.constate >= 0)
                end
                continue
            elseif all(ex -> isequal(ex, coeffs[1]), coeffs)
                slackvar = @variable(m, lowerbound = 0)
                @constraint(m, slackvar >= ccexpr.coeffs[1])
                @constraint(m, slackvar >= -ccexpr.coeffs[1])
                sumvar = sum([getvariance(ccexpr.vars[i]) for i in 1:nterms])
                sqrtsumvar = sqrt(sumvar)
                if cc.sense == :(<=)
                    @constraint(m, sum(getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i=1:nterms) + nu*sqrtsumvar*slackvar + ccexpr.constant <= 0)
                else
                    @assert cc.sense == :(>=)
                    @constraint(m, sum(getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i=1:nterms) - nu*sqrtsumvar*slackvar + ccexpr.constant >= 0)
                end
                continue
            end
            # add auxiliary variables for variance of each term
            varterm = @variable(m, [1:nterms])
            defvar = @constraint(m, [i=1:nterms], varterm[i] == getstdev(ccexpr.vars[i])*ccexpr.coeffs[i])
            slackvar = @variable(m, lowerbound = 0)
            # conic constraint
            @constraint(m, sum( varterm[i]^2 for i in 1:nterms ) <= slackvar^2)
            if cc.sense == :(<=)
                @constraint(m, sum(getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i=1:nterms) + nu*slackvar + ccexpr.constant <= 0)
            else
                @assert cc.sense == :(>=)
                @constraint(m, sum(getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i=1:nterms) - nu*slackvar + ccexpr.constant >= 0)
            end
        end

        lbvar = @variable(m, [1:length(twosidechance)])
        ubvar = @variable(m, [1:length(twosidechance)])
        t = @variable(m, [1:length(twosidechance)], lowerbound = 0)
        for k in 1:length(twosidechance)
            cc = twosidechance[k]
            ccexpr = cc.ccexpr
            cc.lb -= ccexpr.constant
            cc.ub -= ccexpr.constant
            ccexpr.constant = AffExpr()
            nterms = length(ccexpr.vars)
            coeffs = ccexpr.coeffs
            if nterms == 0
                # actually not a chance constraint
                @constraint(m, cc.lb ≤ 0)
                @constraint(m, cc.ub ≥ 0)
            elseif all(ex -> isequal(ex, coeffs[1]), coeffs)
                @constraint(m, t[k] >= ccexpr.coeffs[1])
                @constraint(m, t[k] >= -ccexpr.coeffs[1])
                sumvar = sum([getvariance(ccexpr.vars[i]) for i in 1:nterms])
                sqrtsumvar = sqrt(sumvar)
                ϵ = 1-cc.with_probability
                @constraint(m, lbvar[k] == cc.lb - sum( getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i = 1:nterms))
                @constraint(m, ubvar[k] == cc.ub - sum( getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i = 1:nterms))
                @constraint(m, lbvar[k] <= Φinv(ϵ)*sqrtsumvar*t[k])
                @constraint(m, ubvar[k] >= Φinv(1-ϵ)*sqrtsumvar*t[k])
                if cc.approx == "1.25"
                    @constraint(m, ubvar[k] - lbvar[k] >= -2*Φinv(ϵ/2)*sqrtsumvar*t[k])
                end
                continue
            end
            
            # add auxiliary variables for variance of each term
            varterm = @variable(m, [1:nterms])
            defvar = @constraint(m, [i=1:nterms], varterm[i] == getstdev(ccexpr.vars[i])*ccexpr.coeffs[i])
            # conic constraint
            @constraint(m, sum( varterm[i]^2 for i in 1:nterms ) <= t[k]^2)
            ϵ = 1-cc.with_probability
            @constraint(m, lbvar[k] == cc.lb - sum( getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i = 1:nterms))
            @constraint(m, ubvar[k] == cc.ub - sum( getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i = 1:nterms))
            # (lbvar, ubvar, t) ∈ \bar S_ϵ
            # lbvar/t ≤ Φ^{-1}(ϵ)
            # ubvar/t ≥ Φ^{-1}(1-ϵ)
            # ubvar/t - lbvar/t ≥ -2Φ^{-1}(ϵ/2)
            @constraint(m, lbvar[k] ≤ Φinv(ϵ)*t[k])
            @constraint(m, ubvar[k] ≥ Φinv(1-ϵ)*t[k])
            if cc.approx == "1.25"
                @constraint(m, ubvar[k] - lbvar[k] ≥ -2*Φinv(ϵ/2)*t[k])
            end
        end
        #println(m)

        status = solve(m,suppress_warnings=suppress_warnings, ignore_solve_hook=true)

        #=
        # debug two-sided constraints
        X = Float64[]
        Y = Float64[]
        for k in 1:length(twosidechance)
            l = getvalue(lbvar[k])/getvalue(t[k])
            u = getvalue(ubvar[k])/getvalue(t[k])
            push!(X,l)
            push!(Y,u)
            println("$k $l $u $(Φ(u) - Φ(l))")
        end
        println("done")
        =#

        return status

    else
        #has_twoside && error("Two-sided chance constraints are not currently supported with method = :Cuts. Use method = :Reformuate instead.")
        # check that we have pure chance constraints
        if no_uncertains
            solvecc_cuts(m, suppress_warnings, probability_tolerance, linearize_objective, debug, iteration_limit, objective_linearization_tolerance, lazy_constraints, silent)
        else
            solverobustcc_cuts(m, suppress_warnings, probability_tolerance, linearize_objective, debug, iteration_limit, objective_linearization_tolerance, lazy_constraints, silent)
        end
    end



end

function solvecc_cuts(m::Model, suppress_warnings::Bool, probability_tolerance::Float64, linearize_objective::Bool, debug::Bool, iteration_limit::Int, objective_linearization_tolerance::Float64, lazy_constraints::Bool, silent::Bool)

    ccdata = getCCData(m)

    has_integers = any(c-> c != :Cont, m.colCat)

    # set up slack variables for both one sided and two sided chance constraints
    nchanceconstr = length(ccdata.chanceconstr)
    ntwosidechanceconstr = length(ccdata.twosidechanceconstr)
    slackvar = @variable(m, [1:nchanceconstr], lowerbound = 0)
    t = @variable(m, [1:ntwosidechanceconstr], lowerbound = 0)
    chancevarterm = Dict()
    twosidechancevarterm = Dict()
    linearizechanceconstr = Int[]
    linearizetwosidechanceconstr = Int[]
    
    # one sided chance constraints case
    for i in 1:nchanceconstr
        cc = ccdata.chanceconstr[i]
        ccexpr = cc.ccexpr
        nterms = length(ccexpr.vars)
        nu = quantile(Normal(0,1),cc.with_probability)
        coeffs = ccexpr.coeffs
        if all(ex -> isequal(ex, coeffs[1]), coeffs)
            @constraint(m, slackvar[i] >= ccexpr.coeffs[1])
            @constraint(m, slackvar[i] >= -ccexpr.coeffs[1])
            sumvar = sum([getvariance(ccexpr.vars[i]) for i in 1:nterms])
            sqrtsumvar = sqrt(sumvar)
            if cc.sense == :(<=)
                @constraint(m, sum(getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i=1:nterms) + nu*sqrtsumvar*slackvar[i] + ccexpr.constant <= 0)
            else
                @assert cc.sense == :(>=)
                @constraint(m, sum(getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i=1:nterms) - nu*sqrtsumvar*slackvar[i] + ccexpr.constant >= 0)
            end
        else
            push!(linearizechanceconstr, i)
            if cc.sense == :(<=)
                @constraint(m, sum(getmean(ccexpr.vars[k])*ccexpr.coeffs[k] for k=1:nterms) + nu*slackvar[i] + ccexpr.constant <= 0)
            else
                @constraint(m, sum(getmean(ccexpr.vars[k])*ccexpr.coeffs[k] for k=1:nterms) - nu*slackvar[i] + ccexpr.constant >= 0)
            end
            # auxiliary variables
            chancevarterm[i] = @variable(m, [1:nterms])
            defvar = @constraint(m, [k=1:nterms], chancevarterm[i][k] == getstdev(ccexpr.vars[k])*ccexpr.coeffs[k])
        end
    end
    
    # Two sided chance constraints with :Cuts option
    lbvar = @variable(m, [1:ntwosidechanceconstr])
    ubvar = @variable(m, [1:ntwosidechanceconstr])

    for k in 1:ntwosidechanceconstr
        cc = ccdata.twosidechanceconstr[k]
        ccexpr = cc.ccexpr
        cc.lb -= ccexpr.constant
        cc.ub -= ccexpr.constant
        ccexpr.constant = AffExpr()
        nterms = length(ccexpr.vars)
        coeffs = ccexpr.coeffs
        if all(ex -> isequal(ex, coeffs[1]), coeffs)
            @constraint(m, t[k] >= ccexpr.coeffs[1])
            @constraint(m, t[k] >= -ccexpr.coeffs[1])
            sumvar = sum([getvariance(ccexpr.vars[i]) for i in 1:nterms])
            sqrtsumvar = sqrt(sumvar)
            ϵ = 1-cc.with_probability
            @constraint(m, lbvar[k] == cc.lb - sum( getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i = 1:nterms))
            @constraint(m, ubvar[k] == cc.ub - sum( getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i = 1:nterms))
            @constraint(m, lbvar[k] <= Φinv(ϵ)*sqrtsumvar*t[k])
            @constraint(m, ubvar[k] >= Φinv(1-ϵ)*sqrtsumvar*t[k])
            if cc.approx == "1.25"
                @constraint(m, ubvar[k] - lbvar[k] >= -2*Φinv(ϵ/2)*sqrtsumvar*t[k])
            end
        else
            push!(linearizetwosidechanceconstr, k)
            twosidechancevarterm[k] = @variable(m, [1:nterms])
            defvar = @constraint(m, [i=1:nterms], twosidechancevarterm[k][i] == getstdev(ccexpr.vars[i])*ccexpr.coeffs[i])
            ϵ = 1-cc.with_probability
            @constraint(m, lbvar[k] == cc.lb - sum( getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i = 1:nterms))
            @constraint(m, ubvar[k] == cc.ub - sum( getmean(ccexpr.vars[i])*ccexpr.coeffs[i] for i = 1:nterms))
            @constraint(m, lbvar[k] ≤ Φinv(ϵ)*t[k])
            @constraint(m, ubvar[k] ≥ Φinv(1-ϵ)*t[k])
            if cc.approx == "1.25"
                @constraint(m, ubvar[k] - lbvar[k] ≥ -2*Φinv(ϵ/2)*t[k])
            end
        end
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
        qlinterm = @variable(m, [1:qterms], lowerbound = 0)
        # it would help to add some initial linearizations
        @objective(m, m.objSense, quadobj.aff + sum(qlinterm[i] for i=1:qterms))
    end

    function addcuts(cb)
        in_callback = isa(cb, MathProgBase.MathProgCallbackData)

        nviol = 0
        nviol_obj = 0
        
        # check violated chance constraints
        for i in linearizechanceconstr
            cc::ChanceConstr = ccdata.chanceconstr[i]
            mean = 0.0
            var = 0.0
            ccexpr = cc.ccexpr
            nterms = length(ccexpr.vars)
            for k in 1:nterms
                exprval = getvalue(ccexpr.coeffs[k])
                mean += getmean(ccexpr.vars[k])*exprval
                var += getvariance(ccexpr.vars[k])*exprval^2
            end
            mean += getvalue(ccexpr.constant)
            if var < 1e-10 # corner case, need to handle carefully
                if cc.sense == :(<=)
                    satisfied_prob = (mean <= 1e-7) ? 1.0 : 0.0
                else
                    satisfied_prob = (mean >= -1e-7) ? 1.0 : 0.0
                end
            else
                if cc.sense == :(<=)
                    satisfied_prob = cdf(Normal(mean,sqrt(var)),0.0)
                else
                    satisfied_prob = 1-cdf(Normal(mean,sqrt(var)),0.0)
                end
            end
            if satisfied_prob >= cc.with_probability - probability_tolerance # feasibility tolerance
                # constraint is okay!
                continue
            else
                # check violation of quadratic constraint
                violation = var - getvalue(slackvar[i])^2
                debug && println("Violated: $cc")
                debug && println("$satisfied_prob $mean $var")
                debug && println("VIOL $violation")
                nviol += 1
                # add a linearization
                if in_callback
                    @lazyconstraint(cb, sum( getvalue(chancevarterm[i][k])*chancevarterm[i][k] for k in 1:nterms) <= sqrt(var)*slackvar[i])
                else
                    @constraint(m, sum( getvalue(chancevarterm[i][k])*chancevarterm[i][k] for k in 1:nterms) <= sqrt(var)*slackvar[i])
                end
            end
        end
        
        # check violated two sided chance constraint
        for i in linearizetwosidechanceconstr
            cc::TwoSideChanceConstr = ccdata.twosidechanceconstr[i]
            mean = 0.0
            var = 0.0
            ccexpr = cc.ccexpr
            meanub = 0.0
            meanlb = 0.0
            nterms = length(ccexpr.vars)
            for k in 1:nterms
                exprval = getvalue(ccexpr.coeffs[k])
                mean += getmean(ccexpr.vars[k])*exprval
                var += getvariance(ccexpr.vars[k])*exprval^2
            end
            mean += getvalue(ccexpr.constant)
            meanlb = mean - getvalue(cc.lb)
            meanub = mean - getvalue(cc.ub)
            if var < 1e-10 # corner case, need to handle carefully
                satisfied_prob_ub = (meanub <= 1e-7) ? 1.0 : 0.0
                satisfied_prob_lb = (meanlb >= -1e-7) ? 1.0 : 0.0
            else
                satisfied_prob_ub = cdf(Normal(meanub,sqrt(var)),0.0)
                satisfied_prob_lb = 1-cdf(Normal(meanlb,sqrt(var)),0.0)
            end
            if satisfied_prob_ub >= cc.with_probability - probability_tolerance  && satisfied_prob_lb >= cc.with_probability - probability_tolerance
                # constraint is okay!
                continue
            else
                # check violation of quadratic constraint
                violation = var - getvalue(t[i])^2
                debug && println("Violated: $cc")
                debug && println("$satisfied_prob $mean $var")
                debug && println("VIOL $violation")
                nviol += 1
                # add a linearization
                if in_callback
                    @lazyconstraint(cb, sum( getvalue(twosidechancevarterm[i][k])*twosidechancevarterm[i][k] for k in 1:nterms) <= sqrt(var)*t[i])
                else
                    @constraint(m, sum( getvalue(twosidechancevarterm[i][k])*twosidechancevarterm[i][k] for k in 1:nterms) <= sqrt(var)*t[i])
                end
            end
        end


        if linearize_objective
            # check violated objective linearizations
            for i in 1:qterms
                qval = quadobj.qcoeffs[i]*getvalue(quadobj.qvars1[i])^2
                if getvalue(qlinterm[i]) <= qval - objective_linearization_tolerance
                    debug && println("Linerization violation: ", qval - getvalue(qlinterm[i]))
                    # add another linearization
                    nviol_obj += 1
                    if in_callback
                        @lazyconstraint(cb, qlinterm[i] >= -qval + 2*quadobj.qcoeffs[i]*getvalue(quadobj.qvars1[i])*quadobj.qvars1[i])
                    else
                        @constraint(m, qlinterm[i] >= -qval + 2*quadobj.qcoeffs[i]*getvalue(quadobj.qvars1[i])*quadobj.qvars1[i])
                    end
                end
            end
        end
        return nviol, nviol_obj
    end
    
    do_lazy = has_integers && lazy_constraints
    (debug && !do_lazy) && println("Solving deterministic model")

    if do_lazy
        addlazycallback(m, addcuts, fractional=true)
    end
    #tic()
    status = solve(m, suppress_warnings=suppress_warnings, ignore_solve_hook=true)

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
            silent || println("Done after $niter iterations")
            #toc()
            return :Optimal
        else
            silent || println("Iteration $niter: $nviol constraint violations, $nviol_obj objective linearization violations")
        end
        status = solve(m, suppress_warnings=suppress_warnings, ignore_solve_hook=true)
        if status != :Optimal
            return status
        end
        niter += 1
    end

    return :UserLimit # hit iteration limit


end

function solverobustcc_cuts(m::Model, suppress_warnings::Bool, probability_tolerance::Float64, linearize_objective::Bool, debug::Bool, iteration_limit::Int, objective_linearization_tolerance::Float64, lazy_constraints::Bool, silent::Bool)

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
        qlinterm = @variable(m, [1:qterms], lowerbound = 0)
        # it would help to add some initial linearizations
        @objective(m, m.objSense, quadobj.aff + sum(qlinterm[i] for i=1:qterms))
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

    status = solve(m, suppress_warnings=suppress_warnings, ignore_solve_hook=true)
    
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
            nu = quantile(Normal(0,1),cc.with_probability)
            # sort to determine worst case
            meanvals = zeros(nterms)
            varvals = zeros(nterms)
            deviation_sign = zeros(nterms)
            nominal_mean = getvalue(ccexpr.constant)
            nominal_var = 0.0
            for k in 1:nterms
                exprval = getvalue(ccexpr.coeffs[k])
                #debug && println(ccexpr.coeffs[k], " => ", exprval)
                idx = ccexpr.vars[k].idx
                if cc.sense == :(>=) && exprval < 0.0
                    deviation_sign[k] = 1.0
                elseif cc.sense == :(>=) && exprval >= 0.0
                    deviation_sign[k] = -1.0
                elseif cc.sense == :(<=) && exprval < 0.0
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
            if cc.sense == :(>=)
                worst_mean = nominal_mean - sum(meanvals[sorted_mean_idx[1:cc.uncertainty_budget_mean]])
            else
                worst_mean = nominal_mean + sum(meanvals[sorted_mean_idx[1:cc.uncertainty_budget_mean]])
            end
            worst_var = nominal_var + sum(varvals[sorted_var_idx[1:cc.uncertainty_budget_variance]])
            if worst_var < 1e-10 # corner case, need to handle carefully
                if cc.sense == :(<=) # var = 0
                    satisfied_prob = (worst_mean <= 1e-7) ? 1.0 : 0.0
                else
                    satisfied_prob = (worst_mean >= -1e-7) ? 1.0 : 0.0
                end
            else
                if cc.sense == :(<=)
                    satisfied_prob = cdf(Normal(worst_mean,sqrt(worst_var)),0.0)
                else
                    satisfied_prob = 1-cdf(Normal(worst_mean,sqrt(worst_var)),0.0)
                end
            end
            if satisfied_prob >= cc.with_probability
                # constraint is okay!
                continue
            else
                if satisfied_prob <= cc.with_probability - probability_tolerance
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
                    debug && println(ccexpr.constant, " ===> ", getvalue(ccexpr.constant))
                    debug && println(ccexpr.coeffs[1], " ===> ", getvalue(ccexpr.coeffs[1]))
                    debug && println("VIOL ", 100*(cc.with_probability - satisfied_prob), "%")
                    nviol += 1
                    var_coeffs = [vars_nominal[ccexpr.vars[k].idx] for k in 1:nterms]
                    for k in 1:cc.uncertainty_budget_variance
                        var_coeffs[sorted_var_idx[k]] += vars_deviation[ccexpr.vars[sorted_var_idx[k]].idx]
                    end
                    debug && println("var_coeffs = $var_coeffs")


                    # add a linearization
                    # f(x') + f'(x')(x-x') <= 0
                    if cc.sense == :(<=)
                        @constraint(m, ccexpr.constant + nu*sqrt(worst_var) +
                            sum(means_nominal[ccexpr.vars[k].idx]*(ccexpr.coeffs[k]) for k=1:nterms) +
                            sum(sign(getvalue(ccexpr.coeffs[sorted_mean_idx[r]]))*means_deviation[ccexpr.vars[sorted_mean_idx[r]].idx]*(ccexpr.coeffs[sorted_mean_idx[r]]) for r=1:cc.uncertainty_budget_mean) +
                            (nu/sqrt(worst_var))*sum(var_coeffs[k]*getvalue(ccexpr.coeffs[k])*(ccexpr.coeffs[k]-getvalue(ccexpr.coeffs[k])) for k=1:nterms)  <= 0)
                        debug && println("ADDED: ", m.linconstr[end])
                    else
                        @constraint(m, ccexpr.constant - nu*sqrt(worst_var) +
                            sum(means_nominal[ccexpr.vars[k].idx]*(ccexpr.coeffs[k]) for k=1:nterms) -
                            sum(sign(getvalue(ccexpr.coeffs[sorted_mean_idx[r]]))*means_deviation[ccexpr.vars[sorted_mean_idx[r]].idx]*(ccexpr.coeffs[sorted_mean_idx[r]]) for r=1:cc.uncertainty_budget_mean) -
                            (nu/sqrt(worst_var))*sum(var_coeffs[k]*getvalue(ccexpr.coeffs[k])*(ccexpr.coeffs[k]-getvalue(ccexpr.coeffs[k])) for k=1:nterms)  >= 0)
                        debug && println("ADDED: ", m.linconstr[end])

                    end

                end
            end
        end

        # check violated objective linearizations
        nviol_obj = 0
        if linearize_objective
            for i in 1:qterms
                qval = quadobj.qcoeffs[i]*getvalue(quadobj.qvars1[i])^2
                if getvalue(qlinterm[i]) <= qval - objective_linearization_tolerance
                    # add another linearization
                    nviol_obj += 1
                    @constraint(m, qlinterm[i] >= -qval + 2*quadobj.qcoeffs[i]*getvalue(quadobj.qvars1[i])*quadobj.qvars1[i])
                end
            end
        end

        if nviol == 0 && nviol_obj == 0
            silent || println("Done after $niter iterations")
            return :Optimal
        else
            silent || println("Iteration $niter: $nviol constraint violations, $nviol_obj objective linearization violations")
        end
        status = solve(m, suppress_warnings=suppress_warnings, ignore_solve_hook=true)
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

function addelt!(v::IndexedVector{AffExpr},i::Integer,val::AffExpr)
    if v.empty[i]  # new index
        v.elts[i] = copy(val)
        v.nzidx[v.nnz += 1] = i
        v.empty[i] = false
    else
        append!(v.elts[i], val)
    end
    return nothing
end

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
function merge_duplicates{CoefType <: JuMP.GenericAffExpr}(aff::JuMP.GenericAffExpr{CoefType,IndepNormal}, v::IndexedVector{CoefType}, counts::IndexedVector{Int}, m::Model)
    ccdata = getCCData(m)

    resize!(v, ccdata.numRVs)
    resize!(counts, ccdata.numRVs)
    # do a first pass to precompute sizes
    for ind in 1:length(aff.coeffs)
        var = aff.vars[ind]
        var.m === m || error("Variable does not belong to this model")
        addelt!(counts, aff.vars[ind].idx, length(aff.coeffs[ind].vars))
    end
    for k in 1:counts.nnz
        ind = counts.nzidx[k]
        nterms = counts.elts[k]
        sizehint!(v.elts[ind].vars, nterms)
        sizehint!(v.elts[ind].coeffs, nterms)
    end

    for ind in 1:length(aff.coeffs)
        addelt!(v, aff.vars[ind].idx, aff.coeffs[ind])
    end
    vars = Array{IndepNormal}(v.nnz)
    coeffs = Array{CoefType}(v.nnz)
    for i in 1:v.nnz
        idx = v.nzidx[i]
        vars[i] = IndepNormal(m,idx)
        coeffs[i] = v.elts[idx]
    end
    empty!(v)
    empty!(counts)

    return JuMP.GenericAffExpr(vars, coeffs, aff.constant) # do we need to eliminate duplicates in the constant part also?

end

function merge_duplicates{CoefType <: Number}(aff::JuMP.GenericAffExpr{CoefType,IndepNormal}, v::IndexedVector{CoefType}, m::Model)
    ccdata = getCCData(m)

    resize!(v, ccdata.numRVs)
    for ind in 1:length(aff.coeffs)
        var = aff.vars[ind]
        is(var.m, m) || error("Variable does not belong to this model")
        addelt!(v, aff.vars[ind].idx, aff.coeffs[ind])
    end
    vars = Array{IndepNormal}(v.nnz)
    coeffs = Array{CoefType}(v.nnz)
    for i in 1:v.nnz
        idx = v.nzidx[i]
        vars[i] = IndepNormal(m,idx)
        coeffs[i] = v.elts[idx]
    end
    empty!(v)

    return JuMP.GenericAffExpr(vars, coeffs, aff.constant)

end
