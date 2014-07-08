using Distributions # this takes a while to load

function solvecc(m::Model;method=:Refomulate,debug::Bool = true)
    @assert method == :Reformulate || method == :Cuts

    ccdata = getCCData(m)
    no_uncertains = all([isa(x,Real) for x in ccdata.RVmeans]) && all([isa(x,Real) for x in ccdata.RVvars])

    if method == :Reformulate
        # check that we have pure chance constraints
        no_uncertains || error("Cannot solve using reformulation, uncertain data are present")
        #display(ccdata.chanceconstr)

        for cc::ChanceConstr in ccdata.chanceconstr
            ccexpr = cc.ccexpr
            nu = quantile(Normal(0,1),1-cc.with_probability)
            nterms = length(ccexpr.vars)
            # add auxiliary variables for variance of each term
            # TODO: merge terms for duplicate r.v.'s
            @defVar(m, varterm[1:nterms])
            @addConstraint(m, defvar[i=1:nterms], varterm[i] == getStdev(ccexpr.vars[i])*ccexpr.coeffs[i])
            @defVar(m, slackvar >= 0)
            # conic constraint
            addConstraint(m, sum([varterm[i]^2 for i in 1:nterms]) <= slackvar^2)
            if cc.sense == :(>=)
                @addConstraint(m, sum{getMean(ccexpr.vars[i])*ccexpr.coeffs[i], i=1:nterms} + nu*slackvar + ccexpr.constant <= 0)
            else
                @assert cc.sense == :(<=)
                @addConstraint(m, sum{getMean(ccexpr.vars[i])*ccexpr.coeffs[i], i=1:nterms} - nu*slackvar + ccexpr.constant >= 0)
            end
        end
        println(m)

        return solve(m)

    else
        # check that we have pure chance constraints
        if no_uncertains
            solvecc_cuts(m, debug=debug)
        else
            solverobustcc_cuts(m, debug=debug)
        end
    end



end


function solvecc_cuts(m::Model; debug=true)


    ccdata = getCCData(m)

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

    # TODO: special handling for quadratic objectives
    
    debug && println("Solving deterministic model")

    status = solve(m)

    @assert status == :Optimal

    const MAXITER = 40 # make a parameter
    niter = 0
    while niter < MAXITER

        nviol = 0
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
            var += 1e-13 # avoid numerical issues with var == 0.0
            if cc.sense == :(<=)
                satisfied_prob = cdf(Normal(mean,sqrt(var)),0.0)
            else
                satisfied_prob = 1-cdf(Normal(mean,sqrt(var)),0.0)
            end
            debug && println("$satisfied_prob $mean $var")
            if satisfied_prob <= cc.with_probability
                # constraint is okay!
                continue
            else
                # check violation of quadratic constraint
                violation = var - getValue(slackvar[i])^2
                debug && println("VIOL $violation")
                @assert violation > -1e-13
                if satisfied_prob >= cc.with_probability + 0.001 || violation > 1e-6
                    nviol += 1
                    # add a linearization
                    @addConstraint(m, sum{ getValue(varterm[i][k])*varterm[i][k], k in 1:nterms} <= sqrt(var)*slackvar[i])
                end
            end
        end

        if nviol == 0
            debug && println("Done after $niter iterations")
            return :Optimal
        else
            debug && println("Iteration $niter: $nviol violations")
        end
        status = solve(m)
        @assert status == :Optimal
        niter += 1
    end

    return :UserLimit # hit iteration limit


end

function solverobustcc_cuts(m::Model; debug=true)


    ccdata = getCCData(m)

    # set up slack variables
    nconstr = length(ccdata.chanceconstr)
    #@defVar(m, slackvar[1:nconstr] >= 0)
    #=
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
    end=#

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
            @assert ub >= lb
            means_nominal[i] = (lb+ub)/2
            means_deviation[i] = (ub-lb)/2
        end
        if isa(rvvar,Real)
            vars_nominal[i] = rvvar
        else
            lb,ub = rvvar
            @assert ub >= lb
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

    @assert status == :Optimal

    const MAXITER = 40 # make a parameter
    niter = 0
    while niter < MAXITER

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
            nominal_mean = getValue(ccexpr.constant)
            nominal_var = 0.0
            for k in 1:nterms
                exprval = getValue(ccexpr.coeffs[k])
                debug && println(ccexpr.coeffs[k], " => ", exprval)
                idx = ccexpr.vars[k].idx
                meanvals[k] = means_deviation[idx]*abs(exprval)
                varvals[k] = vars_deviation[idx]*exprval^2
                nominal_mean += means_nominal[idx]*exprval
                nominal_var += vars_nominal[idx]*exprval^2
            end
            debug && println("nominal_mean = $nominal_mean")
            debug && println("nominal_var = $nominal_var")
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
            worst_var += 1e-13 # avoid numerical issues with var == 0.0
            debug && println("worst_mean = $worst_mean")
            debug && println("worst_var = $worst_var")
            if cc.sense == :(<=)
                satisfied_prob = cdf(Normal(worst_mean,sqrt(worst_var)),0.0)
            else
                satisfied_prob = 1-cdf(Normal(worst_mean,sqrt(worst_var)),0.0)
            end
            debug && println("$satisfied_prob $mean $var")
            if satisfied_prob <= cc.with_probability
                # constraint is okay!
                continue
            else
                debug && println("VIOL ", 100*(satisfied_prob - cc.with_probability), "%")
                if satisfied_prob >= cc.with_probability + 0.001
                    nviol += 1
                    var_coeffs = [vars_nominal[ccexpr.vars[k].idx] for k in 1:nterms]
                    for k in 1:cc.uncertainty_budget_variance
                        var_coeffs[sorted_var_idx[k]] += vars_deviation[ccexpr.vars[sorted_var_idx[k]].idx]
                    end


                    # add a linearization
                    # f(x') + f'(x')(x-x') <= 0
                    if cc.sense == :(>=)
                        @addConstraint(m, ccexpr.constant + nu*sqrt(worst_var) + 
                            sum{means_nominal[ccexpr.vars[k].idx]*(ccexpr.coeffs[k]), k=1:nterms} + 
                            sum{sign(getValue(ccexpr.coeffs[sorted_mean_idx[r]]))*means_deviation[sorted_mean_idx[r]]*(ccexpr.coeffs[sorted_mean_idx[r]]), r=1:cc.uncertainty_budget_mean} + 
                            (nu/sqrt(worst_var))*sum{var_coeffs[k]*getValue(ccexpr.coeffs[k])*(ccexpr.coeffs[k]-getValue(ccexpr.coeffs[k])),k=1:nterms}  <= 0)
                    else
                        @addConstraint(m, ccexpr.constant - nu*sqrt(worst_var) + 
                            sum{means_nominal[ccexpr.vars[k].idx]*(ccexpr.coeffs[k]), k=1:nterms} - 
                            sum{sign(getValue(ccexpr.coeffs[sorted_mean_idx[r]]))*means_deviation[sorted_mean_idx[r]]*(ccexpr.coeffs[sorted_mean_idx[r]]), r=1:cc.uncertainty_budget_mean} - 
                            (nu/sqrt(worst_var))*sum{var_coeffs[k]*getValue(ccexpr.coeffs[k])*(ccexpr.coeffs[k]-getValue(ccexpr.coeffs[k])), k=1:nterms}  >= 0)
                        debug && println("ADDED: ", m.linconstr[end])

                    end

                end
            end
        end

        if nviol == 0
            debug && println("Done after $niter iterations")
            return :Optimal
        else
            debug && println("Iteration $niter: $nviol violations")
        end
        status = solve(m)
        @assert status == :Optimal
        niter += 1
    end

    return :UserLimit # hit iteration limit


end
