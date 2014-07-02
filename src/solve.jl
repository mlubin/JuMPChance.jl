using Distributions # this takes a while to load

function solvecc(m::Model;method=:Refomulate,debug::Bool = true)
    @assert method == :Reformulate || method == :Cuts

    ccdata = getCCData(m)

    if method == :Reformulate
        # check that we have pure chance constraints
        @assert all([isa(x,Real) for x in ccdata.RVmeans])
        @assert all([isa(x,Real) for x in ccdata.RVvars])

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
        @assert all([isa(x,Real) for x in ccdata.RVmeans])
        @assert all([isa(x,Real) for x in ccdata.RVvars])

        # set up slack variables and linear constraints
        nconstr = length(ccdata.chanceconstr)
        @defVar(m, slackvar[1:nconstr] >= 0)
        for i in 1:nconstr
            cc = ccdata.chanceconstr[i]
            ccexpr = cc.ccexpr
            nterms = length(ccexpr.vars)
            nu = quantile(Normal(0,1),1-cc.with_probability)
            if cc.sense == :(>=)
                @addConstraint(m, sum{getMean(ccexpr.vars[i])*ccexpr.coeffs[i], i=1:nterms} + nu*slackvar[i] + ccexpr.constant <= 0)
            else
                @addConstraint(m, sum{getMean(ccexpr.vars[i])*ccexpr.coeffs[i], i=1:nterms} - nu*slackvar[i] + ccexpr.constant >= 0)
            end
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
                    if violation > 1e-6
                        nviol += 1
                        # add a linearization
                        @addConstraint(m, sum{ getVar(ccexpr.vars[k])*getValue(ccexpr.coeffs[k])*ccexpr.coeffs[k], k in 1:nterms} <= sqrt(var)*slackvar[i])
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



end
