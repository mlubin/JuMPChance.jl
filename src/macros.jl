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
        varname = esc(x.args[1])
        idxvars = {}
        idxsets = {}
        refcall = Expr(:ref,varname)
        for s in x.args[2:end]
            if isa(s,Expr) && s.head == :(=)
                idxvar = esc(s.args[1])
                idxset = esc(s.args[2])
            else
                idxvar = gensym()
                idxset = esc(s)
            end
            push!(idxvars, idxvar)
            push!(idxsets, idxset)
            push!(refcall.args, idxvar)
        end
        varstr = :(string($(string(varname.args[1])),"["))
        for idxvar in idxvars
            push!(varstr.args,:(string($idxvar)))
            push!(varstr.args,",")
        end
        deleteat!(varstr.args,length(varstr.args))
        push!(varstr.args,"]")
        code = :( $(refcall) = IndepNormal($m, $mean, $var, $varstr) )
        for (idxvar, idxset) in zip(reverse(idxvars),reverse(idxsets))
            code = quote
                for $idxvar in $idxset
                    $code
                end
            end
        end
        
        mac = Expr(:macrocall,symbol("@gendict"),varname,:IndepNormal,idxsets...)
        code = quote 
            $mac
            $code
            nothing
        end
        return code
    end
end
