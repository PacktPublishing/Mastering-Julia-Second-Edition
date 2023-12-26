function traverse!(ex, symbols) end

function traverse!(ex::Symbol, symbols)
    push!(symbols, ex)
end

function traverse!(ex::Expr, symbols)
    if ex.head == :call  # function call
        for arg in ex.args[2:end]
            traverse!(arg, symbols)  # recursive
        end
    else
        for arg in ex.args
            traverse!(arg, symbols)  # recursive
        end
    end
end

function traverse(ex::Expr)
    symbols = Symbol[]
    traverse!(ex, symbols)
    return unique(symbols)  # Don't output duplicate
end
