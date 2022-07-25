triangularnumber(n::Int) = (n * (n + 1)) ÷ 2
idiag(i::Int, j::Int) = i + j
fdiag(i::Int, j::Int) = triangularnumber(idiag(i, j)) + i + 1

function plotpadetable(t::Matrix{NTuple{2,T}}, padecolors; ifticks=false) where {T}
    c = [padecolors[fdiag(i...)+1] for i in t]
    plt = plot(c')
    xl, yl = xlims(), ylims()
    n = xl[2] - xl[1]
    m = yl[2] - yl[1]
    m = Int(m) - 1
    n = Int(n) - 1
    for i = 1:m+2
        plot!((i - 0.5) * ones(2), collect(yl), c=:black, lw=1, label=:none)
    end
    for j = 1:n+2
        plot!(collect(xl), (j - 0.5) * ones(2), c=:black, lw=1, label=:none)
    end
    if ifticks
        plot!(
            xticks=((1:m+1), 0:m),
            yticks=((1:n+1), 0:n),
            xguide="m",
            yguide="n",
            xmirror=true,
            xtickfontsize=5,
            ytickfontsize=5,
            xlabelfontsize=8,
            ylabelfontsize=8
        )
    else
        plot!(axis=false, xticks=false, yticks=false)
    end
    plot!(
        xlims=xl,
        ylims=yl,
    )
    return plt
end
