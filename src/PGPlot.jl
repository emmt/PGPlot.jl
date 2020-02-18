module PGPlot

export
    Figure,
    colorbar,
    figure,
    heatmap!,
    heatmap,
    hist!,
    hist,
    palette,
    plot!,
    plot,
    scatter!,
    scatter

include("bindings.jl")

include("colormaps.jl")
using .Colormaps

include("interface.jl")
using .Plotting

end # module
