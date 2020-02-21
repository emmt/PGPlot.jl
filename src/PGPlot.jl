module PGPlot

export
    Figure,
    colorbar,
    curve!,
    curve,
    figure,
    heatmap!,
    heatmap,
    hist!,
    hist,
    palette,
    scatter!,
    scatter

include("bindings.jl")

include("colormaps.jl")
using .Colormaps

include("interface.jl")
using .Plotting

end # module
