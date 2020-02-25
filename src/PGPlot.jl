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
    mouse,
    palette,
    scatter!,
    scatter

include("bindings.jl")

include("colormaps.jl")
using .Colormaps

include("plotting.jl")
using .Plotting

end # module
