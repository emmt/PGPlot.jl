#
# interface.jl --
#
# Implements high-level interface for plotting with the PGPlot library.
#

module Plotting

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

using Colors, ArrayTools, TwoDimensional
using ..PGPlot.Bindings
using ..PGPlot.Colormaps
using .Bindings: pgarray, DEFAULT_XFORM

# FIXME: type piracy
Base.Tuple(A::AffineTransform) = (A.x, A.xx, A.xy, A.y, A.yx, A.yy)


# Structure used to wrap an integer so that some basic function can be extended
# while avoiding type-piracy.
struct Figure
    num::Int   # Figure number.
end
number(fig::Figure) = fig.num

const ExtentOption = Union{Nothing,NTuple{4,<:Union{Nothing,Real}}}
const AxisStyleOption = Union{Nothing,Symbol,Integer}
const ColorOption = Union{Nothing,AbstractString,Symbol,Integer,Colorant,RGBVec}
const FontOption = Union{Nothing,AbstractString,Symbol,Integer}
const LabelOption = Union{Nothing,AbstractString,Symbol}
const XFormOption = Union{Nothing, AffineTransform,NTuple{6,Real},
                          AbstractVector{<:Real},AbstractMatrix{<:Real}}
const HeightOption = Union{Nothing,Real}
const FigureOption = Union{Nothing,Figure,Integer}

const Marker = Union{Integer,Char} # FIXME: Symbol?
const Markers = Union{Marker,AbstractVector{<:Marker}}

# Color index to use for color set by value/name.
const USER_COLOR_INDEX = PGInt(9)

# Constants for available fonts.
const NORMAL_FONT = PGInt(1)
const ROMAN_FONT  = PGInt(2)
const ITALIC_FONT = PGInt(3)
const SCRIPT_FONT = PGInt(4)

# Transparent text background color.
const TRANSPARENT_BACKGROUND = PGInt(-1)
const BACKGROUND_COLOR_INDEX = PGInt(0)
const FOREGROUND_COLOR_INDEX = PGInt(1)
const DEFAULT_FRAME_COLOR = FOREGROUND_COLOR_INDEX
const DEFAULT_TITLE_COLOR = FOREGROUND_COLOR_INDEX
const DEFAULT_DRAWING_COLOR = FOREGROUND_COLOR_INDEX
const DEFAULT_TEXT_BACKGROUND = TRANSPARENT_BACKGROUND
const DEFAULT_TEXT_COLOR = DEFAULT_DRAWING_COLOR # for annotations
const DEFAULT_TEXT_FONT = NORMAL_FONT
const DEFAULT_TEXT_HEIGHT = PGFloat(1.0)

function check(A::AbstractArray)
    flag = true
    @inbounds @simd for i in eachindex(A)
        flag &= isfinite(A[i])
    end
    flag || throw(ArgumentError("non-finite value(s) in argument"))
end

function check(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) ||
        throw(DimensionMismatch("x and y must have the same length"))
    check(x)
    check(y)
end

"""

```julia
confirm(flag)
```

specifies whether or not interactive devices should ask for user
confirmation before advancing to next plot.

"""
confirm(flag::Bool) = pgask(flag)

# PGPlot has a rather small maximum number (GRIMAX = 8) of open devices.
# We keep the list of open devices/figures in a vector.
const DEFAULT_DEVICE = "/XSERVE"
const DEVICES = Vector{PGInt}(undef, 0)
const BAD_DEVICE = eltype(DEVICES)(-1)
const FIGCNT = Ref{Int}(0)
const MAXFIGCNT = 100 # just used to detect wrong figure numbers

"""

```julia
device(fig)
```

yields the PGPlot device number associated with figure `fig`.  The result is a
positive integer if there is an open PGPlot device for that figure, -1
otherwise.

"""
device(fig::Figure) = device(number(fig))
device(fig::Integer) = (1 ≤ fig ≤ length(DEVICES) ? DEVICES[fig] : BAD_DEVICE)

"""

```julia
figure() -> fig
figure(nothing) -> fig
```

create a new figure using the default plotting device and selects it as the
current plotting device.  The figure is returned.

```julia
figure(dev::AbstractString) -> fig
```

opens the plotting device `dev`, selects it as the current plotting device and
associates a new figure with it.  The figure is returned.

```julia
figure(num::Integer) -> fig
```

if figure number `num` exists, makes it associated plotting device the current
one; otherwise opens the default plotting device, selects it as the current
plotting device and associates it with the figure number `num`.  The figure is
returned.

```julia
figure(fig::Figure) -> fig
```

checks that figure `fig` is associated with an open plotting device and makes
it the current one.  The figure is returned.  An exception is thrown if `fig`
is not associated with an open plotting device.

""" figure

# Create a new figure using the default device.
figure() = figure(DEFAULT_DEVICE)
figure(::Nothing) = figure()

# Create a new figure for the given device string.
figure(dev::AbstractString) = register_device(open_device(dev))

function figure(num::Integer)
    1 ≤ num ≤ MAXFIGCNT || throw_invalid_figure_number(num)
    if num ≤ length(DEVICES)
        dev = DEVICES[num]
        if dev ≥ 1
            pgslct(dev)
            return Figure(num)
        end
    end
    while length(DEVICES) < num
        push!(DEVICES, BAD_DEVICE)
    end
    DEVICES[num] = open_device(DEFAULT_DEVICE)
    return Figure(num)
end

# Use an existing figure.
function figure(fig::Figure)
    dev = device(fig)
    if dev > 0
        pgslct(dev)
        return fig
    end
    1 ≤ number(fig) ≤ length(DEVICES) && throw_device_closed(fig)
    throw_nonexisting_figure(fig)
end

Base.isopen(fig::Figure) = device(fig) > 0

function Base.close(fig::Figure)
    dev = device(fig)
    if dev > 0
        DEVICES[dev] = BAD_DEVICE
        pgclos(dev)
        return nothing
    end
    1 ≤ number(fig) ≤ length(DEVICES) && throw_device_closed(fig)
    throw_nonexisting_figure(fig)
end

@noinline throw_invalid_figure_number(num::Integer) =
    throw(ArgumentError(string("invalid figure number: ", num)))

@noinline throw_invalid_device(dev::Integer) =
    throw(ArgumentError(string("invalid PGPlot device number: ", dev)))

@noinline throw_open_device_failure(dev::AbstractString) =
    error(string("cannot open PGPlot device: ", dev))

throw_nonexisting_figure(fig::Figure) = throw_nonexisting_figure(number(fig))
@noinline throw_nonexisting_figure(num::Integer) =
    throw(ArgumentError(string("Figure #", num, " does not exist")))

throw_device_closed(fig::Figure) = throw_device_closed(number(fig))
@noinline throw_device_closed(num::Integer) =
    throw(ArgumentError(string("Plotting device of Figure #", num,
                               " has been closed")))

function current_figure()
    dev = pgqid()
    dev > 0 || error("there is no open plotting devices")
    return register_device(dev)
end

select_figure(fig::Figure) = figure(fig)
select_figure(num::Integer) = figure(num)
select_figure(::Nothing) = select_figure()
select_figure() = (if pgqid() < 1; figure(); end; nothing)

function open_device(name::AbstractString)
    dev = pgopen(name)
    dev > 0 || throw_open_device_failure(name)
    pgask(false) # do not wait for user input
    if pgqinf("HARDCOPY") == "YES"
        pgscr(0, 1,1,1) # set background color to white
        pgscr(1, 0,0,0) # set foreground color to black
    else
        pgscr(0, 0,0,0) # set background color to black
        pgscr(1, 1,1,1) # set foreground color to white
    end
    return dev
end

function register_device(dev::Integer)
    dev > 0 || throw_invalid_device(dev)
    j = 0
    @inbounds for i in 1:length(DEVICES)
        if DEVICES[i] == dev
            return Figure(i)
        elseif j == 0 && DEVICES[i] < 1
            j = i
        end
    end
    if j != 0
        DEVICES[j] = dev
        return Figure(j)
    end
    if length(DEVICES) < MAXFIGCNT
        push!(DEVICES, dev)
        return Figure(length(DEVICES))
    end
    error("too many registered figures")
end

function forget_device(dev::Integer)
    if dev > 0
        @inbounds for i in 1:length(DEVICES)
            if DEVICES[i] == dev
                DEVICES[i] = BAD_DEVICE
                break
            end
        end
    end
    nothing
end

"""

```julia
frame(...)
```

draws the environment of subsequent plots.

"""
function frame(xmin::Real, xmax::Real, ymin::Real, ymax::Real;
               fig::FigureOption = nothing,
               wait::Bool = false,
               just::Bool = false,
               style::AxisStyleOption = nothing,
               xlabel::LabelOption = nothing,
               ylabel::LabelOption = nothing,
               title::LabelOption = nothing,
               textcolor::ColorOption = nothing,
               framecolor::ColorOption = nothing)
    select_figure(fig)
    pgsci(get_color(framecolor, DEFAULT_FRAME_COLOR))
    pgenv(xmin, xmax, ymin, ymax, just, get_axis_style(style))
    if xlabel !== nothing || ylabel !== nothing || title !== nothing
        # deal with colors, size, text style
        pgsci(get_color(textcolor, DEFAULT_TITLE_COLOR))
        pglab(get_label(xlabel), get_label(ylabel), get_label(title))
    end
end

function frame(x::AbstractVector, y::AbstractVector;
               extent::ExtentOption = nothing,
               kwds...)
    check(x, y)
    xmin, xmax = get_vrange(PGFloat, x, extent, (1,2))
    ymin, ymax = get_vrange(PGFloat, y, extent, (3,4))
    frame(xmin, xmax, ymin, ymax; kwds...)
end

"""

```julia
curve([x,] y; kwds...)
```

draws a frame box (with [`frame`](@ref)) and a curve (as a polyline) joining
the points defined by abscissae `x` and ordinates `y` (both being vectors in
that case).  If `x` is missing, the vector `y` is plotted against its indices.

To plot functions, do one of:

```julia
curve(x, fy; kwds...)      # plot fy(x) against x
curve(fx, y; kwds...)      # plot y against fx(y)
curve(t, fx, fy; kwds...)  # plot fy(t) against fx(t) (parametric plot)
```

To add a curve to an existing plot, call [`curve!`](@ref) instead.

Possible keywords:

- `color` specifies the color of the lines.

Other keywods are transmitted to [`frame`](@ref) used to draw the frame box
(e.g., `fig`, `xlabel`, `ylabel`, `title`, `extent`, ...).

"""
function curve(x::AbstractVector, y::AbstractVector;
              color::ColorOption = nothing,
              kwds...)
    frame(x, y; kwds...)
    _curve(x, y, color)
end

"""

```julia
curve!(args...; kwds...)
```

draws a curve in the current frame preserving existing contents.  Arguments are
the same as for [`curve`](@ref).

Possible keywords: `color`, `fig`.

- `color` specifies the color of the lines.
- `fig` specifies the figure where to draw.

"""
function curve!(x::AbstractVector, y::AbstractVector;
               fig::FigureOption = nothing,
               color::ColorOption = nothing)
    select_figure(fig)
    check(x, y)
    _curve(x, y, color)
end

# Helper for `curve` and `curve!` methods.
function _curve(x::AbstractVector, y::AbstractVector,
                color::ColorOption)
    pgsci(get_color(color))
    pgline(x, y)
end

"""

```julia
hist([x,] y; kwds...)
```

draws a histogram of `y` against `x`, that is a stair-style curve.

"""
function hist(x::AbstractVector, y::AbstractVector;
              center::Bool = true,
              color::ColorOption = nothing,
              kwds...)
    frame(x, y; kwds...)
    _hist(x, y, center, color)
end

function hist!(x::AbstractVector, y::AbstractVector;
               fig::FigureOption = nothing,
               center::Bool = true,
               color::ColorOption = nothing)
    select_figure(fig)
    check(x, y)
    _hist(x, y, center, color)
end

# Helper for `hist` and `hist!` methods.
function _hist(x, y, center::Bool, color::ColorOption)
    pgsci(get_color(color))
    pgbin(x, y, center)
end

"""

```julia
scatter([x,] y, sym; kwds...)
```

draws points with symbol `sym` at coordinates `(x,y)`.

"""
function scatter(x::AbstractVector, y::AbstractVector,
                 sym::Markers;
                 color::ColorOption = nothing,
                 kwds...)
    frame(x, y; kwds...)
    _scatter(x, y, sym, color)
end

function scatter!(x::AbstractVector, y::AbstractVector,
                  sym::Markers;
                  fig::FigureOption = nothing,
                  color::ColorOption = nothing)
    select_figure(fig)
    check(x, y)
    _scatter(x, y, sym, color)
end

# Helper for `scatter` and `scatter!` methods.
function _scatter(x::AbstractVector, y::AbstractVector,
                  sym::Union{Integer,AbstractVector{<:Integer}},
                  color::ColorOption)
    pgsci(get_color(color))
    pgpnts(x, y, sym)
end

# Extend curve-like methods.
for func in (:curve, :curve!, :hist, :hist!)
    @eval begin

        $func(y::AbstractVector; kwds...) =
            $func(axes(y,1), y; kwds...)

        $func(x::AbstractVector, fy; kwds...) =
            $func(x, PGFloat[fy(x[i]) for i in eachindex(x)]; kwds...)

        $func(fx, y::AbstractVector; kwds...) =
            $func(PGFloat[fx(y[i]) for i in eachindex(y)], y; kwds...)

        $func(t::AbstractVector, fx, fy; kwds...) =
            $func(PGFloat[fx(t[i]) for i in eachindex(t)],
                  PGFloat[fy(t[i]) for i in eachindex(t)]; kwds...)
    end
end

for func in (:scatter, :scatter!)
    @eval begin

        $func(y::AbstractVector, sym::Markers; kwds...) =
            $func(axes(y,1), y, sym; kwds...)

        $func(x::AbstractVector, fy, sym::Markers; kwds...) =
            $func(x, PGFloat[fy(x[i]) for i in eachindex(x)], sym; kwds...)

        $func(fx, y::AbstractVector, sym::Markers; kwds...) =
            $func(PGFloat[fx(y[i]) for i in eachindex(y)], y, sym; kwds...)

        $func(t::AbstractVector, fx, fy, sym::Markers; kwds...) =
            $func(PGFloat[fx(t[i]) for i in eachindex(t)],
                  PGFloat[fy(t[i]) for i in eachindex(t)], sym; kwds...)
    end
end

"""

```julia
heatmap(A, tr=[0,1,0, 0,0,1]; kwds...)
```

draws a *heat-map* , that is a map whose cells have pseudo-colors computed from
the values in the 2-dimensional array `A`.  By default, the first and second
dimensions of `A` are assumed to correspond to the abscissa and oridinate axes
respectively.

Keywords (other keywords are passed to [`frame`](@ref)):

- `vmin` and `vmax` specify the range of values to consider.  It is possible to
  have `vmin > vmax` to reverse the order of the colors.  If unspecified, the
  extreme values of `A` are considered.

- `cbar` specifies where to put a color-bar: `cbar = 'L'` `'R'`, `'T'` or `'B'`
  to draw a color-bar on the Left, Right, Top or Bottom side of the frame;
  `cbar = nothing` to draw no color-bar.  By default, `cbar = 'R'`.  Keyword
  `zlabel` can be used to specify a label for the color-bar.  For more options,
  use `cbar = nothing` and call [`colorbar`](@ref) directly after drawing the
  heat-map (the range of values are automatically saved by `heatmap` and need
  not be specified to `colorbar` in this case).

- `just` specifies whether both axis should have the same scale or not.  By
  default, `just = true`.

To add another heat-map to an existing plot, call [`heatmap!`](@ref) instead.

"""
function heatmap(A::AbstractMatrix{<:Real},
                 xform::XFormOption = DEFAULT_XFORM;
                 kwds...)
    heatmap(pgarray(PGFloat, A), get_xform(xform); kwds...)
end

# @btime Plotting.heatmap($z) for a 81×81 image
#   8.112 ms (2 allocations: 25.77 KiB)
function heatmap(A::DenseMatrix{PGFloat}, tr::DenseArray{PGFloat};
                 vmin::Union{Nothing,Real} = nothing,
                 vmax::Union{Nothing,Real} = nothing,
                 cbar::Union{Nothing,Char} = 'R',
                 zlabel::AbstractString = "",
                 cmap::Union{Nothing,AbstractString} = nothing,
                 just::Bool = true,
                 kwds...)
    # Get range of values to plot.
    check(A)
    bg, fg = get_zrange(A, vmin, vmax)

    # Get maximal extent of coordinates.
    I, J = axes(A)
    xmin, xmax, ymin, ymax = get_extent(PGFloat, I, J, tr)

    # Draw viewport.
    frame(xmin, xmax, ymin, ymax; just = just, kwds...)

    # Draw map.
    if cmap !== nothing
        palette(cmap)
    end
    pgimag(A, I, J, bg, fg, tr)

    # Optional color bar.
    if cbar !== nothing
        colorbar(bg, fg; side = cbar, label = zlabel)
    end
end

"""

```julia
heatmap!(A, tr=[0,1,0, 0,0,1]; kwds...)
```

draws a *heat-map* in the current frame preserving existing contents.
Arguments are the same as for [`heatmap`](@ref).  The only supported keywords
are `vmin`, `vmax` and `fig`.

"""
function heatmap!(A::AbstractMatrix{<:Real},
                  xform::XFormOption = DEFAULT_XFORM;
                  kwds...)
    heatmap!(pgarray(PGFloat, A), pgarray(PGFloat, xform); kwds...)
end

function heatmap!(A::DenseMatrix{PGFloat}, tr::DenseArray{PGFloat};
                  fig::FigureOption = nothing,
                  vmin::Union{Nothing,Real} = nothing,
                  vmax::Union{Nothing,Real} = nothing)
    select_figure(fig)
    check(A)
    bg, fg = get_zrange(A, vmin, vmax)
    pgimag(A, bg, fg, tr)
end

# The following are to remember image settings for the color-bar.
const last_bg = Ref{PGFloat}(0)
const last_fg = Ref{PGFloat}(1)

# Like get_vrange but fix equal endpoints and remember values
# for the color-bar.
function get_zrange(A::AbstractMatrix{PGFloat}, vmin, vmax)
    bg, fg = get_vrange(A, vmin, vmax)
    if bg == fg
        # Quick fix.
        bg -= tiny(bg)
        fg += tiny(fg)
    end
    last_bg[] = bg
    last_fg[] = fg
    return (bg, fg)
end

tiny(val::AbstractFloat) =
    (del = oftype(val, 0.001);
     val == zero(val) ? del : abs(val)*del)

"""

```julia
colorbar(vmin, vmax)
```

draws a color-bar for values in the range `vmin` to `vmax`.  It is possible to
have `vmin > vmax` to reverse the order of the colors.  If unspecified, the
last values set by [`heatmap`](@ref) or [`heatmap!`](@ref) are considered.

Keywords:

- `side` specifies where to draw the color-bar: `side = 'L'` `'R'`, `'T'` or
  `'B'` to draw a color-bar on the Left, Right, Top or Bottom side of the
  frame. The default is to draw the color-bar on the right.

- `label` specifies a label for the color-bar.

- `spacing` specifies a distance breween the color-bar and the frame border in
  units of the current character size.  The value may be negative to draw the
  color-bar inside the frame.

- `width` specifies the thickness of the color-bar in units of the current
  character size.

"""
function colorbar(vmin::Real = last_bg[],
                  vmax::Real = last_fg[];
                  side::Char = 'R',
                  label::AbstractString = "",
                  spacing::Real = PGFloat(2),
                  width::Real = PGFloat(5))
    pgwedg(side*'I', spacing, width, vmin, vmax, label)
end

"""

```julia
get_xform(arg) -> tr
```

yields coordinate transform `tr` directly usable by PGPlot for a variety of
arguments: `nothing` to get the default *identity* transform, a 2×3 matrix, a
6-element vector, a 6-tuple or an `AffineTransform` from the `TwoDimensional`
package.

"""
get_xform(::Nothing) = DEFAULT_XFORM
get_xform(A::AffineTransform) = PGFloat[A.x, A.xx, A.xy, A.y, A.yx, A.yy]
get_xform(A::NTuple{6,Real}) = PGFloat[A...]
get_xform(A::AbstractVector) = (@assert size(A) == (6,); pgarray(PGFloat, A))
get_xform(A::AbstractMatrix) = (@assert size(A) == (2,3); pgarray(PGFloat, A))

"""

```julia
get_extent(I, J, tr) -> xmin, xmax, ymin, ymax
```

yields the extent of world coordinates after applying coordinates transform
`tr` to the range of indices `I = imin:imax` and `J = jmin:jmax`.  Ranges of
indices can be specified by their endpoints:

```julia
get_extent(i1, i2, j1, j2, tr) -> xmin, xmax, ymin, ymax
```

In this latter case, the order of `i1` and `i2 (and of `j1` and `j2`) can be
reversed.

"""
function get_extent(i1::Integer, i2::Integer,
                    j1::Integer, j2::Integer,
                    tr::AbstractArray{<:Real})
    get_extent(PGFloat, i1, i2, j1, j2, tr)
end

function get_extent(I::AbstractUnitRange{<:Integer},
                    J::AbstractUnitRange{<:Integer},
                    tr::AbstractArray{<:Real})
    get_extent(first(I), last(I), first(J), last(J), tr)
end

function get_extent(::Type{T},
                    I::AbstractUnitRange{<:Integer},
                    J::AbstractUnitRange{<:Integer},
                    tr::AbstractArray{<:Real}) where {T<:AbstractFloat}
    get_extent(T, first(I), last(I), first(J), last(J), tr)
end

function get_extent(::Type{T},
                    i1::Integer, i2::Integer,
                    j1::Integer, j2::Integer,
                    tr::AbstractArray{<:Real}) where {T<:AbstractFloat}
    @assert length(tr) == 6
    tr1, tr2, tr3 = T(tr[1]), T(tr[2]), T(tr[3])
    tr4, tr5, tr6 = T(tr[4]), T(tr[5]), T(tr[6])
    s = T(0.5)
    imin = min(T(i1), T(i2)) - s
    imax = max(T(i1), T(i2)) + s
    jmin = min(T(j1), T(j2)) - s
    jmax = max(T(j1), T(j2)) + s
    x1 = tr1 + tr2*imin + tr3*jmin
    y1 = tr4 + tr5*imin + tr6*jmin
    x2 = tr1 + tr2*imax + tr3*jmin
    y2 = tr4 + tr5*imax + tr6*jmin
    x3 = tr1 + tr2*imin + tr3*jmax
    y3 = tr4 + tr5*imin + tr6*jmax
    x4 = tr1 + tr2*imax + tr3*jmax
    y4 = tr4 + tr5*imax + tr6*jmax
    xmin = min(x1, x2, x3, x4)
    xmax = max(x1, x2, x3, x4)
    ymin = min(y1, y2, y3, y4)
    ymax = max(y1, y2, y3, y4)
    return (xmin, xmax, ymin, ymax)
end

"""

```julia
get_vrange([T=eltype(A),] A, lo, hi) -> vmin, vmax
```

yields the range of values in array `A` given optional low and high clipping
levels `lo` and `hi` which both can be either a value or `nothing` to indicate
that the corrsponding extreme (minimum or maximum) value in `A` should be
taken.  The result is a 2-tuple of values of type `T` which is, by default,
that of the elements of `A`.

"""
get_vrange(A::AbstractArray{T}, lo, hi) where {T} =
    get_vrange(T, A, lo, hi)
get_vrange(::Type{T}, A::AbstractArray, lo::Nothing, hi::Nothing) where {T} =
    extrema(A)
get_vrange(::Type{T}, A::AbstractArray, lo::Real, hi::Nothing) where {T} =
    (T(lo), maximum(A))
get_vrange(::Type{T}, A::AbstractArray, lo::Nothing, hi::Real) where {T} =
    (minimum(A), T(hi))
get_vrange(::Type{T}, A::AbstractArray, lo::Real, hi::Real) where {T} =
    (T(lo), T(hi))

function get_vrange(A::AbstractArray{T},
                    opt::Union{Nothing,Tuple},
                    I::NTuple{2,Integer}) where {T}
    get_vrange(T, A, opt, I)
end

function get_vrange(::Type{T}, A::AbstractArray, opt::Nothing,
                    I::NTuple{2,Integer}) where {T}
    get_vrange(T, A, nothing, nothing)
end

function get_vrange(::Type{T}, A::AbstractArray, opt::Tuple,
                    I::NTuple{2,Integer}) where {T}
    get_vrange(T, A, opt[I[1]], opt[I[2]])
end

nicerange(x::Tuple{Real,Real}) = nicerange(x[1], x[2])
nicerange(x1::Real, x2::Real) = nicerange(PGFloat(x1), PGFloat(x2))
function nicerange(x1::PGFloat, x2::PGFloat)
    @assert isfinite(x1) && isfinite(x2)
    if x1 != x2
        return pgrnge(x1, x2)
    elseif x1 == 0
        return (PGFloat(-0.001), PGFloat(+0.001))
    else
        d = PGFloat(10^(round(log10(abs(x1))) - 1))
        return (-d, +d)
    end
end

get_label(::Nothing) = ""
get_label(val::AbstractString) = val
get_label(val::Symbol) = string(val)

get_axis_style(::Nothing) = PGInt(1) # default axis style
get_axis_style(val::Integer) =
    (-2 ≤ val ≤ 2 || val == 10 || val == 20 || val == 30 ? PGInt(val) :
     throw_bad_axis_style(val))
get_axis_style(val::Symbol) =
    # Draw no box, axes or labels:
    (val == :none ? PGInt(-2) :
     # Draw box only:
     val == :box ? PGInt(-1) :
     # Draw box and label it with coordinates:
     val == :label ? PGInt(0) :
     # Same as AXIS=0, but also draw the coordinate axes (X=0, Y=0):
     val == :axis ? PGInt(1) :
     # Same as AXIS=1, but also draw grid lines at major increments of the
     # coordinates:
     val == :grid ? PGInt(2) :
     # Draw box and label X-axis logarithmically:
     val == :loglin ? PGInt(10) :
     # Draw box and label Y-axis logarithmically:
     val == :linlog ? PGInt(20) :
     # Draw box and label both axes logarithmically:
     val == :loglog ? PGInt(30) :
     throw_bad_axis_style(val))

@noinline throw_bad_axis_style(val) =
    throw(ArgumentError(string("invalid axis style: ", val)))

get_color(::Nothing) = FOREGROUND_COLOR_INDEX
get_color(val, def) = get_color(val)
get_color(::Nothing, def) = get_color(def)
get_color(ci::Integer) = PGInt(ci)
get_color(val::UInt32) = get_color(RGB24(val))
get_color(val::Symbol) = get_color(String(val))
get_color(val::AbstractString) = get_color(parse(RGB{PGFloat}, val))
get_color(val::RGBVec) = get_color(convert(RGB{PGFloat}, val))
function get_color(val::Colorant)
    pgscr(USER_COLOR_INDEX, val)
    return USER_COLOR_INDEX
end

get_text_height(val::Real) =
    (@assert isfinite(val) && val > 0; PGFloat(val))

get_font(val::AbstractString) =
  (val == "normal" ? NORMAL_FONT :
   val == "roman"  ? ROMAN_FONT  :
   val == "italic" ? ITALIC_FONT :
   val == "script" ? SCRIPT_FONT : throw_bad_text_font(val))
get_font(val::Symbol) =
  (val == :normal ? NORMAL_FONT :
   val == :roman  ? ROMAN_FONT  :
   val == :italic ? ITALIC_FONT :
   val == :script ? SCRIPT_FONT : throw_bad_text_font(val))
get_font(val::Integer) =
  (val == 1 ? NORMAL_FONT  :
   val == 2 ? ROMAN_FONT   :
   val == 3 ? ITALIC_FONT  :
   val == 4 ? SCRIPT_FONT  : throw_bad_text_font(val))

@noinline throw_bad_text_font(val) =
    throw(ArgumentError(string("invalid text font: ", val)))

end # module
