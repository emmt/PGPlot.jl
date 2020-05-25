#
# plotting.jl --
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
    mouse,
    palette,
    scatter!,
    scatter

using Colors, ArrayTools, TwoDimensional
using ..PGPlot.Bindings
using ..PGPlot.Colormaps
using .Bindings: pgarray, DEFAULT_XFORM, get_own_xform, throw_bad_xform_size

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
const XFormOption = Union{Nothing,AffineTransform,NTuple{6,Real},
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

"""

`Box{T}` represents a rectangular region whose side are aligned with the
coordinate axes.

"""
struct Box{T<:Real}
    x1::T
    x2::T
    y1::T
    y2::T
end

Box(A::NTuple{4,Real}) = Box(A[1], A[2], A[3], A[4])
Box(x1, x2, y1, y2) = Box(promote(x1, x2, y1, y2)...)
Box(B::Box) = B
Box{T}(B::Box{T}) where {T<:Real} = B
Box{T}(B::Box) where {T<:Real} = Box{T}(B.x1, B.x2, B.y1, B.y2)
Base.Tuple(B::Box) = (B.x1, B.x2, B.y1, B.y2)

const BAD_COLOR = colorant"orange"

abstract type AbstractScaling{T<:Integer} end

struct LinearScaling{T} <: AbstractScaling{T}
    lo::T   # first color index
    hi::T   # last color index
    bad::T  # color index for bad values
end

Base.first(scl::AbstractScaling) = scl.lo
Base.last(scl::AbstractScaling) = scl.hi
bad(scl::AbstractScaling) = scl.bad

LinearScaling(lo::Integer, hi::Integer) =
    LinearScaling(lo, hi, get_color(BAD_COLOR))

LinearScaling{T}(lo::Integer, hi::Integer) where {T<:Integer} =
    LinearScaling{T}(lo, hi, get_color(BAD_COLOR))

LinearScaling(; badcolor::ColorOption = BAD_COLOR) =
    LinearScaling{PGInt}(pgqcir()..., get_color(badcolor))

(scl::LinearScaling)(A::AbstractArray, vmin, vmax) =
    rescale(scl, A, vmin, vmax)


have_bad_values(A::AbstractArray{<:Integer}) = false

function have_bad_values(A::AbstractArray{<:Real})
    flag = true
    @inbounds @simd for i in eachindex(A)
        flag &= isfinite(A[i])
    end
    return !flag
end

check(A::AbstractArray) =
    (have_bad_values(A)
     && throw(ArgumentError("non-finite value(s) in argument")))

function check(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) ||
        throw(DimensionMismatch("x and y must have the same length"))
    check(x)
    check(y)
end

"""

```julia
twice(x) -> (x, x)
```

yeilds a 2-tuple with the same value repeated twice.

"""
twice(x) = (x, x)

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
newpage()
```

starts a new page with the standard viewport.  It is possible to specify
the viewport bounds in normalized device coordinate units (NDC):

```julia
newpage(B::Box)
newpage(x1, x2, y1, y2)
```

"""
newpage() = (pgpage(); setviewport())
newpage(B::Box) = (pgpage(); setviewport(B))
newpage(x1::Real, x2::Real, y1::Real, y2::Real) =
    (pgpage(); setviewport(x1, x2, y1, y2))

"""

```julia
getviewport(units=0) -> x1, x2, y1, y2
```

yields the viewport limits in requested units, normalized device coordinate
units (NDC) by default.  To retrieve the viewport as a box:

```julia
getviewport(Box{T}, units=0) -> box
```

"""
getviewport(units::Integer=0) = pgqvp(units)
getviewport(::Type{T}, units::Integer=0) where {T<:Box} = T(getviewport(units))

"""

```julia
setviewport()
```

sets the standard viewport.  It is possible to specify the viewport
bounds in normalized device coordinate units (NDC):

```julia
setviewport(B::Box)
setviewport(x1, x2, y1, y2)
```

"""
setviewport() = pgvstd()
setviewport(B::Box) = setviewport(B.x1, B.x2, B.y1, B.y2)
setviewport(x1::Real, x2::Real, y1::Real, y2::Real) =
    pgsvp(x1, x2, y1, y2)

"""

```julia
getworld() -> x1, x2, y1, y2
```

yields the current world coordinate limits.  To retrieve these limits as a
box:

```julia
getworld(Box{T}) -> box
```

"""
getworld() = Box(pgqwin()...)
getworld(::Type{T}) where {T<:Box} = T(getworld())

"""

```julia
setworld(B::Box, just=false)
setworld(x1, x2, y1, y2, just=false)
```

both set the world coordinate limits.  If `just` is true, the viewport is
adjusted so that both axes have the same scaling.

"""
setworld(B::Box, just::Bool=false) = setworld(B.x1, B.x2, B.y1, B.y2, just)
setworld(x1::Real, x2::Real, y1::Real, y2::Real, just::Bool=false) =
    (just ? pgwnad(x1, x2, y1, y2) : pgswin(x1, x2, y1, y2))


"""

```julia
frame(..., adv=false)
```

draws the environment of subsequent plots.

If `adv` is true, the graphic page is erased and a new page is created
befroe drawing the frame.

"""
frame(B::Box, adv::Bool = false; kwds...) =
    frame(B.x1, B.x2, B.y1, B.y2, adv; kwds...)

function frame(x1::Real, x2::Real, y1::Real, y2::Real,
               adv::Bool = false;
               fig::FigureOption = nothing,
               wait::Bool = false,
               just::Bool = false,
               style::AxisStyleOption = nothing,
               xlabel::LabelOption = nothing,
               ylabel::LabelOption = nothing,
               title::LabelOption = nothing,
               textcolor::ColorOption = nothing,
               framecolor::ColorOption = nothing)

    xopts, yopts = get_axes_style(style)
    select_figure(fig)
    pgsci(get_color(framecolor, DEFAULT_FRAME_COLOR))
    adv && newpage()
    setworld(x1, x2, y1, y2, just)
    pgbox(xopts, 0.0, 0, yopts, 0.0, 0)
    #pgenv(x1, x2, y1, y2, just, get_axis_style(style))
    if xlabel !== nothing || ylabel !== nothing || title !== nothing
        # deal with colors, size, text style
        pgsci(get_color(textcolor, DEFAULT_TITLE_COLOR))
        pglab(get_label(xlabel), get_label(ylabel), get_label(title))
    end
end

function frame(x::AbstractVector, y::AbstractVector, adv::Bool = false;
               extent::ExtentOption = nothing,
               kwds...)
    check(x, y)
    x1, x2 = get_vrange(PGFloat, x, extent, (1,2))
    y1, y2 = get_vrange(PGFloat, y, extent, (3,4))
    frame(x1, x2, y1, y2, adv; kwds...)
end

get_axes_style(opt) = get_axes_style(get_axis_style(opt))

get_axes_style(opt::Integer) =
    (opt == -2 ? twice(" ") :
     opt == -1 ? twice("BC") :
     opt ==  0 ? twice("BCNST") :
     opt ==  1 ? twice("ABCNST") :
     opt ==  2 ? twice("ABCGNST") :
     opt == 10 ? ("BCNSTL", "BCNST") :
     opt == 20 ? ("BCNST",  "BCNSTL") :
     opt == 30 ? ("BCNSTL", "BCNSTL") :
     throw(ArgumentError("illegal AXIS argument")))

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
    frame(x, y, true; kwds...)
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
    frame(x, y, true; kwds...)
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
    frame(x, y, true; kwds...)
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
heatmap(A; kwds...)
```

draws a *heat-map*, that is a map whose cells have pseudo-colors computed
from the values in the 2-dimensional array `A`.  By default, the first and
second dimensions of `A` are assumed to correspond to the abscissa and
ordinate axes respectively and the world coordinates are assumed equal to
the indices of `A` (see below for other possibilities).

Keywords (other keywords are passed to [`frame`](@ref)):

- `vmin` and `vmax` specify the range of values to consider.  It is possible to
  have `vmin > vmax` to reverse the order of the colors.  If unspecified, the
  extreme values of `A` are considered.

- `badcolor` is the color to attribute to bad values (like NaN, *not-a-number*)
  in `A`.

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

To change the default coordiante mapping, the world coordinates `x1`, `x2`,
`y1` and `y2` of the corners of the rectangular region covered by `A` may
be specified:

```julia
heatmap(A, x1, x2, y1, y2; kwds...)
```

The rectangular region covered by `A` has its edges aligned with the world
coordinate axes.  To draw a heat-map of `A` in a region of a different
shape, call:

```julia
heatmap(A, tr; kwds...)
```

with `tr` an affine coordinate transform to map index `(i,j)` to world
coordinates `(x,y)` as follows:

```
x = tr[1] + tr[2]*i + tr[3]*j
y = tr[4] + tr[5]*i + tr[6]*j
```

The coordinate transform `tr` is a 6-element vector or a 2×3 matrix.  Hence
identity is given by: `tr = [0,1,0, 0,0,1]`.

"""
function heatmap(A::AbstractMatrix{T}, args...;
                 fig::FigureOption = nothing,
                 just::Bool = true,
                 cbar::Union{Nothing,Char} = 'R',
                 zlabel::AbstractString = "",
                 cmap::Union{Nothing,AbstractString} = nothing,
                 vmin::Union{Nothing,Real} = nothing,
                 vmax::Union{Nothing,Real} = nothing,
                 badcolor::ColorOption = nothing,
                 kwds...) where {T<:Real}
    # Get size of region covered by the map.
    box = get_extent(Box{PGFloat}, A, args...)

    # Select figure, advance page, set default viewport and set world.
    select_figure(fig)
    newpage()
    setworld(box, just)

    # Select colormap.
    cmap === nothing || palette(cmap)

    # Draw map.
    _heatmap(A, args...; vmin=vmin, vmax=vmax, badcolor=badcolor)

    # Draw the frame.
    frame(box, false; just = just, kwds...)

    # Draw optional color bar.
    cbar === nothing || colorbar(; side = cbar, label = zlabel)
end

function _heatmap(A::AbstractMatrix{<:Real}; kwds...)
    _heatmap(pgarray(PGFloat, A), DEFAULT_XFORM; kwds...)
end

function _heatmap(A::AbstractMatrix{<:Real},
                  tr::AbstractArray{<:Real}; kwds...)
    _heatmap(pgarray(PGFloat, A), get_own_xform(tr); kwds...)
end

function _heatmap(A::DenseMatrix{PGFloat}, tr::DenseArray{PGFloat};
                  vmin::Union{Nothing,Real} = nothing,
                  vmax::Union{Nothing,Real} = nothing,
                  badcolor::ColorOption = nothing)
    bad, lo, hi = rangeofvalues(PGFloat, A, vmin, vmax)
    bad && error("cannot plot image with bad value, call with (x1,x2,y1,y2)")
    memorize_zrange(lo, hi)
    pgimag(A, lo, hi, tr)
end

function _heatmap(A::AbstractMatrix{<:Real},
                  x1::Real, x2::Real, y1::Real, y2::Real; kwds...)
    _heatmap(A, PGFloat(x1), PGFloat(x2), PGFloat(y1), PGFloat(y2); kwds...)
end

function _heatmap(A::AbstractMatrix{<:Integer},
                  x1::PGFloat, x2::PGFloat, y1::PGFloat, y2::PGFloat;
                  vmin::Union{Nothing,Real} = nothing,
                  vmax::Union{Nothing,Real} = nothing,
                  badcolor::ColorOption = BAD_COLOR)
    x1, x2 = adjust_interval(x1, x2)
    y1, y2 = adjust_interval(y1, y2)
    # FIXME: `pgpixl` is much slower than `pgimag` so avoid calling `pgpixl`.
    tr = simple_xform!(get_own_xform(), A, x1, x2, y1, y2)
    bad, lo, hi = rangeofvalues(PGFloat, A, vmin, vmax)
    memorize_zrange(lo, hi)
    @assert bad == false
    pgimag(A, lo, hi, tr)
end

function _heatmap(A::AbstractMatrix{<:AbstractFloat},
                  x1::PGFloat, x2::PGFloat, y1::PGFloat, y2::PGFloat;
                  vmin::Union{Nothing,Real} = nothing,
                  vmax::Union{Nothing,Real} = nothing,
                  badcolor::ColorOption = BAD_COLOR)
    x1, x2 = adjust_interval(x1, x2)
    y1, y2 = adjust_interval(y1, y2)
    select_figure(fig)
    # FIXME: `pgpixl` is much slower than `pgimag` so only call `pgpixl` if
    #        there are bad values.
    bad, lo, hi = rangeofvalues(eltype(A), A, vmin, vmax)
    memorize_zrange(lo, hi)
    if bad
        scaling = LinearScaling(badcolor=badcolor)
        pgpixl(rescale(scaling, A, lo, hi), x1, x2, y1, y2)
    else
        tr = simple_xform!(get_own_xform(), A, x1, x2, y1, y2)
        pgimag(A, lo, hi, tr)
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
function heatmap!(A::AbstractMatrix{<:Real}, args...;
                  fig::FigureOption = nothing, kwds...)
    select_figure(fig)
    _heatmap(A, args...; kwds...)
end

# The following are to remember image settings for the color-bar.
const last_bg = Ref{PGFloat}(0)
const last_fg = Ref{PGFloat}(1)

memorize_zrange(bg::Real, fg::Real) =
    memorize_zrange(PGFloat(bg), PGFloat(fg))

function memorize_zrange(bg::PGFloat, fg::PGFloat)
    if bg == fg
        # Quick fix.
        bg -= tiny(bg)
        fg += tiny(fg)
    end
    last_bg[] = bg
    last_fg[] = fg
    return (bg, fg)
end

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

adjust_interval(x1, x2) =
    (x1 == x2 ? (x1 - tiny(x1), x2 + tiny(x2)) : (x1, x2))

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
mouse([mode = :none, [x0, y0]]; kwds...) -> (x, y, c)
```

yields cursor position `(x,y)` and character `c` typed by the user.  The
position is returned in world coordinates. If the device has no cursor or if
some other error occurs, the value `Char(0)` (ASCII NUL character) is returned
for `c`.  For an X-window device, clicking the 1st mouse button is reported as
an `A`, clicking the 2nd mouse button is reported as a `D` and clicking other
mouse buttons is reported as an `X`.

Argument `mode` can be used to specify the feedback for the user:

- If `mode = :none` (the default) no specific feedback is applied.

- If `mode = :line`, a straight line is drawn joining the anchor point and the
  cursor position.

- If `mode = :box`, a hollow rectangle is extended as the cursor is moved, with
  one vertex at the anchor point `(x0,y0)` and the opposite vertex at the
  current cursor position; the edges of the rectangle are horizontal and
  vertical.

- If `mode = :vrange`, two horizontal lines are extended across the width of
  the display, one drawn at ordinate `y0` of the anchor point and the other
  through the moving cursor position. This could be used to select a Y-axis
  range when one end of the range is known.

- If `mode = :hrange`, two vertical lines are extended over the height of the
  display, one drawn at abscissa `x0` of the anchor point and the other through
  the moving cursor position. This could be used to select an X-axis range when
  one end of the range is known.

- If `mode = :hline`, a horizontal line is extended through the cursor position
  over the width of the display. This could be used to select an X-axis value
  such as the start of an X-axis range. The anchor point is ignored.  :vline

- If `mode = :vline`, a vertical line is extended through the cursor position
  over the height of the display. This could be used to select a Y-axis value
  such as the start of a Y-axis range. The anchor point is ignored.  :hline

- If `mode = :cross`, a cross-hair, centered on the cursor, is extended over
  the width and height of the display. The anchor point is ignored.

Argument `(x0,y0)` specifies the position of the anchor point in world
coordinates.  The anchor point may be omitted if `mode` is not `:line`, `:box`,
`:vrange` or `:hrange`.

Keyword `fig` can be used to specify the figure to consider.  The default
is the current one.

Keyword `color` can be used to specify the color of the line(s) drawn for the
feedback.

The returned cursor position can be a `Point` if the first positional argument
is `Point` or if the anchor is specified as a `Point`.  For instance:

```julia
mouse(Point) -> Point(x, y), c
mouse(Point(x0,y0)) -> Point(x, y), c
```

"""
function mouse(mode::Symbol = :none;
               fig::FigureOption = nothing,
               color::ColorOption = nothing)
    md = (mode == :none  ?  PGInt(0) :
          mode == :hline ?  PGInt(5) :
          mode == :vline ?  PGInt(6) :
          mode == :cross ?  PGInt(7) :
          throw_invalid_mouse_mode(mode))
    select_figure(fig)
    mode == :none || pgsci(get_color(color))
    # Start with the (unused) anchor point at the center of the "window".
    xmin, xmax, ymin, ymax = pgqwin()
    pgband(md, false, PGFloat((xmin + xmax)/2), PGFloat((ymin + ymax)/2))
end

function mouse(mode::Symbol, x0::Real, y0::Real;
               fig::FigureOption = nothing,
               color::ColorOption = nothing)
    md = (mode == :none   ?  PGInt(0) :
          mode == :line   ?  PGInt(1) :
          mode == :box    ?  PGInt(2) :
          mode == :vrange ?  PGInt(3) :
          mode == :hrange ?  PGInt(4) :
          mode == :hline  ?  PGInt(5) :
          mode == :vline  ?  PGInt(6) :
          mode == :cross  ?  PGInt(7) :
          throw_invalid_mouse_mode(mode))
    select_figure(fig)
    mode == :none || pgsci(get_color(color))
    pgband(md, false, PGFloat(x0), PGFloat(y0))
end

function mouse(::Type{Point}, args...; kwds...)
    x, y, c = mouse(args...; kwds...)
    return Point(x, y), c
end

function mouse(P::Point; kwds...)
    x, y, c = mouse(P.x, P.y; kwds...)
    return Point(x, y), c
end

@noinline throw_invalid_mouse_mode(val) =
    throw(ArgumentError(val == :line || val == :box ||
                        val == :vrange || val == :hrange ?
                        string("invalid mouse mode ", val,
                               " or specify anchor point") :
                        string("invalid mouse mode ", val)))

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
splat_xform(T, tr) -> Ax,Axx,Axy, Ay,Ayx,Ayy
```

converts the coordinate transform `tr` into a 6-tuple of coefficients of
type `T`.

"""
function splat_xform(::Type{T},
                     tr::AbstractVector{<:Real}) where {T<:AbstractFloat}
    axes(tr,1) == 1:6 || throw_bad_xform_size()
    return (T(tr[1]), T(tr[2]), T(tr[3]),
            T(tr[4]), T(tr[5]), T(tr[6]))
end

function splat_xform(::Type{T},
                     tr::AbstractMatrix{<:Real}) where {T<:AbstractFloat}
    axes(tr) == (1:2, 1:3) || throw_bad_xform_size()
    return (T(tr[1,1]), T(tr[1,2]), T(tr[1,3]),
            T(tr[2,1]), T(tr[2,2]), T(tr[2,3]))
end


"""

```julia
get_extent(T, I, J, tr) -> x1, x2, y1, y2
```

yields the extent of world coordinates after applying coordinates transform
`tr` to the range of indices `I = imin:imax` and `J = jmin:jmax`.  Argument
`T` is the floating-point type of the result.  It can be replace by a
`Box{T}` type to get the resilt as a box instead of a 4-tuplke.

Ranges of indices can be specified by their endpoints:

```julia
get_extent(i1, i2, j1, j2, tr) -> x1, x2, y1, y2
```

In this latter case, the order of `i1` and `i2 (and of `j1` and `j2`) can be
reversed.

Other possibilities:

```julia
get_extent(T, A, tr)
get_extent(T, I, J, tr)
get_extent(T, i1, i2, j1, j2, tr)
```

or

```julia
get_extent(T, A, x1, x2, y1, y2)
get_extent(T, A, box)
```

"""
get_extent(::Type{Box{T}}, args...) where {T<:AbstractFloat} =
    Box{T}(get_extent(T, args...)...)

get_extent(::Type{T}, A::AbstractMatrix) where {T<:AbstractFloat} =
    get_extent(T, axes(A)...)

function get_extent(::Type{T},
                    A::AbstractMatrix,
                    B::Box) where {T<:AbstractFloat}
    return get_extent(T, A, B.x1, B.x2, B.y1, B.y2)
end

function get_extent(::Type{T},
                    A::AbstractMatrix,
                    x1::Real, x2::Real,
                    y1::Real, y2::Real) where {T<:AbstractFloat}
    return get_extent(T, A, T(x1), T(x2), T(y1), T(y2))
end

function get_extent(::Type{T},
                    A::AbstractMatrix,
                    x1::T, x2::T,
                    y1::T, y2::T) where {T<:AbstractFloat}
    return (adjust_interval(x1, x2)..., adjust_interval(y1, y2)...)
end

function get_extent(::Type{T},
                    A::AbstractMatrix,
                    tr::AbstractArray) where {T<:AbstractFloat}
    return get_extent(T, axes(A)..., tr)
end

function get_extent(::Type{T},
                    i1::Integer, i2::Integer,
                    j1::Integer, j2::Integer) where {T<:AbstractFloat}
    return (_get_extent(T, i1, i2)..., _get_extent(T, j1, j2)...)
end

function get_extent(::Type{T},
                    i1::Integer, i2::Integer,
                    j1::Integer, j2::Integer,
                    tr::AbstractArray{<:Real}) where {T<:AbstractFloat}
    _get_extent(T, _get_extent(T, i1, i2)..., _get_extent(T, j1, j2)..., tr)
end

function get_extent(::Type{T},
                    I::AbstractRange{<:Integer},
                    J::AbstractRange{<:Integer}) where {T<:AbstractFloat}
    return (_get_extent(T, I)..., _get_extent(T, J)...)
end

function get_extent(::Type{T},
                    I::AbstractRange{<:Integer},
                    J::AbstractRange{<:Integer},
                    tr::AbstractArray{<:Real}) where {T<:AbstractFloat}
    _get_extent(T, _get_extent(T, I)..., _get_extent(T, J)..., tr)
end

function _get_extent(::Type{T},
                     i1::Integer, i2::Integer) where {T<:AbstractFloat}
    s = (i1 ≤ i2 ? one(T) : -one(T))/T(2)
    return (T(i1) - s, T(I2) + s)
end

function _get_extent(::Type{T},
                     I::AbstractRange{<:Integer}) where {T<:AbstractFloat}
    s = T(step(I))/T(2)
    return (T(first(I)) - s, T(last(I)) + s)
end

function _get_extent(::Type{T},
                     x1::T, x2::T,
                     y1::T, y2::T,
                     tr::AbstractArray{<:Real}) where {T<:AbstractFloat}
    Ax,Axx,Axy, Ay,Ayx,Ayy = splat_xform(T, tr)
    x1p = Ax + Axx*x1 + Axy*y1
    y1p = Ay + Ayx*x1 + Ayy*y1
    x2p = Ax + Axx*x2 + Axy*y1
    y2p = Ay + Ayx*x2 + Ayy*y1
    x3p = Ax + Axx*x1 + Axy*y2
    y3p = Ay + Ayx*x1 + Ayy*y2
    x4p = Ax + Axx*x2 + Axy*y2
    y4p = Ay + Ayx*x2 + Ayy*y2
    xmin = min(x1p, x2p, x3p, x4p)
    xmax = max(x1p, x2p, x3p, x4p)
    ymin = min(y1p, y2p, y3p, y4p)
    ymax = max(y1p, y2p, y3p, y4p)
    return (xmin, xmax, ymin, ymax)
end

"""

```julia
simple_xform!(tr, A, x1p, x2p, y1p, y2p) -> tr
```

overwrites `tr` with the coefficients of the coordinate transform which maps
the indices of the 2-dimensional array `A` to the rectangular region defined by
`x1p`, `x2p`, `y1p` and `y2p`.

Another possibility is to specify the ranges of indices `I` and `J`:

```julia
simple_xform!(tr, I, J, x1p, x2p, y1p, y2p) -> tr
```

Yet another possibility is to specify the limits `x1`, `x2`, `y1` and `y2` of
the rectangular region to map to the output region:

```julia
simple_xform!(tr, x1, x2, y1, y2, x1p, x2p, y1p, y2p) -> tr
```

"""
function simple_xform!(tr::AbstractArray{T},
                       A::AbstractMatrix,
                       x1p::Real, x2p::Real,
                       y1p::Real, y2p::Real) where {T<:AbstractFloat}
    simple_xform!(tr, axes(A)..., x1p, x2p, y1p, y2p)
end

function simple_xform!(tr::AbstractArray{T},
                       I::AbstractRange{<:Integer},
                       J::AbstractRange{<:Integer},
                       x1p::Real, x2p::Real,
                       y1p::Real, y2p::Real) where {T<:AbstractFloat}
    simple_xform!(tr, I, J, T(x1p), T(x2p), T(y1p), T(y2p))
end

function simple_xform!(tr::AbstractArray{T},
                       I::AbstractRange{<:Integer},
                       J::AbstractRange{<:Integer},
                       x1p::T, x2p::T,
                       y1p::T, y2p::T) where {T<:AbstractFloat}
    di = T(step(I)/2)
    x1 = T(first(I)) - di
    x2 = T(last(I)) + di
    dj = T(step(J)/2)
    y1 = T(first(J)) - dj
    y2 = T(last(J)) + dj
    simple_xform!(tr, x1, x2, y1, y2, x1p, x2p, y1p, y2p)
end

function simple_xform!(tr::DenseArray{T},
                       x1::Real, x2::Real,
                       y1::Real, y2::Real,
                       x1p::Real, x2p::Real,
                       y1p::Real, y2p::Real) where {T<:AbstractFloat}
    simple_xform!(tr, T(x1), T(x2), T(y1), T(y2),
                  T(x1p), T(x2p), T(y1p), T(y2p))
end

function simple_xform!(tr::DenseArray{T},
                       x1::T, x2::T,
                       y1::T, y2::T,
                       x1p::T, x2p::T,
                       y1p::T, y2p::T) where {T<:AbstractFloat}
    @assert length(tr) == 6
    ax = (x2p - x1p)/(x2 - x1)
    bx = (x1p*x2 - x2p*x1)/(x2 - x1)
    @assert isfinite(ax) && isfinite(bx)
    ay = (y2p - y1p)/(y2 - y1)
    by = (y1p*y2 - y2p*y1)/(y2 - y1)
    @assert isfinite(ay) && isfinite(by)
    tr[1] = bx
    tr[2] = ax
    tr[3] = zero(T)
    tr[4] = by
    tr[5] = zero(T)
    tr[6] = ay
    return tr
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

"""

```julia
extremefinitevalues(A) -> (bad, amin, amax)
```

yields the extreme finite values in array `A` and a boolean indicating
whether `A` has bad (i.e., non-finite values).  The returned values `amin`
and `amax` are such that `amin ≤ amax` unless there are no finite values in
`A` in which case `(true,typemax(T),typemin(T))` is returned, with `T` the
type of the elements of `A`.

"""
function extremefinitevalues(A::AbstractArray{T}) where {T<:Real}
    amin = typemax(T)
    amax = typemin(T)
    flag = true
    @inbounds @simd for i in eachindex(A)
        val = A[i]
        flag &= isfinite(val)
        amin = ifelse(isfinite(val) & (val < amin), val, amin)
        amax = ifelse(isfinite(val) & (val > amax), val, amax)
    end
    return (!flag, amin, amax)
end

function extremefinitevalues(A::AbstractArray{T}) where {T<:Integer}
    amin = typemax(T)
    amax = typemin(T)
    @inbounds @simd for i in eachindex(A)
        val = A[i]
        amin = ifelse(val < amin, val, amin)
        amax = ifelse(val > amax, val, amax)
    end
    return (false, amin, amax)
end

"""

```julia
minimumfinitevalue(A) -> (bad, amin)
```

yields the smallest finite value in array `A` and a boolean indicating
whether `A` has bad (i.e., non-finite values).  If the elements of `A` are
of a floating-point type `T`, then `amin` is such that `amin < Inf` unless
there are no finite values in `A` in which case `(true,typemax(T))`, that
is `(true,+T(Inf))`, is returned.

"""
function minimumfinitevalue(A::AbstractArray{T}) where {T<:Real}
    amin = typemax(T)
    flag = true
    @inbounds @simd for i in eachindex(A)
        val = A[i]
        flag &= isfinite(val)
        amin = ifelse(isfinite(val) & (val < amin), val, amin)
    end
    return (!flag, amin)
end

function minimumfinitevalue(A::AbstractArray{T}) where {T<:Integer}
    amin = typemax(T)
    @inbounds @simd for i in eachindex(A)
        val = A[i]
        amin = ifelse(val < amin, val, amin)
    end
    return (false, amin)
end

"""

```julia
maximumfinitevalue(A) -> (bad, amax)
```

yields the largest finite value in array `A` and a boolean indicating
whether `A` has bad (i.e., non-finite values).  If the elements of `A` are
of a floating-point type `T`, then `amax` is such that `amax >-Inf` unless
there are no finite values in `A` in which case `(true,typemin(T))`, that
is `(true,-T(Inf))`, is returned.

"""
function maximumfinitevalue(A::AbstractArray{T}) where {T<:Real}
    amax = typemin(T)
    flag = true
    @inbounds @simd for i in eachindex(A)
        val = A[i]
        flag &= isfinite(val)
        amax = ifelse(isfinite(val) & (val > amax), val, amax)
    end
    return (!flag, amax)
end

function maximumfinitevalue(A::AbstractArray{T}) where {T<:Integer}
    amax = typemin(T)
    @inbounds @simd for i in eachindex(A)
        val = A[i]
        amax = ifelse(val > amax, val, amax)
    end
    return (false, amax)
end

"""

```julia
rangeofvalues(T, A, vmin, vmax) -> (bad, amin, amax)
```

"""
function rangeofvalues(::Type{T}, A::AbstractArray,
                       ::Nothing, ::Nothing) where {T}
    bad, amin, amax = extremefinitevalues(A)
    return (bad, T(amin), T(amax))
end

function rangeofvalues(::Type{T}, A::AbstractArray,
                       vmin::Real, ::Nothing) where {T}
    if isfinite(vmin)
        bad, vmax = maximumfinitevalue(A)
        return (bad, T(vmin), T(vmax))
    else
        isnan(vmin) && throw(ArgumentError("vmin must not be a `NaN`"))
        bad, amin, amax = extremefinitevalues(A)
        if vmin > 0
            # Exchange min. and max.
            amin, amax = amax, amin
        end
        return (bad, T(amin), T(amax))
    end
end

function rangeofvalues(::Type{T}, A::AbstractArray,
                       ::Nothing, vmax::Real) where {T}
    if isfinite(vmax)
        bad, vmin = minimumfinitevalue(A)
        return (bad, T(vmin), T(vmax))
    else
        isnan(vmax) && throw(ArgumentError("vmax must not be a `NaN`"))
        bad, amin, amax = extremefinitevalues(A)
        if vmax < 0
            # Exchange min. and max.
            amin, amax = amax, amin
        end
        return (bad, T(amin), T(amax))
    end
end

function rangeofvalues(::Type{T}, A::AbstractArray,
                       vmin::Real, vmax::Real) where {T}
    choice = (isfinite(vmin) ? 1 : 0) | (isfinite(vmax) ? 2 : 0)
    if choice == 3
        return (have_bad_values(A), T(vmin), T(vmax))
    elseif choice == 2
        # vmin is not finite
        isnan(vmin) && throw(ArgumentError("vmin must not be a `NaN`"))
        if vmin ≤ vmax
            bad, amin = minimumfinitevalue(A)
            return (bad, T(amin), T(vmax))
        else
            bad, amax = maximumfinitevalue(A)
            return (bad, T(amax), T(vmax))
        end
    elseif choice == 1
        # vmax is not finite
        isnan(vmax) && throw(ArgumentError("vmax must not be a `NaN`"))
        if vmin ≤ vmax
            bad, amax = maximumfinitevalue(A)
            return (bad, T(vmin), T(amax))
        else
            bad, amin = minimumfinitevalue(A)
            return (bad, T(vmin), T(amin))
        end
    else
        isnan(vmin) && throw(ArgumentError("vmin must not be a `NaN`"))
        isnan(vmax) && throw(ArgumentError("vmax must not be a `NaN`"))
        bad, amin, amax = extremefinitevalues(A)
        if vmax < vmin
            # Exchange min. and max.
            amin, amax = amax, amin
        end
        return (bad, T(amin), T(amax))
    end
end

"""

```julia
rescale(scl, A, vmin, vmax)
```

yields an array of element type `T` (must be an integer) which maps values of
array `A` from the range `vmin:vmax` to `first(scl):last(scl)`.  It is possible
to have `vmin > vmax` or `first(scl) > last(scl)` to reverse the scaling.

If `A` is of floating point type, `bad(scl)` is the value given to bad numbers
(i.e., *NaN*) and, if `vmin` and `vmax` are not both finite, the result is
filled with `bad(scl)`.

See [`rescale!`](@ref) for providing the destination array and
[`rangeofvalues`](@ref) as a mean to determine `vmin` and `vmax`.

"""
function rescale(scl::AbstractScaling{Ti},
                 A::AbstractArray{<:Integer,N},
                 vmin::Real, vmax::Real) where {Ti<:Integer,N}
    Tf = float(promote_type(typeof(vmin), typeof(vmax)))
    rescale!(Array{Ti,N}(undef, size(A)), scl, A, Tf(vmin), Tf(vmax))
end

function rescale(scl::AbstractScaling{Ti},
                 A::AbstractArray{Tf,N},
                 vmin::Real, vmax::Real) where {Ti<:Integer,
                                                Tf<:AbstractFloat,N}
    rescale!(Array{Ti,N}(undef, size(A)), scl, A, Tf(vmin), Tf(vmax))
end

"""

```julia
rescale!(dst, scl, A, vmin, vmax) -> dst
```

overwrites destination array `dst` with the result or rescaling the values in
array `A` from the range `vmin:vmax` and according to `scl`.  See
[`rescale`](@ref) for details.

"""
function rescale!(dst::DenseArray{Ti,N},
                  scl::AbstractScaling{Ti},
                  A::AbstractArray{<:Integer,N},
                  vmin::Real, vmax::Real) where {Ti<:Integer,N}
    Tf = float(promote_type(typeof(vmin), typeof(vmax)))
    rescale!(dst, scl, A, Tf(vmin), Tf(vmax))
end

function rescale!(dst::DenseArray{Ti,N},
                  scl::AbstractScaling{Ti},
                  A::AbstractArray{Tf,N},
                  vmin::Real, vmax::Real) where {Ti<:Integer,
                                                 Tf<:AbstractFloat,N}
    rescale!(dst, scl, A, Tf(vmin), Tf(vmax))
end

function rescale!(dst::DenseArray{Ti,N},
                  scl::AbstractScaling{Ti},
                  A::AbstractArray{<:Integer,N},
                  vmin::Tf, vmax::Tf) where {Ti<:Integer,
                                             Tf<:AbstractFloat,N}
    throw_not_implemented(typeof(scl))
end

function rescale!(dst::DenseArray{Ti,N},
                  scl::AbstractScaling{Ti},
                  A::AbstractArray{Tf,N},
                  vmin::Tf, vmax::Tf) where {Ti<:Integer,
                                             Tf<:AbstractFloat,N}
    throw_not_implemented(typeof(scl))
end

@noinline throw_not_implemented(::Type{T}) where {T<:AbstractScaling} =
    error("scaling of type $T is not implemented")

function rescale!(dst::DenseArray{Ti,N},
                  scl::LinearScaling{Ti},
                  A::AbstractArray{<:Integer,N},
                  vmin::Tf, vmax::Tf) where {Ti<:Integer,
                                             Tf<:AbstractFloat,N}
    cmin, cmax = first(scl), last(scl)
    if vmin == vmax
        cmid = ((cmin + cmax) ÷ 2)
        @inbounds for i in safe_indices(dst, A)
            val = A[i]
            dst[i] = (val > vmax ? cmax :
                      val < vmin ? cmin : cmid)
        end
    else
        a = (Tf(cmax) - Tf(cmin))/(vmax - vmin)
        b = Tf(cmin) - a*vmin
        @inbounds for i in safe_indices(dst, A)
            val = Tf(A[i])
            dst[i] = (val ≥ vmax ? cmax :
                      val ≤ vmin ? cmin : round(Ti, a*val + b))
        end
    end
    return dst
end

function rescale!(dst::DenseArray{Ti,N},
                  scl::LinearScaling{Ti},
                  A::AbstractArray{Tf,N},
                  vmin::Tf, vmax::Tf) where {Ti<:Integer,
                                             Tf<:AbstractFloat,N}
    cmin, cmax, cbad = first(scl), last(scl), bad(scl)
    if isfinite(vmin) && isfinite(vmin)
        if vmin == vmax
            cmid = ((cmin + cmax) ÷ 2)
            @inbounds for i in safe_indices(dst, A)
                val = A[i]
                dst[i] = (isnan(val) ? cbad :
                          val > vmax ? cmax :
                          val < vmin ? cmin : cmid)
            end
        else
            a = (Tf(cmax) - Tf(cmin))/(vmax - vmin)
            b = Tf(cmin) - a*vmin
            @inbounds for i in safe_indices(dst, A)
                val = A[i]
                dst[i] = (isnan(val) ? cbad :
                          val ≥ vmax ? cmax :
                          val ≤ vmin ? cmin : round(Ti, a*val + b))
            end
        end
    else
        # There are no finite values in `A`.
        fill!(dst, cbad)
    end
    return dst
end

end # module
