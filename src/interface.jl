#
# interface.jl --
#
# Implements high-level interface for plotting with the PGPlot library.
#

module Plotting

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

using ArrayTools, TwoDimensional
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
const ColorOption = Union{Nothing,AbstractString,Symbol,Integer}
const FontOption = Union{Nothing,AbstractString,Symbol,Integer}
const LabelOption = Union{Nothing,AbstractString,Symbol}
const XFormOption = Union{Nothing, AffineTransform,NTuple{6,Real},
                          AbstractVector{<:Real},AbstractMatrix{<:Real}}
const HeightOption = Union{Nothing,Real}
const FigureOption = Union{Nothing,Figure,Integer}

#=
PGQCI -- inquire color index
PGQCIR -- inquire color index range
PGQCOL -- inquire color capability
PGQCR -- inquire color representation
PGQTBG -- inquire text background color index
PGSCI -- set color index
PGSCIR -- set color index range
PGSCR -- set color representation
PGSCRN -- set color representation by name
PGSHLS -- set color representation using HLS system
PGSTBG -- set text background color index

PGETXT -- erase text from graphics display
PGMTXT -- write text at position relative to viewport
PGPTXT -- write text at arbitrary position and angle
PGQTXT -- find bounding box of text string
PGTEXT -- write text (horizontal, left-justified)

Character height:
PGQCH -- inquire character height
PGQCS -- inquire character height in a variety of units
PGSCH -- set character height

Text font:
PGQCF -- inquire character font
PGSCF -- set character font
=#

#=
Set the Character Font for subsequent text plotting. Four different
fonts are available:
  1: (default) a simple single-stroke font ("normal" font)
  2: roman font
  3: italic font
  4: script font
This call determines which font is in effect at the beginning of
each text string. The font can be changed (temporarily) within a text
string by using the escape sequences \fn, \fr, \fi, and \fs for fonts
1, 2, 3, and 4, respectively.
=#
const NORMAL_FONT = PGInt(1)
const ROMAN_FONT  = PGInt(2)
const ITALIC_FONT = PGInt(3)
const SCRIPT_FONT = PGInt(4)

function rgb_color(col::UInt32)
    mask = 0x000000FF
    mult = 1/PGFloat(255)
    return (((col >> 16) & mask)*mult,
            ((col >>  8) & mask)*mult,
            ( col        & mask)*mult)
end

const ALT_RED      = rgb_color(0x00c5000b)
const ALT_GREEN    = rgb_color(0x00579d1c)
const ALT_BLUE     = rgb_color(0x000084d1)
const LIGHT_BLUE   = rgb_color(0x0083caff)
const ALT_VIOLET   = rgb_color(0x004b1f6f)
const DARK_BLUE    = rgb_color(0x00004586)
const ALT_ORANGE   = rgb_color(0x00ff420e)
const ALT_YELLOW   = rgb_color(0x00ffd320)
const DARK_VIOLET  = rgb_color(0x007e0021)
const DARK_GREEN   = rgb_color(0x00314004)
const LIGHT_GREEN  = rgb_color(0x00aecf00)
const GOLDEN       = rgb_color(0x00ff950e)

const BLACK        = rgb_color(0x00000000)
const WHITE        = rgb_color(0x00ffffff)
const RED          = rgb_color(0x00ff0000)
const GREEN        = rgb_color(0x0000ff00)
const BLUE         = rgb_color(0x000000ff)
const CYAN         = rgb_color(0x0000ffff)
const MAGENTA      = rgb_color(0x00ff00ff)
const YELLOW       = rgb_color(0x00ffff00)
const ORANGE       = rgb_color(0x00ff8000)
const YELLOW_GREEN = rgb_color(0x0080ff00)
const BLUE_GREEN   = rgb_color(0x0000ff80)
const DEEP_SKY     = rgb_color(0x000080ff)
const BLUE_VIOLET  = rgb_color(0x008000ff)
const DEEP_PINK    = rgb_color(0x00ff0080)
const DARK_GRAY    = rgb_color(0x00555555)
const LIGHT_GRAY   = rgb_color(0x00aaaaaa)

function set_cmap0(dark::Bool = true)
    if dark
        pgscr( 0, BLACK...)
        pgscr( 1, WHITE...)
        pgscr( 2, RED...)
        pgscr( 3, GREEN...)
        pgscr( 4, BLUE...)
        pgscr( 5, CYAN...)
        pgscr( 6, MAGENTA...)
        pgscr( 7, YELLOW...)
        pgscr( 8, ORANGE...)
        pgscr( 9, YELLOW_GREEN...)
        pgscr(10, BLUE_GREEN...)
        pgscr(11, DEEP_SKY...)
        pgscr(12, BLUE_VIOLET...)
        pgscr(13, DEEP_PINK...)
    else
        pgscr( 0, WHITE...)
        pgscr( 1, BLACK...)
        pgscr( 2, ALT_RED...)
        pgscr( 3, ALT_GREEN...)
        pgscr( 4, ALT_BLUE...)
        pgscr( 5, LIGHT_BLUE...)
        pgscr( 6, MAGENTA...)
        pgscr( 7, ALT_YELLOW...)
        pgscr( 8, ALT_ORANGE...)
        pgscr( 9, GOLDEN...)
        pgscr(10, BLUE_GREEN...)
        pgscr(11, DEEP_SKY...)
        pgscr(12, BLUE_VIOLET...)
        pgscr(13, DEEP_PINK...)
    end
    pgscr(14, DARK_GRAY...)
    pgscr(15, LIGHT_GRAY...)
end

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
        pgscr(0, 0,0,0) # set background color to white
        pgscr(1, 1,1,1) # set foreground color to black
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

get_color(ci::Integer) = PGInt(ci)
get_color(::Nothing) = FOREGROUND_COLOR_INDEX
get_color(ci::Integer, def::Integer) = get_color(ci)
get_color(::Nothing, def::Integer) = get_color(def)


function plot(x::AbstractVector, y::AbstractVector;
              fig::FigureOption = nothing,
              color::ColorOption = nothing,
              kwds...)
    select_figure(fig)
    frame(x, y; kwds...)
    _plot(x, y, color)
end

function plot!(x::AbstractVector, y::AbstractVector;
               fig::FigureOption = nothing,
               color::ColorOption = nothing)
    select_figure(fig)
    check(x, y)
    _plot(x, y, color)
end

function _plot(x::AbstractVector, y::AbstractVector,
                color::ColorOption)
    pgsci(get_color(color))
    pgline(x, y)
end

function hist(x::AbstractVector, y::AbstractVector;
              fig::FigureOption = nothing,
              center::Bool = true,
              color::ColorOption = nothing,
              kwds...)
    select_figure(fig)
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

function _hist(x, y, center::Bool, color::ColorOption)
    pgsci(get_color(color))
    pgbin(x, y, center)
end

function scatter(x::AbstractVector, y::AbstractVector,
                 sym::Union{Integer,AbstractVector{<:Integer}};
                 fig::FigureOption = nothing,
                 color::ColorOption = nothing,
                 kwds...)
    select_figure(fig)
    frame(x, y; kwds...)
    _scatter(x, y, sym, color)
end

function scatter!(x::AbstractVector, y::AbstractVector,
                  sym::Union{Integer,AbstractVector{<:Integer}};
                  fig::FigureOption = nothing,
                  color::ColorOption = nothing)
    select_figure(fig)
    check(x, y)
    _scatter(x, y, sym, color)
end

function _scatter(x::AbstractVector, y::AbstractVector,
                  sym::Union{Integer,AbstractVector{<:Integer}},
                  color::ColorOption)
    pgsci(get_color(color))
    pgpnts(x, y, sym)
end


# @btime Plotting.heatmap($z) for a 81×81 image
#   8.112 ms (2 allocations: 25.77 KiB)

function heatmap(A::AbstractMatrix{<:Real},
                 xform::XFormOption = DEFAULT_XFORM;
                 kwds...)
    heatmap(pgarray(PGFloat, A), get_xform(xform); kwds...)
end


function heatmap!(A::AbstractMatrix{<:Real},
                  xform::XFormOption = DEFAULT_XFORM;
                  kwds...)
    heatmap!(pgarray(PGFloat, A), pgarray(PGFloat, xform); kwds...)
end

function heatmap(A::DenseMatrix{PGFloat}, tr::DenseArray{PGFloat};
                 fig::FigureOption = nothing,
                 vmin::Union{Nothing,Real} = nothing,
                 vmax::Union{Nothing,Real} = nothing,
                 cbar::Union{Nothing,Char} = 'R',
                 zlabel::AbstractString = "",
                 cmap::Union{Nothing,AbstractString} = nothing,
                 just::Bool = true,
                 kwds...)
    # Select figure.
    select_figure(fig)

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
const last_bg = Ref{PGFloat}(1)
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

function colorbar(vmin::Real = last_bg[],
                  vmax::Real = last_fg[];
                  side::Char = 'R',
                  grayscale::Bool = false,
                  label::AbstractString = "",
                  spacing::Real = PGFloat(2),
                  width::Real = PGFloat(5))
    pgwedg(side*(grayscale ? 'G' : 'I'), spacing, width, vmin, vmax, label)
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
frame(...)
```

draws the environment of subsequent plots.

"""
function frame(xmin::Real, xmax::Real, ymin::Real, ymax::Real;
               wait::Bool = false,
               just::Bool = false,
               style::AxisStyleOption = nothing,
               xlabel::LabelOption = nothing,
               ylabel::LabelOption = nothing,
               title::LabelOption = nothing,
               textcolor::ColorOption = nothing,
               framecolor::ColorOption = nothing)
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

get_text_height(val::Real) =
    (@assert isfinite(val) && val > 0; PGFloat(val))

get_font(val::AbstractString) = get_text_font(Symbol(val))
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
