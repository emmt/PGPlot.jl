module Plotting

using ArrayTools, TwoDimensional
using PGPlot.Bindings
using PGPlot.Bindings: pgarray, DEFAULT_XFORM

# FIXME: type piracy
Base.Tuple(A::AffineTransform) = (A.x, A.xx, A.xy, A.y, A.yx, A.yy)

# Structure used to wrap an integer so that some basic function can be extended
# while avoiding type-piracy.
struct Figure
    num::Int   # Figure number
end
struct Device
    id::PGInt # PGPlot device identifier
end
number(fig::Figure) = fig.num
identifier(dev::Device) = dev.id

const ExtentOption = Union{Nothing,NTuple{4,<:Union{Nothing,Real}}}
const AxisStyleOption = Union{Nothing,Symbol,Integer}
const ColorOption = Union{Nothing,AbstractString,Symbol,Integer}
const FontOption = Union{Nothing,AbstractString,Symbol,Integer}
const LabelOption = Union{Nothing,AbstractString,Symbol}
const XFormOption = Union{Nothing, AffineTransform,NTuple{6,Real},
                          AbstractVector{<:Real},AbstractMatrix{<:Real}}
const HeightOption = Union{Nothing,Real}

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


const DEFAULT_DEVICE = "/XSERVE"
const FIGURES = Dict{Figure,Device}()
const FIGCNT = Ref{Int}(0)

checkfigure() =
    pgqid() > 0 || error("there is no figure open/selected")

"""

```julia
figure() -> id
```

yields the identifier of the current figure, or of the first available one
or creates a new one (with the default device).

FIXME: restrict to interactive figures

""" figure

# Create a new figure using the default device.
figure() = figure(DEFAULT_DEVICE)

# Create a new figure for the given device string.
function figure(dev::AbstractString)
    id = pgopen(dev)
    id ≥ 1 || error(string("cannot open PGPlot device: ", dev))
    pgask(false)
    return register(Device(id))
end

figure(num::Integer) = figure(Figure(num))

function figure(fig::Figure)
    dev = get(FIGURES, fig, Device(-1))
    identifier(dev) ≥ 1 || throw_nonexisting_figure(fig)
    pgslct(identifier(dev))
end

@noinline throw_nonexisting_figure(fig::Figure) =
    throw(ArgumentError(string("Figure #", number(fig), " does not exist")))

function current_figure()
    id = pgqid()
    for (fig, dev) in FIGURES
        if identifier(dev) == id
            return fig
        end
    end
    return Figure(-1)
end

function register(dev::Device)
    if identifier(dev) > 0
        FIGCNT[] += 1
        fig = Figure(FIGCNT[])
        FIGURES[fig] = dev
    else
        fig = Figure(-1)
    end
    return fig
end

Base.isopen(fig::Figure) = number(fig) ≥ 1

function Base.close(fig::Figure)
    dev = get(FIGURES, fig, Device(-1))
    if identifier(dev) < 1
         println(stderr, "Warning: Figure #", number(fig), " does not exist")
    else
         delete!(FIGURES, fig)
         pgclos(identifier(dev))
    end
    nothing
end

get_color(ci::Integer) = PGInt(ci)
get_color(::Nothing) = FOREGROUND_COLOR_INDEX
get_color(ci::Integer, def::Integer) = get_color(ci)
get_color(::Nothing, def::Integer) = get_color(def)


function plot(x::AbstractVector, y::AbstractVector;
              color::ColorOption = nothing,
              kwds...)
    frame(x, y; kwds...)
    _plot(x, y, color)
end

function plot!(x::AbstractVector, y::AbstractVector;
               color::ColorOption = nothing)
    check(x, y)
    _plot(x, y, color)
end

function _plot(x::AbstractVector, y::AbstractVector,
                color::ColorOption)
    pgsci(get_color(color))
    pgline(x, y)
end

function hist(x::AbstractVector, y::AbstractVector;
              center::Bool = true,
              color::ColorOption = nothing,
              kwds...)
    frame(x, y; kwds...)
    _hist(x, y, center, color)
end

function hist!(x::AbstractVector, y::AbstractVector;
               center::Bool = true,
               color::ColorOption = nothing)
    check(x, y)
    _hist(x, y, center, color)
end

function _hist(x, y, center::Bool, color::ColorOption)
    pgsci(get_color(color))
    pgbin(x, y, center)
end

function scatter(x::AbstractVector, y::AbstractVector,
                 sym::Union{Integer,AbstractVector{<:Integer}};
                 color::ColorOption = nothing,
                 kwds...)
    frame(x, y; kwds...)
    _scatter(x, y, sym, color)
end

function scatter!(x::AbstractVector, y::AbstractVector,
                  sym::Union{Integer,AbstractVector{<:Integer}};
                  color::ColorOption = nothing)
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
                 grayscale::Bool = false,
                 vmin::Union{Nothing,Real} = nothing,
                 vmax::Union{Nothing,Real} = nothing,
                 cbar::Union{Nothing,Char} = 'R',
                 zlabel::AbstractString = "",
                 kwds...)
    # Get range of values to plot.
    check(A)
    bg, fg = get_vrange(A, vmin, vmax)

    # Get maximal extent of coordinates.
    I, J = axes(A)
    xmin, xmax, ymin, ymax = get_extent(PGFloat, I, J, tr)

    # Draw viewport.
    frame(xmin, xmax, ymin, ymax; kwds...)

    # Draw map.
    if grayscale
        pggray(A, I, J, fg, bg, tr)
    else
        pgimag(A, I, J, fg, bg, tr)
    end

    # Optional color bar.
    if cbar !== nothing
        colorbar(bg, fg; side = cbar, grayscale = grayscale, label = zlabel)
    end
end

function heatmap!(A::DenseMatrix{PGFloat}, tr::DenseArray{PGFloat};
                  grayscale::Bool = false,
                  vmin::Union{Nothing,Real} = nothing,
                  vmax::Union{Nothing,Real} = nothing)
    check(A)
    bg, fg = get_vrange(A, vmin, vmax)
    I, J = axes(A)
    if grayscale
        pggray(A, I, J, fg, bg, tr)
    else
        pgimag(A, I, J, fg, bg, tr)
    end
end

# FIXME: remember settings of last color/gray image
function colorbar(vmin::Real = PGFloat(0), vmax::Real = PGFloat(1);
                  side::Char = 'R',
                  grayscale::Bool = false,
                  label::AbstractString = "",
                  spacing::Real = PGFloat(2),
                  width::Real = PGFloat(5))
    pgwedg(side*(grayscale ? 'G' : 'I'), spacing, width, vmax, vmin, label)
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

#=
function get_xrange(extent::Union{Nothing,Tuple{Nothing,Nothing,Any,Any}},
                    x::AbstractVector{<:Real})
    nicerange(extrema(x))
end

function get_yrange(extent::Union{Nothing,Tuple{Any,Any,Nothing,Nothing}},
                    y::AbstractVector{<:Real})
    nicerange(extrema(y))
end

get_xrange(extent::NTuple{4,Union{Nothing,Real}}, x::AbstractVector{<:Real}) =
    nicerange(get_xmin(extent, x), get_xmax(extent, x))
get_yrange(extent::NTuple{4,Union{Nothing,Real}}, y::AbstractVector{<:Real}) =
    nicerange(get_ymin(extent, y), get_ymax(extent, y))

get_xmin(extent::Tuple{Nothing,Any,Any,Any}, x::AbstractVector{<:Real}) = minimum(x)
get_xmin(extent::Tuple{Real,Any,Any,Any}, x::AbstractVector{<:Real}) = extent[1]
get_xmax(extent::Tuple{Any,Nothing,Any,Any}, x::AbstractVector{<:Real}) = maximum(x)
get_xmax(extent::Tuple{Any,Real,Any,Any}, x::AbstractVector{<:Real}) = extent[2]

get_ymin(extent::Tuple{Any,Any,Nothing,Any}, y::AbstractVector{<:Real}) = minimum(y)
get_ymin(extent::Tuple{Real,Any,Any,Any}, y::AbstractVector{<:Real}) = extent[3]
get_ymax(extent::Tuple{Any,Any,Any,Nothing}, y::AbstractVector{<:Real}) = maximum(y)
get_ymax(extent::Tuple{Any,Any,Any,Real}, y::AbstractVector{<:Real}) = extent[4]
=#

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
