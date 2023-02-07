#
# colormaps.jl --
#
# Implements management of colors and colormaps for using with the PGPlot
# library.
#

module Colormaps

export
    RGBVec,
    palette

using Colors
using PGPlot.Bindings
import PGPlot.Bindings: pgqcr, pgscr

const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))

"""

`RGBVec{T}(r,g,b)` represents an RGB color whose components have type `T`.

When `T` is a floating-point, the color triplet can be manipulated as a
3-element *vector*, that is `α*v + β*v` for any reals `α` and `β` and
colors `u` and `v` of type `RGBVec{<:AbstractFloat}` yields a color
of type `RGBVec{<:AbstractFloat}`.  This is useful to interpolate colors
and build colormaps.

Having `T = UInt8` can be used to parse/convert colors.

"""
struct RGBVec{T}
    r::T
    g::T
    b::T
end

Base.eltype(::RGBVec{T}) where {T} = T
RGBVec{T}(col::RGBVec{T}) where {T} = col

function RGBVec{T}(col::RGBVec{UInt8}) where {T<:AbstractFloat}
    a = one(T)/T(255)
    return RGBVec{T}(a*col.r, a*col.g, a*col.b)
end

function RGBVec{UInt8}(col::RGBVec{T}) where {T<:AbstractFloat}
    a = T(255)
    return RGBVec{T}(round(UInt8, a*clamp(col.r, zero(T), one(T))),
                     round(UInt8, a*clamp(col.g, zero(T), one(T))),
                     round(UInt8, a*clamp(col.b, zero(T), one(T))))
end

RGBVec{T}(rgb::NTuple{3,Real}) where {T} = RGB{T}(rgb...)
RGBVec(rgb::NTuple{3,Real}) = RGB(rgb...)

RGBVec{T}(col::Colorant) where {T} = RGBVec{T}(RGB(col))
RGBVec(col::Colorant) = RGBVec(RGB(col))

RGBVec{T}(col::RGB) where {T} = RGBVec{T}(col.r, col.g, col.b)
RGBVec(col::RGB) = RGBVec(col.r, col.g, col.b)

RGBVec{T}(col::Symbol) where {T} = RGBVec{T}(String(col))
RGBVec(col::Symbol) = RGBVec(String(col))

RGBVec{T}(col::AbstractString) where {T} = RGBVec{T}(parse(RGB, col))
RGBVec(col::AbstractString) = RGBVec(parse(RGB, col))

Base.convert(::Type{T}, arg::T) where {T<:RGBVec} = arg
Base.convert(::Type{T}, arg::Colorant) where {T<:RGBVec} = T(arg)
Base.convert(::Type{T}, arg::NTuple{3,Real}) where {T<:RGBVec} = T(arg)
Base.convert(::Type{T}, arg::RGBVec) where {T<:RGB} = T(arg)

black(::Type{RGBVec{T}}) where {T} = RGBVec{T}(zero(T), zero(T), zero(T))
white(::Type{RGBVec{T}}) where {T<:AbstractFloat} =
    RGBVec{T}(one(T), one(T), one(T))
white(::Type{RGBVec{T}}) where {T<:Unsigned} =
    RGBVec{T}(typemax(T), typemax(T), typemax(T))
gray(::Type{RGBVec{T}}, f::T) where {T<:AbstractFloat} = RGBVec{T}(f, f, f)
gray(::Type{RGBVec{T}}, f::Real) where {T<:AbstractFloat} =
    gray(RGBVec{T}, T(f))
gray(::Type{RGBVec{T}}, f::Real) where {T<:Unsigned} =
    (g = round(T, typemax(T)*clamp(f, oftype(f, 0), oftype(f, 1)));
     RGBVec{T}(g, g, g))

background(::Type{T}) where {T<:RGBVec} = pgqcr(T, 0)
foreground(::Type{T}) where {T<:RGBVec} = pgqcr(T, 1)

# Extend arithmetic operators to allow for simple color operations such as
# linear interpolation.  This is restricted to floating-point colorants.
Base.:( + )(c1::RGBVec{T1}, c2::RGBVec{T2}) where {T1<:Real,T2<:Real} =
    (T = float(promote_type(T1, T2)); RGBVec{T}(c1) + RGBVec{T}(c2))
Base.:( + )(c1::RGBVec{<:AbstractFloat}, c2::RGBVec{<:AbstractFloat}) =
    RGBVec(c1.r + c2.r, c1.g + c2.g, c1.b + c2.b)
Base.:( - )(c1::RGBVec{T1}, c2::RGBVec{T2}) where {T1<:Real,T2<:Real} =
    (T = float(promote_type(T1, T2)); RGBVec{T}(c1) - RGBVec{T}(c2))
Base.:( - )(c1::RGBVec{<:AbstractFloat}, c2::RGBVec{<:AbstractFloat}) =
    RGBVec(c1.r - c2.r, c1.g - c2.g, c1.b - c2.b)
Base.:( * )(col::RGBVec, α::Real) = α*col
Base.:( * )(α::Real, col::RGBVec{T}) where {T<:Real} = α*RGB{float(T)}(col)
Base.:( * )(α::Real, col::RGBVec{<:AbstractFloat}) =
    RGBVec(α*col.r, α*col.g, α*col.b)
Base.:( / )(col::RGBVec{T}, α::Real) where {T<:Real} = RGB{float(T)}(col)/α
Base.:( / )(col::RGBVec{<:AbstractFloat}, α::Real) =
    RGBVec(col.r/α, col.g/α, col.b/α)
Base.:( \ )(α::Real, col::RGBVec) = col/α

Base.Tuple(col::RGBVec) = (col.r, col.g, col.b)
Base.clamp(col::RGBVec{T}) where {T<:AbstractFloat} =
    RGBVec{T}(clamp(col.r, zero(T), one(T)),
              clamp(col.g, zero(T), one(T)),
              clamp(col.b, zero(T), one(T)))

function Base.tryparse(::Type{RGBVec{T}}, str::AbstractString) where {T}
    cols = split(str, ' ', keepempty=false)
    length(cols) == 3 || return nothing
    r = tryparse(T, cols[1])
    r === nothing && return nothing
    g = tryparse(T, cols[2])
    g === nothing && return nothing
    b = tryparse(T, cols[3])
    b === nothing && return nothing
    return RGBVec(r,g,b)
end

# Query indexed color as an RGBVec structure.
pgqcr(::Type{RGBVec}, ci::Integer) =
    RGBVec(pgqcr(ci)...)
pgqcr(::Type{RGBVec{T}}, ci::Integer) where {T<:AbstractFloat} =
    RGBVec{T}(pgqcr(ci)...)
pgqcr(::Type{RGBVec{UInt8}}, ci::Integer) =
    RGBVec{UInt8}(pgqcr(RGBVec, ci))

# Set indexed color with an RGBVec structure.
pgscr(ci::Integer, col::RGBVec{UInt8}) = pgscr(PGInt(ci), RGBVec{PGFloat}(col))
pgscr(ci::Integer, col::RGBVec{<:AbstractFloat}) =
    pgscr(PGInt(ci), col.r, col.g, col.b)

"""

```julia
find_file(name) -> path
```

yields the path to a readable graphics file.

"""
find_file(name::AbstractString) = find_file(convert(String, name))
function find_file(name::String)
    if isfile(name)
        return name
    else
        path = joinpath(DATA_DIR, name)
        if isfile(path)
            return path
        else
            throw_file_not_found(name)
        end
    end
end

@noinline throw_file_not_found(name::AbstractString) =
    throw(ArgumentError(string("file \"", name, "\" not found")))

"""

```julia
load_gist(name) -> lut
```

yields the colormap read in Gist file `name`.

"""
function load_gist(name::AbstractString)
    path = find_file(name)
    lut = Vector{RGBVec{UInt8}}(undef, 0)
    open(path, "r") do io
        load_gist!(lut, io)
    end
    return lut
end

load_gist(io::IO) = load_gist!(Vector{RGBVec{UInt8}}(undef, 0), io::IO)

"""

```julia
load_gist!(lut, name) -> lut
```

overwrites the contents of the colormap `lut` with the contents read in Gist
file `name`.

"""
function load_gist!(lut::AbstractVector{T}, io::IO) where {T<:RGBVec}
    resize!(lut, 0)
    while !eof(io)
        line = readline(io)
        rgb = tryparse(T, line)
        if rgb !== nothing
            push!(lut, rgb)
        end
    end
    return lut
end

"""

```julia
palette(cmap)
```

installs the colormap `cmap` (a name or a look-up table) in the current
plotting device.

The color index range may be specified:

```julia
palette(cmap, cmin, cmax)
```

If `cmin:cmax` is larger than the current index range, an attempt is made to
enlarge it.

Also see [`set_color_ramp`](@ref).

"""
function palette(ident::Union{AbstractString,AbstractVector{RGBVec{UInt8}}},
                 cmin::Union{Nothing,Integer} = nothing,
                 cmax::Union{Nothing,Integer} = nothing)
    palette(ident, get_color_index_range(cmin, cmax)...)
end

function palette(name::AbstractString, cmin::Int, cmax::Int)
    if endswith(name, ".gp")
        lut = load_gist(name)
        palette(lut, cmin, cmax)
    elseif name == "gray" || name == "+gray"
        set_color_ramp(cmin, cmax, 0)
    elseif name == "-gray"
        set_color_ramp(cmin, cmax, 1)
    elseif name == "bg-fg"
        set_color_ramp(cmin, cmax, 2)
    elseif name == "fg-bg"
        set_color_ramp(cmin, cmax, 3)
    elseif name == "red" || name == "+red"
        set_color_ramp(cmin, cmax, black(RGBVec{PGFloat}),
                       RGBVec{PGFloat}(1,0,0))
    elseif name == "-red"
        set_color_ramp(cmin, cmax, RGBVec{PGFloat}(1,0,0),
                       black(RGBVec{PGFloat}))
    elseif name == "green" || name == "+green"
        set_color_ramp(cmin, cmax, black(RGBVec{PGFloat}),
                       RGBVec{PGFloat}(0,1,0))
    elseif name == "-green"
        set_color_ramp(cmin, cmax, RGBVec{PGFloat}(0,1,0),
                       black(RGBVec{PGFloat}))
    elseif name == "blue" || name == "+blue"
        set_color_ramp(cmin, cmax, black(RGBVec{PGFloat}),
                       RGBVec{PGFloat}(0,0,1))
    elseif name == "-blue"
        set_color_ramp(cmin, cmax, RGBVec{PGFloat}(0,0,1),
                       black(RGBVec{PGFloat}))
    else
        throw_unknown_colormap(name)
    end
end

function palette(lut::AbstractVector{RGBVec{UInt8}}, cmin::Int, cmax::Int)
    length(lut) > 0 || error("no colors!")
    f = 1/255
    I = axes(lut, 1)
    imin, imax = Int(first(I)), Int(last(I))
    if cmin != cmax
        a = (imax - imin)/(cmax - cmin)
        for c in min(cmin,cmax):max(cmin,cmax)
            t = (c - cmin)*a + imin
            i0 = floor(Int, t)
            i1 = min(i0 + 1, imax)
            a1 = t - i0
            a0 = one(a1) - a1
            r = a0*lut[i0].r + a1*lut[i1].r
            g = a0*lut[i0].g + a1*lut[i1].g
            b = a0*lut[i0].b + a1*lut[i1].b
            pgscr(c, f*r, f*g, f*b)
        end
    else
        i = ((imax + imin + 1) >> 1)
        pgscr(c, f*lut[i].r, f*lut[i].g, f*lut[i].b)
    end
end

@noinline throw_unknown_colormap(name::AbstractString) =
    throw(ArgumentError(string("unknown colormap \"", name, "\"")))

get_color_index_range(::Nothing, ::Nothing) = get_color_index_range()

function get_color_index_range()
    cmin, cmax = pgqcir()
    return (Int(cmin), Int(cmax))
end

function get_color_index_range(cmin::Union{Nothing,Integer},
                               cmax::Union{Nothing,Integer})
    qmin, qmax = get_color_index_range()
    rmin, rmax = qmin, qmax
    if cmin !== nothing
        rmin = oftype(rmin, cmin)
    end
    if cmax !== nothing
        rmax = oftype(rmax, cmax)
    end
    if min(rmin, rmax) < qmin || max(rmin, rmax) > qmax
        pgscir(min(rmin, rmax), max(rmin, rmax))
    end
    return (rmin, rmax)
end

"""

```julia
set_color_ramp([cmin::Integer, cmax::Integer,]
               [flag=0 | lo::RGBVec, hi::RGBVec])
```

sets the current colormap with a linear ramp of shades of grays or of colors
interpolated between the background and the foreground color or between two
given colors.

Optional arguments `cmin` and `cmax` are to specify the range for the color
indices to set.  If unspecified, the full range of indices used for images
(cmap1) is modified.  Note that `cmin > cmax` is allowed to reverse the order
of colors.

Optional argument `flag` is an integer.  If the least significant bit of `flag`
is set.  Then the colors are reversed; if the second least significant bit of
`flag` is set, background and the foreground colors are interpolated;
otherwise, black and white colors are interpolated.

Two RGBVec colors, `lo` and `hi`, can be specified instead of `flag` to
interpolate between these two colors.

"""
set_color_ramp(flag::Integer = 0) = set_color_ramp(pgqcir()..., flag)

set_color_ramp(cmin::Integer, cmax::Integer, flag::Integer = 0) =
    set_color_ramp(Int(cmin), Int(cmax), Int(flag))

function set_color_ramp(cmin::Int, cmax::Int, flag::Int = 0)
    if (flag&2) == 1
        # Use background and foreground colors.
        col0 = background(RGBVec{PGFloat})
        col1 = foreground(RGBVec{PGFloat})
    else
        # Force black and white.
        col0 = black(RGBVec{PGFloat})
        col1 = white(RGBVec{PGFloat})
    end
    if (flag&1) == 1
        # Reverse colors.
        col0, col1 = col1, col0
    end
    set_color_ramp(col0, col1, cmin, cmax)
end

set_color_ramp(lo::RGBVec, hi::RGBVec) = set_color_ramp(pgqcir()..., lo, hi)

set_color_ramp(lo::RGBVec, hi::RGBVec, cmin::Integer, cmax::Integer) =
    set_color_ramp(cmin, cmax, lo, hi)

set_color_ramp(cmin::Integer, cmax::Integer, lo::RGBVec, hi::RGBVec) =
    set_color_ramp(Int(cmin), Int(cmax),
                   RGBVec{PGFloat}(lo), RGBVec{PGFloat}(hi))

function set_color_ramp(cmin::Int, cmax::Int,
                        lo::RGBVec{PGFloat}, hi::RGBVec{PGFloat})
    lo = clamp(lo)
    hi = clamp(hi)
    if cmin == cmax
        # Set all color indices to the mean level.
        col = (lo + hi)/2
        for c in min(cmin,cmax):max(cmin,cmax)
            pgscr(c, col)
        end
    else
        # Interpolate the given colors.
        f = one(PGFloat)/PGFloat(cmax - cmin)
        for c in min(cmin,cmax):max(cmin,cmax)
            a1 = (c - cmin)*f
            a0 = one(a1) - a1
            pgscr(c, a0*lo + a1*hi)
        end
    end
end

function set_standard_colors()
    pgscr(0, 0.0,0.0,0.0)
    pgscr(1, 1.0,1.0,1.0)
    pgscr(2, 1.0,0.0,0.0)
    pgscr(3, 0.0,1.0,0.0)
    pgscr(4, 0.0,0.0,1.0)
    pgscr(5, 0.0,1.0,1.0)
    pgscr(6, 1.0,0.0,1.0)
    pgscr(7, 1.0,1.0,0.0)
    pgscr(8, 1.0,0.5,0.0)
    pgscr(9, 0.5,1.0,0.0)
    pgscr(10, 0.0,1.0,0.5)
    pgscr(11, 0.0,0.5,1.0)
    pgscr(12, 0.5,0.0,1.0)
    pgscr(13, 1.0,0.0,0.5)
    pgscr(14, 0.333,0.333,0.333)
    pgscr(15, 0.667,0.667,0.667)
end

end # module
