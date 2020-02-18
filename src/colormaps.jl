#
# colormaps.jl --
#
# Implements management of colors and colormaps for using with the PGPlot
# library.
#

module Colormaps

export
    RGB,
    palette

using PGPlot.Bindings
import PGPlot.Bindings: pgqcr, pgscr

const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))

struct RGB{T}
    r::T
    g::T
    b::T
end

RGB{T}(col::RGB{T}) where {T} = col
function RGB{T}(col::RGB{UInt8}) where {T<:AbstractFloat}
    a = one(T)/T(255)
    return RGB{T}(a*col.r, a*col.g, a*col.b)
end

function RGB{UInt8}(col::RGB{T}) where {T<:AbstractFloat}
    a = T(255)
    return RGB{T}(round(UInt8, a*clamp(col.r, zero(T), one(T))),
                  round(UInt8, a*clamp(col.g, zero(T), one(T))),
                  round(UInt8, a*clamp(col.b, zero(T), one(T))))
end

black(::Type{RGB{T}}) where {T} = RGB{T}(zero(T), zero(T), zero(T))
white(::Type{RGB{T}}) where {T<:AbstractFloat} = RGB{T}(one(T), one(T), one(T))
white(::Type{RGB{T}}) where {T<:Unsigned} =
    RGB{T}(typemax(T), typemax(T), typemax(T))
gray(::Type{RGB{T}}, f::T) where {T<:AbstractFloat} = RGB{T}(f, f, f)
gray(::Type{RGB{T}}, f::Real) where {T<:AbstractFloat} = gray(RGB{T}, T(f))
gray(::Type{RGB{T}}, f::Real) where {T<:Unsigned} =
    (g = round(T, typemax(T)*clamp(f, oftype(f, 0), oftype(f, 1)));
     RGB{T}(g, g, g))

background(::Type{T}) where {T<:RGB} = pgqcr(T, 0)
foreground(::Type{T}) where {T<:RGB} = pgqcr(T, 1)

# Extend arithmetic operators to allow for simple color operations such as
# linear interpolation.  This is restricted to floating-point colorants.
Base.:( + )(c1::RGB{<:AbstractFloat}, c2::RGB{<:AbstractFloat}) =
    RGB(c1.r + c2.r, c1.g + c2.g, c1.b + c2.b)
Base.:( - )(c1::RGB{<:AbstractFloat}, c2::RGB{<:AbstractFloat}) =
    RGB(c1.r - c2.r, c1.g - c2.g, c1.b - c2.b)
Base.:( * )(α::Real, col::RGB{<:AbstractFloat}) = RGB(α*col.r, α*col.g, α*col.b)
Base.:( * )(col::RGB{<:AbstractFloat}, α::Real) = α*col
Base.:( / )(col::RGB{<:AbstractFloat}, α::Real) = RGB(col.r/α, col.g/α, col.b/α)
Base.:( \ )(α::Real, col::RGB{<:AbstractFloat}) = col/α

Base.Tuple(col::RGB) = (col.r, col.g, col.b)
Base.clamp(col::RGB{T}) where {T<:AbstractFloat} =
    RGB{T}(clamp(col.r, zero(T), one(T)),
           clamp(col.g, zero(T), one(T)),
           clamp(col.b, zero(T), one(T)))

function Base.tryparse(::Type{RGB{T}}, str::AbstractString) where {T}
    cols = split(str, ' ', keepempty=false)
    length(cols) == 3 || return nothing
    r = tryparse(T, cols[1])
    r === nothing && return nothing
    g = tryparse(T, cols[2])
    g === nothing && return nothing
    b = tryparse(T, cols[3])
    b === nothing && return nothing
    return RGB(r,g,b)
end

# Query indexed color as an RGB structure.
pgqcr(::Type{RGB}, ci::Integer) where {T<:AbstractFloat} =
    RGB(pgqcr(ci)...)
pgqcr(::Type{RGB{T}}, ci::Integer) where {T<:AbstractFloat} =
    RGB{T}(pgqcr(ci)...)
pgqcr(::Type{RGB{UInt8}}, ci::Integer) where {T<:AbstractFloat} =
    RGB{UInt8}(pgqcr(RGB, ci))

# Set indexed color with an RGB structure.
pgscr(ci::Integer, col::RGB{UInt8}) = pgscr(PGInt(ci), RGB{PGFloat}(col))
       pgscr(ci::Integer, col::RGB{<:AbstractFloat}) =
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
    lut = Vector{RGB{UInt8}}(undef, 0)
    open(path, "r") do io
        load_gist!(lut, io)
    end
    return lut
end

load_gist(io::IO) = load_gist!(Vector{RGB{UInt8}}(undef, 0), io::IO)

"""

```julia
load_gist!(lut, name) -> lut
```

overwrites the contents of the colormap `lut` with the contents read in Gist
file `name`.

"""
function load_gist!(lut::AbstractVector{T}, io::IO) where {T<:RGB}
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
function palette(ident::Union{AbstractString,AbstractVector{RGB{UInt8}}},
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
        set_color_ramp(cmin, cmax, black(RGB{PGFloat}), RGB{PGFloat}(1,0,0))
    elseif name == "-red"
        set_color_ramp(cmin, cmax, RGB{PGFloat}(1,0,0), black(RGB{PGFloat}))
    elseif name == "green" || name == "+green"
        set_color_ramp(cmin, cmax, black(RGB{PGFloat}), RGB{PGFloat}(0,1,0))
    elseif name == "-green"
        set_color_ramp(cmin, cmax, RGB{PGFloat}(0,1,0), black(RGB{PGFloat}))
    elseif name == "blue" || name == "+blue"
        set_color_ramp(cmin, cmax, black(RGB{PGFloat}), RGB{PGFloat}(0,0,1))
    elseif name == "-blue"
        set_color_ramp(cmin, cmax, RGB{PGFloat}(0,0,1), black(RGB{PGFloat}))
    else
        throw_unknown_colormap(name)
    end
end

function palette(lut::AbstractVector{RGB{UInt8}}, cmin::Int, cmax::Int)
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
set_color_ramp([cmin::Integer, cmax::Integer,] [flag=0 | lo::RGB, hi::RGB])
```

sets the current colormap with a linear ramp of shades of grays or of colors
interpolated between the background and the foreground color or between two
given colors.

Optional arguments `cmin` and `cmax` are to specify the range for the color
indices to set.  If unspecified, the full range of indices used for images
(cmap1) is modified.  Note that `cmin > cmax` is allowed to reverse the
order of colors.

Optional argument `flag` is an integer.  If the least significant bit of `flag`
is set.  Then the colors are reversed; if the second least significant bit of
`flag` is set, background and the foreground colors are interpolated; otherwise,
black and white colors are interpolated.

Two RGB colors, `lo` and `hi`, can be specified instead of `flag` to
interpolate between these two colors.

"""
set_color_ramp(flag::Integer = 0) = set_color_ramp(pgqcir()..., flag)

set_color_ramp(cmin::Integer, cmax::Integer, flag::Integer = 0) =
    set_color_ramp(Int(cmin), Int(cmax), Int(flag))

function set_color_ramp(cmin::Int, cmax::Int, flag::Int = 0)
    if (flag&2) == 1
        # Use background and foreground colors.
        col0 = background(RGB{PGFloat})
        col1 = foreground(RGB{PGFloat})
    else
        # Force black and white.
        col0 = black(RGB{PGFloat})
        col1 = white(RGB{PGFloat})
    end
    if (flag&1) == 1
        # Reverse colors.
        col0, col1 = col1, col0
    end
    set_color_ramp(col0, col1, cmin, cmax)
end

set_color_ramp(lo::RGB, hi::RGB) = set_color_ramp(pgqcir()..., lo, hi)

set_color_ramp(lo::RGB, hi::RGB, cmin::Integer, cmax::Integer) =
    set_color_ramp(cmin, cmax, lo, hi)

set_color_ramp(cmin::Integer, cmax::Integer, lo::RGB, hi::RGB) =
    set_color_ramp(Int(cmin), Int(cmax), RGB{PGFloat}(lo), RGB{PGFloat}(hi))

function set_color_ramp(cmin::Int, cmax::Int,
                        lo::RGB{PGFloat}, hi::RGB{PGFloat})
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
