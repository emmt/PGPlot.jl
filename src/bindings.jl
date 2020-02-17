#
# bindings.jl --
#
# Implements low-level bindings to PGPlot library.
#

module Bindings

# --------------------------------
# PGPlot                 PGPlot.jl
# --------------------------------
# pgarro                 pgarro
# pgenv                  pgenv
# pglab                  pglab
# pgline                 pgline
# pgpt, pgpt1, pgpnts    pgpt
# pgpage                 pgpage
# pgsvp                  pgsvp
# pgswin                 pgswin
# pgbox                  pgbox
# pgmtxt                 pgmtxt
# pgmove                 pgmove
# pgbbuf / pgebuf        pgbbuf / pgebuf
# pgsls                  pgsls
# pgslw                  pgslw
# pgcont                 pgcont
# pgfunt                 pgfunct
# pgfunx                 pgfuncx
# pgfuny                 pgfuncy
# --------------------------------

export
    PGArray,
    PGBool,
    PGFloat,
    PGInt,
    PGMatrix,
    PGPlotError,
    PGSymbol,
    PGVector,
    pgarro,
    pgask,
    pgaxis,
    pgband,
    pgbbuf,
    pgbeg,
    pgbin,
    pgbox,
    pgcirc,
    pgclos,
    pgconb,
    pgconf,
    pgconl,
    pgcons,
    pgcont,
    pgctab,
    pgcurs,
    pgdraw,
    pgebuf,
    pgend,
    pgenv,
    pgeras,
    pgerr1,
    pgerrb,
    pgerrx,
    pgerry,
    pgetxt,
    pgfunt,
    pgfunx,
    pgfuny,
    pggray,
    pghi2d,
    pghist,
    pgiden,
    pgimag,
    pglab,
    pglcur!,
    pgldev,
    pglen,
    pgline,
    pgmove,
    pgmtxt,
    pgncur,
    pgnumb,
    pgolin,
    pgopen,
    pgpage,
    pgpanl,
    pgpap,
    pgpixl,
    pgpnts,
    pgpoly,
    pgpt,
    pgptxt,
    pgqah,
    pgqcf,
    pgqch,
    pgqci,
    pgqcir,
    pgqclp,
    pgqcol,
    pgqcr,
    pgqcs,
    pgqdt,
    pgqfs,
    pgqhs,
    pgqid,
    pgqinf,
    pgqitf,
    pgqls,
    pgqlw,
    pgqndt,
    pgqpos,
    pgqtbg,
    pgqtxt,
    pgqvp,
    pgqvsz,
    pgqwin,
    pgrect,
    pgrnd,
    pgrnge,
    pgsah,
    pgsave,
    pgscf,
    pgsch,
    pgsci,
    pgscir,
    pgsclp,
    pgscr,
    pgscrl,
    pgscrn,
    pgsfs,
    pgshls,
    pgshs,
    pgsitf,
    pgslct,
    pgsls,
    pgslw,
    pgstbg,
    pgsubp,
    pgsvp,
    pgswin,
    pgtbox,
    pgtext,
    pgtick,
    pgunsa,
    pgupdt,
    pgvect,
    pgvsiz,
    pgvstd,
    pgwedg,
    pgwnad

isfile(joinpath(@__DIR__, "..", "deps", "deps.jl")) ||
    error("PGPlot not properly installed.  Please run Pkg.build(\"PGPlot\")")
include(joinpath("..", "deps", "deps.jl"))

# Basic types used by PGPlot library.  These definitions are to make changes easy.
const PGBool     = Cint
const PGInt      = Cint
const PGFloat    = Cfloat
const PGChar     = Cchar
const PGString   = Cstring
const PGIntegers = Union{Integer,Char} # used for symbols, can be converted to PGInt

# Arrays passed to PGPlot must have contiguous elements in colum-major order.
const PGArray{T,N}  = DenseArray{T,N}
const PGVector{T}   = PGArray{T,1}
const PGMatrix{T}   = PGArray{T,2}
const PGVecOrRef{T} = Union{Ref{T},PGVector{T}}

# Custom exception to report errors.
struct PGPlotError <: Exception
    func::Symbol
    code::PGInt
end

pgerror(func::Symbol, code::Integer = PGInt(0)) =
    throw(PGPlotError(func, code))

# Convert returned values.
output(x::AbstractFloat) = Float64(x)
output(x::Integer) = Int(x)
output(x::PGChar) = Char(x)
output(x::Ref) = output(x[])
output(x1, x2) = (output(x1), output(x2))
output(x1, x2, x3) = (output(x1), output(x2), output(x3))
output(x1, x2, x3, x4) = (output(x1), output(x2), output(x3), output(x4))
output(::Type{Bool}, x::PGInt) = (x != 0)
output(::Type{Bool}, x::Ref{PGInt}) = output(Bool, x[])

"""

```julia
pgarray(T, x) -> A
```

if `x` is an array, yields an array `A` with the same values as `x`, element
type `T` and suitable for PGPlot routines; if `x` is a scalar, yields a
reference (`Ref{T}`) with the value of `x`.  `T` should be `PGFloat` or
`PGInt`.

If `x` is already a dense array with element type `T`, `x` is returned.

"""
pgarray(::Type{T}, x::PGIntegers) where {T<:Integer} = Ref{T}(x)
pgarray(::Type{T}, x::Real) where {T<:AbstractFloat} = Ref{T}(x)
pgarray(::Type{T}, x::PGArray{T,N}) where {T<:Integer,N} = x
pgarray(::Type{T}, x::PGArray{T,N}) where {T<:AbstractFloat,N} = x
pgarray(::Type{T}, x::AbstractArray{<:PGIntegers,N}) where {T<:Integer,N} =
    convert(Array{T,N}, x)
pgarray(::Type{T}, x::AbstractArray{<:Real,N}) where {T<:AbstractFloat,N} =
    convert(Array{T,N}, x)

"""

`pgintarr(x)` is just an alias for `pgarray(PGInt,x)`.

"""
pgintarr(x) = pgarray(PGInt, x)

"""

`pgfltarr(x)` is just an alias for `pgarray(PGFloat,x)`.

"""
pgfltarr(x) = pgarray(PGFloat, x)

# Default coordinate transform.
const DEFAULT_XFORM = PGFloat[0, 1, 0,
                              0, 0, 1]
default_xform() = DEFAULT_XFORM

# Workspace arrays.
const XFORM_WORKSPACE = Vector{Vector{PGFloat}}(undef, 0)

function __init__()
    n = Threads.nthreads()
    resize!(XFORM_WORKSPACE, n)
    for i in 1:n
        XFORM_WORKSPACE[i] = zeros(PGFloat, 6)
    end
end

"""

```julia
get_own_xform()
```

yields a workspace array owned by the calling thread and suitable to store the
coordinate transform for some PGPlot routines.

```julia
get_own_xform(tr)
```

checks that `tr` is a valid coordinate transform and copy it in a workspace
array onwed by the caller and which is returned.  The returned array can be
modified by the caller.

"""
get_own_xform() = XFORM_WORKSPACE[Threads.threadid()]

function get_own_xform(tr::AbstractArray{<:Real, N}) where {N}
    (N == 1 && length(tr) == 6) || (N == 2 && size(tr) == (2,3)) ||
        throw_bad_xform_size()
    wtr = get_own_xform()
    k = 0
    @inbounds for i in eachindex(tr)
        k += 1
        wtr[k] = tr[i]
    end
    return wtr
end

@noinline throw_bad_xform_size() =
    throw(DimensionMismatch("coordinate transform must be a 3-element vector or a 2×3 matrix"))

"""

```julia
@check_vectors n x y z ...
```

checks that vectors `x`, `y`, `z`, etc. have the same length and store this
length in local variable `n`.  A `DimensionMismatch` exception is thrown with
appropriate message if the lengths are not the same.  There must be at least
2 vectors.

"""
macro check_vectors(args::Symbol...)
    esc(_check_vectors(args))
end

function _check_vectors(args::NTuple{N,Symbol}) where {N}
    N ≥ 3 || error("expecting at least 2 vectors to check")
    io = IOBuffer()
    print(io, "vectors")
    for i in 2:N
        print(io, (i == 2 ? " `" : i == N ? " and `" : ", `"), args[i], "`")
    end
    print(io, " must have the same length")
    mesg = String(take!(io))
    #println(mesg)
    expr = Expr(:comparison, :(($(args[1]) = length($(args[2])))))
    for i in 3:N
        push!(expr.args, :(==), :(length($(args[i]))))
    end
    :(($expr) || throw(DimensionMismatch($mesg)))
end

"""

```julia
submatrix(A, i1, i2, j1, j2, tr) -> A′, idim′, jdim′, i1′, i2′, j1′, j2′, tr′
```

checks that intervals `i1:i2` and `j1:j2` and coordinate transform `tr` are
valid for matrix `A` and returns 8-tuple suitable for PGPlot routines that deal
with sub-matrices.

Arguments `i1, `i2, `j1` and `j2` can be omitted to consider the full matrix
`A` or replaced by intervals `I = i1:i2` and `J = j1:j2`.

"""
function submatrix(A::AbstractMatrix,
                   tr::AbstractArray{<:Real})
    I, J = axes(A)
    return submatrix(A, I, J, tr)
end

function submatrix(A::AbstractMatrix,
                   I::AbstractUnitRange{<:Integer},
                   J::AbstractUnitRange{<:Integer},
                   tr::AbstractArray{<:Real})
    return submatrix(A, first(I), last(I), first(J), last(J), tr)
end

# Array is already in a form directly usable by PGPlot.
function submatrix(A::PGMatrix{PGFloat},
                   i1::Integer, i2::Integer,
                   j1::Integer, j2::Integer,
                   tr::AbstractArray{<:Real})
    wtr = get_own_xform(tr)
    idim, jdim = size(A)
    1 ≤ i1 ≤ i2 ≤ idim ||
        throw(ArgumentError("out of range indices `i1:i2`"))
    1 ≤ j1 ≤ j2 ≤ jdim ||
        throw(ArgumentError("out of range indices `j1:j2`"))
    return (A, PGInt(idim), PGInt(jdim),
            PGInt(i1), PGInt(i2), PGInt(j1), PGInt(j2), wtr)
end

# Array must be converted in a form usable by PGPlot.  Only the useful part is
# converted. FIXME:  add a ±1 pixel margin?
function submatrix(A::AbstractMatrix{<:Real},
                   i1::Integer, i2::Integer,
                   j1::Integer, j2::Integer,
                   tr::AbstractArray{<:Real})
    wtr = get_own_xform(tr)
    I, J = axes(A)
    imin, imax = axis_limits(I)
    jmin, jmax = axis_limits(J)
    imin ≤ i1 ≤ i2 ≤ imax ||
        throw(ArgumentError("out of range indices `i1:i2`"))
    jmin ≤ j1 ≤ j2 ≤ jmax ||
        throw(ArgumentError("out of range indices `j1:j2`"))
    ioff = i1 - 1
    joff = j1 - 1
    if ioff != 0 || joff != 0
        # Account for offsets in coordinate transform.
        di = PGFloat(ioff)
        dj = PGFloat(joff)
        wtr[1] += wtr[2]*di + xtr[3]*dj
        wtr[4] += wtr[5]*di + xtr[6]*dj
    end
    idim = i2 - ioff
    jdim = j2 - joff
    B = Array{PGFloat,2}(undef, idim, jdim)
    k = 0
    @inbounds for j in Int(j1):Int(j2), i in Int(i1):Int(i2)
        k += 1
        B[k] = A[i,j]
    end
    _one, _idim, _jdim = one(PGInt), PGInt(idim), PGInt(jdim)
    return (B, _idim, _jdim, _one, _idim, _one, _jdim, wtr)
end

# pgconb, pgconf, pgconl, pgcons, pgcont,
# pggray, pgimag,
# pgpixl

function submatrix(::Type{T},
                   A::AbstractMatrix) where {T}
    idim, jdim = size(A)
    _one, _idim, _jdim = one(PGInt), PGInt(idim), PGInt(jdim)
    return (pgarray(T, A), _idim, _jdim, _one, _idim, _one, _jdim)
end

function submatrix(::Type{T},
                   A::AbstractMatrix,
                   I::AbstractUnitRange{<:Integer},
                   J::AbstractUnitRange{<:Integer}) where {T}
    return submatrix(T, A, first(I), last(I), first(J), last(J))
end

# No conversion is needed, just check the arguments.
function submatrix(::Type{T},
                   A::PGMatrix{T},
                   i1::Integer, i2::Integer,
                   j1::Integer, j2::Integer,
                   tr::AbstractArray{<:Real}) where {T}
    idim, jdim = size(A)
    1 ≤ i1 ≤ i2 ≤ idim ||
        throw(ArgumentError("out of range indices `i1:i2`"))
    1 ≤ j1 ≤ j2 ≤ jdim ||
        throw(ArgumentError("out of range indices `j1:j2`"))
    return (A, PGInt(idim), PGInt(jdim),
            PGInt(i1), PGInt(i2), PGInt(j1), PGInt(j2))
end

# A conversion is needed, the smallest possible part is converted.
function submatrix(::Type{T},
                   A::AbstractMatrix{<:Real},
                   i1::Integer, i2::Integer,
                   j1::Integer, j2::Integer) where {T}
    I, J = axes(A)
    imin, imax = axis_limits(I)
    jmin, jmax = axis_limits(J)
    imin ≤ i1 ≤ i2 ≤ imax ||
        throw(ArgumentError("out of range indices `i1:i2`"))
    jmin ≤ j1 ≤ j2 ≤ jmax ||
        throw(ArgumentError("out of range indices `j1:j2`"))
    idim = i2 + 1 - i1
    jdim = j2 + 1 - j1
    B = Array{T,2}(undef, idim, jdim)
    k = 0
    @inbounds for j in Int(j1):Int(j2)
        for i in Int(i1):Int(i2)
            k += 1
            B[k] = A[i,j]
        end
    end
    _one, _idim, _jdim = one(PGInt), PGInt(idim), PGInt(jdim)
    return (B, _idim, _jdim, _one, _idim, _one, _jdim)
end

autorange(t::AbstractVector) = autorange(extrema(t)...)

autorange(tmin::Real, tmax::Real) = autorange(PGFloat(tmin), PGFloat(tmax))

function autorange(tmin::PGFloat, tmax::PGFloat)
    @assert tmin ≤ tmax
    dt = PGFloat(0.05)*(tmax - tmin)
    if dt ≤ zero(dt)
        dt = one(dt)
    end
    return (tmin - dt, tmax + dt)
end

# Helper to query a single parameter.
function query(::Type{T}, func::Function) where {T<:Union{PGInt,PGFloat}}
    val = Ref{T}()
    func(val)
    return output(val)
end

function query(::Type{Bool}, func::Function)
    val = Ref{PGInt}()
    func(val)
    return output(Bool, val)
end

# Helper to query two parameters.
function query(::Type{T1}, ::Type{T2},
               func::Function) where {T1<:Union{PGInt,PGFloat},
                                      T2<:Union{PGInt,PGFloat}}
    val1 = Ref{T1}()
    val2 = Ref{T2}()
    func(val1, val2)
    return output(val1, val2)
end

# Helper to query three parameters.
function query(::Type{T1}, ::Type{T2}, ::Type{T3},
               func::Function) where {T1<:Union{PGInt,PGFloat},
                                      T2<:Union{PGInt,PGFloat},
                                      T3<:Union{PGInt,PGFloat}}
    val1 = Ref{T1}()
    val2 = Ref{T2}()
    val3 = Ref{T3}()
    func(val1, val2, val3)
    return output(val1, val2, val3)
end

# Helper to query four parameters given an argument.
function query(::Type{T1}, ::Type{T2}, ::Type{T3}, ::Type{T4},
               func::Function) where {T1<:Union{PGInt,PGFloat},
                                           T2<:Union{PGInt,PGFloat},
                                           T3<:Union{PGInt,PGFloat},
                                           T4<:Union{PGInt,PGFloat}}
    val1 = Ref{T1}()
    val2 = Ref{T2}()
    val3 = Ref{T3}()
    val4 = Ref{T4}()
    func(val1, val2, val3, val4)
    return output(val1, val2, val3, val4)
end

# Helper to query four parameters given an argument.
function query(::Type{T1}, ::Type{T2}, ::Type{T3}, ::Type{T4},
               func::Function, arg) where {T1<:Union{PGInt,PGFloat},
                                           T2<:Union{PGInt,PGFloat},
                                           T3<:Union{PGInt,PGFloat},
                                           T4<:Union{PGInt,PGFloat}}
    val1 = Ref{T1}()
    val2 = Ref{T2}()
    val3 = Ref{T3}()
    val4 = Ref{T4}()
    func(arg, val1, val2, val3, val4)
    return output(val1, val2, val3, val4)
end

function to_string(buf::PGVector{UInt8}, len::Integer)
    @assert 0 ≤ len < length(buf)
    buf[len+1] = 0
    unsafe_string(pointer(buf), len)
end

"""

```julia
axis_limits(I) = (i0,i1)
```

yields the limits `i0` and `i1` of index range `I` as a 2-tuple of `Int`'s and
such that `i0:i1` represents the same indices as `I` (although not in the same
order if `step(I) < 0`).  If `step(I)` is not equal to ±1, an `ArgumentError`
exception is thrown.

"""
axis_limits(I::AbstractUnitRange{<:Integer}) =
    (Int(first(I)), Int(last(I)))
axis_limits(I::AbstractRange{<:Integer}) =
    ((i0, i1, s) = (Int(first(I)), Int(last(I)), step(I));
     (s == +1 ? (i0,i1) :
      s == -1 ? (i1,i0) : throw_invalid_range_step()))

@noinline throw_invalid_range_step() =
    throw(ArgumentError("expecting a range with a step equal to ±1"))

"""

```julia
pgunits(units)
```

yields the integer code corresponding to the symbolic `units` which can be one
of `:ndc` for *normalized device coordinates*, `:in` for *inches*, `:mm` for
*millimeters* or `:pix` for *device units* (usually *pixels*).

""" pgunits

const NDC_UNITS   = PGInt(0)
const INCHES      = PGInt(1)
const MILLIMETERS = PGInt(2)
const PIXELS      = PGInt(3)
const WORLD       = PGInt(4)

pgunits(val::Symbol) =
    (val == :ndc ? NDC_UNITS   :
     val == :in  ? INCHES      :
     val == :mm  ? MILLIMETERS :
     val == :pix ? PIXELS      : throw_bad_units(val))

pgunits(val::Integer) =
    (0 ≤ val ≤ 3 ? PGInt(val) : throw_bad_units(val))

@noinline throw_bad_units(val) =
    throw(ArgumentError(string("unknown units: ", val)))

#------------------------------------------------------------------------------
# PGPlot routines.

"""

```julia
pgarro(x1, y1, x2, y2)
```

Draw an arrow from the point with world-coordinates `(x1,y1)` to
`(x2,y2)`. The size of the arrowhead at `(x2,y2)` is determined by the
current character size set by routine [`pgsch`](@ref). The default size is
1/40th of the smaller of the width or height of the view surface.  The
appearance of the arrowhead (shape and solid or open) is controlled by
routine [`pgsah`][@ref).

"""
pgarro(x1::Real, y1::Real, x2::Real, y2::Real) =
    ccall((:cpgarro, pgplotlib), Cvoid,
          (PGFloat, PGFloat, PGFloat, PGFloat), x1, y1, x2, y2)

"""

```julia
pgask(flag)
```

Change the *prompt state* of PGPlot. If the prompt state is *on*,
[`pgpage`](@ref) will type

> Type RETURN for next page:

and will wait for the user to type a carriage-return before starting a new
page.  The initial prompt state (after the device has been opened) is *on* for
interactive devices. Prompt state is always *off* for non-interactive devices.

"""
pgask(flag::Bool) = ccall((:cpgask, pgplotlib), Cvoid, (PGBool,), flag)

"""

```julia
pgaxis(opt, x1, y1, x2, y2, v1, v2, step, nsub, dmajl, dmajr, fmin, disp,
       orient)
```

Draw a labelled graph axis from world-coordinate position `(x1,y1)` to
`(x2,y2)`.

Normally, this routine draws a standard *linear* axis with equal
subdivisions.  The quantity described by the axis runs from `v1` to `v2`;
this may be, but need not be, the same as `x` or `y`.

If the `'L'` option is specified, the routine draws a *logarithmic* axis.
In this case, the quantity described by the axis runs from `10^v1` to
`10^v2`. A logarithmic axis always has major, labeled, tick marks spaced by
one or more decades. If the major tick marks are spaced by one decade (as
specified by the `step` argument), then minor tick marks are placed at 2,
3, .., 9 times each power of 10; otherwise minor tick marks are spaced by
one decade. If the axis spans less than two decades, numeric labels are
placed at 1, 2, and 5 times each power of ten.

If the axis spans less than one decade, or if it spans many decades, it is
preferable to use a linear axis labeled with the logarithm of the quantity
of interest.

Arguments:

- `opt`: a string containing single-letter codes for various options. The
  options currently recognized are:

  - `L`: draw a logarithmic axis;

  - `N`: write numeric labels;

  - `1`: force decimal labelling, instead of automatic choice (see
    [`pgnumb`](@ref));

  - `2`: force exponential labelling, instead of automatic.

- `(x1,y1)`: world coordinates of one endpoint of the axis.

- `(x2,y2)`: world coordinates of the other endpoint of the axis.

- `(v1,v2)`: axis value at first abd second endpoints.

- `step`: major tick marks are drawn at axis value 0.0 plus or minus
  integer multiples of `step`. If `step = 0`, a value is chosen
  automatically.

- `nsub`: minor tick marks are drawn to divide the major divisions into
  `nsub` equal subdivisions (ignored if `step = 0`). If `nsub ≤ 1`, no minor
  tick marks are drawn. `nsub` is ignored for a logarithmic axis.

- `dmajl`: length of major tick marks drawn to left of axis (as seen
  looking from first endpoint to second), in units of the character height.

- `dmajr`: length of major tick marks drawn to right of axis, in units of
  the character height.

- `fmin`: length of minor tick marks, as fraction of major.

- `disp`: displacement of baseline of tick labels to right of axis, in
  units of the character height.

- `orient`: orientation of label text, in degrees; angle between baseline
  of text and direction of axis (0-360°).

"""
function pgaxis(opt::AbstractString, x1::Real, y1::Real, x2::Real, y2::Real,
                v1::Real, v2::Real, step::Real, nsub::Integer, dmajl::Real,
                dmajr::Real, fmin::Real, disp::Real, orient::Real)
    ccall((:cpgaxis, pgplotlib), Cvoid,
          (PGString, PGFloat, PGFloat, PGFloat, PGFloat, PGFloat, PGFloat,
           PGFloat, PGInt, PGFloat, PGFloat, PGFloat, PGFloat, PGFloat),
          opt, x1, y1, x2, y2, v1, v2, step, nsub, dmajl, dmajr, fmin,
          disp, orient)
end

"""

```julia
pgband(mode, posn, x0=0, y0=0, x1=x0, y1=y0) -> (x, y, c)
```

yields the cursor position and a character typed by the user.  The position
is returned in world coordinates.  `pgband` positions the cursor at the
position specified (if `posn=1`), allows the user to move the cursor using
the mouse or arrow keys or whatever is available on the device. When he has
positioned the cursor, the user types a single character on the keyboard;
[`pgband`](@ref) then returns this character and the new cursor position
(in world coordinates).

Some interactive devices offer a selection of cursor types, implemented as
thin lines that move with the cursor, but without erasing underlying
graphics. Of these types, some extend between a stationary anchor-point at
`(x0,y0)`, and the position of the cursor, while others simply follow
the cursor without changing shape or size. The cursor type is specified
with one of the following `mode` values. Cursor types that are not
supported by a given device, are treated as `mode=0`.

- If `mode=0`, the anchor point is ignored and the routine behaves
  like [`pgcurs`](@ref).

- If `mode=1`, a straight line is drawn joining the anchor point
  and the cursor position.

- If `mode=2`, a hollow rectangle is extended as the cursor is moved, with
  one vertex at the anchor point and the opposite vertex at the current
  cursor position; the edges of the rectangle are horizontal and vertical.

- If `mode=3`, two horizontal lines are extended across the width of the
  display, one drawn through the anchor point and the other through the
  moving cursor position. This could be used to select a Y-axis range when
  one end of the range is known.

- If `mode=4`, two vertical lines are extended over the height of the
  display, one drawn through the anchor point and the other through the
  moving cursor position. This could be used to select an X-axis range when
  one end of the range is known.

- If `mode=5`, a horizontal line is extended through the cursor position
  over the width of the display. This could be used to select an X-axis
  value such as the start of an X-axis range. The anchor point is ignored.

- If `mode=6`, a vertical line is extended through the cursor position over
  the height of the display. This could be used to select a Y-axis value
  such as the start of a Y-axis range. The anchor point is ignored.

- If `mode=7`, a cross-hair, centered on the cursor, is extended over the
  width and height of the display. The anchor point is ignored.

Arguments:

- `mode`: display mode (0, 1, ..7: see above).

- `posn`: if true, `pgband` attempts to place the cursor at point
  `(x1,y1)`; otherwise, it leaves the cursor at its current
  position. (On some devices this request may be ignored.)

- `(x0,y0)`: the world coordinates of the anchor point.

- `(x1,y1)`: the initial world coordinates of the cursor if `posn` is true.

Outputs:

- `(x,y)`: the world coordinates of the cursor when a character was
  typed. The returned cursor coordinates `(x,y)` may be different from
  `(x1,y1)` even if the device has no cursor or if the user does not
  move the cursor.  Under these circumstances, the position returned in
  `(x,y)` is that of the pixel nearest to the requested position.

- `c`: the character typed by the user; if the device has no cursor or if
  some other error occurs, the value `Char(0)` (ASCII NUL character) is
  returned.  For an X-window device, clicking the 1st mouse button is
  reported as an `A`, clicking the 2nd mouse button is reported as a `D`
  and clicking other mouse buttons is reported as an `X`.

"""
function pgband(mode::Integer, posn::Bool,
                x0::Real = PGFloat(0), y0::Real = PGFloat(0),
                x1::Real = x0, y1::Real = y0)
    @assert isfinite(x0) && isfinite(y0) && isfinite(x1) && isfinite(y1)
    x = Ref{PGFloat}(x1)
    y = Ref{PGFloat}(y1)
    c = Ref{PGChar}(0)
    code = ccall((:cpgband, pgplotlib), PGInt,
                 (PGInt, PGInt, PGFloat, PGFloat, Ref{PGFloat},
                  Ref{PGFloat}, Ref{PGChar}),
                 mode, posn, x0, y0, x, y, c)
    # PGBAND returns 1 if the call was successful; 0 if the device has no
    # cursor or some other error occurs.
    code == 0 && error("the device has no cursor or some other error occurs")
    return output(x, y, c)
end

"""

```julia
pgbbuf()
```

Begin saving graphical output commands in an internal buffer; the commands are
held until a matching [`pgebuf`](@ref) call (or until the buffer is emptied by
[`pgupdt`](@ref)). This can greatly improve the efficiency of PGPlot.  `pgbbuf`
increments an internal counter, while [`pgebuf`](@ref) decrements this counter
and flushes the buffer to the output device when the counter drops to zero.
`pgbbuf` and [`pgebuf`](@ref) calls should always be paired.

"""
pgbbuf() = ccall((:cpgbbuf, pgplotlib), Cvoid, (),)

"""

```julia
pgbeg([unit=0,] file, nxsub, nysub) -> id
```

opens a graphical device or file and prepares it for subsequent plotting.  A
device must be opened with [`pgbeg`](@ref) or [`pgopen`](@ref) before any other
calls to PGPlot subroutines for the device.

If any device is already open for PGPlot output, it is closed before the new
device is opened.

Note: new programs should use [`pgopen`](@ref) rather than [`pgbeg`](@ref).
[`pgbeg`](@ref) is retained for compatibility with existing programs.  Unlike
[`pgopen`](@ref), [`pgbeg`](@ref) closes any graphics devices that are already
open, so it cannot be used to open devices to be used in parallel.


Returns an integer status.  A value of 1 indicates successful completion, any
other value indicates an error.  In the event of error a message is written on
the standard error output.  To test the return value, call [`pgbeg`](@ref) as a
function, *e.g.* `ier = pgbeg(...)`.

Arguments:

- `unit` : this argument is ignored by `pgbeg` (use zero).

- `file` : the "device specification" for the plot device.  (For explanation,
  see description of `pgopen`.)

- `nxsub` : the number of subdivisions of the view surface in X (`nxsub > 0` or
  `nxsub < 0`).

- `nysub` : the number of subdivisions of the view surface in Y (`nysub > 0`).

PGPlot puts `nxsub × nysub` graphs on each plot page or screen; when the view
surface is sub-divided in this way, [`pgpage`](@ref) moves to the next panel,
not the next physical page.  If `nxsub > 0`, PGPlot uses the panels in row
order; if `nxsub < 0`, PGPlot uses them in column order.

"""
pgbeg(unit::Integer, file::AbstractString, nxsub::Integer, nysub::Integer) =
    ccall((:cpgbeg, pgplotlib), PGInt,
          (PGInt, PGString, PGInt, PGInt),
          unit, file, nxsub, nysub)

pgbeg(file::AbstractString, nxsub::Integer, nysub::Integer) =
    pgbeg(0, file, nxsub, nysub)

"""

```julia
pgbin(x, y, center=true)
```

plots a histogram of values with abscissae `x` and ordinates `y`. Bin width is
spacing between `x` values.  If `center` is `true`, the `x` values denote the
center of the bins; otherwise, the `x` values denote the left edge of the bins.

"""
function pgbin(x::AbstractVector{<:Real},
               y::AbstractVector{<:Real},
               cen::Bool = true)
    @check_vectors n x y
    _pgbin(n, pgfltarr(x), pgfltarr(y), cen)
end

_pgbin(n::Int, x::PGVector{PGFloat}, y::PGVector{PGFloat}, cen::Bool) =
    ccall((:cpgbin, pgplotlib), Cvoid,
          (PGInt, Ptr{PGFloat}, Ptr{PGFloat}, PGBool), n, x, y, cen)

"""

```julia
pgbox(xopt, xtick, nxsub, yopt, ytick, nysub)
```

annotates the viewport with frame, axes, numeric labels, etc.  `pgbox` is
called by on the user's behalf by [`pgenv`](@ref), but may also be called
explicitly.

Arguments:

- `xopt`, `yopt`: strings of options for X (horizontal) and Y (vertical) axes
  of plot. Options are single letters, and may be in any order (see below).

- `xtick`, `ytick`: world coordinate intervals between major tick marks on X
  and Y axes. If 0, the interval is chosen by `pgbox`, so that there will be at
  least 3 major tick marks along the axis.

- `nxsub`, `nysub`: the numbers of subintervals to divide the major coordinate
  interval into.  If the tick interval or number of subintervals is zero for an
  axis, the number of subintervals is chosen by `pgbox`.

Options (for arguments `xopt` and `yopt`):

- `A`: draw Axis (X axis is horizontal line `y = 0`, Y axis is vertical line `X
  = 0`).

- `B` : draw bottom (X) or left (Y) edge of frame.

- `C`: draw top (X) or right (Y) edge of frame.

- `G`: draw Grid of vertical (X) or horizontal (Y) lines.

- `I`: Invert the tick marks; ie draw them outside the viewport instead of
  inside.

- `L`: label axis Logarithmically (see below).

- `N`: write Numeric labels in the conventional location below the viewport (X)
  or to the left of the viewport (Y).

- `P`: extend ("Project") major tick marks outside the box (ignored if option I
  is specified).

- `M`: write numeric labels in the unconventional location above the viewport
  (X) or to the right of the viewport (Y).

- `T`: draw major Tick marks at the major coordinate interval.

- `S`: draw minor tick marks (Subticks).

- `V`: orient numeric labels Vertically. This is only applicable to Y.  The
  default is to write Y-labels parallel to the axis.

- `1`: force decimal labelling, instead of automatic choice (see
  [`pgnumb`](@ref)).

- `2`: force exponential labelling, instead of automatic.

To get a complete frame, specify `BC` in both `xopt` and `yopt`.  Tick marks,
if requested, are drawn on the axes or frame or both, depending which are
requested. If none of `ABC` is specified, tick marks will not be drawn. When
[`pgenv`](@ref) calls `pgbox`, it sets both `xopt` and `yopt` according to the
value of its parameter `axis`: -1: 'BC', 0: 'BCNST', 1: 'ABCNST', 2: 'ABCGNST'.

For a logarithmic axis, the major tick interval is always 1.0. The numeric
label is `10^x` where `x` is the world coordinate at the tick mark. If subticks
are requested, 8 subticks are drawn between each major tick at equal
logarithmic intervals.

To label an axis with time (days, hours, minutes, seconds) or angle (degrees,
arcmin, arcsec), use routine [`pgtbox`](@ref).

"""
function pgbox(xopt::AbstractString, xtick::Real, nxsub::Integer,
               yopt::AbstractString, ytick::Real, nysub::Integer)
    ccall((:cpgbox, pgplotlib), Cvoid,
          (PGString, PGFloat, PGInt, PGString, PGFloat, PGInt),
          xopt, xtick, nxsub, yopt, ytick, nysub)
end

"""

```julia
pgcirc(x, y, r)
```

draws a circle of radius `r` centered at world coordinates '(x,y)` . The action
of this routine depends on the setting of the Fill-Area Style attribute. If
Fill-Area Style is *solid* (the default), the interior of the circle is
solid-filled using the current Color Index. If Fill-Area Style is *hollow*, the
outline of the circle is drawn using the current line attributes (color index,
line-style, and line-width).

"""
pgcirc(x::Real, y::Real, r::Real) =
    ccall((:cpgcirc, pgplotlib), Cvoid, (PGFloat, PGFloat, PGFloat), x, y, r)

"""

```julia
pgclos()
```

Close the currently selected graphics device. After the device has been closed,
either another open device must be selected with [`pgslct`](@ref) or another
device must be opened with [`pgopen`](@ref) before any further plotting can be
done. If the call to `pgclos` is omitted, some or all of the plot may be lost.

"""
pgclos() = ccall((:cpgclos, pgplotlib), Cvoid, (), )

"""

```julia
pgconb(A, [i1, i2, j1, j2,] c, tr=[0,1,0,0,0,1], blank=NaN)
```

Draw a contour map of sub-array `A[i1:i2,j1:j2]` at level(s) specified by
`c`. This routine is the same as [`pgcons`](@ref), except that array elements
that have the *magic value* defined by argument `blank` are ignored, making
gaps in the contour map. The routine may be useful for data measured on most
but not all of the points of a grid.

Sub-array indices `i1`, `i2`, `j1` and `j2` may be specified as ranges
`I = i1:i2`, `J = j1:j2` or may be unspecified to plot the whole array `A`.

Argument `tr` specifies a transformation matrix used to calculate the world
coordinates of the center of the *cell* that represents each array element. The
world coordinates of the center of the cell corresponding to array element
`A[i,j]` are given by:

    x = tr(1) + tr(2)*i + tr(3)*j
    y = tr(4) + tr(5)*i + tr(6)*j

Hence `tr` must be 6-element vector or a 2×3 matrix.

"""
function pgconb(A::AbstractMatrix{<:Real},
                c::Union{Real,AbstractVector{<:Real}},
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM,
                blank::Real = PGFloat(NaN))
    _pgconb(submatrix(A, tr)..., pgfltarr(c), length(c), blank)
end

function pgconb(A::AbstractMatrix{<:Real},
                I::AbstractUnitRange{<:Integer},
                J::AbstractUnitRange{<:Integer},
                c::Union{Real,AbstractVector{<:Real}},
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM,
                blank::Real = PGFloat(NaN))
    _pgconb(submatrix(A, I, J, tr)..., pgfltarr(c), length(c), blank)
end

function pgconb(A::AbstractMatrix{<:Real},
                i1::Integer, i2::Integer,
                j1::Integer, j2::Integer,
                c::Union{Real,AbstractVector{<:Real}},
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM,
                blank::Real = PGFloat(NaN))
    _pgconb(submatrix(A, i1, i2, j1, j2, tr)..., pgfltarr(c), length(c), blank)
end

# Argument `tr` moved to directly use the output of `submatrix`.
_pgconb(a, idim, jdim, i1, i2, j1, j2, tr, c, nc, blank) =
    ccall((:cpgconb, pgplotlib), Cvoid,
          (Ptr{PGFloat}, PGInt, PGInt, PGInt, PGInt, PGInt, PGInt,
           Ptr{PGFloat}, PGInt, Ptr{PGFloat}, PGFloat),
          a, idim, jdim, i1, i2, j1, j2, c, nc, tr, blank)

"""

```julia
pgconf(A, [i1, i2, j1, j2,] c1, c2, tr=[0,1,0,0,0,1])
```
Shade the region between two contour levels of a function defined on
the nodes of a rectangular grid. The routine uses the current fill
attributes, hatching style (if appropriate), and color index.

If you want to both shade between contours and draw the contour
lines, call this routine first (once for each pair of levels) and
then CALL PGCONT (or PGCONS) to draw the contour lines on top of the
shading.

Note 1: This routine is not very efficient: it generates a polygon
fill command for each cell of the mesh that intersects the desired
area, rather than consolidating adjacent cells into a single polygon.

Note 2: If both contours intersect all four edges of a particular
mesh cell, the program behaves badly and may consider some parts
of the cell to lie in more than one contour range.

Note 3: If a contour crosses all four edges of a cell, this
routine may not generate the same contours as PGCONT or PGCONS
(these two routines may not agree either). Such cases are always
ambiguous and the routines use different approaches to resolving
the ambiguity.

"""
function pgconf(A::AbstractMatrix{<:Real},
                c1::Real, c2::Real,
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgconf(submatrix(A, tr)..., c1, c2)
end

function pgconf(A::AbstractMatrix{<:Real},
                I::AbstractUnitRange{<:Integer},
                J::AbstractUnitRange{<:Integer},
                c1::Real, c2::Real,
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgconf(submatrix(A, I, J, tr)..., c1, c2)
end

function pgconf(A::AbstractMatrix{<:Real},
                i1::Integer, i2::Integer,
                j1::Integer, j2::Integer,
                c::Union{Real,AbstractVector{<:Real}},
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM,
                blank::Real = PGFloat(NaN))
    _pgconf(submatrix(A, i1, i2, j1, j2, tr)..., c1, c2)
end

# Argument `tr` moved to directly use the output of `submatrix`.
_pgconf(a, idim, jdim, i1, i2, j1, j2, tr, c1, c2) =
    ccall((:cpgconf, pgplotlib), Cvoid,
          (Ptr{PGFloat}, PGInt, PGInt, PGInt, PGInt, PGInt, PGInt,
           PGFloat, PGFloat, Ptr{PGFloat}),
          a, idim, jdim, i1, i2, j1, j2, c1, c2, tr)

"""

```julia
pgconl(A, [i1, i2, j1, j2,] c, lebel, intval, minint, tr=[0,1,0,0,0,1])
```

Label a contour map drawn with routine PGCONT. Routine PGCONT should
be called first to draw the contour lines, then this routine should be
called to add the labels. Labels are written at intervals along the
contour lines, centered on the contour lines with lettering aligned
in the up-hill direction. Labels are opaque, so a part of the under-
lying contour line is obscured by the label. Labels use the current
attributes (character height, line width, color index, character
font).

The first 9 arguments are the same as those supplied to PGCONT, and
should normally be identical to those used with PGCONT. Note that
only one contour level can be specified; tolabel more contours, call
PGCONL for each level.

The Label is supplied as a character string in argument LABEL.

The spacing of labels along the contour is specified by parameters
INTVAL and MININT. The routine follows the contour through the
array, counting the number of cells that the contour crosses. The
first label will be written in the MININT'th cell, and additional
labels will be written every INTVAL cells thereafter. A contour
that crosses less than MININT cells will not be labelled. Some
experimentation may be needed to get satisfactory results; a good
place to start is INTVAL=20, MININT=10.

"""
function pgconl(A::AbstractMatrix{<:Real},
                c::Real,
                label::AbstractString,
                intval::Integer,
                minint::Integer,
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgconl(submatrix(A, tr)..., c, label, intval, minint)
end

function pgconl(A::AbstractMatrix{<:Real},
                I::AbstractUnitRange{<:Integer},
                J::AbstractUnitRange{<:Integer},
                c::Real,
                label::AbstractString,
                intval::Integer,
                minint::Integer,
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgconl(submatrix(A, I, J, tr)..., c, label, intval, minint)
end

function pgconl(A::AbstractMatrix{<:Real},
                i1::Integer, i2::Integer,
                j1::Integer, j2::Integer,
                c::Real,
                label::AbstractString,
                intval::Integer,
                minint::Integer,
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgconl(submatrix(A, i1, i2, j1, j2, tr)..., c, label, intval, minint)
end

# Argument `tr` moved to directly use the output of `submatrix`.
_pgconl(a, idim, jdim, i1, i2, j1, j2, tr, c, label, intval, minint) =
    ccall((:cpgconl, pgplotlib), Cvoid,
          (Ptr{PGFloat}, PGInt, PGInt, PGInt, PGInt, PGInt, PGInt,
           PGFloat, Ptr{PGFloat}, PGString, PGInt, PGInt),
          a, idim, jdim, i1, i2, j1, j2, c, tr, label, intval, minint)


# FIXME: fix doc. below.
"""
Contour map of a 2D data array (contour-following)
==================================================

`pgcont(a, c)` draws a contour map of the 2-D array `a` at level(s) `c`.
The map is truncated if necessary at the boundaries of the viewport.  Each
contour line is drawn with the current line attributes (color index, style,
and width); except that the line style is set to 1 (solid) for positive
contours or 2 (dashed) for negative contours.

Keyword `samestyle` can be set true to force all contours to have the same
line style.

Keyword `transform` can be used to specify transformation between the
`(i,j)` grid of the array and the world coordinates.  The argument is a
6-element array `tr` and the world coordinates of the array point `a[i,j]`
are given by:

    x = tr[1] + tr[2]*i + tr[3]*j
    y = tr[4] + tr[5]*i + tr[6]*j

Usually `tr[3]` and `tr[5]` are zero - unless the coordinate transformation
involves a rotation or shear.  The default coordinate transform is:

    tr = [0, 1, 0, 0, 0, 1]

The following call:

    pgcont(a, i1:i2, j1:j2, c)

can be used to restrict the contours to the region corresponding to
`a[i1:i2,j1:j2]` without the needs to fix the coordinate transform.


  Each contour line
is drawn with the current line attributes (color index, style, and
width); except that if argument NC is positive (see below), the line
style is set by [`pgcont`](@ref) to 1 (solid) for positive contours or 2
(dashed) for negative contours.

A lower level method is provided by:

    pgcont(a, i1, i2, j1, j2, c, nc, tr)

which draws a contour map of the 2-D sub-array `a[i1:i2,j1:j2]` at level(s)
`c`.  `abs(nc)` is the number of contour levels to draw (must be at most
`length(c)`).  If `nc` is positive, it is the number of contour levels, and
the line-style is chosen automatically as described above.  If `nc` is
negative, it is minus the number of contour levels, and the current setting
of line-style is used for all the contours.

""" pgcont

@doc @doc(pgcont) pgcons

function pgcons(A::AbstractMatrix{<:Real},
                c::Union{Real,AbstractVector{<:Real}},
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgcons(submatrix(A, tr)..., pgfltarr(c), length(c))
end

function pgcons(A::AbstractMatrix{<:Real},
                I::AbstractUnitRange{<:Integer},
                J::AbstractUnitRange{<:Integer},
                c::Union{Real,AbstractVector{<:Real}},
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgcons(submatrix(A, I, J, tr)..., pgfltarr(c), length(c))
end

function pgcons(A::AbstractMatrix{<:Real},
                i1::Integer, i2::Integer,
                j1::Integer, j2::Integer,
                c::Union{Real,AbstractVector{<:Real}},
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgcons(submatrix(A, i1, i2, j1, j2, tr)..., pgfltarr(c), length(c))
end

# Argument `tr` moved to directly use the output of `submatrix`.
_pgcons(a, idim, jdim, i1, i2, j1, j2, tr, c, nc) =
    ccall((:cpgcons, pgplotlib), Cvoid,
          (Ptr{PGFloat}, PGInt, PGInt, PGInt, PGInt, PGInt, PGInt,
           Ptr{PGFloat}, PGInt, Ptr{PGFloat}),
          a, idim, jdim, i1, i2, j1, j2, c, nc, tr)

function pgcont(A::AbstractMatrix{<:Real},
                c::Union{Real,AbstractVector{<:Real}},
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgcont(submatrix(A, tr)..., pgfltarr(c), length(c))
end

function pgcont(A::AbstractMatrix{<:Real},
                I::AbstractUnitRange{<:Integer},
                J::AbstractUnitRange{<:Integer},
                c::Union{Real,AbstractVector{<:Real}},
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgcont(submatrix(A, I, J, tr)..., pgfltarr(c), length(c))
end

function pgcont(A::AbstractMatrix{<:Real},
                i1::Integer, i2::Integer,
                j1::Integer, j2::Integer,
                c::Union{Real,AbstractVector{<:Real}},
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgcont(submatrix(A, i1, i2, j1, j2, tr)..., pgfltarr(c), length(c))
end

# Argument `tr` moved to directly use the output of `submatrix`.
_pgcont(a, idim, jdim, i1, i2, j1, j2, tr, c, nc) =
    ccall((:cpgcont, pgplotlib), Cvoid,
          (Ptr{PGFloat}, PGInt, PGInt, PGInt, PGInt, PGInt, PGInt,
           Ptr{PGFloat}, PGInt, Ptr{PGFloat}),
          a, idim, jdim, i1, i2, j1, j2, c, nc, tr)

"""

```julia
pgctab(l, r, g, b, contra=1.0, bright=0.5)
```

Use the given color table to change the color representations of all color
indexes marked for use by [`pgimag`](@ref).  To change which color indexes are
thus marked, call [`pgscir`](@ref) before calling `pgctab` or [`pgimag`](@ref).
On devices that can change the color representations of previously plotted
graphics, `pgctab` will also change the colors of existing graphics that were
plotted with the marked color indexes.  This feature can then be combined with
[`pgband`](@ref) to interactively manipulate the displayed colors of data
previously plotted with [`pgimag`](@ref).

In the following, the term *color table* refers to the input L,R,G,B arrays,
whereas *color ramp* refers to the resulting ramp of colors that would be seen
with [`pgwedg`](@ref). Arguments:

- `L` is an array of normalized ramp-intensity levels corresponding to the RGB
   primary color intensities in `r`, `g` and `b` arrays.  Colors on the ramp
   are linearly interpolated from neighbouring levels.  Levels must be sorted
   in increasing order.  0.0 places a color at the beginning of the ramp.  1.0
   places a color at the end of the ramp.  Colors outside these limits are
   legal, but will not be visible if `contra=1.0` and `bright=0.5`.

- `r`, `g` and `b` are arrays of normalized red, green and blue intensities.

- `contra` is the contrast of the color ramp (normally 1.0).  Negative values
   reverse the direction of the ramp.

- `bright` is the brightness of the color ramp. This is normally 0.5, but can
  sensibly hold any value between 0.0 and 1.0. Values at or beyond the latter
  two extremes, saturate the color ramp with the colors of the respective end
  of the color table.

To reverse the sense of a color table, change the chosen contrast and
brightness to `-contra` and `1-bright`.

Limitations:

1. Some devices do not propagate color representation changes to previously
   drawn graphics.

2. Some devices ignore requests to change color representations.

3. The appearance of specific color representations on grey-scale devices is
   device-dependent.

"""
function pgctab(l::AbstractVector{<:Real},
                r::AbstractVector{<:Real},
                g::AbstractVector{<:Real},
                b::AbstractVector{<:Real},
                contra::Real = PGFloat(1.0),
                bright::Real = PGFloat(0.5))
    @check_vectors nc l r g b
    _pgctab(pgfltarr(l), pgfltarr(r), pgfltarr(g), pgfltarr(b), nc,
            contra, bright)
end

_pgctab(l, r, g, b, nc, contra, bright) =
    ccall((:cpgctab, pgplotlib), Cvoid,
          (Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat},
           PGInt, PGFloat, PGFloat),
          l, r, g, b, nc, contra, bright)

"""

```julia
pgcurs(x0,y0) -> (x, y, c)
```

yields the cursor position and a character typed by the user.  The position is
returned in world coordinates.  `pgcurs` positions the cursor at the position
specified `(x0,y0)`, allows the user to move the cursor using the joystick or
arrow keys or whatever is available on the device. When he has positioned the
cursor, the user types a single character on the keyboard; `pgcurs` then
returns this character and the new cursor position (in world coordinates).  See
[`pgband`](@ref) for a descrition of the output `(x,y,c)`.

"""
pgcurs(x0::Real, y0::Real) = pgcurs(PGFloat(x0), PGFloat(y0))

function pgcurs(x0::PGFloat, y0::PGFloat)
    @assert isfinite(x0) && isfinite(y0)
    x = Ref{PGFloat}(x0)
    y = Ref{PGFloat}(y0)
    c = Ref{PGChar}(0)
    code = ccall((:cpgcurs, pgplotlib), PGInt,
                 (Ref{PGFloat}, Ref{PGFloat}, Ref{PGChar}), x, y, c)
    # PGCURS returns 1 if the call was successful; 0 if the device has no
    # cursor or some other error occurs.
    code == 0 && error("the device has no cursor or some other error occurs")
    return output(x, y, c)
end

"""

```julia
pgdraw(x,y)
```

Draw a line from the current pen position to the point with world-coordinates
`(x,y)`. The line is clipped at the edge of the current window. The new pen
position is `(x,y)` in world coordinates.

"""
pgdraw(x::Real, y::Real) = pgdraw(PGFloat(x), PGFloat(y))
function pgdraw(x::PGFloat, y::PGFloat)
    # FIXME: @assert isfinite(x) && isfinite(y)
    ccall((:cpgdraw, pgplotlib), Cvoid, (PGFloat, PGFloat), x, y)
end

"""

```julia
pgebuf()
```

marks the end of a batch of graphical output begun with the last call of
[`pgbbuf`](@ref).  [`pgbbuf`](@ref) and `pgebuf` calls should always be
paired. Each call to [`pgbbuf`](@ref) increments a counter, while each call to
`pgebuf` decrements the counter. When the counter reaches 0, the batch of
output is written on the output device.

"""
pgebuf() = ccall((:cpgebuf, pgplotlib), Cvoid, (), )

"""

```julia
pgend()
```

Close and release any open graphics devices. All devices must be closed by
calling either [`pgclos`](@ref) (for each device) or `pgend` before the program
terminates. If a device is not closed properly, some or all of the graphical
output may be lost.

"""
pgend() = ccall((:cpgend, pgplotlib), Cvoid, (), )

"""

```julia
pgenv(xmin, xmax, ymin, ymax, just, axis)
```

Set PGPlot "Plotter Environment".  `pgenv` establishes the scaling for
subsequent calls to [`pgpt`](@ref), [`pgline`](@ref), etc.  The plotter is
advanced to a new page or panel, clearing the screen if necessary.  If the
"prompt state" is ON (see [`pgask`](@ref)), confirmation is requested from the
user before clearing the screen.  If requested, a box, axes, labels, etc. are
drawn according to the setting of argument `axis`.

Arguments:

- `xmin`: the world x-coordinate at the bottom left corner of the viewport.

- `xmax`: the world x-coordinate at the top rightthat ` corner` of the viewport (note
  that `xmax` may be less than `xmin`).

- `ymin`: the world y-coordinate at the bottom left corner of the viewport.

- `ymax`: the world y-coordinate at the top right corner of the viewport (note
  that `ymaxthat `` may` be less than `ym`in`)`.


- `just`: if true, the scales of the x and y axes (in world coordinates per
  inch) will be equal, otherwise they will be scaled independently.

- `axis`: control s the plotting of axes, tick marks, etc:

  - `axis = -2`: draw no box, axes or labels;
  - `axis = -1`: draw box only;
  - `axis =  0`: draw box and label it with coordinates;
  - `axis =  1`: same as AXIS=0, but also draw the coordinate axes (X=0, Y=0);
  - `axis =  2`: same as AXIS=1, but also draw grid lines at major increments
     of the coordinates;
  - `axis = 10`: draw box and label X-axis logarithmically;
  - `axis = 20`: draw box and label Y-axis logarithmically;
  - `axis = 30`: draw box and label both axes logarithmically.

For other axis options, use routine [`pgbox`](@ref). `pgenv` can be persuaded to
call [`pgbox`](@ref) with additional axis options by defining an environment
parameter `PGPLOT_ENVOPT` containing the required option codes. Examples:

    PGPLOT_ENVOPT=P      # draw Projecting tick marks
    PGPLOT_ENVOPT=I      # Invert the tick marks
    PGPLOT_ENVOPT=IV     # Invert tick marks and label y Vertically

"""
function pgenv(xmin::Real, xmax::Real, ymin::Real, ymax::Real,
               just::Bool, axis::Integer)
    ccall((:cpgenv, pgplotlib), Cvoid,
          (PGFloat, PGFloat, PGFloat, PGFloat, PGInt, PGInt),
          xmin, xmax, ymin, ymax, just, axis)
end

"""

```julia
pgeras()
```

Erase all graphics from the current page (or current panel, if the view surface
has been divided into panels with [`pgsubp`](@ref)).

"""
pgeras() = ccall((:cpgeras, pgplotlib), Cvoid, (), )

"""

```julia
pgerrb(dir, x, y, e, t)
```

Plot error bar(s) in the direction specified by `dir`.  This routine only draws
error bar(s); to mark the data point at the start of the error bar(s), an
additional call to [`pgpt`](@ref) is required.

Arguments:

- `dir`: direction to plot the error bar relative to the data point.  One-sided
  error bar, dir` is:

  - 1 for `+x` (`x` to `x+e`);
  - 2 for `+y` (`y` to `y+e`);
  - 3 for `-x` (`x` to `x-e`);
  - 4 for `-y` (`y` to `y-e`).

  Two-sided error bar, dir` is:

  - 5 for `+/-x` (`x-e` to `x+e`);
  - 6 for `+/-y` (`y-e` to `y+e`).

- `x`: world x-coordinate of the data.

- `y`: world y-coordinate of the data.

- `e`: value of error bar distance to be added to the data position in world
   coordinates.

- `t`: length of terminals to be drawn at the ends of the error bar, as a
  multiple of the default length; if `T = 0.0`, no terminals will be drawn.

"""
pgerrb(dir::Integer, x::Real, y::Real, e::Real, t::Real) =
    ccall((:cpgerr1, pgplotlib), Cvoid,
          (PGInt, PGFloat, PGFloat, PGFloat, PGFloat),
          dir, x, y, e, t)

const pgerr1 = pgerrb

function pgerrb(dir::Integer,
                x::AbstractVector{<:Real},
                y::AbstractVector{<:Real},
                e::AbstractVector,
                t::Real)
    pgerrb(dir, pgfltarr(x), pgfltarr(y), pgfltarr(e), t)
end

function pgerrb(dir::Integer,
                x::PGVecOrRef{PGFloat},
                y::PGVecOrRef{PGFloat},
                e::PGVecOrRef{PGFloat},
                t::Real)
    @check_vectors n x y e
    ccall((:cpgerrb, pgplotlib), Cvoid,
          (PGInt, PGInt, Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat}, PGFloat),
          dir, n, x, y, e, t)
end

"""

```julia
pgerrx(x1, x2, y, t)
```

Plot horizontal error bar(s).  This routine only draws error bar(s); to mark
the data point in the middle of the error bar(s), an additional call to
[`pgpt`](@ref) or [`pgerry`](@ref) is required.

Arguments:

- `x1`: world x-coordinates of lower end of the error bars.
- `x2`: world x-coordinates of upper end of the error bars.
- `y`: world y-coordinates of the data.
- `t`: length of terminals to be drawn at the ends of the error bar, as a
  multiple of the default length; if `t = 0.0`, no terminals will be drawn.

"""
pgerrx(x1::Real, x2::Real, y::Real, t::Real) =
    pgerrx(Ref{PGFloat}(x1), Ref{PGFloat}(x2), Ref{PGFloat}(y), t)

function pgerrx(x1::AbstractVector{<:Real},
                x2::AbstractVector{<:Real},
                y::AbstractVector{<:Real},
                t::Real)
    pgerrx(pgfltarr(x1), pgfltarr(x2), pgfltarr(y), t)
end

function pgerrx(x1::PGVecOrRef{PGFloat},
                x2::PGVecOrRef{PGFloat},
                y::PGVecOrRef{PGFloat},
                t::Real)
    @check_vectors n x1 x2 y
    ccall((:cpgerrx, pgplotlib), Cvoid,
          (PGInt, Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat}, PGFloat),
          n, x1, x2, y, t)
end

"""

```julia
pgerry(x, y1, y2, t)
```

Plot vertical error bar(s).  This routine only draws error bar(s); to mark the
data point in the middle of the error bar(s), an additional call to
[`pgpt`](@ref) or [`pgerrx`](@ref) is required.

Arguments:

- `x`: world x-coordinates of the data.
- `y1`: world y-coordinates of top end of the error bars.
- `y2`: world y-coordinates of of bottom end of the error bars.
- `t`: length of terminals to be drawn at the ends of the error bar, as a
  multiple of the default length; if `t = 0.0`, no terminals will be drawn.

"""
pgerry(x::Real, y1::Real, y2::Real, t::Real) =
    pgerry(Ref{PGFloat}(x), Ref{PGFloat}(y1), Ref{PGFloat}(y2), t)

function pgerry(x::AbstractVector{<:Real},
                y1::AbstractVector{<:Real},
                y2::AbstractVector{<:Real},
                t::Real)
    pgerry(pgfltarr(x), pgfltarr(y1), pgfltarr(y2), t)
end

function pgerry(x::PGVecOrRef{PGFloat},
                y1::PGVecOrRef{PGFloat},
                y2::PGVecOrRef{PGFloat},
                t::Real)
    @check_vectors n x y1 y2
    ccall((:cpgerry, pgplotlib), Cvoid,
          (PGInt, Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat}, PGFloat),
          n, x, y1, y2, t)
end

"""

```julia
pgetxt()
```

Some graphics terminals display text (the normal interactive dialog) on the
same screen as graphics. This routine erases the text from the view surface
without affecting the graphics. It does nothing on devices which do not display
text on the graphics screen, and on devices which do not have this capability.

"""
pgetxt() = ccall((:cpgetxt, pgplotlib), Cvoid, (), )


"""

```julia
pgfunt(fx, fy, t, flag = false)
```

Draw a curve defined by parametric equations `x = fx(t)`, `y = fy(t)` for all
values in `t`.  Arguments `fx` and `fy` must be callable objects.  If `flag` is
true, the curve is plotted in the current window and viewport; otherwise,
[`pgenv`](@ref) is called automatically by [`pgfunt](@ref) to start a new plot
with automatic scaling.

"""
function pgfunt(fx, fy, t::AbstractVector, flag::Bool = false)
    (n = length(t)) ≥ 2 ||
        throw(ArgumentError("`t` must have at least 2 values"))
    pgbbuf()
    if flag
        # Directly draw function in current viewport.
        for val in t
            x = PGFloat(fx(val))
            y = PGFloat(fy(val))
            if flag
                flag = false
                pgmove(x, y)
            else
                pgdraw(x, y)
            end
        end
    else
        # Define viewport according to extreme values.
        xmin = ymin = typemax(PGFloat)
        xmax = ymax = typemin(PGFloat)
        wrk = Array{PGFloat}(undef, 2, n)
        i = 0
        @inbounds for val in t
            i += 1; @assert i ≤ n
            xi = PGFloat(fx(val))
            yi = PGFloat(fy(val))
            xmin = min(xmin, xi)
            xmax = max(xmax, xi)
            ymin = min(ymin, yi)
            ymax = max(ymax, yi)
            wrk[1,i] = xi
            wrk[2,i] = yi
        end
        @assert i == n
        pgenv(autorange(xmin, xmax)..., autorange(ymin, ymax)..., false, 0)
        pgmove(wrk[1,1], wrk[2,1])
        @inbounds for i in 2:n
            pgdraw(wrk[1,i], wrk[2,i])
        end
    end
    pgebuf()
end

"""

```julia
pgfunx(x, fy, flag = false)
```

Draw a curve defined by `y = fy(x)` for all values in `x`.  Argument `fy` must
be a callable object.  If `flag` is true, the curve is plotted in the current
window and viewport; otherwise, [`pgenv`](@ref) is called automatically by
[`pgfunx](@ref) to start a new plot with automatic scaling.

"""
function pgfunx(x::AbstractVector, fy, flag::Bool = false)
    (n = length(x)) ≥ 2 ||
        throw(ArgumentError("`x` must have at least 2 values"))
    pgbbuf()
    if flag
        # Directly draw function in current viewport.
        for val in x
            y = fy(val)
            if flag
                flag = false
                pgmove(val, y)
            else
                pgdraw(val, y)
            end
        end
    else
        # Define viewport according to extreme values.
        xmin = ymin = typemax(PGFloat)
        xmax = ymax = typemin(PGFloat)
        wrk = Array{PGFloat}(undef, 2, n)
        i = 0
        @inbounds for val in x
            i += 1; @assert i ≤ n
            xi = PGFloat(val)
            yi = PGFloat(fy(val))
            xmin = min(xmin, xi)
            xmax = max(xmax, xi)
            ymin = min(ymin, yi)
            ymax = max(ymax, yi)
            wrk[1,i] = xi
            wrk[2,i] = yi
        end
        @assert i == n
        pgenv(autorange(xmin, xmax)..., autorange(ymin, ymax)..., false, 0)
        pgmove(wrk[1,1], wrk[2,1])
        @inbounds for i in 2:n
            pgdraw(wrk[1,i], wrk[2,i])
        end
    end
    pgebuf()
end

"""

```julia
pgfuny(fx, y, flag = false)
```

Draw a curve defined by `x = fx(y)` for all values in `y`.  Argument `fx` must
be a callable object.  If `flag` is true, the curve is plotted in the current
window and viewport; otherwise, [`pgenv`](@ref) is called automatically by
[`pgfuny](@ref) to start a new plot with automatic scaling.

"""
function pgfuny(fx, y::AbstractVector, flag::Bool = false)
    (n = length(y)) ≥ 2 ||
        throw(ArgumentError("`y` must have at least 2 values"))
    pgbbuf()
    if flag
        # Directly draw function in current viewport.
        for val in y
            x = fx(val)
            if flag
                flag = false
                pgmove(x, val)
            else
                pgdraw(x, val)
            end
        end
    else
        # Define viewport according to extreme values.
        xmin = ymin = typemax(PGFloat)
        xmax = ymax = typemin(PGFloat)
        wrk = Array{PGFloat}(undef, 2, n)
        i = 0
        @inbounds for val in y
            i += 1; @assert i ≤ n
            xi = PGFloat(fx(val))
            yi = PGFloat(val)
            xmin = min(xmin, xi)
            xmax = max(xmax, xi)
            ymin = min(ymin, yi)
            ymax = max(ymax, yi)
            wrk[1,i] = xi
            wrk[2,i] = yi
        end
        @assert i == n
        pgenv(autorange(xmin, xmax)..., autorange(ymin, ymax)..., false, 0)
        pgmove(wrk[1,1], wrk[2,1])
        @inbounds for i in 2:n
            pgdraw(wrk[1,i], wrk[2,i])
        end
    end
    pgebuf()
end

"""

```julia
pggray(A, [i1, i2, j1, j2,] fg, bg, tr=[0,1,0, 0,0,1])
```

draws a gray-scale map of an array in current window. The subsection of the
array `A` defined by indices `(i1:i2, j1:j2)` is mapped onto the view surface
world-coordinate system by the transformation matrix `tr`. The resulting
quadrilateral region is clipped at the edge of the window and shaded with the
shade at each point determined by the corresponding array value.  The shade is
a number in the range 0 to 1 obtained by linear interpolation between the
background level (`bg`) and the foreground level (`fg`), i.e.,

    shade = (A[i,j] - bg)/(fg - bg)

The background level `bg` can be either less than or greater than the
foreground level `fg`.  Points in the array that are outside the range `bg` to
`fg` are assigned shade 0 or 1 as appropriate.

`pggray` uses two different algorithms, depending how many color indices are
available in the color index range specified for images.  (This range is set
with routine [`pgscir`](@ref), and the current or default range can be queried
by calling routine [`pgqcir`](@ref)).

If 16 or more color indices are available, `pggray` first assigns color
representations to these color indices to give a linear ramp between the
background color (color index 0) and the foreground color (color index 1), and
then calls [`pgimag`](@ref) to draw the image using these color indices. In
this mode, the shaded region is "opaque": every pixel is assigned a color.

If less than 16 color indices are available, `pggray` uses only color index 1,
and uses a *dithering* algorithm to fill in pixels, with the shade (computed as
above) determining the faction of pixels that are filled. In this mode the
shaded region is "transparent" and allows previously-drawn graphics to show
through.

The transformation matrix `tr` is used to calculate the world coordinates of
the center of the *cell* that represents each array element. The world
coordinates of the center of the cell corresponding to array element `A[i,j]`
are given by:

    x = tr[1] + tr[2]*i + tr[3]*j
    y = tr[4] + tr[5]*i + tr[6]*j

Usually `tr[3]` and `tr[5]` are zero -- unless the coordinate transformation
involves a rotation or shear.  The corners of the quadrilateral region that is
shaded by `pggray` are given by applying this transformation to
`(i1-0.5,j1-0.5)`, `(i2+0.5, j2+0.5)`.

Intervals `I = i1:i2` and `J = j1:j2` may be used to specify index ranges.
If `i1`, `i2`, `j1` and `j2` are omitted, the full array `A` is considered.

"""
function pggray(A::AbstractMatrix{<:Real},
                fg::Real, bg::Real,
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = default_xform())
    _pggray(submatrix(A, tr)..., fg, bg)
end

function pggray(A::AbstractMatrix{<:Real},
                I::AbstractUnitRange{<:Integer},
                J::AbstractUnitRange{<:Integer},
                fg::Real, bg::Real,
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = default_xform())
    _pggray(submatrix(A, I, J, tr)..., fg, bg)
end

function pggray(A::AbstractMatrix{<:Real},
                i1::Integer, i2::Integer,
                j1::Integer, j2::Integer,
                fg::Real, bg::Real,
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = default_xform())
    _pggray(submatrix(A, i1, i2, j1, j2, tr)..., fg, bg)
end

# Argument `tr` moved to directly use the output of `submatrix`.
_pggray(a, i1, i2, j1, j2, tr, fg, bg) =
    ccall((:cpggray, pgplotlib), Cvoid,
          (Ptr{PGFloat}, PGInt, PGInt, PGInt, PGInt,
           PGInt, PGInt, PGFloat, PGFloat, Ptr{PGFloat}),
          a, idim, jdim, i1, i2, j1, j2, fg, bg, tr)

"""

```julia
pghi2d(A, [i1, i2, j1, j2,] x, off, bias, center,
       wrk = Vector{Cfloat}(undef, length(x)))
```

plots a series of cross-sections through a 2D data array `A`.  Each
cross-section is plotted as a hidden line histogram.  The plot can be slanted
to give a pseudo-3D effect -- if this is done, the call to [`pgenv`](@ref) may
have to be changed to allow for the increased X range that will be needed.

`pghi2d` plots a subset of the input array DATA.  This subset is delimited in
the first (x) dimension by `i1` and `i2` and the 2nd (y) by `j1` and `j2`,
inclusively. Note: `j2 < j1` is permitted, resulting in a plot with the
cross-sections plotted in reverse Y order.  However, `i2 ≥ i1` is mandatory.

Intervals `I = i1:i2` and `J = j1:j2` may be used to specify index ranges.
If `i1`, `i2`, `j1` and `j2` are omitted, the full array `A` is considered.

Argument `x` specifies the abscissae of the bins to be plotted. That is, `x[1]`
should be the X value for `A[i1,j1]`, and `x` should have `i2-i1+1` elements.
The program assumes that the X value for `A(x,y)` is the same for all y.

`off` is an offset in array elements applied to successive cross-sections to
produce a slanted effect.  A plot with `off >` 0 slants to the right, one with
`off < 0` slants left.

`bias` is a bias value applied to each successive cross-section in order to
raise it above the previous cross-section.  This is in the same units as the
values in `A`.

If `center` is `true`, the X values denote the center of the bins; otherwise
the X values denote the lower edges (in X) of the bins.

`wrk` is an optional workspace.  Should be an array of at least `i2-i1+1`
elements.

"""
function pghi2d(A::AbstractMatrix{<:Real},
                x::AbstractVector{<:Real},
                off::Integer,
                bias::Real,
                center::Bool,
                wrk::PGVector{PGFloat} = Vector{PGFloat}(undef, length(x)))
    I, J = axes(A)
    pghi2d(A, I, J, x, off, bias, center, wrk)
end

function pghi2d(A::AbstractMatrix{<:Real},
                I::AbstractUnitRange{<:Integer},
                J::AbstractUnitRange{<:Integer},
                x::AbstractVector{<:Real},
                off::Integer,
                bias::Real,
                center::Bool,
                wrk::PGVector{PGFloat} = Vector{PGFloat}(undef, length(x)))
    pghi2d(A, first(I), last(I), first(J), last(J), x, off, bias, center, wrk)
end

function pghi2d(A::AbstractMatrix{<:Real},
                i1::Integer, i2::Integer,
                j1::Integer, j2::Integer,
                x::AbstractVector{<:Real},
                off::Integer,
                bias::Real,
                center::Bool,
                wrk::PGVector{PGFloat} = Vector{PGFloat}(undef, length(x)))
    i1 ≤ i2 || error("invalid range `i1:i2`")
    length(x) == i2 + 1 - i1 || error("`x` must have `i2 + 1 - i1` values")
    length(wrk) ≥ length(x) || error("`wrk` must have at least `i2 + 1 - i1` values")
    _pghi2d(submatrix(PGFloat, i1, i2, j1, j2)..., pgfltarr(x), ioff, bias, center, wrk)
end

_pghi2d(a, idim, jdim, i1, i2, j1, j2, x, ioff, bias, center, wrk) =
    ccall((:cpghi2d, pgplotlib), Cvoid,
          (Ptr{PGFloat}, PGInt, PGInt, PGInt, PGInt,
           PGInt, PGInt, Ptr{PGFloat}, PGInt, PGFloat,
           PGBool, Ptr{PGFloat}),
          a, idim, jdim, i1, i2, j1, j2, x, ioff, bias, center, wrk)

"""

```julia
```

Draw a histogram of N values of a variable in array
DATA(1...N) in the range DATMIN to DATMAX using NBIN bins.  Note
that array elements which fall exactly on the boundary between
two bins will be counted in the higher bin rather than the
lower one; and array elements whose value is less than DATMIN or
greater than or equal to DATMAX will not be counted at all.

Arguments:
 N      (input)  : the number of data values.
 DATA   (input)  : the data values. Note: the dimension of array
                   DATA must be greater than or equal to N. The
                   first N elements of the array are used.
 DATMIN (input)  : the minimum data value for the histogram.
 DATMAX (input)  : the maximum data value for the histogram.
 NBIN   (input)  : the number of bins to use: the range DATMIN to
                   DATMAX is divided into NBIN equal bins and
                   the number of DATA values in each bin is
                   determined by PGHIST.  NBIN may not exceed 200.
 PGFLAG (input)  : if PGFLAG = 1, the histogram is plotted in the
                   current window and viewport; if PGFLAG = 0,
                   PGENV is called automatically by PGHIST to start
                   a new plot (the x-limits of the window will be
                   DATMIN and DATMAX; the y-limits will be chosen
                   automatically.
                   IF PGFLAG = 2,3 the histogram will be in the same
                   window and viewport but with a filled area style.
                   If pgflag=4,5 as for pgflag = 0,1, but simple
                   line drawn as for PGBIN

"""
function pghist(data::AbstractVector, datmin::Real, datmax::Real,
                nbin::Integer, flag::Integer)
    _cpghist(length(data), pgfltarr(data), PGFloat(datmin), PGFloat(datmax),
             PGInt(nbin), PGInt(flag))
end

_cpghist(n, data, datmin, datmax, nbin, pgflag) =
    ccall((:cpghist, pgplotlib), Cvoid,
          (PGInt, Ptr{PGFloat}, PGFloat, PGFloat, PGInt, PGInt),
          n, data, datmin, datmax, nbin, pgflag)

"""

```julia
pgiden()
```

writes username, date, and time at bottom of plot.

"""
pgiden() = ccall((:cpgiden, pgplotlib), Cvoid, (), )

"""

```julia
pgimag()
```

"""
function pgimag(A::AbstractMatrix{<:Real},
                a1::Real, a2::Real,
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgimag(submatrix(A, tr)..., a1, a2)
end

function pgimag(A::AbstractMatrix{<:Real},
                I::AbstractUnitRange{<:Integer},
                J::AbstractUnitRange{<:Integer},
                a1::Real, a2::Real,
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgimag(submatrix(A, I, J, tr)..., a1, a2)
end

function pgimag(A::AbstractMatrix{<:Real},
                i1::Integer, i2::Integer,
                j1::Integer, j2::Integer,
                a1::Real, a2::Real,
                tr::Union{AbstractVector{<:Real},
                          AbstractMatrix{<:Real}} = DEFAULT_XFORM)
    _pgimag(submatrix(A, i1, i2, j1, j2, tr)..., a1, a2)
end

# Argument `tr` moved to directly use the output of `submatrix`.
_pgimag(a, idim, jdim, i1, i2, j1, j2, tr, a1, a2) =
    ccall((:cpgimag, pgplotlib), Cvoid,
          (Ptr{PGFloat}, PGInt, PGInt, PGInt, PGInt, PGInt, PGInt,
           PGFloat, PGFloat, Ptr{PGFloat}),
          a, idim, jdim, i1, i2, j1, j2, a1, a2, tr)

"""

```julia
pglab()
```

"""
pglab(xlbl::AbstractString, ylbl::AbstractString, toplbl::AbstractString) =
    ccall((:cpglab, pgplotlib), Cvoid,
          (PGString, PGString, PGString),
          xlbl, ylbl, toplbl)

"""

```julia
pglcur!(x, y, n=0) -> npt
```

Interactive routine for user to enter a polyline by use of the cursor.  Routine
allows user to Add and Delete vertices; vertices are joined by straight-line
segments.

Arguments `x` and `y` are vectors of single precision floating point values
which are overwritten by the coordinates of the selected points.

Optional argument `n` is the number of valid points already defined in `x` and
`y`, these points will be plotted first.

The returned value is the total number of defined points, that `n` plus the
selected points during the call of `pglcur!`.

Notes:

1. On return from the program, cursor points are returned in the order they
   were entered. Routine may be (re-)called with points already defined in `x`
   and `y` (with `n > 0`), and they will be plotted first, before editing.

2. User commands: the user types single-character commands after positioning
   the cursor: the following are accepted:

  - `A` / mouse button 1 (Add) to add point at current cursor location.
  - `D` / mouse button 2 (Delete) to delete the last point entered.
  - `X` / any other mouse button (eXit) to leave subroutine.

"""
function pglcur!(x::PGVector{PGFloat}, y::PGVector{PGFloat}, n::Integer = 0)
    @check_vectors maxpt x y
    0 ≤ n ≤ maxpt || throw(ArgumentError("bad number of previous points"))
    npt = Ref{PGInt}(n)
    ccall((:cpglcur, pgplotlib), Cvoid,
          (PGInt, Ref{PGInt}, Ptr{PGFloat}, Ptr{PGFloat}),
          maxpt, npt, x, y)
    return output(npt)
end

_pglcur(maxpt, npt, x, y) =
    ccall((:cpglcur, pgplotlib), Cvoid,
          (PGInt, Ref{PGInt}, Ptr{PGFloat}, Ptr{PGFloat}),
          maxpt, npt, x, y)

"""

```julia
pgldev()
```

"""
pgldev() = ccall((:cpgldev, pgplotlib), Cvoid, (), )

# Cvoid cpglen(PGInt units, PGString string, Ptr{PGFloat} xl, Ptr{PGFloat} yl);
"""

```julia
pglen()
```

"""
function pglen(units::Integer, string::AbstractString, xl::PGVector{PGFloat}, yl::PGVector{PGFloat})
    ccall((:cpglen, pgplotlib), Cvoid,
          (PGInt, PGString, Ptr{PGFloat}, Ptr{PGFloat}),
          units, string, xl, yl)
end

"""

```julia
pgline(x, y)
```

draws a polyline as connected straight-line segments whose end-points are
`(x,y)`.  The polyline is drawn using the current setting of attributes
color-index, line-style, and line-width. The polyline is clipped at the
edge of the window.

"""
function pgline(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    @check_vectors n x y
    _cpgline(n, pgfltarr(x), pgfltarr(y))
end

_cpgline(n::Int, x::PGVector{PGFloat}, y::PGVector{PGFloat}) =
    ccall((:cpgline, pgplotlib), Cvoid,
          (PGInt, Ptr{PGFloat}, Ptr{PGFloat}), n, x, y)

"""
Move pen (change current pen position)
======================================

    pgmove(x, y)

moves the *pen* to the point with world coordinates `(x,y)`.  No line is
drawn.

"""
pgmove(x::Real, y::Real) =
    ccall((:cpgmove, pgplotlib), Cvoid, (PGFloat, PGFloat), x, y)


"""
Write text at position relative to viewport
===========================================

    cpgmtxt(side, disp, coord, just, text)

writes text at a position specified relative to the viewport (outside or
inside).  This routine is useful for annotating graphs. It is used by routine
`pglab`.  The text is written using the current values of attributes
color-index, line-width, character-height, and character-font.

Arguments:

* `side` : must include one of the characters 'B', 'L', 'T', or 'R' signifying
  the Bottom, Left, Top, or Right margin of the viewport. If it includes 'LV'
  or 'RV', the string is written perpendicular to the frame rather than
  parallel to it.

* `disp` : the displacement of the character string from the specified edge of
  the viewport, measured outwards from the viewport in units of the character
  height.  Use a negative value to write inside the viewport, a positive value
  to write outside.

* `coord` : the location of the character string along the specified edge of
  the viewport, as a fraction of the length of the edge.

* `just` : controls justification of the string parallel to the specified edge
  of the viewport. If `just = 0`, the left-hand end of the string will be
  placed at `coord`; if `just` = 0.5, the center of the string will be placed
  at `coord`; if `just = 1`, the right-hand end of the string will be placed at
  at `coord`.  Other values between 0 and 1 give intermediate placing, but they
  are not very useful.

* `text` : the text string to be plotted.  Trailing spaces are ignored when
  justifying the string, but leading spaces are significant.

"""
function pgmtxt(side::AbstractString, disp::Real, coord::Real, just::Real,
                text::AbstractString)
    ccall((:cpgmtxt, pgplotlib), Cvoid,
          (PGString, PGFloat, PGFloat, PGFloat, PGString),
          side, disp, coord, just, text)
end

"""

```julia
pgncur!(x, y, sym, n=0) -> npt
```

Interactive routine for user to enter data points by use of the cursor.
Routine allows user to Add and Delete points.  The points are returned in order
of increasing x-coordinate, not in the order they were entered (unlike
[`pgolin!`](@ref)).

Arguments `x` and `y` are vectors of single precision floating point values
which are overwritten by the coordinates of the selected points.

Argument `sym` is the code number of symbol to use for marking entered points
(see [`pgpt`](@ref)).

Optional argument `n` is the number of valid points already defined in `x` and
`y`, these points will be plotted first.

The returned value is the total number of defined points, that `n` plus the
selected points during the call of `pgncur!`.

Notes:

1. On return from the program, cursor points are returned in the order they
   were entered.  Routine may be (re-)called with points already defined in `x`
   and `y` (with `n > 0`), and they will be plotted first, before editing.

2. User commands: the user types single-character commands after positioning
   the cursor: the following are accepted:

  - `A` / mouse button 1 (Add) to add point at current cursor location.
  - `D` / mouse button 2 (Delete) to delete the last point entered.
  - `X` / any other mouse button (eXit) to leave subroutine.


"""
function pgncur!(x::PGVector{PGFloat}, y::PGVector{PGFloat},
                 sym::Integer, n::Integer = 0)
    @check_vectors maxpt x y
    0 ≤ n ≤ maxpt || throw(ArgumentError("`n` is out of range"))
    npt = Ref{PGInt}(n)
    _pgncur(maxpt, npt, x, y, sym)
    return output(npt)
end

_pgncur(maxpt, npt, x, y, sym) =
    ccall((:cpgncur, pgplotlib), Cvoid,
          (PGInt, Ref{PGInt}, Ptr{PGFloat}, Ptr{PGFloat}, PGInt),
          maxpt, npt, x, y, sym)

# Cvoid cpgnumb(PGInt mm, PGInt pp, PGInt form, PGString string, Ptr{PGInt} string_length);
"""

```julia
pgnumb()
```

"""
function pgnumb(mm::Integer, pp::Integer, form::Integer, string::AbstractString, string_length::PGVector{PGInt})
    ccall((:cpgnumb, pgplotlib), Cvoid,
          (PGInt, PGInt, PGInt, PGString, Ptr{PGInt}),
          mm, pp, form, string, string_length)
end

"""

```julia
pgpolin!(x, y, sym, n=0) -> npt
```

Interactive routine for user to enter data points by use of the cursor.
Routine allows user to Add and Delete points.  The points are returned in the
order that they were entered (unlike [`pgncur!`](@ref)).

Arguments `x` and `y` are vectors of single precision floating point values
which are overwritten by the coordinates of the selected points.

Argument `sym` is the code number of symbol to use for marking entered points
(see [`pgpt`](@ref)).

Optional argument `n` is the number of valid points already defined in `x` and
`y`, these points will be plotted first.

The returned value is the total number of defined points, that `n` plus the
selected points during the call of `pgolin!`.

Notes:

1. On return from the program, cursor points are returned in the order they
   were entered.  Routine may be (re-)called with points already defined in `x`
   and `y` (with `n > 0`), and they will be plotted first, before editing.

2. User commands: the user types single-character commands after positioning
   the cursor: the following are accepted:

  - `A` / mouse button 1 (Add) to add point at current cursor location.
  - `D` / mouse button 2 (Delete) to delete the last point entered.
  - `X` / any other mouse button (eXit) to leave subroutine.

"""
function pgolin!(x::PGVector{PGFloat}, y::PGVector{PGFloat},
                 sym::Integer, n::Integer = 0)
    @check_vectors maxpt x y
    0 ≤ n ≤ maxpt || throw(ArgumentError("bad number of previous points"))
    npt = Ref{PGInt}(n)
    _pgolin(maxpt, npt, x, y, sym)
    return output(npt)
end

_pgolin(maxpt, npt, x, y, sym) =
    ccall((:cpgolin, pgplotlib), Cvoid,
          (PGInt, Ref{PGInt}, Ptr{PGFloat}, Ptr{PGFloat}, PGInt),
          maxpt, npt, x, y, sym)

"""

```julia
pgopen(dev="?") -> id
```

Opens a graphics device for PGPlot output. If the device is opened successfully,
it becomes the selected device to which graphics output is directed until
another device is selected with [`pgslct`](@ref) or the device is closed with
[`pgclos`](@ref).

Returns a positive value, the identifier of the graphics device for use with
[`pgslct`](@ref). In the event of error a message is written on the standard
error output and a `PGPlotError` exception is thrown.

The `dev` argument is a text string; its value should be one of the following:

1. A complete device specification of the form `device/type` or `file/type`,
   where `type` is one of the allowed PGPLOT device types
   (installation-dependent) and `device` or `file` is the name of a graphics
   device or disk file appropriate for this type.  The `device` or `file` may
   contain `/` characters; the final `/` delimits the `type`. If necessary to
   avoid ambiguity, the `device` part of the string may be enclosed in double
   quotation marks.

2. A device specification of the form `/type`, where `type` is one
    of the allowed PGPLOT device types. PGPLOT supplies a default
    file or device name appropriate for this device type.

3. A device specification with `/type` omitted; in this case the type is taken
   from the environment variable `PGPLOT_TYPE`, if defined (e.g., `setenv
   PGPLOT_TYPE PS`). Because of possible confusion with `/` in file-names,
   omitting the device type in this way is not recommended.

4. A blank string (` `); in this case, [`pgopen`](@ref) will use the value of
   environment variable `PGPLOT_DEV` as the device specification, or `/NULL` if
   the environment variable is undefined.

5. A single question mark, with optional trailing spaces (`?`); in this case,
   PGPlot will prompt the user to supply the device specification, with a
   prompt string of the form:

   > Graphics device/type (? to see list, default XXX):

   where `XXX` is the default (value of environment variable `PGPLOT_DEV`).

6. A non-blank string in which the first character is a question mark (e.g.,
   `?Device: `); in this case, PGPlot will prompt the user to supply the device
   specification, using the supplied string as the prompt (without the leading
   question mark but including any trailing spaces).

In cases (5) and (6), the device specification is read from the standard
input. The user should respond to the prompt with a device specification of the
form (1), (2), or (3). If the user types a question-mark in response to the
prompt, a list of available device types is displayed and the prompt is
re-issued. If the user supplies an invalid device specification, the prompt is
re-issued. If the user responds with an end-of-file character, e.g., ctrl-D in
UNIX, program execution is aborted; this avoids the possibility of an infinite
prompting loop.  A programmer should avoid use of PGPlot-prompting if this
behavior is not desirable.

The device type is case-insensitive (e.g., `/ps` and `/PS` are equivalent). The
device or file name may be case-sensitive in some operating systems.

Examples of valid `dev` arguments:

1. `"plot.ps/ps"`, `"dir/plot.ps/ps"`, `"\"dir/plot.ps\"/ps"`,
   `"user:[tjp.plots]plot.ps/PS"`.

2. `"/ps"` (PGPlot interprets this as `pgplot.ps/ps`).

3. `"plot.ps"` (if `PGPLOT_TYPE` is defined as `ps`, PGPlot interprets this as
   `"plot.ps/ps"`).

4. `"   "` (if `PGPLOT_DEV` is defined).

5. `"?  "`

6.  `"?Device specification for PGPLOT: "`.

"""
function pgopen(dev::AbstractString = "?")
    id = ccall((:cpgopen, pgplotlib), PGInt, (PGString,), dev)
    id > 0 || pgerror(:pgopen, id)
    return id
end

"""

```julia
pgpage()
```

Advance plotter to a new page or panel, clearing the screen if
necessary. If the "prompt state" is `on` (see [`pgask`](@ref)),
confirmation is requested from the user before clearing the screen. If the
view surface has been subdivided into panels with [`pgbeg`](@ref) or
[`pgsubp`](@ref), then [`pgpage`](@ref) advances to the next panel, and if
the current panel is the last on the page, [`pgpage`](@ref) clears the
screen or starts a new sheet of paper.  [`pgpage`](@ref) does not change
the PGPlot window or the viewport (in normalized device coordinates); but
note that if the size of the view-surface is changed externally (e.g., by a
workstation window manager) the size of the viewport is changed in
proportion.

"""
pgpage() = ccall((:cpgpage, pgplotlib), Cvoid, (), )

"""

```julia
pgpanl()
```

"""
pgpanl(nxc::Integer, nyc::Integer) = ccall((:cpgpanl, pgplotlib), Cvoid,
                                           (PGInt, PGInt), nxc, nyc)

"""

```julia
pgpap()
```

"""
pgpap(width::Real, aspect::Real) = ccall((:cpgpap, pgplotlib), Cvoid,
                                         (PGFloat, PGFloat), width, aspect)

"""

```julia
pgpixl(A, [i1, i2, j1, j2,], x1=i1-0.5, x2=i2+0.5, y1, y2)
```

"""

function pgpixl(A::AbstractMatrix{<:Integer},
                x1::Real, x2::Real, y1::Real, y2::Real)
    _pgpixl(submatrix(PGInt, A)..., x1, x2, y1, y2)
end

function pgpixl(A::AbstractMatrix{<:Integer},
                I::AbstractUnitRange{<:Integer},
                J::AbstractUnitRange{<:Integer},
                x1::Real, x2::Real, y1::Real, y2::Real)
    _pgpixl(submatrix(PGInt, A, I, J)..., x1, x2, y1, y2)
end

function pgpixl(A::AbstractMatrix{<:Integer},
                i1::Integer, i2::Integer,
                j1::Integer, j2::Integer,
                x1::Real, x2::Real, y1::Real, y2::Real)
    _pgpixl(submatrix(PGInt, A, i1, i2, j1, j2)..., x1, x2, y1, y2)
end

_pgpixl(a, idim, jdim, i1, i2, j1, j2, x1, x2, y1, y2) =
    ccall((:cpgpixl, pgplotlib), Cvoid,
          (Ptr{PGInt}, PGInt, PGInt, PGInt, PGInt, PGInt, PGInt,
           PGFloat, PGFloat, PGFloat, PGFloat),
          a, idim, jdim, i1, i2, j1, j2, x1, x2, y1, y2)

"""

```julia
pgpnts(x, y, sym)
```

Draw Graph Markers at world coordinates `(x,y)` with symbols `sym`.  Unlike
PGPT, this routine can draw a different symbol at each point. The markers are
drawn using the current values of attributes color-index, line-width, and
character-height (character-font applies if the symbol number is >31).  If the
point to be marked lies outside the window, no marker is drawn.  The *pen
position* is changed to the last point (if any).

Argument `sym` specifies the code number of the symbol(s) to be plotted at each
point (see PGPT).  The first symbol is re-used for all `i`-th data such that `i
≥ length(sym)`.

"""
pgpnts(x::Real, y::Real, sym::PGIntegers) = pgpt(x, y, sym)

function pgpnts(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                sym::PGIntegers)
    @check_vectors n x y
    _cpgpnts(n, pgfltarr(x), pgfltarr(y), Ref{PGInt}(sym), 1)
end

function pgpnts(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                sym::AbstractVector{<:PGIntegers})
    @check_vectors n x y
    (ns = length(sym)) ≥ 1 ||
        throw(DimensionMismatch("`sym` must have at least one element"))
    _cpgpnts(n, pgfltarr(x), pgfltarr(y), pgintarr(sym), ns)
end

function _cpgpnts(n::Integer, x::PGVecOrRef{PGFloat}, y::PGVecOrRef{PGFloat},
                  sym::PGVecOrRef{PGInt}, ns::Integer)
    ccall((:cpgpnts, pgplotlib), Cvoid,
          (PGInt, Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGInt}, PGInt),
          n, x, y, sym, ns)
end

"""

```julia
pgpoly(x, y)
```

draws a polygon defined by the vertices `(x[i],y[i])` in world coordinates for
`i = 1:n` (with `n` the length of `x` and `y`) and using current fill-area
attributes (see [`pgsfs`](@ref)).  The polygon is clipped at the edge of the
window.  If the polygon is not convex, a point is assumed to lie inside the
polygon if a straight line drawn to infinity intersects and odd number of the
polygon's edges.  The pen position is changed to `(x[1],y[1])` in world
coordinates (if `n > 1`).

"""
function pgpoly(x::AbstractVector, y::AbstractVector)
    @check_vectors n x y
    _cpgpoly(n, pgfltarr(x), pgfltarr(y))
end

_cpgpoly(n, x, y) = ccall((:cpgpoly, pgplotlib), Cvoid,
                          (PGInt, Ptr{PGFloat}, Ptr{PGFloat}), n, x, y)

"""

```julia
pgpt(x, y, sym)
```

draws graph marker(s) at position(s) `(x,y)` with symbol `sym`.  If `sym`
is a vector, the first markers are drawn with `sym[1]`, `sym[2]`, *etc.*
and all the remaining markers are drawn with `sym[end]`.

The following symbols are available:
* `-1`, `-2`  : a single dot (diameter = current line width).
* `-3..-31` : a regular polygon with `abs(sym)` edges (style set by
  current fill style).
* `0..31` : standard marker symbols.
* `32..127` : ASCII characters (in current font); *e.g.*, to use letter F
  as a marker, let `sym = 'F'`.
* `sym > 127` : a Hershey symbol number.

"""
pgpt(x::Real, y::Real, sym::Union{Char,Integer}) = _cpgpt1(x, y, PGInt(sym))

function pgpt(x::AbstractVector, y::AbstractVector, sym::Union{Char,Integer})
    @check_vectors n x y
    _cpgpt(n, pgfltarr(x), pgfltarr(y), PGInt(sym))
end

_cpgpt(n, x, y, sym) = ccall((:cpgpt, pgplotlib), Cvoid,
                             (PGInt, Ptr{PGFloat}, Ptr{PGFloat}, PGInt),
                             n, x, y, sym)

_cpgpt1(x, y, sym) = ccall((:cpgpt1, pgplotlib), Cvoid,
                           (PGFloat, PGFloat, PGInt),
                           x, y, sym)

"""

```julia
pgptxt(x, y, angle, just, text)
```

"""
pgptxt(x::Real, y::Real, angle::Real, just::Real, text::AbstractString) =
    ccall((:cpgptxt, pgplotlib), Cvoid,
          (PGFloat, PGFloat, PGFloat, PGFloat, PGString),
          x, y, angle, just, text)


"""

```julia
pgqah() -> fs, angle, barb
```

yields arrow-head style: `fs` is 1 for filled, 2 for outline, `angle` is the
acute angle of the arrow point, in degrees and `barb` is the fraction of the
triangular arrow-head that is cut away from the back.

"""
pgqah() = query(PGInt, PGFloat, PGFloat, _pgqah)

_pgqah(fs, angle, barb) = ccall((:cpgqah, pgplotlib), Cvoid,
                                (Ptr{PGInt}, Ptr{PGFloat}, Ptr{PGFloat}),
                                fs, angle, barb)

"""

```julia
pgqcf() -> font
```

yields the current character font (set by routine [`pgscf`](@ref)).

"""
pgqcf() = query(PGInt, _pgqcf)

_pgqcf(font) = ccall((:cpgqcf, pgplotlib), Cvoid, (Ptr{PGInt},), font)

"""

```julia
pgqch() -> size
```

yields the current character size attribute (set by routine [`pgsch`](@ref)).

"""
pgqch() = query(PGFloat, _pgqch)

_pgqch(size) = ccall((:cpgqch, pgplotlib), Cvoid, (Ptr{PGFloat},), size)

"""

```julia
pgqci() -> ci
```

yields the current color index set by [`pgsci`](@ref).

"""
pgqci() = query(PGInt, _pgqci) # 7.297 ns (0 allocations: 0 bytes)

_pgqci(ci) = ccall((:cpgqci, pgplotlib), Cvoid, (Ptr{PGInt},), ci)

"""

```julia
pgqcir() -> (lo, hi)
```

yields color index range for plotting images with [`pggray`](@ref) or
[`pgimag`](@ref), as set by routine [`pgscir`](@ref) or by device default.

"""
pgqcir() = query(PGInt, PGInt, _pgqcir)

_pgqcir(lo, hi) = ccall((:cpgqcir, pgplotlib), Cvoid,
                        (Ptr{PGInt}, Ptr{PGInt}), lo, hi)

"""

```julia
pgqclp() -> state
```

yields whether clipping is active (set by routine [`pgsclp`](@ref)).

"""
pgqclp() = query(Bool, _pgqclp)

_pgqclp(state) = ccall((:cpgqclp, pgplotlib), Cvoid, (Ptr{PGInt},), state)

"""

```julia
pgqcol() -> (cmin, cmax)
```

yields the range of color indices available on the current device.  `cmin` is
the minimum available color index. This will be either 0 if the device can
write in the background color, or 1 if not.  `cmax` is the maximum available
color index. This will be 1 if the device has no color capability, or a larger
number (e.g., 3, 7, 15, 255).

"""
pgqcol() = query(PGInt, PGInt, _pgqcol)

_pgqcol(lo, hi) = ccall((:cpgqcol, pgplotlib), Cvoid,
                        (Ptr{PGInt}, Ptr{PGInt}), lo, hi)


"""

```julia
pgqcr(ci) -> (r,g,b)
```

yields the RGB colors associated with color index `ci`.  The result is a
3-tuple of values in the range 0.0 to 1.0.

"""
function pgqcr(ci::Integer)
    cr = Ref{PGFloat}()
    cg = Ref{PGFloat}()
    cb = Ref{PGFloat}()
    _pgqcr(ci, cr, cg, cb)
    output(cr, cg, cb)
end

_pgqcr(ci, cr, cg, cb) = ccall((:cpgqcr, pgplotlib), Cvoid,
                               (PGInt, Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat}),
                               ci, cr, cg, cb)

"""

```julia
pgqcs(units) -> (xch, ych)
```

yields the current PGPlot character height in a variety of units: `xch` is the
character height for text written with a vertical baseline while `ych` is the
character height for text written with a horizontal baseline (the usual case).
Argument `units` is one of:

- 0: normalized device coordinates;
- 1: inches;
- 2: millimeters;
- 3: pixels;
- 4: world coordinates.

Other values give an error message, and are treated as 0.

If `units=1` or `units=2`, `xch` and `ych` both have the same value.

If `units=3`, `xch` is the height in horizontal pixel units, and `ych` is the
height in vertical pixel units; on devices for which the pixels are not square,
`xch` and `ych` will be different.

If `units=4`, `xch` is the height in horizontal world coordinates (as used for
the x-axis), and `ych` is the height in vertical world coordinates (as used for
the y-axis). Unless special care has been taken to achieve equal
world-coordinate scales on both axes, the values of XCH and YCH will be
different.

If `units=0`, `xch` is the character height as a fraction of the horizontal
dimension of the view surface, and `ych` receives the character height as a
fraction of the vertical dimension of the view surface.

This routine provides facilities that are not available via [`pgqch`](@ref).
Use `pgqcs` if the character height is required in units other than those used
in [`pgsch`](@ref).

The [`pgplot`](@ref) *character height* is a dimension that scales with the
size of the view surface and with the scale-factor specified with routine
[`pgsch`](@ref). The default value is 1/40th of the height or width of the view
surface (whichever is less); this value is then multiplied by the scale-factor
supplied with [`pgsch`](@ref). Note that it is a nominal height only; the
actual character size depends on the font and is usually somewhat smaller.

"""
function pgqcs(units::Integer)
    xch = Ref{PGFloat}()
    ych = Ref{PGFloat}()
    _pgqcs(units, xch, ych)
    output(xch, ych)
end

_pgqcs(unit, xch, ych) = ccall((:cpgqcs, pgplotlib), Cvoid,
                               (PGInt, Ptr{PGFloat}, Ptr{PGFloat}),
                               units, xch, ych)

"""

```julia
pgqdt(n) -> (type, descr, inter)
```

returns informations about the `n`-th device type: the device type with a
description and whether the device is interactive.  The number of graphics
devices can be retrieved by calling [`pgqndt`](@ref).

"""
function pgqdt(id::Integer)
    type_buffer = zeros(UInt8, 9)
    type_length = Ref{PGInt}(0)
    descr_buffer = zeros(UInt8, 65)
    descr_length = Ref{PGInt}(0)
    inter = Ref{PGInt}(0)
    _pgqdt(id, type_buffer, type_length, descr_buffer, descr_length, inter)
    return (to_string(type_buffer, type_length),
            to_string(descr_buffer, descr_length),
            output(Bool, inter))
end

_pgqdt(id, type_buffer, type_length, descr_buffer, descr_length, inter) =
    ccall((:cpgqdt, pgplotlib), Cvoid,
          (PGInt, Ptr{UInt8}, Ptr{PGInt}, Ptr{UInt8}, Ptr{PGInt}, Ptr{PGInt}),
          id, type_buffer, type_length, descr_buffer, descr_length, inter)

# void cpgqfs(int *fs);
"""

```julia
pgqfs() -> fs
```

yields the current Fill-Area Style attribute (set by routine [`pgsfs`](@ref)):
`fs = 1` for *solid* (default), `fs = 2` for *outline*, `fs = 3` for *hatched*
and `fs = 4` for *cross-hatched*.

"""
pgqfs() = query(PGInt, _pgqfs)

_pgqfs(fs) = ccall((:cpgqfs, pgplotlib), Cvoid, (Ptr{PGInt},), fs)

"""

```julia
pgqhs() -> (angle, sepn, phase)
```

yields the style to be used hatching (fill area with fill-style 3):

- `angle` is the angle the hatch lines make with the horizontal, in degrees,
  increasing counterclockwise (this is an angle on the view surface, not in
  world-coordinate space);

- `sepn` is the spacing of the hatch lines. The unit spacing is 1 percent of
  the smaller of the height or width of the view surface;

- `phase` is a real number between 0 and 1; the hatch lines are displaced by
  this fraction of `sepn` from a fixed reference.  Adjacent regions hatched with
  the same `phase` have contiguous hatch lines.

"""
pgqhs() = query(PGFloat, PGFloat, PGFloat, _pgqhs)

_pgqhs(angle, sepn, phase) = ccall((:cpgqhs, pgplotlib), Cvoid,
                                   (Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat}),
                                   angle, sepn, phase)

"""

```julia
pgqid() -> id
```

returns the identifier of the currently selected device, or 0 if no device is
selected.  The identifier is assigned when [`pgopen`](@ref) is called to open
the device, and may be used as an argument to [`pgslct`](@ref).  Each open
device has a different identifier.

"""
pgqid() = query(PGInt, _pgqid)

_pgqid(id) = ccall((:cpgqid, pgplotlib), Cvoid, (Ptr{PGInt},), id)

"""

```julia
pgqinf(item) -> str
```

This routine can be used to obtain miscellaneous information about the PGPlot
environment.  Argument `item` is a character string defining the information
required, and output is a character string containing the requested
information.

The following item codes are accepted (note that the strings must match
exactly, except for case, but only the first 8 characters are significant). For
items marked `*`, PGPlot must be in the *open* state for the inquiry to
succeed. If the inquiry is unsuccessful, either because the item code is not
recognized or because the information is not available, a question mark (`"?"`)
is returned.

- `"VERSION"`: version of PGPLOT software in use.

- `"STATE"`: status of PGPLOT (`"OPEN"` if a graphics device is open for
  output, `"CLOSED"` otherwise).

- `"USER"`: the username associated with the calling program.

- `"NOW"`: current date and time (e.g., `"17-FEB-1986 10:04"`).

- `"DEVICE"` (*): current PGPLOT device or file.

- `"FILE"` (*): current PGPLOT device or file.

- `"TYPE"`(*): device-type of the current PGPLOT device.

- 'DEV/TYPE'(*): current PGPLOT device and type, in a form which is acceptable
  as an argument for PGBEG.

- `"HARDCOPY"`(*): is the current device a hardcopy device? (`"YES"` or
  `"NO"`).

- `"TERMINAL"` (*): is the current device the user's interactive terminal?
  (`"YES"` or `"NO"`).

- `"CURSOR"` (*): does the current device have a graphics cursor? (`"YES"` or
  `"NO"`).

- `"SCROLL"`(*): does current device have rectangle-scroll capability (`"YES"`
  or `"NO"`); see PGSCRL.

"""
function pgqinf(item::AbstractString)
    maxlen = 128 # FIXME: 64 should be enough (according to the demo code)
    buf = Array{UInt8}(undef, maxlen+1)
    len = Ref{PGInt}(maxlen)
    ccall((:cpgqinf, pgplotlib), Cvoid, (PGString, Ptr{UInt8}, Ptr{PGInt}),
          item, buf, len)
    len[] ≤ maxlen || error("too long string (unexpected)")
    i = len[]
    while i ≥ 1 && buf[i] == ' '
        i -= 1
    end
    buf[i+1] = 0
    return unsafe_string(pointer(buf))
end

"""

`pgqitf()` returns the Image Transfer Function as set by default or by a
previous call to [`pgsitf`](@ref). The Image Transfer Function is used by
[`pgimag`](@ref), [`pggray`](@ref), and [`pgwedg`](@ref).

"""
pgqitf() = query(PGInt, _pgqitf)

_pgqitf(itf) = ccall((:cpgqitf, pgplotlib), Cvoid, (Ptr{PGInt},), itf)

"""

`pgqls()` returns the current Line Style attribute (set by [`pgsls`](@ref)).

"""
pgqls() = query(PGInt, _pgqls)

_pgqls(ls) = ccall((:cpgqls, pgplotlib), Cvoid, (Ptr{PGInt},), ls)

"""

`pgqlw()` returns the current Line-Width attribute (set by [`pgslw`](@ref)).

"""
pgqlw() = query(PGInt, _pgqlw)

_pgqlw(lw) = ccall((:cpgqlw, pgplotlib), Cvoid, (Ptr{PGInt},), lw)

"""

`pgqndt()` returns the number of available device types. This routine is
usually used in conjunction with [`pgqdt`](@ref) to get a list of the available
device types.

"""
pgqndt() = query(PGInt, _pgqndt)

_pgqndt(n) = ccall((:cpgqndt, pgplotlib), Cvoid, (Ptr{PGInt},), n)

"""

`pgqpos() -> (x,y)` returns the current "pen" position in world coordinates.

"""
pgqpos() = query(PGFloat, PGFloat, _pgqpos)

_pgqpos(x, y) = ccall((:cpgqpos, pgplotlib), Cvoid,
                      (Ptr{PGFloat}, Ptr{PGFloat},), x, y)

"""

`pgqtbg()` returns the current Text Background Color Index (set by
[`pgstbg`](@ref)).

"""
pgqtbg() = query(PGInt, _pgqtbg)

_pgqtbg(tbci) = ccall((:cpgqtbg, pgplotlib), Cvoid, (Ptr{PGInt},), tbci)

"""

```julia
_pgqtxt(x, y, angle, just, text) -> (x1, x2, x3, x4), (y1, y2, y3, y4)
```

yields a bounding box for a text string.  Arguments `x`, `y`, `angle`, `just`
and `text` are the as the corrresponding arguments in [`pgptxt`](@ref) but
instead of drawing the string as routine [`pgptxt`](@ref) does, `pgqtxt`
returns the coordinates of the corners of a rectangle parallel to the string
baseline that just encloses the string. The four corners are in the order:
lower left, upper left, upper right, lower right (where left and right refer to
the first and last characters in the string).

If the string is blank or contains no drawable characters, all four coordinates
of the returned bonding box are assigned the starting point of the string,
`(x,y)`.

"""
function pgqtxt(x::Real, y::Real, angle::Real, just::Real, text::AbstractString)
    xbox = Vector{PGFloat}(undef, 4)
    ybox = Vector{PGFloat}(undef, 4)
    _pgqtxt(x, y, angle, just, text, xbox, ybox)
    return ((output(xbox[1]), output(xbox[2]), output(xbox[3]), output(xbox[4])),
            (output(ybox[1]), output(ybox[2]), output(ybox[3]), output(ybox[4])))
end

_pgqtxt(x, y, angle, just, text, xbox, ybox) =
    ccall((:cpgqtxt, pgplotlib), Cvoid,
          (PGFloat, PGFloat, PGFloat, PGFloat, PGString, Ptr{PGFloat}, Ptr{PGFloat}),
          x, y, angle, just, text, xbox, ybox)

"""

```julia
pgqvp(units) -> x1, x2, y1, y2
```

yields the current viewport settings in the required units: `:ndc` or 0 for
normalized device coordinates, `:in` or 1 for inches, `:mm` or 2 for
millimeters, `:pix` or 3 for pixels or device units.

The result is a 4-tuple of values `(x1,x2,y1,y2)` with `x1` the x-coordinate of
the bottom left corner of the viewport, `x2` the x-coordinate of the top right
corner of the viewport, `y1` the y-coordinate of the bottom left corner of the
viewport and `y2` the y-coordinate of the top right corner of the viewport.

"""
pgqvp(units::Union{Symbol,Integer}) =
    query(PGFloat, PGFloat, PGFloat, PGFloat, _pgqvp, pgunits(units))

_pgqvp(units, x1, x2, y1, y2) =
    ccall((:cpgqvp, pgplotlib), Cvoid,
          (PGInt, Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat}),
          units, x1, x2, y1, y2)

"""

```julia
pgqvsz(units) -> width, height
```

yields the dimensions of the view surface (the maximum plottable area) of the
currently selected graphics device, in the required units: `:ndc` or 0 for
normalized device coordinates, `:in` or 1 for inches, `:mm` or 2 for
millimeters, `:pix` or 3 for pixels or device units.  The result is a 2-tuple
of values: the width and the height of the view surface.

The size of the view surface is device-dependent and is established when the
graphics device is opened. On some devices, it can be changed by calling
[`pgpap`](@ref) before starting a new page with [`pgpage`](@ref). On some
devices, the size can be changed (e.g., by a workstation window manager)
outside PGPlot, and PGPlot detects the change when [`pgpage`](@ref) is
used. Call this routine after [`pgpage`](@ref) to find the current size.

Notes:

1. The width and the height of the view surface in normalized device
   coordinates are both always equal to 1.0.

2. When the device is divided into panels (see [`pgsubp`](@ref)), the view
   surface is a single panel.

"""
function pgqvsz(units::Union{Symbol,Integer})
    x1 = Ref{PGFloat}()
    x2 = Ref{PGFloat}()
    y1 = Ref{PGFloat}()
    y2 = Ref{PGFloat}()
    _pgqvsz(pgunits(units), x1, x2, y1, y2)
    @assert x1[] == 0
    @assert y1[] == 0
    return output(x2, y2)
end

_pgqvsz(units, x1, x2, y1, y2) =
    ccall((:cpgqvsz, pgplotlib), Cvoid,
          (PGInt, Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat}),
          units, x1, x2, y1, y2)

"""

```julia
pgqwin() -> x1, x2, y1, y2
```

yields the world coordinates of the current window boundaries: `x1` is the
X-coordinate of the bottom left corner of the window, `x2` is the X-coordinate
of the top right corner of the window, `y1` is the Y-coordinate of the bottom
left corner of the window and `y2` is the Y-coordinate of the top right corner
of the window.

"""
pgqwin() = query(PGFloat, PGFloat, PGFloat, PGFloat, _pgqwin)

_pgqwin(x1, x2, y1, y2) =
    ccall((:cpgqwin, pgplotlib), Cvoid,
          (Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat}, Ptr{PGFloat}),
          x1, x2, y1, y2)

"""

```julia
pgrect(x1, x2, y1, y2)
```

draws a rectangle, using fill-area attributes.  The rectangle has vertices at
`(x1,y1)`, `(x1,y2)`, `(x2,y2)`, and `(x2,y1)`.  This routine can be used
instead of [`pgpoly`](@ref) for the special case of drawing a rectangle aligned
with the coordinate axes; only two vertices need be specified instead of four.
On most devices, it is faster to use `pgrect` than [`pgpoly`](@ref) for drawing
rectangles.

"""
pgrect(x1::Real, x2::Real, y1::Real, y2::Real) =
    ccall((:cpgrect, pgplotlib), Cvoid, (PGFloat, PGFloat, PGFloat, PGFloat),
          x1, x2, y1, y2)

"""

```
pgrnd(x) -> xrnd, nsub
```

finds the smallest *round* number `xrnd` larger than `x` (in magnitude), a
*round* number being 1, 2 or 5 times a power of 10.  Second returned value,
`nsub` is a suitable number of subdivisions for subdividing the *nice*
number: 2 or 5.  If `x` is negative, `pgrnd(x) = -pgrnd(abs(x))`. *e.g.*
`pgrnd(8.7)` yields `xrnd = 10.0`, while `pgrnd(-0.4)` yields `xrnd =
-0.5`.  If `x` is zero, the value returned is zero.  This routine is used
by `pgbox` for choosing tick intervals.

See also [`pgrnge`](@ref).

"""
function pgrnd(x::Real)
    nsub = Ref{PGInt}()
    xrnd = ccall((:cpgrnd, pgplotlib), PGFloat,
                 (PGFloat, Ref{PGInt}), x, nsub)
    return output(xrnd, nsub)
end

"""

```julia
pgrnge(x1, x2) -> xlo, xhi
```

yields plotting limits `xlo` and `xhi` which encompass the data range `x1`
to `x2`.  That is such that `xlo < x1 < x2 < xhi` or `xlo > x1 > x2 > xhi`
or `xlo = x1 = x2 = xhi` depending on the ordering of `x1` and `x2`.

Input argument can also be a 2-tuple of values `(x1,x2)` or a vector of
values, say `x`, to consider `x1 = minimum(x)` and `x2 = maximum(x)`.

See also [`pgrnd`](@ref).

"""
pgrnge(x::AbstractVector{<:Real}) = pgrnge(extrema(x))

pgrnge(x::Tuple{Real, Real}) = pgrnge(x...)

function pgrnge(x1::Real, x2::Real)
    xlo = Ref{PGFloat}()
    xhi = Ref{PGFloat}()
    _pgrnge(x1, x2, xlo, xhi)
    return output(xlo, xhi)
end

_pgrnge(x1, x2, xlo, xhi) =
    ccall((:cpgrnge, pgplotlib), Cvoid,
          (PGFloat, PGFloat, Ptr{PGFloat}, Ptr{PGFloat}),
          x1, x2, xlo, xhi)

"""

```julia
pgsah(fs=1, angle=45.0, barb=0.3)
```

sets the style to be used for arrowheads drawn with routine [`pgarro`](@ref):

- `fs` is 1 for filled; 2 for outline.  Other values are treated as 2. Default
  1.

- `angle` is the acute angle of the arrow point, in degrees; angles in the
  range 20.0 to 90.0 give reasonable results. Default 45.0.

- `barb` is the fraction of the triangular arrow-head that is cut away from the
  back. 0.0 gives a triangular wedge arrow-head; 1.0 gives an open >. Values
  0.3 to 0.7 give reasonable results. Default 0.3.

"""
function pgsah(fs::Integer = PGInt(1),
               angle::Real = PGFloat(45),
               barb::Real = PGFloat(0.3))
    ccall((:cpgsah, pgplotlib), Cvoid,
          (PGInt, PGFloat, PGFloat),
          fs, angle, barb)
end

"""

```julia
pgsave()
```

This routine saves the current PGPlot attributes in a private storage
area.  They can be restored by calling [`pgunsa`](@ref) (unsave). Attributes
saved are: character font, character height, color index, fill-area style, line
style, line width, pen position, arrow-head style, hatching style, and clipping
state. Color representation is not saved.

Calls to `pgsave` and [`pgunsa`](@ref) should always be paired. Up to 20 copies
of the attributes may be saved. [`pgunsa`](@ref) always retrieves the
last-saved values (last-in first-out stack).

Note that when multiple devices are in use, [`pgunsa`](@ref) retrieves the
values saved by the last `pgsave` call, even if they were for a different
device.

"""
pgsave() = ccall((:cpgsave, pgplotlib), Cvoid, (), )

"""

```julia
pgunsa()
```

This routine restores the PGPlot attributes saved in the last call to
[`pgsave`](@ref).

"""
pgunsa() = ccall((:cpgunsa, pgplotlib), Cvoid, (), )

# Cvoid cpgscf(PGInt font);
"""

```julia
pgscf(font = 1)
```

sets the character font for subsequent text plotting. Four different fonts are
available:

- 1: (default) a simple single-stroke font ("normal" font);
- 2: roman font;
- 3: italic font;
- 4: script font.

This call determines which font is in effect at the beginning of each text
string. The font can be changed (temporarily) within a text string by using the
escape sequences `\\fn', '\\fr', '\\fi', and '\\fs' for fonts 1, 2, 3, and 4,
respectively.

"""
pgscf(font::Integer = PGInt(1)) =
    ccall((:cpgscf, pgplotlib), Cvoid, (PGInt,), font)

"""

```julia
pgsch(size = 1.0)
```

sets the character size attribute.  The size affects all text and graph markers
drawn later in the program. The default character size is 1.0, corresponding to
a character height about 1/40 the height of the view surface.  Changing the
character size also scales the length of tick marks drawn by [`pgbox`](@ref)
and terminals drawn by [`pgerrx`](@ref) and [`pgerry`](@ref).

"""
pgsch(size::Real = PGFloat(1.0)) =
    ccall((:cpgsch, pgplotlib), Cvoid, (PGFloat,), size)

"""

```julia
cpgsci(ci)
```

sets the Color Index for subsequent plotting, if the output device permits
this.  The default color index is 1, usually white on a black background for
video displays or black on a white background for printer plots.  The color
index is an integer in the range 0 to a device-dependent maximum.  Color index
0 corresponds to the background color; lines may be *erased* by overwriting
them with color index 0 (if the device permits this).

If the requested color index is not available on the selected device, color
index 1 will be substituted.

The assignment of colors to color indices can be changed with subroutine
[`pgscr`](@ref) (set color representation).  Color indices 0-15 have predefined
color representations (see the PGPlot manual), but these may be changed with
[`pgscr`](@ref).  Color indices above 15 have no predefined representations: if
these indices are used, [`pgscr](@ref) must be called to define the
representation.

"""
pgsci(ci::Integer) = ccall((:cpgsci, pgplotlib), Cvoid, (PGInt,), ci)

"""

```julia
pgscir(cmin, cmax)
```

sets the color index range to be used for producing images with
[`pggray`](@ref) or [`pgimag`](@ref).  If the range is not all within the range
supported by the device, a smaller range will be used.  `cmin` is the lowest
color index to use for images while `cmax` is the highest color index to use
for images.  The number of different colors available for images is `cmax -
cmin + 1`.

"""
pgscir(icilo::Integer, icihi::Integer) =
    ccall((:cpgscir, pgplotlib), Cvoid, (PGInt, PGInt), icilo, icihi)

"""

```julia
pgsclp(state)
```

Normally all PGPLOT primitives except text are *clipped* at the edge of the
viewport: parts of the primitives that lie outside the viewport are not drawn.
If clipping is disabled by calling this routine, primitives are visible
wherever they lie on the view surface.  Argument `state` is 0 to disable
clipping, or 1 to enable clipping.  The default (clipping enabled) is
appropriate for almost all applications.

"""
pgsclp(state::Integer) = ccall((:cpgsclp, pgplotlib), Cvoid, (PGInt,), state)


"""

```julia
pgscr(ci, cr, cg, cb)
```

sets color representation: i.e., defines the color to be associated with color
index `ci`.  Arguments `rr`, `cg` and `cb` are thre red, green, and blue
intensities, in the range 0.0 to 1.0.  Ignored for devices which do not support
variable color or intensity.  Color indices 0-15 have predefined color
representations (see the PGPLOT manual), but these may be changed with `pgscr`.
Color indices 16-maximum have no predefined representations: if these indices
are used, `pgscr` must be called to define the representation. On monochrome
output devices (e.g. VT125 terminals with monochrome monitors), the monochrome
intensity is computed from the specified Red, Green, Blue intensities as:

    0.30*R + 0.59*G + 0.11*B

as in US color television systems, NTSC encoding.  Note that most devices do
not have an infinite range of colors or monochrome intensities available; the
nearest available color is used.  Examples: for black, set `cr=cg=cb=0.0`; for
white, set `cr=cg=cb=1.0`; for medium gray, set `cr=cg=cb=0.5`; for medium
yellow, set `cr=cg=0.5`, `cb=0.0`.

"""
pgscr(ci::Integer, cr::Real, cg::Real, cb::Real) =
    ccall((:cpgscr, pgplotlib), Cvoid,
          (PGInt, PGFloat, PGFloat, PGFloat),
          ci, cr, cg, cb)

"""

```julia
pgscrl(dx, dy)
```

This routine moves the window in world-coordinate space while leaving the
viewport unchanged. On devices that have the capability, the pixels within the
viewport are scrolled horizontally, vertically or both in such a way that
graphics previously drawn in the window are shifted so that their world
coordinates are unchanged.

If the old window coordinate range was `(x1,x2,y1,y2)`, the new coordinate
range will be approximately `(x1+dx,x2+dx,y1+dy,y2+dy)`.  The size and scale of
the window are unchanged.

Thee window can only be shifted by a whole number of pixels (device
coordinates). If `dx` and `dy` do not correspond to integral numbers of pixels,
the shift will be slightly different from that requested. The new
window-coordinate range, and hence the exact amount of the shift, can be
determined by calling [`pgqwin`](@ref) after this routine.

Pixels that are moved out of the viewport by this operation are lost
completely; they cannot be recovered by scrolling back.  Pixels that are
*scrolled into* the viewport are filled with the background color (color index
0).

If the absolute value of `dx` is bigger than the width of the window, or the
aboslute value of `dy` is bigger than the height of the window, the effect will
be the same as zeroing all the pixels in the viewport.

Not all devices have the capability to support this routine.  It is only
available on some interactive devices that have discrete pixels. To determine
whether the current device has scroll capability, call [`pgqinf`](@ref).

"""
pgscrl(dx::Real, dy::Real) = ccall((:cpgscrl, pgplotlib), Cvoid,
                                   (PGFloat, PGFloat), dx, dy)

# void cpgscrn(int ci, const char *name, int *ier);

# Cvoid cpgscrn(PGInt ci, PGString name, Ptr{PGInt} ier);
"""

```julia
pgscrn(ci, name)
```

Set color representation: i.e., define the color to be associated with a color
index.  Ignored for devices which do not support variable color or intensity.
This is an alternative to routine [`pgscr`](@ref). The color representation is
defined by name instead of (R,G,B) components.

Color names are defined in an external file which is read the first time that
`pgscrn` is called. The name of the external file is found as follows:

1. if environment variable (logical name) `PGPLOT_RGB` is defined, its value is
   used as the file name;

2. otherwise, if environment variable `PGPLOT_DIR` is defined, a file `rgb.txt`
   in the directory named by this environment variable is used;

3. otherwise, file `rgb.txt` in the current directory is used.

If all of these fail to find a file, an error is reported and the routine does
nothing.

Each line of the file defines one color, with four blank- or tab-separated
fields per line. The first three fields are the R, G, B components, which are
integers in the range 0 (zero intensity) to 255 (maximum intensity). The fourth
field is the color name. The color name may include embedded blanks. Example:

> 255   0   0 red
> 255 105 180 hot pink
> 255 255 255 white
>   0   0   0 black

Arguments:

- `ci` is the color index to be defined, in the range 0-max.  If the color
 index greater than the device maximum is specified, the call is ignored. Color
 index 0 applies to the background color.

- `NAME` is the name of the color to be associated with this color index. This
  name must be in the external file. The names are not case-sensitive.  If the
  color is not listed in the file, the color representation is not changed.

  IER (output) : returns 0 if the routine was successful, 1 if an error
  occurred (either the external file could not be read, or the requested color
  was not defined in the file).


"""
function pgscrn(ci::Integer, name::AbstractString)
    ier = Ref{PGInt}()
    _pgscrn(ci, name, ier)
    ier[] == 0 ||
        error("file \"rgb.txt\" could not be read, or the requested color was not defined in this file")
end

_pgscrn(ci, name, ier) = ccall((:cpgscrn, pgplotlib), Cvoid,
                               (PGInt, PGString, Ptr{PGInt}),
                               ci, name, ier)

"""

```julia
pgsfs(fs)
```

Set the Fill-Area Style attribute for subsequent area-fill by [`pgpoly`](@ref),
[`pgrect`](@ref), or [`pgcirc`](@ref).  Four different styles are available:
solid (fill polygon with solid color of the current color-index), outline (draw
outline of polygon only, using current line attributes), hatched (shade
interior of polygon with parallel lines, using current line attributes), or
cross-hatched. The orientation and spacing of hatch lines can be specified with
routine [`pgshs`](@ref) (set hatch style).

Argument `fs` is the fill-area style to be used for subsequent plotting:

- `fs = 1` for solid (default);
- `fs = 2` for outline;
- `fs = 3` for hatched;
- `fs = 4` for cross-hatched.

Other values give an error message and are treated as 2.

"""
pgsfs(fs::Integer) = ccall((:cpgsfs, pgplotlib), Cvoid, (PGInt,), fs)

"""

```julia
pgshls(ci, ch, cl, cs)
```

Set color representation: i.e., define the color to be associated with a color
index.  This routine is equivalent to [`pgscr`](@ref), but the color is defined
in the Hue-Lightness-Saturation model instead of the Red-Green-Blue model. Hue
is represented by an angle in degrees, with red at 120, green at 240, and blue
at 0 (or 360). Lightness ranges from 0.0 to 1.0, with black at lightness 0.0
and white at lightness 1.0. Saturation ranges from 0.0 (gray) to 1.0 (pure
color). Hue is irrelevant when saturation is 0.0.

| Color        ||  H  |  L  |  S  ||  R     G     B
|:-------------||:---:|:---:|:---:||:---:|:---:|:---:|
| black        || any | 0.0 | 0.0 || 0.0 | 0.0 | 0.0 |
| white        || any | 1.0 | 0.0 || 1.0 | 1.0 | 1.0 |
| medium gray  || any | 0.5 | 0.0 || 0.5 | 0.5 | 0.5 |
| red          || 120 | 0.5 | 1.0 || 1.0 | 0.0 | 0.0 |
| yellow       || 180 | 0.5 | 1.0 || 1.0 | 1.0 | 0.0 |
| pink         || 120 | 0.7 | 0.8 || 0.94| 0.46| 0.46|

Reference: SIGGRAPH Status Report of the Graphic Standards Planning Committee,
Computer Graphics, Vol.13, No.3, Association for Computing Machinery, New York,
NY, 1979. See also: J. D. Foley et al, ``Computer Graphics: Principles and
Practice'', second edition, Addison-Wesley, 1990, section 13.3.5.

Argument:

- `ci` is the color index to be defined, in the range 0-max.  If the color
  index greater than the device maximum is specified, the call is
  ignored. Color index 0 applies to the background color.

- `ch` is the hue, in range 0.0 to 360.0.

- `cl` is the lightness, in range 0.0 to 1.0.

- `cs` is the saturation, in range 0.0 to 1.0.

"""
pgshls(ci::Integer, ch::Real, cl::Real, cs::Real) =
    ccall((:cpgshls, pgplotlib), Cvoid,
          (PGInt, PGFloat, PGFloat, PGFloat),
          ci, ch, cl, cs)

"""

```julia
pgshs(angle, sepn, phase)
```

Set the style to be used for hatching (fill area with fill-style 3).
The default style is ANGLE=45.0, SEPN=1.0, PHASE=0.0.

Arguments:
 ANGLE  (input)  : the angle the hatch lines make with the
                   horizontal, in degrees, increasing
                   counterclockwise (this is an angle on the
                   view surface, not in world-coordinate space).
 SEPN   (input)  : the spacing of the hatch lines. The unit spacing
                   is 1 percent of the smaller of the height or
                   width of the view surface. This should not be
                   zero.
 PHASE  (input)  : a real number between 0 and 1; the hatch lines
                   are displaced by this fraction of SEPN from a
                   fixed reference.  Adjacent regions hatched with the
                   same PHASE have contiguous hatch lines. To hatch
                   a region with alternating lines of two colors,
                   fill the area twice, with PHASE=0.0 for one color
                   and PHASE=0.5 for the other color.

"""
pgshs(angle::Real, sepn::Real, phase::Real) =
    ccall((:cpgshs, pgplotlib), Cvoid,
          (PGFloat, PGFloat, PGFloat),
          angle, sepn, phase)

"""

```julia
pgsitf(itf)
```

sets the Image Transfer Function for subsequent images drawn by
[`pgimag`](@ref), [`pggray`](@ref), or [`pgwedg`](@ref). The Image Transfer
Function is used to map array values into the available range of color indices
specified with routine [`pgscir`](@ref) or (for [`pggray`](@ref) on some
devices) into dot density.  Argument can be: `0` for a linear ITF, `1` for a
logarithmic ITF or `2` for a square-root ITF.

"""
pgsitf(itf::Integer) = ccall((:cpgsitf, pgplotlib), Cvoid, (PGInt,), itf)

"""

```julia
pgslct(id)
```

Select one of the open graphics devices and direct subsequent plotting to
it. The argument is the device identifier returned by [`pgopen`](@ref) when
the device was opened. If the supplied argument is not a valid identifier
of an open graphics device, a warning message is issued and the current
selection is unchanged.

"""
pgslct(id::Integer) = ccall((:cpgslct, pgplotlib), Cvoid,
                            (PGInt,), id)

"""

```julia
pgsls(ls)
```

Set the line style attribute for subsequent plotting. This
attribute affects line primitives only; it does not affect graph
markers, text, or area fill.
Five different line styles are available, with the following codes:
1 (full line), 2 (dashed), 3 (dot-dash-dot-dash), 4 (dotted),
5 (dash-dot-dot-dot). The default is 1 (normal full line).

Argument:
 LS     (input)  : the line-style code for subsequent plotting
                   (in range 1-5).

"""
pgsls(ls::Integer) = ccall((:cpgsls, pgplotlib), Cvoid, (PGInt,), ls)

"""

```julia
pgslw(lw)
```

Set the line-width attribute. This attribute affects lines, graph
markers, and text. The line width is specified in units of 1/200
(0.005) inch (about 0.13 mm) and must be an integer in the range
1-201. On some devices, thick lines are generated by tracing each
line with multiple strokes offset in the direction perpendicular to
the line.

Argument:
 LW     (input)  : width of line, in units of 0.005 inch (0.13 mm)
                   in range 1-201.

"""
pgslw(lw::Integer) = ccall((:cpgslw, pgplotlib), Cvoid, (PGInt,), lw)

"""

```julia
pgstbg(tbci)
```

Set the Text Background Color Index for subsequent text. By default
text does not obscure underlying graphics. If the text background
color index is positive, however, text is opaque: the bounding box
of the text is filled with the color specified by PGSTBG before
drawing the text characters in the current color index set by PGSCI.
Use color index 0 to erase underlying graphics before drawing text.

Argument:
 TBCI   (input)  : the color index to be used for the background
                   for subsequent text plotting:
                     TBCI < 0  => transparent (default)
                     TBCI >= 0 => text will be drawn on an opaque
                   background with color index TBCI.

"""
pgstbg(tbci::Integer) = ccall((:cpgstbg, pgplotlib), Cvoid, (PGInt,), tbci)

"""

```julia
pgsubp(nxsub, nysub)
```

PGPlot divides the physical surface of the plotting device (screen,
window, or sheet of paper) into NXSUB x NYSUB `panels'. When the
view surface is sub-divided in this way, PGPAGE moves to the next
panel, not the next physical page. The initial subdivision of the
view surface is set in the call to PGBEG. When PGSUBP is called,
it forces the next call to PGPAGE to start a new physical page,
subdivided in the manner indicated. No plotting should be done
between a call of PGSUBP and a call of PGPAGE (or PGENV, which calls
PGPAGE).

If NXSUB > 0, PGPLOT uses the panels in row order; if <0,
PGPLOT uses them in column order, e.g.,

 NXSUB=3, NYSUB=2            NXSUB=-3, NYSUB=2

+-----+-----+-----+         +-----+-----+-----+
|  1  |  2  |  3  |         |  1  |  3  |  5  |
+-----+-----+-----+         +-----+-----+-----+
|  4  |  5  |  6  |         |  2  |  4  |  6  |
+-----+-----+-----+         +-----+-----+-----+

PGPLOT advances from one panels to the next when PGPAGE is called,
clearing the screen or starting a new page when the last panel has
been used. It is also possible to jump from one panel to another
in random order by calling PGPANL.

Arguments:
 NXSUB  (input)  : the number of subdivisions of the view surface in
                   X (>0 or <0).
 NYSUB  (input)  : the number of subdivisions of the view surface in
                   Y (>0).

"""
pgsubp(nxsub::Integer, nysub::Integer) =
    ccall((:cpgsubp, pgplotlib), Cvoid, (PGInt, PGInt), nxsub, nysub)


"""

```julia
pgsvp(xleft, xright, ybot, ytop)
```

Change the size and position of the viewport, specifying
the viewport in normalized device coordinates.  Normalized
device coordinates run from 0 to 1 in each dimension. The
viewport is the rectangle on the view surface "through"
which one views the graph.  All the PG routines which plot lines
etc. plot them within the viewport, and lines are truncated at
the edge of the viewport (except for axes, labels etc drawn with
PGBOX or PGLAB).  The region of world space (the coordinate
space of the graph) which is visible through the viewport is
specified by a call to PGSWIN.  It is legal to request a
viewport larger than the view surface; only the part which
appears on the view surface will be plotted.

Arguments:
 XLEFT  (input)  : x-coordinate of left hand edge of viewport, in NDC.
 XRIGHT (input)  : x-coordinate of right hand edge of viewport,
                   in NDC.
 YBOT   (input)  : y-coordinate of bottom edge of viewport, in NDC.
 YTOP   (input)  : y-coordinate of top  edge of viewport, in NDC.

"""
pgsvp(xleft::Real, xright::Real, ybot::Real, ytop::Real) =
    ccall((:cpgsvp, pgplotlib), Cvoid,
          (PGFloat, PGFloat, PGFloat, PGFloat),
          xleft, xright, ybot, ytop)

"""

```julia
pgswin(x1, x2, y1, y2)
```

Change the window in world coordinate space that is to be mapped on
to the viewport.  Usually PGSWIN is called automatically by PGENV,
but it may be called directly by the user.

Arguments:
 X1     (input)  : the x-coordinate of the bottom left corner
                   of the viewport.
 X2     (input)  : the x-coordinate of the top right corner
                   of the viewport (note X2 may be less than X1).
 Y1     (input)  : the y-coordinate of the bottom left corner
                   of the viewport.
 Y2     (input)  : the y-coordinate of the top right corner
                   of the viewport (note Y2 may be less than Y1).

"""
pgswin(x1::Real, x2::Real, y1::Real, y2::Real) =
    ccall((:cpgswin, pgplotlib), Cvoid,
          (PGFloat, PGFloat, PGFloat, PGFloat),
          x1, x2, y1, y2)

"""

```julia
pgtbox(xopt, xtick, nxsub, yopt, ytick, nysub)
```

Draw a box and optionally label one or both axes with (DD) HH MM SS
style numeric labels (useful for time or RA - DEC plots).   If this
style of labelling is desired, then PGSWIN should have been called
previously with the extrema in SECONDS of time.

In the seconds field, you can have at most 3 places after the decimal
point, so that 1 ms is the smallest time interval you can time label.

Large numbers are coped with by fields of 6 characters long.  Thus
you could have times with days or hours as big as 999999.  However,
in practice, you might have trouble with labels overwriting  themselves
with such large numbers unless you a) use a small time INTERVAL,
b) use a small character size or c) choose your own sparse ticks in
the call to PGTBOX.

PGTBOX will attempt, when choosing its own ticks, not to overwrite
the labels, but this algorithm is not very bright and may fail.

Note that small intervals but large absolute times such as
TMIN = 200000.0 s and TMAX=200000.1 s will cause the algorithm
to fail.  This is inherent in PGPLOT's use of single precision
and cannot be avoided.  In such cases, you should use relative
times if possible.

PGTBOX's labelling philosophy is that the left-most or bottom tick of
the axis contains a full label.  Thereafter, only changing fields are
labelled.  Negative fields are given a '-' label, positive fields
have none.   Axes that have the DD (or HH if the day field is not
used) field on each major tick carry the sign on each field.  If the
axis crosses zero, the zero tick will carry a full label and sign.

This labelling style can cause a little confusion with some special
cases, but as long as you know its philosophy, the truth can be divined.
Consider an axis with TMIN=20s, TMAX=-20s.   The labels will look like

       +----------+----------+----------+----------+
    0h0m20s      10s      -0h0m0s      10s        20s

Knowing that the left field always has a full label and that
positive fields are unsigned, informs that time is decreasing
from left to right, not vice versa.   This can become very
unclear if you have used the 'F' option, but that is your problem !

Exceptions to this labelling philosophy are when the finest time
increment being displayed is hours (with option 'Y') or days.
Then all fields carry a label.  For example,

       +----------+----------+----------+----------+
     -10h        -8h        -6h        -4h        -2h


PGTBOX can be used in place of PGBOX; it calls PGBOX and only invokes
time labelling if requested. Other options are passed intact to PGBOX.

Inputs:
 XOPT   :  X-options for PGTBOX.  Same as for PGBOX plus

            'Z' for (DD) HH MM SS.S time labelling
            'Y' means don't include the day field so that labels
                are HH MM SS.S rather than DD HH MM SS.S   The hours
                will accumulate beyond 24 if necessary in this case.
            'X' label the HH field as modulo 24.  Thus, a label
                such as 25h 10m would come out as 1h 10m
            'H' means superscript numbers with d, h, m, & s  symbols
            'D' means superscript numbers with    o, ', & '' symbols
            'F' causes the first label (left- or bottom-most) to
                be omitted. Useful for sub-panels that abut each other.
                Care is needed because first label carries sign as well.
            'O' means omit leading zeros in numbers < 10
                E.g.  3h 3m 1.2s rather than 03h 03m 01.2s  Useful
                to help save space on X-axes. The day field does not
                use this facility.

 YOPT   :  Y-options for PGTBOX.  See above.
 XTICK  :  X-axis major tick increment.  0.0 for default.
 YTICK  :  Y-axis major tick increment.  0.0 for default.
           If the 'Z' option is used then XTICK and/or YTICK must
           be in seconds.
 NXSUB  :  Number of intervals for minor ticks on X-axis. 0 for default
 NYSUB  :  Number of intervals for minor ticks on Y-axis. 0 for default

 The regular XOPT and YOPT axis options for PGBOX are

 A : draw Axis (X axis is horizontal line Y=0, Y axis is vertical
     line X=0).
 B : draw bottom (X) or left (Y) edge of frame.
 C : draw top (X) or right (Y) edge of frame.
 G : draw Grid of vertical (X) or horizontal (Y) lines.
 I : Invert the tick marks; ie draw them outside the viewport
     instead of inside.
 L : label axis Logarithmically (see below).
 N : write Numeric labels in the conventional location below the
     viewport (X) or to the left of the viewport (Y).
 P : extend ("Project") major tick marks outside the box (ignored if
     option I is specified).
 M : write numeric labels in the unconventional location above the
     viewport (X) or to the right of the viewport (Y).
 T : draw major Tick marks at the major coordinate interval.
 S : draw minor tick marks (Subticks).
 V : orient numeric labels Vertically. This is only applicable to Y.
     The default is to write Y-labels parallel to the axis.
 1 : force decimal labelling, instead of automatic choice (see PGNUMB).
 2 : force exponential labelling, instead of automatic.

     The default is to write Y-labels parallel to the axis


       ******************        EXCEPTIONS       *******************

       Note that
         1) PGBOX option 'L' (log labels) is ignored with option 'Z'
         2) The 'O' option will be ignored for the 'V' option as it
            makes it impossible to align the labels nicely
         3) Option 'Y' is forced with option 'D'

       ***************************************************************

"""
function pgtbox(xopt::AbstractString, xtick::Real, nxsub::Integer,
                yopt::AbstractString, ytick::Real, nysub::Integer)
    ccall((:cpgtbox, pgplotlib), Cvoid,
          (PGString, PGFloat, PGInt, PGString, PGFloat, PGInt),
          xopt, xtick, nxsub, yopt, ytick, nysub)
end

"""

```julia
pgtext(x, y, text)
```

Write text. The bottom left corner of the first character is placed
at the specified position, and the text is written horizontally.
This is a simplified interface to the primitive routine PGPTXT.
For non-horizontal text, use PGPTXT.

Arguments:
 X      (input)  : world x-coordinate of start of string.
 Y      (input)  : world y-coordinate of start of string.
 TEXT   (input)  : the character string to be plotted.

"""
pgtext(x::Real, y::Real, text::AbstractString) =
    ccall((:cpgtext, pgplotlib), Cvoid,
          (PGFloat, PGFloat, PGString),
          x, y, text)

"""

```julia

```

Draw and label single tick mark on a graph axis. The tick mark is
a short line perpendicular to the direction of the axis (which is not
drawn by this routine). The optional text label is drawn with its
baseline parallel to the axis and reading in the same direction as
the axis (from point 1 to point 2). Current line and text attributes
are used.

Arguments:
 X1, Y1 (input)  : world coordinates of one endpoint of the axis.
 X2, Y2 (input)  : world coordinates of the other endpoint of the axis.
 V      (input)  : draw the tick mark at fraction V (0<=V<=1) along
                   the line from (X1,Y1) to (X2,Y2).
 TIKL   (input)  : length of tick mark drawn to left of axis
                   (as seen looking from first endpoint to second), in
                   units of the character height.
 TIKR   (input)  : length of major tick marks drawn to right of axis,
                   in units of the character height.
 DISP   (input)  : displacement of label text to
                   right of axis, in units of the character height.
 ORIENT (input)  : orientation of label text, in degrees; angle between
                   baseline of text and direction of axis (0-360°).
 STR    (input)  : text of label (may be blank).

"""
function pgtick(x1::Real, y1::Real, x2::Real, y2::Real,
                v::Real, tikl::Real, tikr::Real, disp::Real,
                orient::Real, str::AbstractString)
    ccall((:cpgtick, pgplotlib), Cvoid,
          (PGFloat, PGFloat, PGFloat, PGFloat, PGFloat, PGFloat,
           PGFloat, PGFloat, PGFloat, PGString),
          x1, y1, x2, y2, v, tikl, tikr, disp, orient, str)
end

"""

```julia
pgupdt()
```

Update the graphics display: flush any pending commands to the
output device. This routine empties the buffer created by PGBBUF,
but it does not alter the PGBBUF/PGEBUF counter. The routine should
be called when it is essential that the display be completely up to
date (before interaction with the user, for example) but it is not
known if output is being buffered.

"""
pgupdt() = ccall((:cpgupdt, pgplotlib), Cvoid, (), )

# void cpgvect(const float *a, const float *b, int idim, int jdim, int i1,
#              int i2, int j1, int j2, float c, int nc, const float *tr,
#              float blank);
"""

```julia
pgvect(A, B, [i1, i2, j1, j2,] c, nc, tr=[0,1,0,0,0,1], blank=NaN)
```

Draw a vector map of two arrays.  This routine is similar to
PGCONB in that array elements that have the "magic value" defined by
the argument BLANK are ignored, making gaps in the vector map.  The
routine may be useful for data measured on most but not all of the
points of a grid. Vectors are displayed as arrows; the style of the
arrowhead can be set with routine PGSAH, and the the size of the
arrowhead is determined by the current character size, set by PGSCH.

Arguments:
 A      (input)  : horizontal component data array.
 B      (input)  : vertical component data array.
 IDIM   (input)  : first dimension of A and B.
 JDIM   (input)  : second dimension of A and B.
 I1,I2  (input)  : range of first index to be mapped (inclusive).
 J1,J2  (input)  : range of second index to be mapped (inclusive).
 C      (input)  : scale factor for vector lengths, if 0.0, C will be
                   set so that the longest vector is equal to the
                   smaller of TR(2)+TR(3) and TR(5)+TR(6).
 NC     (input)  : vector positioning code.
                   <0 vector head positioned on coordinates
                   >0 vector base positioned on coordinates
                   =0 vector centered on the coordinates
 TR     (input)  : array defining a transformation between the I,J
                   grid of the array and the world coordinates. The
                   world coordinates of the array point A(I,J) are
                   given by:
                     X = TR(1) + TR(2)*I + TR(3)*J
                     Y = TR(4) + TR(5)*I + TR(6)*J
                   Usually TR(3) and TR(5) are zero - unless the
                   coordinate transformation involves a rotation
                   or shear.
 BLANK   (input) : elements of arrays A or B that are exactly equal to
                   this value are ignored (blanked).

"""
function pgvect(a::PGVector{PGFloat}, b::PGVector{PGFloat},
                idim::Integer, jdim::Integer,
                i1::Integer, i2::Integer,
                j1::Integer, j2::Integer,
                c::Real, nc::Integer,
                tr::PGVector{PGFloat} = DEFAULT_XFORM,
                blank::Real = PGFloat(NaN))
    error("FIXME: not yet implemented")
end

_pgvect(a, b, idim, jdim, i1, i2, j1, j2, c, nc, tr, blank) =
    ccall((:cpgvect, pgplotlib), Cvoid,
          (Ptr{PGFloat}, Ptr{PGFloat}, PGInt, PGInt, PGInt, PGInt,
           PGInt, PGInt, PGFloat, PGInt, Ptr{PGFloat}, PGFloat),
          a, b, idim, jdim, i1, i2, j1, j2, c, nc, tr, blank)

"""

```julia
pgvsiz(xleft, xright, ybot, ytop)
```

Change the size and position of the viewport, specifying the viewport in
physical device coordinates (inches).  The viewport is the rectangle on the
view surface "through" which one views the graph.  All the PG routines which
plot lines etc. plot them within the viewport, and lines are truncated at the
edge of the viewport (except for axes, labels etc drawn with [`pgbox`](@ref) or
[`pglab`](@ref)).  The region of world space (the coordinate space of the
graph) which is visible through the viewport is specified by a call to
[`pgswin`](@ref).  It is legal to request a viewport larger than the view
surface; only the part which appears on the view surface will be plotted.

Arguments:
 XLEFT  (input)  : x-coordinate of left hand edge of viewport, in
                   inches from left edge of view surface.
 XRIGHT (input)  : x-coordinate of right hand edge of viewport, in
                   inches from left edge of view surface.
 YBOT   (input)  : y-coordinate of bottom edge of viewport, in
                   inches from bottom of view surface.
 YTOP   (input)  : y-coordinate of top  edge of viewport, in inches
                   from bottom of view surface.

"""
pgvsiz(xleft::Real, xright::Real, ybot::Real, ytop::Real) =
    ccall((:cpgvsiz, pgplotlib), Cvoid,
          (PGFloat, PGFloat, PGFloat, PGFloat),
          xleft, xright, ybot, ytop)

"""

```julia
pgvstd()
```

Define the viewport to be the standard viewport.  The standard viewport is the
full area of the view surface (or panel), less a margin of 4 character heights
all round for labelling.  It thus depends on the current character size, set by
[`pgsch`](@ref).

"""
pgvstd() = ccall((:cpgvstd, pgplotlib), Cvoid, (), )

"""

```julia
pgwedg(side, disp, width, fg, bg, label)
```

Plot an annotated grey-scale or color wedge parallel to a given axis of the the
current viewport. This routine is designed to provide a brightness/color scale
for an image drawn with [`pgimag`](@ref) or [`pggray`](@ref).  The wedge will
be drawn with the transfer function set by [`pgsitf`](@ref) and using the color
index range set by [`pgscir`](@ref).

Arguments:

- `side` is The first character must be one of the characters `'B'`, `'L'`,
  `'T'`, or `'R'` signifying the Bottom, Left, Top, or Right edge of the
  viewport.  The second character should be `'I'` to use [`pgimag`](@ref) to
  draw the wedge, or `'G'` to use [`pggray`](@ref).

- `disp` is the displacement of the wedge from the specified edge of the
  viewport, measured outwards from the viewport in units of the character
  height. Use a negative value to write inside the viewport, a positive value
  to write outside.

- `width` is The total width of the wedge including annotation, in units of the
  character height.

- `fg` is The value which is to appear with shade 1 (*foreground*). Use the
  values of `fg` and `bg` that were supplied to [`pggray`](@ref) or
  [`pgimag`](@ref).

- `bg` is the value which is to appear with shade 0 (*background*).

- `label` is Optional units label. If no label is required use `" "`.

"""
function pgwedg(side::AbstractString, disp::Real, width::Real,
                fg::Real, bg::Real, label::AbstractString = " ")
    ccall((:cpgwedg, pgplotlib), Cvoid,
          (PGString, PGFloat, PGFloat, PGFloat, PGFloat, PGString),
          side, disp, width, fg, bg, label)
end

"""

```julia
pgwnad(x1, x2, y1, y2)
```

Change the window in world coordinate space that is to be mapped on to the
viewport, and simultaneously adjust the viewport so that the world-coordinate
scales are equal in `x` and `y`. The new viewport is the largest one that can
fit within the previously set viewport while retaining the required aspect
ratio.

Arguments are: `x1` the x-coordinate of the bottom left corner of the viewport,
`x2` the x-coordinate of the top right corner of the viewport (note `x2` may be
less than `x1`), `y1` the y-coordinate of the bottom left corner of the
viewport and `y2` the y-coordinate of the top right corner of the viewport
(note `y2` may be less than `y1`)

"""
pgwnad(x1::Real, x2::Real, y1::Real, y2::Real) =
    ccall((:cpgwnad, pgplotlib), Cvoid,
          (PGFloat, PGFloat, PGFloat, PGFloat),
          x1, x2, y1, y2)

end # module
