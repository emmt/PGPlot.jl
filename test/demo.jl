module PGPlotDemo

using Printf

using PGPlot.Bindings
using PGPlot.Bindings: DEFAULT_XFORM

const DEFAULT_DEVICE = "/XSERVE"

function setupdevice()
    if pgqinf("HARDCOPY") == "YES"
        pgask(false)    # do not wait for user input
        pgscr(0, 0,0,0) # set background color to white
        pgscr(1, 1,1,1) # set foreground color to black
    else
        pgask(true)     # do wait for user input
        pgscr(0, 0,0,0) # set background color to black
        pgscr(1, 1,1,1) # set foreground color to white
    end
end

# Make sure adevice is currently open.
function opendevice()
    if pgqid() < 1
        pgopen(DEFAULT_DEVICE)
        pgask(false)    # do not wait for user input
        pgscr(0, 0,0,0) # set background color to black
        pgscr(1, 1,1,1) # set foreground color to white
    end
end

#-----------------------------------------------------------------------
# Demonstration program for PGPlot. The main program opens the output
# device and calls a series of subroutines, one for each sample plot.
#-----------------------------------------------------------------------
function pgdem1(dev::AbstractString = DEFAULT_DEVICE)
    # Call PGOPEN to initiate PGPlot and open the output device; PGOPEN
    # will prompt the user to supply the device name and type. Always
    # check the return code from PGOPEN.
    id = pgopen(dev)
    setupdevice()

    # Print information about device.
    pgex0()

    # Call the demonstration subroutines (4,5 are put on one page).
    pgex1()
    pgex2()
    pgex3()
    pgsubp(2,1)
    pgex4()
    pgex5()
    pgsubp(1,1)
    pgex6()
    pgex7()
    pgex8()
    pgex9()
    pgex10()
    pgex11()
    pgex12()
    pgex13()
    pgex14()
    pgex15()

    # Finally, call PGCLOS to terminate things properly.
    pgclos()
end

#-----------------------------------------------------------------------
# This example tests PGQINF and displays the information returned on
# the standard output.
# -----------------------------------------------------------------------
function pgex0()
    opendevice() # make sure a plotting device is open

    # Information available from PGQINF:
    println("Device information (from PGQINF):")
    for item in ("version", "state", "user", "now", "device", "file",
                 "type", "dev/type", "hardcopy", "terminal", "cursor")
        println("  ", item, " = ", pgqinf(item))
    end

    # Get view surface dimensions (in inches):
    w, h = pgqvsz(:in)
    println()
    @printf("Plot dimensions: %.2f × %.2f inches\n", w, h)
    @printf("                 %.2f × %.2f millimeters\n", 25.4w, 25.4h)
end

#-----------------------------------------------------------------------
# This example illustrates the use of PGENV, PGLAB, PGPT, PGLINE.
#-----------------------------------------------------------------------
function pgex1()
    opendevice() # make sure a plotting device is open

    # Compute the function at N points.
    n = 60
    dx = PGFloat(0.1)
    xr = [i*dx for i in 1:n]
    yr = xr.^2

    # Call PGENV to specify the range of the axes and to draw a box, and
    # PGLAB to label it.  The x-axis runs from 0 to 10, and y from 0 to 20.
    pgenv(0.0, 10.0, 0.0, 20.0, false, 1)
    pglab("(x)", "(y)", "PGPlot Example 1:  y = x\\u2\\d")

    # Mark five points using symbol number 9.
    sub = 10:10:n
    pgpt(xr[sub], yr[sub], 9)

    # Use PGLINE to draw the function.
    pgline(xr, yr)
    #xt = 0:2:10
    #pgline(xt, sqrt.(xt))
end

#-----------------------------------------------------------------------
# Repeat the process for another graph. This one is a graph of the
# sinc (sin x over x) function.
#-----------------------------------------------------------------------
function pgex2()
    opendevice() # make sure a plotting device is open
    N = 100
    xr = PGFloat[(i-20)/6 for i in 1:N]
    yr = [xr[i] != 0 ? sin(xr[i])/xr[i] : one(PGFloat) for i in 1:N]
    pgenv(-2.0, 10.0, -0.4, 1.2, false, 1)
    pglab("(x)", "sin(x)/x", "PGPlot Example 2: Sinc Function")
    pgline(xr, yr)
end

#----------------------------------------------------------------------
# This example illustrates the use of PGBOX and attribute routines to
# mix colors and line-styles.
#----------------------------------------------------------------------
function pgex3()
    opendevice() # make sure a plotting device is open
    N = 360

    # Call PGENV to initialize the viewport and window; the AXIS argument
    # is -2, so no frame or labels will be drawn.
    pgenv(0.0, 720.0, -2.0, 2.0, false, -2)
    pgsave()

    # Set the color index for the axes and grid (index 5 = cyan).  Call
    # PGBOX to draw first a grid at low brightness, and then a frame and
    # axes at full brightness. Note that as the x-axis is to represent an
    # angle in degrees, we request an explicit tick interval of 90 deg with
    # subdivisions at 30 deg, as multiples of 3 are a more natural division
    # than the default.
    pgsci(14)
    pgbox("G", 30.0, 0, "G", 0.2, 0)
    pgsci(5)
    pgbox("ABCTSN", 90.0, 3, "ABCTSNV", 0.0, 0)

    # Call PGLAB to label the graph in a different color (3=green).
    pgsci(3)
    pglab("x[degrees]", "f(x)", "PGPlot Example 3")

    # Compute the function to be plotted: a trig function of an angle in
    # degrees, computed every 2 degrees.
    xr = [2i for i in 1:360]
    yr = [(arg = x*(π/180);
           sin(arg) + 0.5*cos(2.0*arg) + 0.5*sin(1.5*arg + π/3)) for x in xr]

    # Change the color (6=magenta), line-style (2=dashed), and line width
    # and draw the function.
    pgsci(6)
    pgsls(2)
    pgslw(3)
    pgline(xr, yr)

    # Restore attributes to defaults.
    pgunsa()
end

#-----------------------------------------------------------------------
# Demonstration program for PGPlot: draw histograms.
#-----------------------------------------------------------------------
function pgex4()
    opendevice() # make sure a plotting device is open

    # Obtain 1000 samples from a normal distribution.
    data = randn(PGFloat, 1000)

    # Draw a histogram of these values.
    pgsave()
    pghist(data, -3.1, 3.1, 31, 0)

    # Samples from another normal distribution.
    data = 1.0 .+ 0.5.*randn(200)

    # Draw another histogram (filled) on same axes.
    pgsci(15)
    pghist(data, -3.1, 3.1, 31, 3)
    pgsci(0)
    pghist(data, -3.1, 3.1, 31, 1)
    pgsci(1)

    # Redraw the box which may have been clobbered by the histogram.
    pgbox("BST", 0.0, 0, "", 0.0, 0)

    # Label the plot.
    pglab("Variate", "", "PGPlot Example 4:  Histograms (Gaussian)")

    # Superimpose the theoretical distribution.
    x = [0.01*(i - 1) - 3.1 for i in 1:621]
    c = 0.2*1000/sqrt(2π)
    y = [c*exp(-0.5*xi^2) for xi in x]
    pgline(x, y)
    pgunsa()
end

#----------------------------------------------------------------------
# Demonstration program for the PGPlot plotting package.  This example
# illustrates how to draw a log-log plot.
# PGPlot subroutines demonstrated:
#    PGENV, PGERRY, PGLAB, PGLINE, PGPT, PGSCI.
#----------------------------------------------------------------------
function pgex5()
    opendevice() # make sure a plotting device is open
    red = 2
    green = 3
    cyan = 5
    freq = [26.0, 38.0, 80.0, 160.0, 178.0, 318.0, 365.0, 408.0, 750.0,
            1400.0, 2695.0, 2700.0, 5000.0, 10695.0, 14900.0]
    flux = [38.0, 66.4, 89.0, 69.8, 55.9, 37.4, 46.8, 42.4, 27.0, 15.8,
            9.09, 9.17, 5.35, 2.56, 1.73]
    err = [6.0, 6.0, 13.0, 9.1, 2.9, 1.4, 2.7, 3.0, 0.34, 0.8, 0.2,
           0.46, 0.15, 0.08, 0.01]

    # Call PGENV to initialize the viewport and window; the AXIS argument
    # is 30 so both axes will be logarithmic. The X-axis (frequency) runs
    # from 0.01 to 100 GHz, the Y-axis (flux density) runs from 0.3 to 300
    # Jy. Note that it is necessary to specify the logarithms of these
    # quantities in the call to PGENV. We request equal scales in x and y
    # so that slopes will be correct.  Use PGLAB to label the graph.
    pgsave()
    pgsci(cyan)
    pgenv(-2.0, 2.0, -0.5, 2.5, true, 30)
    pglab("Frequency, \\gn (GHz)", "Flux Density, S\\d\\gn\\u (Jy)",
          "PGPlot Example 5:  Log-Log plot")

    # Draw a fit to the spectrum (don't ask how this was chosen). This
    # curve is drawn before the data points, so that the data will write
    # over the curve, rather than vice versa.
    x = [1.3 + i*0.03 for i in 1:100]
    xp = [xi - 3.0 for xi in x]
    yp = [5.18 - 1.15*xi - 7.72*exp(-xi) for xi in x]
    pgline(xp, yp)

    # Plot the measured flux densities: here the data are installed with a
    # DATA statement; in a more general program, they might be read from a
    # file. We first have to take logarithms (the -3.0 converts MHz to GHz).
    xp = log10.(freq) .- 3
    yp = log10.(flux)
    pgsci(green)
    pgpt(xp, yp, 17)

    # Draw +/- 2 sigma error bars: take logs of both limits.
    yhi = log10.(flux .+ 2err)
    ylo = log10.(flux .- 2err)
    pgerry(xp, ylo, yhi, 1.0)
    pgunsa()
end

#----------------------------------------------------------------------
# Demonstration program for the PGPlot plotting package.  This example
# illustrates the use of PGPOLY, PGCIRC, and PGRECT using SOLID,
# OUTLINE, HATCHED, and CROSS-HATCHED fill-area attributes.
#----------------------------------------------------------------------
function pgex6()
    opendevice() # make sure a plotting device is open
    n1 = (3, 4, 5, 5, 6, 8)
    n2 = (1, 1, 1, 2, 1, 3)
    lab = ("Fill style 1 (solid)",
           "Fill style 2 (outline)",
           "Fill style 3 (hatched)",
           "Fill style 4 (cross-hatched)")

    # Initialize the viewport and window.
    pgbbuf()
    pgsave()
    pgpage()
    pgsvp(0.0, 1.0, 0.0, 1.0)
    pgwnad(0.0, 10.0, 0.0, 10.0)

    # Label the graph.
    pgsci(1)
    pgmtxt("T", -2.0, 0.5, 0.5,
           "PGPlot fill area: routines PGPOLY, PGCIRC, PGRECT")

    # Draw assorted polygons.
    for k in 1:4
        pgsci(1)
        y0 = 10.0 - 2.0*k
        pgtext(0.2, y0 + 0.6, lab[k])
        pgsfs(k)
        for i in 1:length(n1)
            c = 2π*n2[i]/n1[i]
            pgsci(i)
            angle = [c*(j - 1) for j in 1:n1[i]]
            x = [i + 0.5*cos(a) for a in angle]
            y = [y0 + 0.5*sin(a) for a in angle]
            pgpoly(x, y)
        end
        pgsci(7)
        pgcirc(7.0, y0, 0.5)
        pgsci(8)
        pgrect(7.8, 9.5, y0-0.5, y0 + 0.5)
    end
    pgunsa()
    pgebuf()
end

#-----------------------------------------------------------------------
# A plot with a large number of symbols; plus test of PGERR1.
#-----------------------------------------------------------------------
function pgex7()
    opendevice() # make sure a plotting device is open

    # Window and axes.
    pgbbuf()
    pgsave()
    pgsci(1)
    pgenv(0.0, 5.0, -0.3, 0.6, false, 1)
    pglab("\\fix", "\\fiy", "PGPlot Example 7: scatter plot")

    # Random data points.
    xs = [5*rand() for i in 1:300]
    ys = [x*exp(-x) + 0.05*randn() for x in xs]
    pgsci(3)
    pgpt(xs, ys, 3)
    sub = 101:200
    pgpt(xs[sub], ys[sub], 17)
    sub = 201:300
    pgpt(xs[sub], ys[sub], 21)

    # Curve defining parent distribution.
    xr = [0.05*(i - 1) for i in 1:101]
    yr = [x*exp(-x) for x in xr]
    pgsci(2)
    pgline(xr, yr)

    # Test of PGERR1/PGPT1.
    xp = xs[101]
    yp = ys[101]
    xsig = 0.2
    ysig = 0.1
    pgsci(5)
    pgsch(3.0)
    pgerr1(5, xp, yp, xsig, 1.0)
    pgerr1(6, xp, yp, ysig, 1.0)
    pgpt(xp, yp, 21)

    # Restore defaults.
    pgunsa()
    pgebuf()
end

#-----------------------------------------------------------------------
# Demonstration program for PGPlot. This program shows some of the
# possibilities for overlapping windows and viewports.
# T. J. Pearson  1986 Nov 28
#-----------------------------------------------------------------------
function pgex8()
    opendevice() # make sure a plotting device is open

    # Color indices:
    black   = 0
    white   = 1
    red     = 2
    green   = 3
    blue    = 4
    cyan    = 5
    magenta = 6
    yellow  = 7

    # Line style:
    full   = 1
    dashed = 2
    dotdsh = 3
    dotted = 4
    fancy  = 5

    # Character font:
    normal = 1
    roman  = 2
    italic = 3
    script = 4

    # Fill-area style:
    solid  = 1
    hollow = 2

    pgpage()
    pgbbuf()
    pgsave()

    # Define the Viewport
    pgsvp(0.1, 0.6, 0.1, 0.6)

    # Define the Window
    pgswin(0.0, 630.0, -2.0, 2.0)

    # Draw a box
    pgsci(cyan)
    pgbox("ABCTS", 90.0, 3, "ABCTSV", 0.0, 0)

    # Draw labels
    pgsci(red)
    pgbox("N",90.0, 3, "VN", 0.0, 0)

    # Draw SIN line
    xr = [2.0*i for i in 1:360]
    yr = [sin(x/57.29577951) for x in xr]

    pgsci(magenta)
    pgsls(dashed)
    pgline(xr, yr)

    # Draw COS line by redefining the window
    pgswin(90.0, 720.0, -2.0, 2.0)
    pgsci(yellow)
    pgsls(dotted)
    pgline(xr, yr)
    pgsls(full)

    # Re-Define the Viewport
    pgsvp(0.45, 0.85, 0.45, 0.85)

    # Define the Window, and erase it
    pgswin(0.0, 180.0, -2.0, 2.0)
    pgsci(0)
    pgrect(0.0, 180.0, -2.0, 2.0)

    # Draw a box
    pgsci(blue)
    pgbox("ABCTSM", 60.0, 3, "VABCTSM", 1.0, 2)

    # Draw SIN line
    pgsci(white)
    pgsls(dashed)
    pgline(xr,yr)

    pgunsa()
    pgebuf()
end

#----------------------------------------------------------------------
# Demonstration program for the PGPlot plotting package.  This example
# illustrates curve drawing with PGFUNT; the parametric curve drawn is
# a simple Lissajous figure.
# T. J. Pearson  1983 Oct 5
#----------------------------------------------------------------------
function pgex9()
    opendevice() # make sure a plotting device is open

    # Call PGFUN to draw the function (autoscaling).
    pgbbuf()
    pgsave()
    pgsci(5)
    pgfunt(t -> sin(5t),
           t -> sin(4t), range(0, stop=2π, length=360), false)

    # Call PGLAB to label the graph in a different color.
    pgsci(3)
    pglab("x","y","PGPlot Example 9:  routine PGFUN")
    pgunsa()
    pgebuf()
end

#----------------------------------------------------------------------
# Demonstration program for the PGPlot plotting package.  This example
# illustrates curve drawing with PGFUNX.
#                              T. J. Pearson  1983 Oct 5
#----------------------------------------------------------------------
function pgex10()
    opendevice() # make sure a plotting device is open

    # The following define mnemonic names for the color indices and
    # linestyle codes.
    black   = 0
    white   = 1
    red     = 2
    green   = 3
    blue    = 4
    cyan    = 5
    magenta = 6
    yellow  = 7
    full    = 1
    dash    = 2
    dotd    = 3

    # Call PGFUN twice to draw two functions (autoscaling the first time).
    pgbbuf()
    pgsave()
    pgsci(yellow)
    pgfunx(range(0, stop=10π, length=500), pgbsj0, false)
    pgsci(red)
    pgsls(dash)
    pgfunx(range(0, stop=10π, length=500), pgbsj1, true)

    # Call PGLAB to label the graph in a different color. Note the use of
    # "\\f" to change font.  Use PGMTXT to write an additional legend
    # inside the viewport.
    pgsci(green)
    pgsls(full)
    pglab("\\fix", "\\fiy", "\\frPGPlot Example 10: routine PGFUNX")
    pgmtxt("T", -4.0, 0.5, 0.5, "\\frBessel Functions")

    # Call PGARRO to label the curves.
    pgarro(8.0, 0.7, 1.0, pgbsj0(1.0))
    pgarro(12.0, 0.5, 9.0, pgbsj1(9.0))
    pgstbg(green)
    pgsci(0)
    pgptxt(8.0, 0.7, 0.0, 0.0, " \\fiy = J\\d0\\u(x)")
    pgptxt(12.0, 0.5, 0.0, 0.0, " \\fiy = J\\d1\\u(x)")
    pgunsa()
    pgebuf()
end

#-----------------------------------------------------------------------
# Bessel function of order 0 (approximate).
# Reference: Abramowitz and Stegun: Handbook of Mathematical Functions.
#-----------------------------------------------------------------------
function pgbsj0(xx)
    x = abs(convert(Float64, xx))
    if x ≤ 3
        t = (x/3)^2
        return 1.0 + t*(-2.2499997 + t*( 1.2656208 + t*(-0.3163866 + t*( 0.0444479 + t*(-0.0039444 + t*( 0.0002100))))))
    else
        t = 3/x
        f0 = 0.79788456 + t*(-0.00000077 + t*(-0.00552740 + t*(-0.00009512 + t*( 0.00137237 + t*(-0.00072805 + t*( 0.00014476))))))
        theta0 = x - 0.78539816 + t*(-0.04166397 + t*(-0.00003954 + t*( 0.00262573 + t*(-0.00054125 + t*(-0.00029333 + t*( 0.00013558))))))
        return f0*cos(theta0)/sqrt(x)
    end
end

#-----------------------------------------------------------------------
# Bessel function of order 1 (approximate).
# Reference: Abramowitz and Stegun: Handbook of Mathematical Functions.
#-----------------------------------------------------------------------
function pgbsj1(xx)
    x = abs(convert(Float64, xx))
    if x ≤ 3
        t = (x/3)^2
        return (0.5 + t*(-0.56249985 + t*( 0.21093573 + t*(-0.03954289 + t*( 0.00443319 + t*(-0.00031761 + t*( 0.00001109)))))))*x
      else
         t = 3/x
         f1 = 0.79788456 + t*( 0.00000156 + t*( 0.01659667 + t*( 0.00017105 + t*(-0.00249511 + t*( 0.00113653 + t*(-0.00020033))))))
         theta1 = x - 2.35619449 + t*( 0.12499612 + t*( 0.00005650 + t*(-0.00637879 + t*( 0.00074348 + t*( 0.00079824 + t*(-0.00029166))))))
        r = f1*cos(theta1)/sqrt(x)
        if xx < 0
            r = -r
        end
        return r
    end
end

"""
     dodecahedron(T = Cdouble)

yields the Cartesian coordinates of the 20 vertices of a dodecahedron.

"""
function dodecahedron(::Type{T} = Cdouble) where {T <: AbstractFloat}
    # Cartesian coordinates of the 20 vertices.
    nvert = 20
    o = zero(T)
    l = one(T)
    t = convert(T, 1.618)
    t1 = l + t
    t2 = -t
    t3 = -t1
    return reshape([t,t,t, t,t,t2, t,t2,t, t,t2,t2, t2,t,t, t2,t,t2, t2,t2,t,
                    t2,t2,t2, t1,l,o, t1,-l,o, t3,l,o, t3,-l,o,
                    o,t1,l, o,t1,-l, o,t3,l, o,t3,-l, l,o,t1,
                    -l,o,t1, l,o,t3, -l,o,t3], 3, 20)
end

#-----------------------------------------------------------------------
# Test routine for PGPlot: draws a skeletal dodecahedron.
#-----------------------------------------------------------------------

function pgex11()
    opendevice() # make sure a plotting device is open

    # Cartesian coordinates of the 20 vertices.
    vert = dodecahedron(PGFloat)
    nvert = size(vert, 2)

    # Initialize the plot (no labels).
    pgbbuf()
    pgsave()
    pgenv(-4.0, 4.0, -4.0, 4.0, true, -2)
    pgsci(2)
    pgsls(1)
    pgslw(1)

    # Write a heading.
    pglab("", "", "PGPlot Example 11:  Dodecahedron")

    # Mark the vertices.
    for i in 1:nvert
        xi, yi, zi = vert[1,i], vert[2,i], vert[3,i]
        pgpt(xi + 0.2*zi, yi + 0.3*zi, 9)
    end

    # Draw the edges - test all vertex pairs to find the edges of the
    # correct length.
    x = Array{PGFloat}(undef, 2)
    y = Array{PGFloat}(undef, 2)
    pgslw(3)
    for i in 2:nvert
        for j in 1:i-1
              r = 0.0
              for k in 1:3
                  r += (vert[k,i] - vert[k,j])^2
              end
              r = sqrt(r)
              if abs(r - 2.0) ≤ 0.1
                  zi = vert[3,i]
                  x[1] = vert[1,i] + 0.2*zi
                  y[1] = vert[2,i] + 0.3*zi
                  zj = vert[3,j]
                  x[2] = vert[1,j] + 0.2*zj
                  y[2] = vert[2,j] + 0.3*zj
                  pgline(x,y)
              end
        end
    end
    pgunsa()
    pgebuf()
end

#-----------------------------------------------------------------------
# Test routine for PGPlot: draw arrows with PGARRO.
#-----------------------------------------------------------------------
function pgex12()
    opendevice() # make sure a plotting device is open

    # Number of arrows.
    nv = 16

    # Select a square viewport.
    pgbbuf()
    pgsave()
    pgsch(0.7)
    pgsci(2)
    pgenv(-1.05, 1.05, -1.05, 1.05, true, -1)
    pglab("", "", "PGPlot Example 12:  PGARRO")
    pgsci(1)

    # Draw the arrows
    k = 1
    d = 360.0/57.29577951/nv
    a = -d
    for i in 1:nv
        a += d
        x = cos(a)
        y = sin(a)
        xt = 0.2*cos(a - d)
        yt = 0.2*sin(a - d)
        pgsah(k, 80.0-3.0*i, 0.5*real(i)/real(nv))
        pgsch(0.25*i)
        pgarro(xt, yt, x, y)
        k = (k == 1 ? 2 : 1)
   end
    pgunsa()
    pgebuf()
end

#----------------------------------------------------------------------
# This example illustrates the use of PGTBOX.
#----------------------------------------------------------------------
function pgex13()
    opendevice() # make sure a plotting device is open

    x1 = [0.0, 0.0, 0.0, 0.0, -8000.0, 100.3, 205.3, -45000.0, 0.0, 0.0]
    x2 = [8000.0, 8000.0, 8000.0, 8000.0, 8000.0, 101.3, 201.1,
          -100000.0, -100000.0, -100000.0]
    xopt = ("BSTN", "BSTNZ", "BSTNZH", "BSTNZD", "BSNTZHFO",
            "BSTNZD", "BSTNZHI", "BSTNZHP", "BSTNZDY", "BSNTZHFOY")
    @assert length(x1) == length(x2) == length(xopt)
    n = length(xopt)

    pgpage()
    pgsave()
    pgbbuf()
    pgsch(0.7)
    for i in 1:n
        pgsvp(0.15, 0.85, (0.7 + (n - i))/n, (0.7 + (n - i + 1))/n)
        pgswin(x1[i], x2[i], 0.0, 1.0)
        pgtbox(xopt[i], 0.0, 0, "", 0.0, 0)
        pglab("Option = $(xopt[i])", "", "")
        if i == 1
            pgmtxt("B", -1.0, 0.5, 0.5, "\\fiAxes drawn with PGTBOX")
        end
    end
    pgebuf()
    pgunsa()
end

#-----------------------------------------------------------------------
# Test routine for PGPlot: polygon fill and color representation.
#-----------------------------------------------------------------------
function pgex14()
    opendevice() # make sure a plotting device is open

    n = 33
    m = 8
    thinc = 2π/n
    xi = zeros(PGFloat, n)
    yi = zeros(PGFloat, n)
    xo = zeros(PGFloat, n)
    yo = zeros(PGFloat, n)
    xt = zeros(PGFloat, 3)
    yt = zeros(PGFloat, 3)
    pgbbuf()
    pgsave()
    pgenv(-1.0, 1.0, -1.0, 1.0, true, -2)
    pglab("", "", "PGPlot Example 14: PGPOLY and PGSCR")
    for j in 1 : m
        r = 1.0
        g = 1.0 - j/m
        b = g
        pgscr(j, r, g, b)
        theta = -j*thinc/2
        r = j/m
        for i in 1 : n
            theta += thinc
            xo[i] = r*cos(theta)
            yo[i] = r*sin(theta)
        end
        for i in 1 : n
            xt[1] = xo[i]
            yt[1] = yo[i]
            xt[2] = xo[mod(i,n) + 1]
            yt[2] = yo[mod(i,n) + 1]
            xt[3] = xi[i]
            yt[3] = yi[i]
            pgsci(j)
            pgsfs(1)
            pgpoly(xt, yt)
            pgsfs(2)
            pgsci(1)
            pgpoly(xt, yt)
        end
        for i in 1 : n
            xi[i] = xo[i]
            yi[i] = yo[i]
        end
    end
    pgunsa()
    pgebuf()
end

#----------------------------------------------------------------------
# This is a line-drawing test; it draws a regular n-gon joining
# each vertex to every other vertex. It is not optimized for pen
# plotters.
#----------------------------------------------------------------------
function pgex15()
    opendevice() # make sure a plotting device is open

    # Set the number of vertices, and compute the coordinates for unit
    # circumradius.
    nv = 17
    c = 360.0/57.29577951/nv
    x = [cos((i - 1)*c) for i in 1 : nv]
    y = [sin((i - 1)*c) for i in 1 : nv]

    # Select a square viewport.
    pgbbuf()
    pgsave()
    pgsch(0.5)
    pgsci(2)
    pgenv(-1.05, 1.05, -1.05, 1.05, true, -1)
    pglab("", "", "PGPlot Example 15: PGMOVE and PGDRAW")
    pgscr(0, 0.2, 0.3, 0.3)
    pgscr(1, 1.0, 0.5, 0.2)
    pgscr(2, 0.2, 0.5, 1.0)
    pgsci(1)

    # Draw the polygon.
    for i in 1 : nv - 1
        for j in i + 1 : nv
            pgmove(x[i], y[i])
            pgdraw(x[j], y[j])
        end
    end

    # Flush the buffer.
    pgunsa()
    pgebuf()
end

#-----------------------------------------------------------------------
# Returns a normally distributed deviate with zero mean and unit
# variance. The routine uses the Box-Muller transformation of uniform
# deviates. For a more efficient implementation of this algorithm,
# see Press et al., Numerical Recipes, Sec. 7.2.
#
# Arguments:
#  ISEED  (in/out) : seed used for PGRAND random-number generator.
#
# Subroutines required:
#  PGRAND -- return a uniform random deviate between 0 and 1.
#
# History:
#  1995 Dec 12 - TJP.
#-----------------------------------------------------------------------
function pgrnrm(seed::Ref{T}) where {T <: Integer}
    while true
        x = 2.0*pgrand(seed) - 1.0
        y = 2.0*pgrand(seed) - 1.0
        r = x^2 + y^2
        if 0.0 < r < 1.0
            return x*sqrt(-2.0*log(r)/r)
        end
    end
end

#-----------------------------------------------------------------------
# Returns a uniform random deviate between 0.0 and 1.0.
#
# NOTE: this is not a good random-number generator; it is only
# intended for exercising the PGPlot routines.
#
# Based on: Park and Miller's "Minimal Standard" random number
#   generator (Comm. ACM, 31, 1192, 1988)
#
# Arguments:
#  ISEED  (in/out) : seed.
#-----------------------------------------------------------------------
function pgrand(seed::Ref{T}) where {T <: Integer}
    im = convert(T, 2_147_483_647)
    ia = convert(T,        16_807)
    iq = convert(T,       127_773)
    ir = convert(T,         2_836)
    k = div(seed[], iq)
    seed[] = ia*(seed[] - k*iq) - ir*k
    if seed[] < 0
        seed[] += im
    end
    return seed[]/im
end

#-----------------------------------------------------------------------
# Demonstration program for PGPlot. The main program opens the output
# device and calls a series of subroutines, one for each sample plot.
#-----------------------------------------------------------------------
function pgdem2(dev::AbstractString = DEFAULT_DEVICE)
    # Call PGOPEN to initiate PGPlot and open the output device; PGOPEN
    # will prompt the user to supply the device name and type.
    id = pgopen(dev)
    setupdevice()

    # Call the demonstration subroutines.
    pgex21()
    pgex22()
    pgex23()
    pgex24()
    pgex25()
    pgex26()

    # Finally, call PGCLOS to terminate things properly.
    pgclos()
end

#-----------------------------------------------------------------------
# Test subroutine for PGPlot: screen alignment and color palette.
#-----------------------------------------------------------------------
function pgex21()
    opendevice() # make sure a plotting device is open

    # Get PGPlot information.
    gver = pgqinf("VERSION")
    gtype = pgqinf("TYPE")
    pgbbuf()

    # Alignment test: clear the screen, and draw a box and grid using
    # three monochrome intensities (color indices 1, 14, and 15).  The
    # plot uses the largest available square viewport and and unit window.
    pgpage()
    pgsvp(0.0, 1.0, 0.0, 1.0)
    pgwnad(0.0, 1.0, 0.0, 1.0)
    pgsci(14)
    pgbox("g", 0.02, 1, "g", 0.02, 1)
    pgsci(15)
    pgbox("g", 0.1, 5, "g", 0.1, 5)
    pgsci(1)
    pgbox("bc", 0.1, 5, "bc", 0.1, 5)

    # Color palette test.
    for i in 0 : 15
        pgsci(i)
        x1 = 0.31 + mod(i,4)*0.1
        y1 = 0.61 - div(i,4)*0.1
        x2 = x1 + 0.08
        y2 = y1 + 0.08
        pgrect(x1, x2, y1, y2)
    end

    # Write the device type on the plot.
    pgsci(0)
    pgrect(0.31, 1.0 - 0.31, 0.85, 0.97)
    pgsci(1)
    pgsfs(2)
    pgrect(0.31, 1.0 - 0.31, 0.85, 0.97)
    pgptxt(0.5, 0.91, 0.0, 0.5, "PGPlot $gver")
    pgptxt(0.5, 0.87, 0.0, 0.5, "Device $gtype")

    pgebuf()
end

#-----------------------------------------------------------------------
# Demonstration program for the PGPlot plotting package.
# Plot a table of the standard PGPlot graph marker symbols. This
# program also illustrates how windows and viewports may be manipulated.
#-----------------------------------------------------------------------
function pgex22()
    opendevice() # make sure a plotting device is open

    # Determine size of view surface.
    # Lower left corner is (X1,Y1), upper right (X2, Y2) [inches].
    pgpage()
    pgsvp(0.0, 1.0, 0.0, 1.0)
    x1, x2, y1, y2 = pgqvp(1)
    x = x2 - x1
    y = y2 - y1

    # Determine device resolution (pixels/inch), and use it to choose
    # line width.
    xpix1, xpix2, ypix1, ypix2 = pgqvp(3)
    res = abs(xpix2 - xpix1)/abs(x)
    lw = (res > 166.0 ? 2 : 1)

    # Choose horizontal or vertical arrangement depending on
    # device aspect ratio.
    if x > y
        nx = 8
        ny = 5
    else
        nx = 5
        ny = 8
    end
    dx = min(x/nx, 0.95*y/ny)
    dy = dx
    ix = nx
    jy = 1
    xoff = x1 + (x-nx*dx)*0.5
    yoff = y1 + (0.95*y-ny*dy)*0.5
    pgbbuf()

    # Each symbol will be drawn in a standard window; the window is moved
    # by manipulating the viewport.
    pgswin(-1.0,1.0,-1.0,1.0)

    # Loop through all known symbols (N=0-31 and -1 to -8).
    for n in 0 : 39
        label = @sprintf("(%2d)", (n ≤ 31 ? n : 31 - n))

        # Define window and viewport. The loop allows the plot to extend over
        # more than one page if necessary; each page is labelled at the top.
        ix += 1
        if ix > nx
            ix = 1
            jy -= 1
        end
        if jy < 1
            jy = ny
            n != 0 && pgpage()
            pgsch(1.2)
            pgvsiz(xoff, xoff + nx*dx, yoff, yoff + ny*dy)
            pgslw(lw)
            pgmtxt("T", 1.0, 0.5, 0.5, "\\fiPGPlot \\frMarker Symbols")
        end
        pgvsiz(xoff + (ix-1)*dx, xoff + ix*dx, yoff + (jy-1)*dy, yoff + jy*dy)

        # Call PGBOX to draw a box and PGMTXT to label it.
        pgslw(1)
        pgbox("BC", 10.0, 0, "BC", 10.0, 0)
        pgsch(0.5)
        pgmtxt("T", -1.5, 0.05, 0.0, label)

        # Call PGPT1 to draw the symbol.
        pgslw(lw)
        pgsch(1.5)
        pgpt(0.0, 0.0, (n ≤ 31 ? n : 31 - n))
    end
    pgebuf()
end

#-----------------------------------------------------------------------
# Demonstration program for the PGPlot plotting package.
#-----------------------------------------------------------------------
function pgex23()
    opendevice() # make sure a plotting device is open

    sample = ("Normal:  \\fnABCDQ efgh 1234 \\ga\\gb\\gg\\gd \\gL\\gH\\gD\\gW",
              "Roman:  \\frABCDQ efgh 1234 \\ga\\gb\\gg\\gd \\gL\\gH\\gD\\gW",
              "Italic:  \\fiABCDQ efgh 1234 \\ga\\gb\\gg\\gd \\gL\\gH\\gD\\gW",
              "Script:  \\fsABCDQ efgh 1234 \\ga\\gb\\gg\\gd \\gL\\gH\\gD\\gW",
              "\\fif\\fr(\\fix\\fr) = \\fix\\fr\\u2\\dcos(2\\gp\\fix\\fr)e\\u\\fix\\fr\\u2",
              "\\fiH\\d0\\u \\fr= 75 \\(2233) 25 km s\\u-1\\d Mpc\\u-1\\d",
              "\\fsL/L\\d\\(2281)\\u\\fr = 5\\.6 \\x 10\\u6\\d (\\gl1216\\A)",
              "Markers: 3=\\m3, 8=\\m8, 12=\\m12, 28=\\m28.",
              "Cyrillic: \\(2830)\\(2912)\\(2906)\\(2911)\\(2919)\\(2917)\\(2915).")
    n = length(sample)

    # Call PGENV to initialize the viewport and window.
    # Call PGLAB to label the graph.
    pgenv(0.0, 20.0, n, 0.0, false, -2)
    pglab("","","\\fiPGPlot \\frFonts")

    # Use PGTEXT to write the sample character strings.
    pgsch(1.6)
    for i in 1 : n
        pgtext(0.0, i - 0.5, sample[i])
    end
    pgsch(1.0)
end

#-----------------------------------------------------------------------
# Demonstration program for the PGPlot plotting package. This example
# illustrates the different line widths.
#                              T. J. Pearson  1982 Dec 28
#----------------------------------------------------------------------
function pgex24()
    opendevice() # make sure a plotting device is open

    x = Array{PGFloat}(undef, 2)
    y = Array{PGFloat}(undef, 2)

    # Call PGENV to initialize the viewport and window.
    pgbbuf()
    pgenv(0.0, 15.0, 0.0, 15.0, false, 0)

    # Call PGLAB to label the graph.
    pglab("Line Width","","\\fiPGPlot \\frLine Widths")

    # Draw 14 oblique lines in different thicknesses.
    for iw in 1 : 14
        x[1] = iw
        y[1] = 0
        x[2] = 0
        y[2] = iw
        pgslw(iw)
        pgline(x, y)
    end

    # Draw another set of lines, dashed instead of solid.
    pgsls(2)
    for iw in 1 : 14
        x[1] = iw
        y[1] = 15
        x[2] = 15
        y[2] = iw
        pgslw(iw)
        pgline(x, y)
    end
    pgsls(1)
    pgslw(1)
    pgebuf()
end

function duck(::Type{T} =  Cdouble) where {T}
    px = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.4, 17.0, 17.3,
          17.8, 18.5, 20.0, 22.0, 24.0, 26.0, 28.0, 29.0, 28.8, 27.2,
          25.0, 23.0, 21.5, 21.1, 21.5, 22.8, 24.1, 25.1, 25.2, 24.2,
          22.1, 20.0, 18.0, 16.0, 14.0, 12.0, 10.0,  8.0,  6.1,  4.2,  3.0,  1.3]
    py = [8.8, 7.6, 7.1, 7.4, 8.0, 8.9, 9.6, 9.9, 9.4, 9.7, 12.0, 14.0,
          16.1, 17.0, 17.0, 16.0, 13.9, 13.1, 13.2, 12.3, 11.5, 11.5, 11.5,
          11.2, 10.5, 9.0, 8.0, 7.0, 5.1, 3.6, 1.9, 1.1, 0.9, 0.7, 0.8, 1.0,
          1.0, 1.2, 1.8, 2.1, 2.9, 4.1, 6.0]
    return (convert(Array{T,1}, px), convert(Array{T,1}, py))
end

#-----------------------------------------------------------------------
# Demonstration program for the PGPlot plotting package. This program
# tests polygon clipping on polygons and circles, and tests that
# markers are clipped correctly. Note that markers exactly on the edge
# of the window are supposed to be visible.
#                              T. J. Pearson  1994 Nov 25
#-----------------------------------------------------------------------
function pgex25()
    opendevice() # make sure a plotting device is open

    px, py = duck(PGFloat)
    sx = [10.0, 10.0, 20.0, 30.0, 15.0]
    sy = [ 0.0, -6.0, -6.0,  5.0, -3.5]
    rx = [26.0, 27.0, 26.0]
    ry = [-4.0, -3.0, -3.0]
    #
    pgpage()
    pgvstd()
    pgbbuf()
    pgsave()

    # Set window.
    pgwnad(5.0, 25.0, -5.0, 15.0)
    pgsci(1)
    pgslw(1)
    pgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0)

    # Test clipping of polygons and circles
    pgsfs(1)
    pgsci(2)
    pgpoly(px, py)
    pgsci(0)
    pgsfs(3)
    pgshs(30.0, 2.0, 0.0)
    pgpoly(px, py)
    pgsci(1)
    pgshs(30.0, 4.0, 0.25)
    pgpoly(px, py)
    pgsci(1)
    pgsfs(2)
    pgpoly(px, py)
    #
    pgsfs(1)
    pgsci(4)
    pgpoly(sx, sy)
    pgsci(0)
    pgsfs(4)
    pgshs(0.0, 1.6, 0.0)
    pgpoly(sx, sy)
    pgsci(1)
    pgsfs(2)
    pgpoly(sx, sy)

    # The next polygon should be invisible.
    pgsfs(1)
    pgsci(4)
    pgpoly(rx, ry)
    pgsci(1)
    pgsfs(2)
    pgpoly(rx, ry)
    #
    pgsfs(1)
    pgsci(3)
    pgcirc(8.0, 12.0, 3.5)
    pgsfs(2)
    pgsci(1)
    pgcirc(8.0, 12.0, 3.5)

    # Test clipping of markers: all should be visible.
    pgsci(1)
    pgslw(1)
    for i in 0 : 5 : 30
        for j in -5 : 5 :25
            pgpt(i,j,9)
        end
    end

    # Draw box.
    pgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0)
    pglab("", "", "PGPlot: clipping polygons and markers")

    pgunsa()
    pgebuf()
end

function pgex26()
    opendevice() # make sure a plotting device is open

    px, py = duck(PGFloat)
    device = pgqinf("DEV/TYPE")
    gver = pgqinf("VERSION")
    gtype = pgqinf("TYPE")
    pgbbuf()

    # Clear the screen; set background and foreground colors.
    pgpage()
    pgscr(0, 0.0, 0.0, 0.35)
    pgscr(1, 1.0, 1.0, 1.0)
    pgeras()

    # Draw a frame at the physical extremities of the plot.
    # Dimensions are X by Y (inches).
      pgsvp(0.0, 1.0, 0.0, 1.0)
      x1, x2, y1, y2 = pgqvp(1)
      x = x2 - x1
      y = y2 - y1
      pgswin(0.0, x, 0.0, y)
      pgsfs(2)
      pgrect(0.0, x, 0.0, y)
      pgmove(0.5*x, 0.0)
      pgdraw(0.5*x, y)
      pgmove(0.0, 0.5*y)
      pgdraw(x, 0.5*y)

    # Draw a circle of diameter 0.5 x min(x,y)
    r = 0.25*min(x,y)
    pgcirc(x*0.5, y*0.5, r)

    # Draw some more circles with different line-styles; this tests
    # the dashing algorithm on curved lines.
    pgsls(2)
    pgcirc(x*0.5, y*0.5, r*1.1)
    pgsls(3)
    pgcirc(x*0.5, y*0.5, r*1.2)
    pgsls(2)
    pgslw(3)
    pgcirc(x*0.5, y*0.5, r*1.3)
    pgsls(1)
    pgslw(1)

    # Demonstrate different line-styles
    for i in 1 : 5
        pgsls(i)
        pgmove(i*(x/20.0), 0.0)
        pgdraw(i*(x/20.0), y)
    end
    pgsls(1)

    # Demonstrate different line-widths
    for i in 1 : 5
        pgslw(i)
        pgmove(0.0, i*(y/20.0))
        pgdraw(x, i*(y/20.0))
    end
    pgslw(1)

    # Demonstrate different line-colors
    pgslw(4)
    for i in 0 : 15
        pgsci(i)
        xi = (i + 20)*(x/40.0)
        pgmove(xi,0.0)
        pgdraw(xi,y)
    end
    pgsci(1)
    pgslw(1)

    # Draw dots in different thicknesses.
    xp = (14 + 20)*(x/40.0)
    for i in 1 : 21
        yp = i*y/22.0
        pgslw(i)
        pgpt(xp, yp, -1)
    end
    pgslw(1)

    # Demonstrate fill area
    for j in 1 : length(px)
        px[j] = (px[j] + 50.0)/100.0*x
        py[j] = (py[j] + 75.0)/100.0*y
    end
    for i in 0 : 3
        pgsci(i)
        pgsfs(1)
        pgpoly(px, py)
        pgsci(1)
        pgsfs(2)
        pgpoly(px, py)
        for j in 1 : length(py)
             py[j] -= 0.25*y
        end
    end

    # Write the device type on the plot.
    pgswin(0.0, 1.0, 0.0, 1.0)
    pgsfs(1)
    pgsci(0)
    pgrect(0.31, 1.0-0.31, 0.85, 0.97)
    pgsci(1)
    pgsfs(2)
    pgrect(0.31, 1.0-0.31, 0.85, 0.97)
    pgptxt(0.5, 0.91, 0.0, 0.5, "PGPlot $gver")
    pgptxt(0.5, 0.87, 0.0, 0.5, "Device $gtype")

    pgebuf()
end

#-----------------------------------------------------------------------
# PGDEM3
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Demonstration program for PGPlot contouring routines.
#-----------------------------------------------------------------------
function pgdem3(dev::AbstractString = DEFAULT_DEVICE)
    println(" Demonstration of PGPlot contouring routines")

    # Call PGBEG to initiate PGPlot and open the output device; PGBEG
    # will prompt the user to supply the device name and type.
    pgbeg(0,dev,1,1)
    setupdevice()

    # Call the demonstration subroutines.
    println(" Routine PGCONT")
    pgex31()
    println(" Routine PGCONS")
    pgex32()
    println(" Routine PGCONT with PGCONL labels")
    pgex36()
    println(" Routine PGCONB")
    pgex33()
    #println(" Routine PGCONX with arrow labels")
    #pgex37()
    #println(" Routine PGCONX")
    #pgex34()
    println(" Routine PGCONT and PGCONL")
    pgex36()

    # Finally, call PGEND to terminate things properly.
    pgend()
end

#-----------------------------------------------------------------------
# Demonstration of contouring routine PGCONT.
#-----------------------------------------------------------------------
function pgex31()
    opendevice() # make sure a plotting device is open

    ni = 40
    nj = 40
    #tr = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0]

    # Compute a suitable function.
    f = [cos(0.3*sqrt(2i) - 0.4*j/3)*cos(0.4*i/3) + (i - j)/40
         for i in 1:ni, j in 1:nj]
    fmin, fmax = extrema(f)

    # Clear the screen. Set up window and viewport.
    pgpage()
    pgsvp(0.05, 0.95, 0.05, 0.95)
    pgswin(1.0, 40.0, 1.0, 40.0)
    pgbox("bcts", 0.0, 0, "bcts", 0.0, 0)
    pgmtxt("t", 1.0, 0.0, 0.0, "Contouring using PGCONT")

    # Draw the map.  PGCONT is called once for each contour, using
    # different line attributes to distinguish contour levels.
    pgbbuf()
    for i in 1 : 21
        lev = fmin + (i - 1)*(fmax - fmin)/20
        if mod(i,5) == 0
            pgslw(3)
        else
            pgslw(1)
        end
        if i < 10
            pgsci(2)
            pgsls(2)
        else
            pgsci(3)
            pgsls(1)
        end
        pgcont(f, lev)
    end
    pgslw(1)
    pgsls(1)
    pgsci(1)
    pgebuf()
end

#-----------------------------------------------------------------------
# Demonstration of contouring routine PGCONS.
#-----------------------------------------------------------------------
function pgex32()
    opendevice() # make sure a plotting device is open

    ni = 40
    nj = 40

    # Compute a suitable function.
    f = [cos(0.3*sqrt(2i) - 0.4*j/3)*cos(0.4*i/3) + (i - j)/40
         for i in 1:ni, j in 1:nj]
    fmin, fmax = extrema(f)

    # Clear the screen. Set up window and viewport.
    pgpage()
    pgsvp(0.05, 0.95, 0.05, 0.95)
    pgswin(1.0, 40.0, 1.0, 40.0)
    pgbox("bcts", 0.0, 0, "bcts", 0.0, 0)
    pgmtxt("t", 1.0, 0.0, 0.0, "Contouring using PGCONS")

    # Draw the map.  PGCONS is called once for each contour, using
    # different line attributes to distinguish contour levels.
    pgbbuf()
    for i in 1 : 21
        lev = fmin + (i - 1)*(fmax - fmin)/20
        if mod(i,5) == 0
            pgslw(3)
        else
            pgslw(1)
        end
        if i < 10
            pgsci(2)
            pgsls(2)
        else
            pgsci(3)
            pgsls(1)
        end
        pgcons(f, lev)
    end
    pgslw(1)
    pgsls(1)
    pgsci(1)
    pgebuf()
end

#-----------------------------------------------------------------------
# Demonstration of contouring routine PGCONB.
#-----------------------------------------------------------------------
function pgex33()
    opendevice() # make sure a plotting device is open

    blank = -1.2e20
    ni = 40
    nj = 40
    tr = DEFAULT_XFORM

    # Compute a suitable function.
    f = [cos(0.3*sqrt(2i) - 0.4*j/3)*cos(0.4*i/3) + (i - j)/40
         for i in 1:ni, j in 1:nj]
    fmin, fmax = extrema(f)

    # "Blank" the data outside an annulus.
    for j in 1 : nj, i in 1 : ni
        r = sqrt((i - 20.5)^2 + (j - 20.5)^2)
        if r > 20.0 || r < 3.0
            f[i,j] = blank
        end
    end

    # Clear the screen. Set up window and viewport.
    pgpage()
    pgsvp(0.05, 0.95, 0.05, 0.95)
    pgswin(1.0, 40.0, 1.0, 40.0)
    pgbox("bcts", 0.0, 0, "bcts", 0.0, 0)
    pgmtxt("t", 1.0, 0.0, 0.0,"Contouring using PGCONB")
    pgbbuf()
    for i in 1 : 21
        level = fmin + (i - 1)*(fmax - fmin)/20
        if mod(i,5) == 0
            pgslw(3)
        else
            pgslw(1)
        end
        if i < 10
            pgsci(2)
            pgsls(2)
        else
            pgsci(3)
            pgsls(1)
        end
        pgconb(f, level, tr, blank)
    end
    pgebuf()

    # Mark the blanked points for easy identification.
    pgbbuf()
    pgsci(1)
    for j in 1 : nj, i in 1 : ni
        if f[i,j] == blank
            x = tr[1] + i*tr[2] + j*tr[3]
            y = tr[4] + i*tr[5] + j*tr[6]
            pgpt(x, y, -1)
        end
    end
    pgebuf()
end

#-----------------------------------------------------------------------
# Demonstration of contouring routine PGCONT and PGCONL.
#-----------------------------------------------------------------------
function pgex36()
    opendevice() # make sure a plotting device is open

    ni = 40
    nj = 40
    tr = DEFAULT_XFORM

    # Compute a suitable function.
    f = [cos(0.3*sqrt(2i) - 0.4*j/3)*cos(0.4*i/3) + (i - j)/40
         for i in 1:ni, j in 1:nj]
    fmin, fmax = extrema(f)

    # Clear the screen. Set up window and viewport.
    pgpage()
    pgsvp(0.05, 0.95, 0.05, 0.95)
    pgswin(1.0, 40.0, 1.0, 40.0)
    pgbox("bcts", 0.0, 0, "bcts", 0.0, 0)
    pgmtxt("t", 1.0, 0.0, 0.0, "Contouring using PGCONT and PGCONL labels")

    # Draw the map.  PGCONT is called once for each contour, using
    # different line attributes to distinguish contour levels.
    pgbbuf()
    for i in 1 : 21
        level = fmin + (i - 1)*(fmax-fmin)/20
        if mod(i,5) == 0
            pgslw(3)
        else
            pgslw(1)
        end
        if i < 10
            pgsci(2)
            pgsls(2)
        else
            pgsci(3)
            pgsls(1)
        end
        pgcont(f, level)
    end
    pgslw(1)
    pgsls(1)
    pgebuf()

    # Label the contours with PGCONL. Only even-numbered contours are
    # labelled.
    pgbbuf()
    for i in 2 : 2 : 21
        level = fmin + (i - 1)*(fmax - fmin)/20
        label = repr(i)
        # label = @sprintf("%.2f", level)
        if i < 10
            pgsci(2)
        else
            pgsci(3)
        end
        pgconl(f, level, label, 16, 8, tr)
    end
    pgsci(1)
    pgebuf()
end

end # module
