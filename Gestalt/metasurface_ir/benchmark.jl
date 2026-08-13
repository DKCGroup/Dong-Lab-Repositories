# benchmark.jl — timing vs key parameters for batch training-set generation.
#
# Run:  julia --project=. -t auto metasurface_ir/benchmark.jl
#
# Measures the quantities that dominate per-run wall time:
#   1) single RCWA solve time  vs truncation order N  (the O(N^6) bottleneck);
#   2) full spectral sweep time vs n_lambda (linear);
#   3) pattern-layer count via staircase slices nslices (linear);
#   4) homogeneous-only stack (block-diagonal fast path) vs one patterned layer.

using GESTALT
using LinearAlgebra
include(joinpath(@__DIR__, "metasurface_ir.jl"))
BLAS.set_num_threads(1)               # sweep uses Threads.@threads; keep BLAS single

const NTHREADS = Threads.nthreads()
println("threads = ", NTHREADS)

vo2 = vo2_metallic_sampled()
al  = aluminum_material()
p   = 5e-6
img = shape_image("circle"; N=64)
f   = c₀ / 10e-6                       # single wavelength 10 µm

# min-of-reps wall time for one TE+TM pair; returns per-solve seconds.
function time_solve_pair(grid, layers; reps=3)
    solve(layers, grid, f; pol_TE=1.0, pol_TM=0.0)
    solve(layers, grid, f; pol_TE=0.0, pol_TM=1.0)
    tmin = Inf
    for _ in 1:reps
        t = @elapsed begin
            solve(layers, grid, f; pol_TE=1.0, pol_TM=0.0)
            solve(layers, grid, f; pol_TE=0.0, pol_TM=1.0)
        end
        tmin = min(tmin, t)
    end
    tmin / 2
end

# ---- 1) single solve time vs N ------------------------------------------------
println("\n=== 1) single solve (TE or TM) time vs N  [circle pattern, one patterned layer] ===")
for N in 2:8
    grid, layers = build_stack(img, vo2, al; p=p, N=N)
    t = time_solve_pair(grid, layers)
    println("N=", N, "  NH=", lpad((2N+1)^2, 3), "  t_solve=", round(t*1e3, digits=2), " ms")
end

# ---- 2) sweep time vs n_lambda (fixed N=4) ------------------------------------
println("\n=== 2) sweep time vs n_lambda (N=4, threads=", NTHREADS, ") ===")
grid, layers = build_stack(img, vo2, al; p=p, N=4)
for nλ in (20, 50, 100, 200, 400)
    λs = range(8e-6, 14e-6; length=nλ)
    sweep_spectrum(layers, grid, λs)                      # warmup
    t = @elapsed sweep_spectrum(layers, grid, λs)
    println("n_lambda=", lpad(nλ, 3), "  t_sweep=", round(t, digits=3),
            " s   (per-λ ", round(t/nλ*1e3, digits=2), " ms)")
end

# ---- 3) staircase slices (pattern-layer count) vs time ------------------------
println("\n=== 3) multi-height grayscale: nslices -> pattern layers -> time (N=6) ===")
gmax = 255
for nslices in (1, 2, 4, 8, 16)
    imgm = zeros(UInt8, 64, 64)
    # radial ramp: many heights => exercises the staircase path
    for j in 1:64, i in 1:64
        x = (i - 0.5) / 64; y = (j - 0.5) / 64
        r = sqrt((x - 0.5)^2 + (y - 0.5)^2)
        imgm[i, j] = r < 0.34 ? round(UInt8, clamp(255 * (1 - r / 0.34) + 50, 0, 255)) : 0x00
    end
    grid, layers = build_stack(imgm, vo2, al; p=p, N=6, nslices=nslices)
    npat = count(l -> l isa PeriodicLayer, layers)
    t = time_solve_pair(grid, layers)
    println("nslices=", lpad(nslices, 2), "  patterned_layers=", lpad(npat, 2),
            "  t_solve=", round(t*1e3, digits=1), " ms")
end

# ---- 4) homogeneous-only vs one patterned layer -------------------------------
println("\n=== 4) homogeneous multilayer (block-diagonal) vs one patterned layer (N=6) ===")
# 4a) pure homogeneous stack: air | 10x Ge spacer | Al | air
air = HomogeneousLayer(NKMaterial(1.0), Inf)
ge  = HomogeneousLayer(NKMaterial(4.0), 500e-9)
alh = HomogeneousLayer(al, 200e-9)
hom_layers = AbstractLayer[air, [ge for _ in 1:10]..., alh, air]
ghom = RCWAGrid((2π/p, 0.0), (0.0, 2π/p); Nx=6, Ny=6)
t_hom = time_solve_pair(ghom, hom_layers)
# 4b) patterned: air | VO2 circle | Al | air
gpat, lpat = build_stack(img, vo2, al; p=p, N=6)
t_pat = time_solve_pair(gpat, lpat)
println("10 homogeneous layers : t_solve=", round(t_hom*1e3, digits=3), " ms")
println("1  patterned layer    : t_solve=", round(t_pat*1e3, digits=2), " ms")
println("speedup (pattern vs hom stack) = ", round(t_pat / t_hom, digits=1), "x  (homogeneous is much cheaper)")

println("\nDone.")
