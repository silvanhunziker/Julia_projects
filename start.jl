#=
test
=#

using LinearAlgebra
using DifferentialEquations
using GLMakie
GLMakie.activate!()

function isdSim_getAccelLorentz_averaged(simtime, state, q, m, beta)
    # calculate acceleration of dust particles based and gravity, radiation pressure and Lorentz force.
    state_pos = state[1:3]
    state_vel = state[4:6]

    # Define constants
    VSWmag = 400.e3  # 400km / s
    # B0 = 3.d - 9  #3nT
    B0 = 2300.e-9  # 2300 nT based on 5nT near Earth
    rs = 5.96e8  # solar radius in meter (from t.cravens, physics of solar system plasmas p.160)
    r0 = 10. * rs  # base for parker spiral: 10 Rs
    gm_sol = 1.3271244e+20
    #r0 = isdSimPar_getAU()  # 1 AU
    omegasun = (2. * pi) / (27. * 24. * 60. * 60.)  # 27d period

    # t0etsolcycle : 1976 = solar MIN,
    # tilt angle epsilon = 0 degrees, N+, S- , 0 < epsilon < 360
    # t0etsolcycle = -7.5742555d+08
    # t0etsolcycle : 1974.5 = solar MIN,
    # corresponding to Markus' 1991 solar max
    t0etsolcycle = -8.0480515e+08

    # get particle vel. w.r.t. solar wind, assuming SW to
    # go radially outward at constant speed

    # = np.linalg.norm(state_pos)
    r = sqrt(dot(state_pos, state_pos))
    #rtemp = np.linalg.norm(state_pos[0:2])
    rtemp = sqrt(dot(state_pos[1:2], state_pos[1:2]))
    sincolat = rtemp / r
    coscolat = state_pos[3] / r
    sinlon = state_pos[2] / rtemp
    coslon = state_pos[1] / rtemp
    tanlat = coscolat / sincolat

    vp_rel = state_vel - VSWmag * state_pos / r

    epsilonR = ((0. + 180. / (11. * 365.25 * 3600. * 24.) * (simtime - (r - r0) / VSWmag - t0etsolcycle)) % 360.) * pi / 180.

    taneps = tan(epsilonR)

    p = sign(cos(epsilonR) * coscolat)
    zeta = p

    if abs(tanlat) < abs(taneps)
        zeta = abs(2. / pi * asin(tanlat / taneps))*p
    end

    Br = zeta * B0 * (r0 / r) ^ 2
    Bphi = zeta * B0 * omegasun / VSWmag * (r - r0) * (r0 / r) ^ 2 * sincolat

    Bx = Br * sincolat * coslon + Bphi * sinlon
    By = Br * sincolat * sinlon - Bphi * coslon
    Bz = Br * coscolat

    lor = cross(vp_rel, [Bx, By, Bz]).*q./m

    accel = lor + state_pos * (gm_sol / r^3 * (beta - 1.))
    return accel
end


function define_diff_equation(u, p, t)
    # function that puts together the vectorized differential equation for the ODE solver
    accel = isdSim_getAccelLorentz_averaged(t, u, p[1], p[2], p[3])
    return [u[4], u[5], u[6], accel[1], accel[2], accel[3]]
end


function get_init_state_rand()
    # function that produces initial state vector for dust particles
    AU = 149597870700.0
    x = rand(Float64, 2).*100. .- 50.
    return [-50*AU, x[1]*AU, x[2]*AU, 26000, 0, 0]
end


u0 = get_init_state_rand()  # get initial state vector
tspan = (0.0,500000000.0)  # define starting time and end time for the integrator
p = [20.5, 1., 0.5]  # define dust parameters q, m and beta
alg = Tsit5()  # specify algorithm that should be used in the solver (e.g. Euler, RK, ...)
prob = ODEProblem(define_diff_equation, u0, tspan, p)  # define problem that should be solved
sol = solve(prob, alg, reltol=1e-8, abstol=1e-8)  # calculate trajectory

@time solve(prob, alg, reltol=1e-8, abstol=1e-8)

#@time begin
#for i = 1:1000
#    u_tmp = get_init_state_rand()
#    sol_tmp = solve(prob, alg, reltol=1e-8, abstol=1e-8)
#end
#end

# make 3D plot of 200 trajectories with GLMakie
AU = 149597870700.0
fig = Figure(resolution=(1600, 1600))
ax = Axis3(fig[1, 1]; aspect=(1, 1, 1), perspectiveness=0.5)
for i = 1:200
    u0_tmp = get_init_state_rand()
    prob_tmp = ODEProblem(define_diff_equation, u0_tmp, tspan, p)
    sol_tmp = solve(prob_tmp, alg, reltol=1e-8, abstol=1e-8)
    lines!(ax, sol_tmp[1,:]/AU, sol_tmp[2,:]/AU, sol_tmp[3,:]/AU; linewidth=2, color="black")
end
scatter!(ax, [0], [0], [0]; markersize=50, color="yellow")
xlims!(-50, 50)
ylims!(-50, 50)
zlims!(-50, 50)
fig
