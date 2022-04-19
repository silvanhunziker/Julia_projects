#=
Program that calculates and plots ISD trajectories. The dust particle beta and Q/m parameters can be changed on the fly by sliders in the plot, and the trajectories are re-drawn in real time.
=#

using LinearAlgebra
using DifferentialEquations
using GLMakie
GLMakie.activate!()

#=
Function that calculates dust particle trajectories.
=#
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

#=
Function that creates the proper differential equation that can be integrated by the Julia ODE solver.
=#
function define_diff_equation(u, p, t)
    # function that puts together the vectorized differential equation for the ODE solver
    accel = isdSim_getAccelLorentz_averaged(t, u, p[1], p[2], p[3])
    return [u[4], u[5], u[6], accel[1], accel[2], accel[3]]
end

#=
Function that creates random starting conditions (positions and velocities) for dust particles on a plane at x=-50AU.
=#
function get_init_state_rand()
    # function that produces initial state vector for dust particles
    AU = 149597870700.0
    x = rand(Float64, 2).*100. .- 50.
    return [-50*AU, x[1]*AU, x[2]*AU, 26000, 0, 0]
end

u0 = get_init_state_rand()  # get initial state vector
tspan = (0.0,500000000.0)  # define starting time and end time for the integrator
alg = Tsit5()  # specify algorithm that should be used in the solver (e.g. Euler, RK, ...)

t_range = LinRange(tspan[1], tspan[2], 200)


AU = 149597870700.0
traj_nr = 50  # how many trajectories should be drawn at the same time


# make figure
fig = Figure()
display(fig)
ax = Axis3(fig[1, 1]; aspect=(1, 1, 1), perspectiveness=0.5)  # draw 3d axis with equal aspect ratios
scatter!(ax, [0], [0], [0]; markersize=50, color="yellow")  # draw the sun at the center
xlims!(-50, 50) # set axis limits
ylims!(-50, 50)
zlims!(-50, 50)

# create lists for trajectories
trajectories_x_obs = []
trajectories_y_obs = []
trajectories_z_obs = []

# create sliders for the 3 parameters
lsgrid = labelslidergrid!(
    fig,
    ["q", "m", "beta"],
    [0:0.1:20, 0:0.1:2, 0:0.1:4];
    formats = [x -> "$(round(x, digits = 1))$s" for s in ["C", "kg", ""]],
    width = 550,
    tellheight = true)
fig[2, 1] = lsgrid.layout

# define initial position of the sliders
set_close_to!(lsgrid.sliders[1], 0.)
set_close_to!(lsgrid.sliders[2], 1.)
set_close_to!(lsgrid.sliders[3], 1.)

# calculate the first set of 50 trajectories and store them in: trajectories_x_obs, trajectories_y_obs, trajectories_z_obs
for i in 1:traj_nr
    sliderobservables = [s.value.val for s in lsgrid.sliders]
    u0_new = get_init_state_rand()
    prob_new = ODEProblem(define_diff_equation, u0_new, tspan, sliderobservables)
    sol_new = solve(prob_new, alg, reltol=1e-8, abstol=1e-8)
    x_traj = Observable(sol_new(t_range)[1,:]/AU)
    y_traj = Observable(sol_new(t_range)[2,:]/AU)
    z_traj = Observable(sol_new(t_range)[3,:]/AU)
    push!(trajectories_x_obs, x_traj)
    push!(trajectories_y_obs, y_traj)
    push!(trajectories_z_obs, z_traj)
end

# draw all trajectories stored in : trajectories_x_obs, trajectories_y_obs, trajectories_z_obs
function draw_traj()
    for i in 1:traj_nr
        lines!(ax, trajectories_x_obs[i], trajectories_y_obs[i], trajectories_z_obs[i])
    end
end


# refresh trajectories stored in: trajectories_x_obs, trajectories_y_obs, trajectories_z_obs, based on the current sliders values (re-draws the plot automatically because the stored values are observables).
function refresh_traj_observables()
    for i in 1:traj_nr
        sliderobservables = [s.value.val for s in lsgrid.sliders]
        u0_new = get_init_state_rand()
        prob_new = ODEProblem(define_diff_equation, u0_new, tspan, sliderobservables)
        sol_new = solve(prob_new, alg, reltol=1e-8, abstol=1e-8)
        trajectories_x_obs[i][] = sol_new(t_range)[1,:]/AU
        trajectories_y_obs[i][] = sol_new(t_range)[2,:]/AU
        trajectories_z_obs[i][] = sol_new(t_range)[3,:]/AU
    end
end


# loop 
draw_traj()
while true
    sleep(0.5)
    refresh_traj_observables()
end

