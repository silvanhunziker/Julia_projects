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
#p = [20.5, 1., 0.5]  # define dust parameters q, m and beta
q = 20.5
m = 1.
beta = 0.5
#p=[q, m, beta]
alg = Tsit5()  # specify algorithm that should be used in the solver (e.g. Euler, RK, ...)
prob = ODEProblem(define_diff_equation, u0, tspan, [q, m, beta])  # define problem that should be solved
sol = solve(prob, alg, reltol=1e-8, abstol=1e-8)  # calculate trajectory

t_range = LinRange(tspan[1], tspan[2], 100)


AU = 149597870700.0

#x_tmp = Observable(sol(t_range)[1,:]/AU)
#y_tmp = Observable(sol(t_range)[2,:]/AU)
#z_tmp = Observable(sol(t_range)[3,:]/AU)
x_tmp = sol(t_range)[1,:]/AU
y_tmp = sol(t_range)[2,:]/AU
z_tmp = sol(t_range)[3,:]/AU
fig = Figure()
display(fig)
ax = Axis3(fig[1, 1]; aspect=(1, 1, 1), perspectiveness=0.5)
#lines!(ax, x_tmp, y_tmp, z_tmp)
scatter!(ax, [0], [0], [0]; markersize=50, color="yellow")
xlims!(-50, 50)
ylims!(-50, 50)
zlims!(-50, 50)

#sl_q = Slider(fig[2, 1], range = 0:0.01:20, startvalue = 0)
#sl_m = Slider(fig[3, 1], range = 0:0.01:20, startvalue = 1)
#sl_beta = Slider(fig[4, 1], range = 0:0.01:4, startvalue = 1)

lsgrid = labelslidergrid!(
    fig,
    ["q", "m", "beta"],
    [0:0.1:20, 0:0.1:2, 0:0.1:4];
    formats = [x -> "$(round(x, digits = 1))$s" for s in ["C", "kg", ""]],
    width = 550,
    tellheight = true)

fig[2, 1] = lsgrid.layout

set_close_to!(lsgrid.sliders[1], 0.)
set_close_to!(lsgrid.sliders[2], 1.)
set_close_to!(lsgrid.sliders[3], 1.)

#sliderobservables = [s.value.val for s in lsgrid.sliders]
#params = lift(sliderobservables...) do slvalues...
#    [slvalues...]
#end

#println(params[1])
#println(params[1].val)

function update_fig()
    sliderobservables = [s.value.val for s in lsgrid.sliders]
    u0_new = get_init_state_rand()
    prob_new = ODEProblem(define_diff_equation, u0_new, tspan, sliderobservables)
    sol_new = solve(prob_new, alg, reltol=1e-8, abstol=1e-8)

    #x_tmp[] = sol_new(t_range)[1,:]/AU
    #y_tmp[] = sol_new(t_range)[2,:]/AU
    #z_tmp[] = sol_new(t_range)[3,:]/AU
    x_tmp = sol_new(t_range)[1,:]/AU
    y_tmp = sol_new(t_range)[2,:]/AU
    z_tmp = sol_new(t_range)[3,:]/AU
    lines!(ax, x_tmp, y_tmp, z_tmp)
end



for i in 1:1000
    update_fig()
    sleep(1.)
end

