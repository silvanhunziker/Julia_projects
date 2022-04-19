#=
test
=#
#=
#using LinearAlgebra
using GLMakie
GLMakie.activate!()


fig = Figure(resolution=(1600, 1600))
ax = Axis(fig[1, 1], aspect = 1)
hidedecorations!(ax)
scatter!([0, 1], [0, 1]; markersize=50)
fig

scatter!([1, 10], [0, 1]; markersize=50)
fig
=#

using GLMakie
GLMakie.activate!()


fig = Figure(resolution = (1000, 1000))

ax = Axis(fig[1, 1], aspect = 1, limits = (-1, 1, -1, 1))
hidedecorations!(ax)

points = Node(Point2f0[])
colors = Node(Bool[])

function draw!()
    p = rand(Point2f0) .* 2 .- 1
    c = sqrt(p[1] ^ 2 + p[2] ^ 2) <= 1
    push!(points[], p)
    push!(colors[], c)
end

sc = scatter!(points, color = colors, colormap = [:red, :blue],
    markersize = 2, strokewidth = 0)

display(fig)

for i in 1:100
    for i in 1:1000
        draw!()
    end
    notify(points)
    notify(colors)
    sleep(1/30)
end

#=
using Random
using GLMakie
GLMakie.activate!()

function scatters_in_3D()
    Random.seed!(123)
    xyz = randn(10, 3)
    x, y, z = xyz[:, 1], xyz[:, 2], xyz[:, 3]
    fig = Figure(resolution=(1600, 400))
    ax1 = Axis3(fig[1, 1]; aspect=(1, 1, 1), perspectiveness=0.5)
    ax2 = Axis3(fig[1, 2]; aspect=(1, 1, 1), perspectiveness=0.5)
    ax3 = Axis3(fig[1, 3]; aspect=:data, perspectiveness=0.5)
    scatter!(ax1, x, y, z; markersize=50)
    meshscatter!(ax2, x, y, z; markersize=0.25)
    hm = meshscatter!(ax3, x, y, z; markersize=0.25,
        marker=FRect3D(Vec3f(0), Vec3f(1)), color=1:size(xyz)[2],
        colormap=:plasma, transparency=false)
    Colorbar(fig[1, 4], hm, label="values", height=Relative(0.5))
    fig
end

scatters_in_3D()
=#