using FactCheck: facts, @fact
using LatBo: D2Q9, D3Q19

facts("Lattice direction inversion") do
    for lattice_name in [:D2Q9, :D3Q19]
        context("for $(lattice_name)") do
            lattice = eval(lattice_name)
            ncomponents = size(lattice.celerities, 2)

            @fact lattice.inversion .>= 1 => all
            @fact lattice.inversion .<= ncomponents => all
            @fact length(unique(lattice.inversion)) => ncomponents
            @fact lattice.celerities[:, lattice.inversion] => -lattice.celerities
        end
    end
end
