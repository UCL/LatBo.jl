using FactCheck: facts, context, @fact
using LatBo

facts("something") do

    α = 5

    context("some context") do
        @fact LatBo.dummy() => 5
    end
end
