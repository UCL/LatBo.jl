module parallel
using Graphs: adjlist, add_edge!
using LatBo.Indices
using LatBo.LB
using LatBo.Playground
function fluid_graph(indexing::Indices.Indexing, lattice::LB.Lattice,
    playground::Array; fluid=Playground.FLUID)
  fluids = zeros(Int, length(indexing))
  count = 1
  for (i, feature) in enumerate(playground)
    if feature == fluid
      fluids[i] = count
      count += 1
    else
      fluids[i] = 0
    end
  end
  graph = adjlist(filter(x -> playground[x] == fluid, 1:length(indexing)), is_directed=false)
  graph.adjlist = Array{Int64, 1}[Int64[] for i in 1:length(graph.vertices)]

  for (site, index) in enumerate(fluids)
    if index != 0
      for direction in 1:length(lattice)
        const to = neighbor_index(indexing, site, lattice, direction)
        if site < to && fluids[to] != 0
          add_edge!(graph, site, to)
        end
      end
    end
  end
  graph
end
end
