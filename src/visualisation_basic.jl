####function for 2d plot of pipeline####
function visualise_basic(playpen)

#declare use of Gadfly and DataFrames
using Gadfly
using DataFrames
using LatBo: plot_frame

#create DataFrame
df = plot_frame(play_pen)

#create simple 2D plot using discrete colours
plot(df,x="A",y="B",color="C",Geom.point,Scale.discrete_color())

end