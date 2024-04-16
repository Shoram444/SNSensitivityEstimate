# This file should contain site-specific commands to be executed on Julia startup;
# Users may store their own personal commands in `~/.julia/config/startup.jl`.

import Pkg
let
    pkgs = ["Revise", "AbbreviatedStackTraces"]
    for pkg in pkgs
    if Base.find_package(pkg) === nothing
        Pkg.add(pkg)
    end
    end
end

using Revise
using AbbreviatedStackTraces

printstyled("using Revise \nusing AbbreviatedStackTraces"; color=:light_blue)

const PLOTS_DEFAULTS = Dict(
	:theme => :dao,
	:size  => (1200, 800),
	:legend => :topleft,
	:guidefontsize  => 16,
	:tickfontsize   => 12,
	:titlefontsize  => 16,
	:legendfontsize => 12,
	:dpi            => 200,
	:colorbar_titlefontsize => 20,
	:thickness_scaling => 1.4,
	:widen => :false,
	:markerstrokewidth => 1,
	:markerstrokecolor => :black,
	:palette => ["#00a0f9", "#ba3030", "#22ac74", "#707070", "#9452bd", "#80ff00", "#ffcc00", "#ff00ff", "#00ffff", "#cc9900"],
	:linewidth => 3,
	:fontfamily => "Serif"
)
