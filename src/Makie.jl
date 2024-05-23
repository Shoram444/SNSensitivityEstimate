function dao(basewidth = 2, spinewidth = basewidth * 0.75) 
	dao_attr = MakieCore.Attributes(
        linewidth = basewidth,
		
        Axis = MakieCore.Attributes(
            spinewidth = spinewidth,
            bottomspinevisible = true,
            topspinevisible = true,
            leftspinevisible = true,
            rightspinevisible = true,
            xgridwidth = basewidth,
            ygridwidth = basewidth,
            xgridstyle = :dot,
            ygridstyle = :dot,
			minorgrid= true,
			xminorgridstyle= :dash,
			yminorgridstyle= :dash,
			xautolimitmargin = (0.0f0, 0.0f0),
			yautolimitmargin = (0.0f0, 0.05f0),
        ),
        Legend = MakieCore.Attributes(
            backgroundcolor = :transparent,
            framevisible = true,
        ),
        Scatter = MakieCore.Attributes(
			linecolor = ["#65ADC2","#111111",  "#233B43", "#E84646", "#C29365", "#362C21", "#316675", "#168E7F", "#109B37"],
            markersize = 15,
            strokewidth = 0,
        ),
		Heatmap = MakieCore.Attributes(
			colormap = :viridis,
		),
		
		)
    # return MakieCore.merge(MakieCore.theme_latexfonts(),dao_attr)
    return dao_attr
end

function Makie.plot(ch::Chains)
	fig = Figure(size = (1200, 400), fontsize= 22, fonts = (; regular = "TeX"))
	for (ind, param) in enumerate(ch.name_map.parameters)
		ax = Axis(fig[1, ind], title=string(param), xticklabelrotation = 45, xlabel ="sample value", ylabel = "counts" )
		for (ind2, datavec) in enumerate(eachcol(getindex(ch, param).data))
			# Get current default colorpalette
			colors = Makie.current_default_theme().attributes[:palette][][:color][]
			Makie.stephist!(ax, datavec, color=(colors[ind2], 1.0), bins = 30,
				linewidth = 2,
			)
			Makie.xlims!(ax, extrema(datavec))
		end
	end
	return fig
end