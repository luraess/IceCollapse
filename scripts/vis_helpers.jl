@views function update_error_plots!(ax,err_plots,iter_evo,errs_evo)
    for iplt in eachindex(err_plots.lines)
        err_plots.lines[iplt][1]   = Point2.(iter_evo,errs_evo[iplt,:])
        err_plots.scatter[iplt][1] = Point2.(iter_evo,errs_evo[iplt,:])
    end
    autolimits!(ax)
    return
end

@views function update_heatmaps!(hmaps,fields)
    for prop in propertynames(hmaps)
        hm    = getproperty(hmaps,prop)
        if prop == :η
            hm[3] = log10.(getproperty(fields,prop))
        elseif prop == :phase
            hm[3] = getproperty(fields,prop)
        else
            hm[3] = getproperty(fields,prop).*fields.phase
        end
    end
    return
end