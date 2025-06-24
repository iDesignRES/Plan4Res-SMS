module iDesignRES_p4r

    # Importer les libs
    using CSV
    using DataFrames
    using Dates
    using GLMakie
    #using Statistics
    
    # Export some function
    export read_base_data
    export make_generator_data
    export make_generator_baselists

    export plot_redispatch_results
    export plot_redispatch_evol
    export compute_redispatch_cost
    export total_prod_changes
    export total_tech_generation
    export use_of_id_markets

    include("../input_handling/import_base_data.jl")
    include("../output_visualisation/examine_redispatch.jl")

    function my_main()
        #plot_redispatch_results()
        zz = 5
    end



    function plotVolumetrics()
        fig = Figure(size=(700, 700))
        ga  = fig[1, 1]
        nbT = 24

        axis = Axis(ga[1,1], ylabel="Deviation (MW)", xlabel="Time Steps", xticks=1:nbT)
        
        cVector=["burlywood1", "sandybrown", "chocolate1", "brown1", "burlywood4", "sienna4", "grey84", "grey80", "grey75", "grey70", "grey66", "grey61", "grey56", "grey51","grey47","grey42","grey38","grey33","grey28","grey23","grey18","grey13", "grey8", "grey3"]

        devP = total_prod_changes("id_g1_pomatwo", "srdacopf_id_g1")
        scatterlines!(axis, 1:nbT, devP, color=:blue, label="Δ(g1(r)-g1(m)")
        devP = total_prod_changes("id_g2_pomatwo", "srdacopf_id_g2")
        scatterlines!(axis, 1:nbT, devP, color=:lightblue, label="Δ(g2(r)-g2(m))")
        
        #devP = total_prod_changes("srdacopf", "srdacopf_ij_1")
        #scatterlines!(axis, 1:nbT, devP, color=:cyan, label="Δ(c00-DA)")

        #devP = total_prod_changes("srdacopf", "srdacopf_ij_2")
        #scatterlines!(axis, 1:nbT, devP, color=:violet, label="Δ(c01-DA)")

        #devP = total_prod_changes("srdacopf", "srdacopf_ij_5")
        #scatterlines!(axis, 1:nbT, devP, color=:orange, label="Δ(c04-DA)")
        #devP = total_prod_changes("srdacopf", "srdacopf_ij_23")
        #scatterlines!(axis, 1:nbT, devP, color=:red, label="Δ(c22-DA)")


        for i=1:24
            devP = total_prod_changes(string("ij_",i,"_pomatwo"), string("srdacopf_ij_",i))
            scatterlines!(axis, i:nbT, devP[i:end], color=cVector[i], label=string("Δ(c",lpad(i-1,2,"0"),"(r) - c",lpad(i-1,2,"0"),"(m))") )
        end
        #else
        #    devP = total_prod_changes("da_pomatwo", "id_g1_pomatwo")
        #    scatterlines!(axis, 1:nbT, devP, color=:blue, label="Δ(g1-DA)")
        #    devP = total_prod_changes("da_pomatwo", "id_g2_pomatwo")
        #    scatterlines!(axis, 1:nbT, devP, color=:lightblue, label="Δ(g2-DA)")
    
        #    for i=1:24
        #        devP = total_prod_changes("da_pomatwo", string("ij_",i, "_pomatwo"))
        #        scatterlines!(axis, i:nbT, devP[i:end], color=:orange, label=string("Δ(c",lpad(i-1,2,"0"),"-DA)") )
        #    end
        #end
        fig[1, 2] = Legend(fig, axis, "Deviations (Redispatch)", framevisible = false)

        display(fig)

        fig_name = string("../../", "redispatch_volumes2", ".png")
        save(fig_name, fig)

    end
    
end