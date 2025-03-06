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
    export compute_redispatch_cost

    include("../input_handling/import_base_data.jl")
    include("../output_visualisation/examine_redispatch.jl")

    function my_main()
        #plot_redispatch_results()
        zz = 5
    end
    
end