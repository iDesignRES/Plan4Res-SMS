#
#
# Check some plants and the Redispatch

function plot_redispatch_results(sumP=nothing)

    if isnothing(sumP) 
        sumP = true
    end
    #println("the flag is ", sumP)

    # Import the needed data
    #
    github_local_d = string(@__DIR__, "/../..")
    github_smsppin = "/smspp_in"
    github_smsppout = "/smspp_out"
    #
    # File names as put together by NTNU
    #
    ntnu_gen_file = string(github_local_d, "/input_data/", "generation.csv") #"Generators.csv"

    # Current Results Name
    res_type = "srdacopf"
    res_name = string("results_", res_type)

    with_selfpomatwo = true #Consider the p4r computation as "POMATWO" -> consistency for Hydro mostly

    # Load the generator data
    gen_data = CSV.read(ntnu_gen_file, DataFrame; delim=',')
    # Do some cleaning on the names
    # if The generator data has a unit_id field, then we can replace the names with that one
    if ("unit_id" in names(gen_data))
        sbgd = filter(row -> (ismissing(row.name)), gen_data)
        for uid in sbgd.unit_id
            #println(uid)
            i0 = findall(gen_data.unit_id .== uid)
            gen_data.name[i0[1]] = uid
        end
    end
    #

    # Import Pomatwo stuff
    if (with_selfpomatwo)
        self_pom_folder = "results_da_pomatwo"
        r_dir = string(github_local_d, github_smsppout, "/nutsx/", self_pom_folder)
        outfile_ext = "OUT"
        res_pomatwo = CSV.read(string(r_dir, "/ActivePower/ActivePower", outfile_ext, ".csv"), DataFrame; delim=',')
    else
        # Refer to the real pomatwo
        pomatwo_extra_extension = "_inflow_restriction"
        pomatwo_dir = string(github_local_d, "/../POMATWO/results/dayahead")
        pomatwo_res = string(pomatwo_dir, "/pomatwo_DA_results_GEN", pomatwo_extra_extension, ".csv")
        res_pomatwo_t = CSV.read(pomatwo_res, DataFrame; delim=',')
        res_pomatwo = unstack(res_pomatwo_t, :Time, :index, :GEN)
    end

    # Import output
    r_dir = string(github_local_d, github_smsppout, "/nutsx/", res_name)
    outfile_ext = "OUT"
    ts_prod = CSV.read(string(r_dir, "/ActivePower/ActivePower", outfile_ext, ".csv"), DataFrame; delim=',')

    nbT = size(ts_prod)[1]
    dt = collect(1:nbT)
    nbperLine = 7

    #Makie.inline!(true)
    fig = Figure(size=(700, 700))
    ga = fig[1, 1]

    global line_number = 1
    Tlist = unique(gen_data.technology)

    for t in Tlist
        I = findall(gen_data.technology .== t)
        genIds = gen_data.unit_id[I]

        dev_prod = Vector{Float64}(undef, nbT)
        dev_prod .= 0.0

        for id in genIds
            if (sumP)
                dev_prod += (ts_prod[:, id] - res_pomatwo[:, id])
            else
                dev_prod += abs.(ts_prod[:, id] - res_pomatwo[:, id])
            end
        end
        if (sumP)
            dev_prod = abs.(dev_prod)
        end
        println(string("techno", t, maximum(dev_prod)))
        axis = Axis(ga[mod(line_number - 1, nbperLine)+1, div(line_number - 1, nbperLine)+1], ylabel=string(t), xlabel="Time Steps", xticks=1:nbT)
        scatterlines!(axis, 1:nbT, dev_prod)
        global line_number += 1
    end
    display(fig)

    if (sumP)
        fig_name = string(r_dir, "/", "redispatch_abssum_", res_type, ".png")
    else
        fig_name = string(r_dir, "/", "redispatch_sumabs_", res_type, ".png")
    end
    save(fig_name, fig)

end

"""
    Plot redispatch production evolution

"""
function plot_redispatch_evol()

    cVector=["blue", "lightblue", "cyan", "burlywood1", "sandybrown", "chocolate1", "brown1", "burlywood4", "sienna4", "grey84", "grey80", "grey75", "grey70", "grey66", "grey61", "grey56", "grey51","grey47","grey42","grey38","grey33","grey28","grey23","grey18","grey13", "grey8", "grey3"]

    # Import the needed data
    #
    github_local_d = string(@__DIR__, "/../..")
    github_smsppin = "/smspp_in"
    github_smsppout = "/smspp_out"
    #
    # File names as put together by NTNU
    #
    ntnu_gen_file = string(github_local_d, "/input_data/", "generation.csv") #"Generators.csv"

    res_list = ["srdacopf", "srdacopf_id_g1", "srdacopf_id_g2"]
    f_idx = [1,1,1]
    naam_list = ["DA", "id_g1", "id_g2"]
    for i=1:24
        append!(res_list, [string("srdacopf_ij_",i)])
        append!(f_idx, [i])
        append!(naam_list, [string("c",lpad(i-1,2,"0"))])
    end

    # Load the generator data
    gen_data = CSV.read(ntnu_gen_file, DataFrame; delim=',')
    # Do some cleaning on the names
    # if The generator data has a unit_id field, then we can replace the names with that one
    if ("unit_id" in names(gen_data))
        sbgd = filter(row -> (ismissing(row.name)), gen_data)
        for uid in sbgd.unit_id
            #println(uid)
            i0 = findall(gen_data.unit_id .== uid)
            gen_data.name[i0[1]] = uid
        end
    end

   nbperLine = 7

    #Makie.inline!(true)
    fig = Figure(size=(1500, 1000))
    ga = fig[1, 1]

    global line_number = 1
    Tlist = unique(gen_data.technology)
    setdiff!(Tlist, ["Offshore floating"])

    global legaxis
    for t in Tlist
        I = findall(gen_data.technology .== t)
        genIds = gen_data.unit_id[I]

        # We will do multiple comparisons
        axis = Axis(ga[mod(line_number - 1, nbperLine)+1, div(line_number - 1, nbperLine)+1], ylabel=string(t), xlabel="Time Steps", xticks=1:24)
        legaxis = axis
        for kk=1:length(res_list)-1
            # Current Results Name
            res_name1 = string("results_", res_list[kk])
            res_name2 = string("results_", res_list[kk+1])

            # Load the results
            # Import output
            r_dir = string(github_local_d, github_smsppout, "/nutsx/", res_name1)
            outfile_ext = "OUT"
            ts_prd1 = CSV.read(string(r_dir, "/ActivePower/ActivePower", outfile_ext, ".csv"), DataFrame; delim=',')

            r_dir = string(github_local_d, github_smsppout, "/nutsx/", res_name2)
            outfile_ext = "OUT"
            ts_prd2 = CSV.read(string(r_dir, "/ActivePower/ActivePower", outfile_ext, ".csv"), DataFrame; delim=',')

            nbT = size(ts_prd1)[1]
            dt = collect(1:nbT)
         
            dev_prod = Vector{Float64}(undef, nbT)
            dev_prod .= 0.0

            for id in genIds
                dev_prod += (ts_prd2[:, id] - ts_prd1[:, id])
            end
            #println(string("techno", t, maximum(dev_prod)))
            
            scatterlines!(axis, f_idx[kk]:nbT, dev_prod[f_idx[kk]:end], color=cVector[kk], label=string("Δ(",naam_list[kk+1]," - ", naam_list[kk],")") )

        end
        global line_number += 1
    end
    fig[1,2] = Legend(fig, legaxis, "Legend", framevisible = false)

    display(fig)
    fig_name = string("../../", "redispatch_by_tech", ".png")
    save(fig_name, fig)

end




"""
Compute the redispatch cost given two files to compare

    base_type : the base folder (POMATWO or self POMATWO result)
    rtype     : the redispatch result
    thf_tab   : thermal unit table

"""
function compute_redispatch_cost(base_type, rtype, thf_tab )
    # Recover the relevant directory structure
    #
    github_local_d = string(@__DIR__, "/../..")
    github_smsppin = "/smspp_in"
    github_smsppout = "/smspp_out"

    # Current Results Name
    base_name = string("results_", base_type)
    # Redispatch Results name
    res_name = string("results_", rtype)

    # Import output
    rb_dir = string(github_local_d, github_smsppout, "/nutsx/", base_name)
    outfile_ext = "OUT"
    tsb_prod = CSV.read(string(rb_dir, "/ActivePower/ActivePower", outfile_ext, ".csv"), DataFrame; delim=',')
   
    r_dir = string(github_local_d, github_smsppout, "/nutsx/", res_name)
    outfile_ext = "OUT"
    tsr_prod = CSV.read(string(r_dir, "/ActivePower/ActivePower", outfile_ext, ".csv"), DataFrame; delim=',')

    if ( size(tsb_prod) != size(tsr_prod ) )
        error("Mismatch in results : no comparison is possible")
    end

    total_slack_vol = 0.0;
    total_thf_rcost = 0.0;

    t1 = 0.0;
    t2 = 0.0;
    for uid in names(tsb_prod[:,2:end])
        #println(uid)
        I = findall( thf_tab.Name .== uid )
        if ( !isempty(I) )
            total_thf_rcost += thf_tab.VariableCost[I[1]]*sum( tsr_prod[:,uid] - tsb_prod[:,uid] )

            t1 += thf_tab.VariableCost[I[1]]*sum( tsb_prod[:,uid] )
            t2 += thf_tab.VariableCost[I[1]]*sum( tsr_prod[:,uid] )
            #println("Unit ", uid, "is thermal")
        end
        if ( occursin("SlackUnit_", uid) )
            #println("Unit ", uid, "is slack")
            total_slack_vol += sum( tsr_prod[:,uid] - tsb_prod[:,uid] )
        end
    end
    println( " costs ", t1, " ", t2 )

    return (total_thf_rcost, total_slack_vol)
   
end

function total_prod_changes(base_type, rtype )
    # Recover the relevant directory structure
    #
    github_local_d = string(@__DIR__, "/../..")
    github_smsppin = "/smspp_in"
    github_smsppout = "/smspp_out"

    # Current Results Name
    base_name = string("results_", base_type)
    # Redispatch Results name
    res_name = string("results_", rtype)

    # Import output
    rb_dir = string(github_local_d, github_smsppout, "/nutsx/", base_name)
    outfile_ext = "OUT"
    tsb_prod = CSV.read(string(rb_dir, "/ActivePower/ActivePower", outfile_ext, ".csv"), DataFrame; delim=',')
   
    r_dir = string(github_local_d, github_smsppout, "/nutsx/", res_name)
    outfile_ext = "OUT"
    tsr_prod = CSV.read(string(r_dir, "/ActivePower/ActivePower", outfile_ext, ".csv"), DataFrame; delim=',')

    if ( size(tsb_prod) != size(tsr_prod ) )
        error("Mismatch in results : no comparison is possible")
    end

    nbT = size(tsb_prod)[1]
    totDelta = Vector{Float64}(undef, nbT)
    totDelta .= 0.0
    for uid in names(tsb_prod[:,2:end])
        if ( ! occursin("SlackUnit_", uid) )
            totDelta += tsr_prod[:,uid] - tsb_prod[:,uid]
        end
    end
    return totDelta
end