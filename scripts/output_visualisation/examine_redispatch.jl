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
    res_type = "srdacopf_d2"
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
        self_pom_folder = "results_d2_da_pomatwo"
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
        println(string("techno: ", t, " ", maximum(dev_prod), " Σ = ", sum(dev_prod)))
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

    # Switch the comparison of the market computations vs. redispatch ones
    res_list = ["srdacopf", "srdacopf_id_g1", "srdacopf_id_g2"]
    #res_list = ["da_pomatwo", "id_g1_pomatwo", "id_g2_pomatwo"]
    f_idx   = [1,1,1]  #first index of result to show
    l_idx   = [24,24,24]
    ref_idx = [-1,1,2] #the index with which to compare the result
    naam_list = ["DA", "id_g1", "id_g2"]
    for i=1:24
        append!(res_list, [string("srdacopf_ij_",i)])
        #append!(res_list, [string("ij_",i, "_pomatwo")])
        append!(f_idx, [i])
        append!(l_idx, [i])
        append!(ref_idx, [3]) #we will always compare the continuous market with g2
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
        for kk=2:length(res_list)
            # Current Results Name
            #res_name1 = string("results_", res_list[kk])
            res_name1 = string("results_", res_list[ref_idx[kk]])
            res_name2 = string("results_", res_list[kk])

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
            
            scatterlines!(axis, f_idx[kk]:l_idx[kk], dev_prod[f_idx[kk]:l_idx[kk]], color=cVector[ref_idx[kk]], label=string("Δ(",naam_list[kk]," - ", naam_list[ref_idx[kk]],")") )

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

function total_tech_generation( base_type )
   # Recover the relevant directory structure
    #
    github_local_d = string(@__DIR__, "/../..")
    github_smsppin = "/smspp_in"
    github_smsppout = "/smspp_out"

    # Current Results Name
    base_name = string("results_", base_type)

    # Import output
    rb_dir = string(github_local_d, github_smsppout, "/nutsx/", base_name)
    outfile_ext = "OUT"
    tsb_prod = CSV.read(string(rb_dir, "/ActivePower/ActivePower", outfile_ext, ".csv"), DataFrame; delim=',')
  
    # File names as put together by NTNU
    #
    ntnu_gen_file = string(github_local_d, "/input_data/", "generation.csv") #"Generators.csv"

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

    nbT = size(tsb_prod)[1]

    Tlist = unique(gen_data.technology)
    setdiff!(Tlist, ["Offshore floating"])

    res_data = DataFrame( TimeStamp=Vector{Int}(undef, nbT) )
    res_data[!,:TimeStamp] .= tsb_prod[:,1]

    for t in Tlist
        I = findall(gen_data.technology .== t)
        genIds = gen_data.unit_id[I]

        dev_prod_p = Vector{Float64}(undef, nbT)
        dev_prod_p .= 0.0

        dev_prod_m = Vector{Float64}(undef, nbT)
        dev_prod_m .= 0.0

        for id in genIds
            dev_prod_p += max.(tsb_prod[:,id],0.0)
            dev_prod_m += min.(tsb_prod[:,id],0.0)
        end
        if ( minimum( dev_prod_m) >= -1e-6 )
            # There is in fact only one thing
            insertcols!(res_data, t=>Vector{Float64}(undef, nbT) )
            res_data[!, t] .= dev_prod_p
        else
            insertcols!(res_data, string(t,"_p")=>Vector{Float64}(undef, nbT) )
            res_data[!, string(t,"_p")] .= dev_prod_p
            insertcols!(res_data, string(t,"_m")=>Vector{Float64}(undef, nbT) )
            res_data[!, string(t,"_m")] .= dev_prod_m
        end        
    end

    return res_data

end

function use_of_id_markets()
    #
    # Compute redispatch balance volumes vs balance after redispatch
    #
    github_local_d = string(@__DIR__, "/../..")
    github_smsppin = "/smspp_in"
    github_smsppout = "/smspp_out"
    outfile_ext = "OUT"

    # DA productions
    da_name = "results_srdacopf"
    r_dir = string(github_local_d, github_smsppout, "/nutsx/", da_name)
    ts_da = CSV.read(string(r_dir, "/ActivePower/ActivePower", outfile_ext, ".csv"), DataFrame; delim=',')

    # Balance Results
    bal_name = "results_srdacopf_bal_vs_da"
    rb_dir = string(github_local_d, github_smsppout, "/nutsx/", bal_name)
    ts_bal = CSV.read(string(rb_dir, "/ActivePower/ActivePower", outfile_ext, ".csv"), DataFrame; delim=',')
    
    nbT = size(ts_da)[1]
    unit_names = names(ts_da[:,2:end])

    t_id_prod = similar(ts_bal,nbT)
    t_bal_pid = similar(ts_bal,nbT)
    # Establish the intra day schedule
    for i=1:24
        base_name = string("results_", string("srdacopf_ij_",i))
        r_dir = string(github_local_d, github_smsppout, "/nutsx/", base_name)

        ts_prod = CSV.read(string(r_dir, "/ActivePower/ActivePower", outfile_ext, ".csv"), DataFrame; delim=',')

        if ( size( ts_prod ) != size(ts_bal) )
            println(string("incorrect size for id ", i))
        end
        for uid in unit_names
            if ( !occursin("IMP", uid) )
                #println(uid)
                #somehow columns are not always in the same order
                t_id_prod[i,uid] = ts_prod[i,uid]
                if ( i > 1 )
                    t_bal_pid[i-1,uid] = ts_prod[i-1,uid] #at the ith hour, the i-1 th hour has passed into the balance market
                end
                # For proxy purposes we let the last hour of the post id balance, be the id
                t_bal_pid[i,uid] = ts_prod[i,uid]
            end
        end
    end
    
    # We can now compare volumes
    totDAbal = Vector{Float64}(undef, nbT)
    totDAbal .= 0.0
    totIDbal = Vector{Float64}(undef, nbT)
    totIDbal .= 0.0
    for uid in names(ts_da[:,2:end])
        if ( ! occursin("SlackUnit_", uid) && !occursin("IMP", uid) )
            totDAbal += abs.(ts_bal[:,uid] - ts_da[:,uid])
            totIDbal += abs.(t_bal_pid[:,uid] - t_id_prod[:,uid])
        end
    end

    fig = Figure(size=(700, 700))
    ga = fig[1, 1]
    axis = Axis(ga, ylabel=string("Volumes"), xlabel="Time Steps", xticks=1:nbT)
    scatterlines!(axis, 1:nbT, totDAbal, color=:red, label="Balance volumes without ID")
    scatterlines!(axis, 1:nbT, totIDbal, color=:blue, label="Balance volumes with ID")
    #fig[2,1] = 
    #Legend(fig, axis, "Legend", framevisible = false)
    axislegend(axis, merge = true, unique = true, position=:lt)

    display(fig)
    
    fig_name = string(github_local_d, github_smsppout, "/nutsx/", "balancing_volumes", ".png")
    save(fig_name, fig)
end