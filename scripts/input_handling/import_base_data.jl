"""
Read the basic data from the ntnu files and output the corresponding DataFrames

Usage read_base_data(bus_file, line_file, gen_file, load_file)
    with 
        bus_file  : file to the bus data
        line_file : file to the data on the power lines
        gen_file  : file to the data regarding the generators
        load_file : file to the data regarding the load

    output : 
        Collection of data frames with the following tuples:
        (busdata=bus_data, linedata=lines_data, gdata=gen_data, loaddata=load_data)
"""
function read_base_data(bus_file, line_file, gen_file, load_file)

    # Load the data of the buses
    #
    bus_data = CSV.read(bus_file, DataFrame; delim=',')

    # Load the data of the lines
    #
    lines_data = CSV.read(line_file, DataFrame; delim=',')

    # Check if we have multiple lines with the same name
    # 
    if (size(lines_data)[1] != length(unique(lines_data.line_id)))
        for linid in lines_data.line_id
            if (length(findall(lines_data.line_id .== linid)) > 1)
                println(string("Line with name ", linid, " has duplicates"))
            end
        end
        error("Duplicate lines found")
    end

    # Load the generator data
    #
    gen_data = CSV.read(gen_file, DataFrame; delim=',')

    # if The generator data has a unit_id field, then we can replace the names with that one
    if ("unit_id" in names(gen_data))
        sbgd = filter(row -> (ismissing(row.name)), gen_data)
        for uid in sbgd.unit_id
            #println(uid)
            i0 = findall(gen_data.unit_id .== uid)
            gen_data.name[i0[1]] = uid
        end
    end

    # Check if the names in the name column are unique
    #if ( size(gen_data)[1] != length(unique(gen_data.name)) )    
    #    for naam in unique(gen_data.name)
    #        if (length( findall( gen_data.name .== naam ))>1 )
    #            println(string("Plant with name ", naam, " has duplicates"))
    #        end
    #    end
    #    error("Duplicate names found")
    #end

    # Do some cleaning on the names
    new_names = replace.(gen_data[:, :name], "Á" => "A", "É" => "E", "Í" => "I", "Ñ" => "N", "Ó" => "O", "Ú" => "U", "Ü" => "U")
    gen_data[:, :name] .= new_names

    # Read the load data
    #
    load_data = CSV.read(load_file, DataFrame; delim=',')

    return (busdata=bus_data, linedata=lines_data, gdata=gen_data, loaddata=load_data)
end

"""
Set up the base lists for various kinds of data
    st_idx   : the index of the relevant stochastic stage
"""
function make_generator_baselists(st_idx)
    ##
    ## Underlying base data
    ##
    Thf_list = ["Coal", "Gas", "Biomass", "Waste", "Nuclear", "Oil", "Oil/Diesel", "Diesel", "Oil/Gas/Diesel", "Oil/Gas"]
    Res_list = ["Wind", "Solar", "Hydro|run_of_river"]
    STS_list = ["Hydro|pumped_storage"]
    SS_list = ["Hydro|reservoir"]

    # List of technology related costs and information (prop cost, fixed cost, InvestmentCost)
    #
    cThf_l = Dict([("Coal", (13.05, 42380, 1695200)), ("Gas", (31.24, 18011.5, 663122.3529)),
        ("Biomass", (64.022, 55242.33, 2511015)), ("Waste", (62.00, 55242.33, 2511015)),
        ("Nuclear", (12.42, 100122.75, 6357000)), ("Oil", (78.92, 6897.345, 688675)),
        ("Oil/Diesel", (85.00, 6500.345, 688675)), ("Diesel", (90.00, 6500.345, 688675)),
        ("Oil/Gas/Diesel", (83.00, 6500.345, 688675)), ("Oil/Gas", (82.00, 6500.345, 688675))])

    # List of technology related information 
    # the last indicates the column of the res load factor => -1 indicates necessary use of the TS file 
    #
    cRes_l = Dict([("Wind", (1059500, "EDF__WindOnshore-LoadFactor-PresentClimate-", "__13082019__13082019__v1.csv", 1.0, 1)),
        ("Solar", (476775, "EDF__PV-LoadFactor-PresentClimate-", "__13082019__13082019__v1.csv", 1.0, 2)),
        ("Hydro", (0, "EDF__RunOfRiver-HourlyCoefficient-PresentClimate-", "__18092019__18092019__v1.csv", 2500.0, -1)),
        ("Hydro|run_of_river", (0, "EDF__RunOfRiver-HourlyCoefficient-PresentClimate-", "__18092019__18092019__v1.csv", 2500.0, -1))])

    # List of technology related information maxVol, turEff, pumpEff
    cSTS_l = Dict([("Hydro|pumped_storage", (100, 1, 0.866)),
        ("Hydro (Pumped storage with natural inflow)", (100, 1, 0.866))])

    # List of technology related stuff
    cSS_l = Dict([("Hydro|reservoir", (100, "EDF__Inflow-HourlyCoefficient-PresentClimate-", "__18092019__18092019__v1.csv"))])

    # We can rely on vol_data to find out the right percentage
    # st_idx indicates the next stochastic stage, so we are in st_idx - 1
    # println(string("vol init FR : ", string(vol_data[st_idx*ssv_step, 2]/16256280.72), " vol init ES : ", string(vol_data[st_idx*ssv_step, 3]/16409036.88) ))
    fr_per = st_idx > 1 ? round(vol_data[(st_idx-1)*ssv_step, 2] / 16256280.72, digits=3) : 0.3
    es_per = st_idx > 1 ? round(vol_data[(st_idx-1)*ssv_step, 3] / 16409036.88, digits=3) : 0.3
    # maxVol, Hydrosystem name, inflows, %initial volume vs maxvol
    ss_d = Dict([("FR", (16256280.72, 0, 96118939.76, fr_per)), ("ES", (16409036.88, 1, 27075319.09, es_per))])

    return (lists=(Thf_list, Res_list, STS_list, SS_list), dicts=(cThf_l,cRes_l,cSTS_l,cSS_l, ss_d))
end


"""
Set up the generator data of various kinds from the basic information

function make_generator_data(bus_data, gen_data, st_idx)
    with
        bus_data : the data on the whereabouts of the various buses
        gen_data : the base data regarding the generators
        st_idx   : the index of the relevant stochastic stage
    
    output:
        Tables with:
            Data on Thermal generating units
            Data on RES units
            Data on Pumped Storage and Small scale storage
            Data on Seasonal Storage Units 
"""
function make_generator_data(bus_data, gen_data, st_idx)
    blist = make_generator_baselists(st_idx)
    (Thf_list, Res_list, STS_list, SS_list) = blist.lists
    (cThf_l,cRes_l,cSTS_l,cSS_l, ss_d) = blist.dicts

    # Type conversion on the column
    gen_data[!, :primary_fuel] = convert.(String, gen_data[:, :primary_fuel])
    # Some fix on the type of hydro
    for ig = 1:size(gen_data)[1]
        if (gen_data.primary_fuel[ig] == "Hydro")
            #println( gen_data.unit_id[ig], "_", gen_data.technology[ig] )
            gen_data.primary_fuel[ig] = string(gen_data.primary_fuel[ig], "|", gen_data.technology[ig])
        end
    end
    #
    println("Found the following technologies: ", unique(gen_data.primary_fuel))

    #
    # The Thermal units
    #
    nb_thf = 0
    (nT, nC) = size(gen_data)
    for i = 1:nT
        if (gen_data.primary_fuel[i] in Thf_list)
            nb_thf += 1
        end
    end
    # Make the thermal unit file from this
    tu_thf_data = DataFrame(Zone=Vector{String}(undef, nb_thf), Name=Vector{String}(undef, nb_thf), NumberUnits=Vector{Int64}(undef, nb_thf),
        MaxPower=Vector{Float64}(undef, nb_thf), MaxPowerProfile=Vector{String}(undef, nb_thf), VariableCost=Vector{Float64}(undef, nb_thf), FixedCost=Vector{Float64}(undef, nb_thf),
        InvestmentCost=Vector{Float64}(undef, nb_thf), Capacity=Vector{Float64}(undef, nb_thf), Energy=Vector{Float64}(undef, nb_thf),
        MaxAddedCapacity=Vector{Float64}(undef, nb_thf))

    tu_thf_data[!, :NumberUnits] .= 1
    tu_thf_data[!, :Energy] .= 0.0
    tu_thf_data[!, :MaxPowerProfile] .= ""

    i_thf = 0
    for i = 1:nT
        if (gen_data.primary_fuel[i] in Thf_list)
            i_thf += 1
            # Add this fellow - first check if it actually exists at some existing bus
            l0 = length(findall(bus_data.bus_id .== gen_data.bus_id[i]))
            if (l0 == 0)
                error("The generator ", string(i), " is situated at some non existing bus:", string(gen_data.bus_id[i]))
            end
            i0 = findall(bus_data.bus_id .== gen_data.bus_id[i])

            tu_thf_data[i_thf, :Zone] = bus_data.bus_id[i0][1] #string.( bus_data.country[i0], "_", bus_data.bus_id[i0] )[1]
            tu_thf_data[i_thf, :Name] = gen_data.unit_id[i] #gen_data.name[i]

            tu_thf_data[i_thf, :MaxPower] = gen_data.capacity_mw[i]

            cost_info = cThf_l[gen_data.primary_fuel[i]]
            tu_thf_data[i_thf, :VariableCost] = cost_info[1]
            tu_thf_data[i_thf, :FixedCost] = cost_info[2]
            tu_thf_data[i_thf, :InvestmentCost] = cost_info[3]
        end
    end
    tu_thf_data[!, :Capacity] .= tu_thf_data[!, :MaxPower]
    tu_thf_data[!, :MaxAddedCapacity] .= ceil.(0.1 * tu_thf_data[!, :MaxPower]) #round.( 0.1*tu_thf_data[!,:MaxPower]; digits=2 )

    #
    # The RES units
    #
    nb_res = 0
    for i = 1:nT
        if (gen_data.primary_fuel[i] in Res_list)
            nb_res += 1
        end
    end

    # res DataFrame
    res_units_data = DataFrame(Name=Vector{String}(undef, nb_res), Zone=Vector{String}(undef, nb_res), NumberUnits=Vector{Int64}(undef, nb_res),
        MaxPower=Vector{Float64}(undef, nb_res), MinPower=Vector{Float64}(undef, nb_res), MaxPowerProfile=Vector{String}(undef, nb_res),
        Energy=Vector{Float64}(undef, nb_res), Kappa=Vector{Float64}(undef, nb_res), Capacity=Vector{Float64}(undef, nb_res),
        MaxAddedCapacity=Vector{Float64}(undef, nb_res), MaxRetCapacity=Vector{Float64}(undef, nb_res), InvestmentCost=Vector{Float64}(undef, nb_res),
        LoadFactorColumn=Vector{Int64}(undef, nb_res))

    res_units_data[!, :NumberUnits] .= 1
    res_units_data[!, :Kappa] .= 1.0
    res_units_data[!, :Energy] .= 0.0

    i_res = 0
    for i = 1:nT
        if (gen_data.primary_fuel[i] in Res_list)
            i_res += 1
            # Add this fellow - first check if it actually exists at some existing bus
            l0 = length(findall(bus_data.bus_id .== gen_data.bus_id[i]))
            if (l0 == 0)
                error("The generator ", string(i), " is situated at some non existing bus:", string(gen_data.bus_id[i]))
            end
            i0 = findall(bus_data.bus_id .== gen_data.bus_id[i])

            res_units_data[i_res, :Zone] = bus_data.bus_id[i0][1] #string.( bus_data.country[i0], "_", bus_data.bus_id[i0] )[1]
            res_units_data[i_res, :Name] = gen_data.unit_id[i] #gen_data.name[i]

            res_units_data[i_res, :Capacity] = gen_data.capacity_mw[i]

            # Get technology related stuff
            t_info = cRes_l[gen_data.primary_fuel[i]]
            res_units_data[i_res, :MaxPowerProfile] = string(t_info[2], bus_data.country[i0][1], t_info[3])
            res_units_data[i_res, :InvestmentCost] = t_info[1]
            # Moving from Pmax to Energy for Hydro related things
            res_units_data[i_res, :MaxPower] = gen_data.capacity_mw[i] * t_info[4]
            # Indicate the res load factor column
            res_units_data[i_res, :LoadFactorColumn] = t_info[5]
        end
    end
    res_units_data[!, :MinPower] .= 0.0
    res_units_data[!, :MaxRetCapacity] .= 0.0

    #res_units_data[!,:Capacity] .= res_units_data[!,:MaxPower]
    res_units_data[!, :MaxAddedCapacity] .= ceil.(0.1 * res_units_data[!, :MaxPower])

    #
    # Now the Pumped Storage file
    #
    nb_sts = 0
    for i = 1:nT
        if (gen_data.primary_fuel[i] in STS_list)
            nb_sts += 1
        end
    end

    sts_data = DataFrame(Name=Vector{String}(undef, nb_sts), Zone=Vector{String}(undef, nb_sts), NumberUnits=Vector{Int64}(undef, nb_sts),
        MaxPower=Vector{Float64}(undef, nb_sts), MaxVolume=Vector{Float64}(undef, nb_sts), TurbineEfficiency=Vector{Float64}(undef, nb_sts),
        PumpingEfficiency=Vector{Float64}(undef, nb_sts), MinPower=Vector{Float64}(undef, nb_sts), MinVolume=Vector{Float64}(undef, nb_sts),
        Energy=Vector{Float64}(undef, nb_sts), Inflows=Vector{Float64}(undef, nb_sts), InitialVolume=Vector{Float64}(undef, nb_sts),
        AddPumpedStorage=Vector{Float64}(undef, nb_sts),
        MaxAddedCapacity=Vector{Float64}(undef, nb_sts), MaxRetCapacity=Vector{Float64}(undef, nb_sts), InvestmentCost=Vector{Float64}(undef, nb_sts))

    sts_data[!, :NumberUnits] .= 1
    sts_data[!, :MinVolume] .= 0.0
    sts_data[!, :Energy] .= 0.0
    sts_data[!, :Inflows] .= 0.0
    sts_data[!, :InitialVolume] .= 0.0
    sts_data[!, :AddPumpedStorage] .= 0.0
    sts_data[!, :MaxAddedCapacity] .= 0.0
    sts_data[!, :MaxRetCapacity] .= 0.0
    sts_data[!, :InvestmentCost] .= 0.0

    i_sts = 0
    for i = 1:nT
        if (gen_data.primary_fuel[i] in STS_list)
            i_sts += 1
            # Add this fellow - first check if it actually exists at some existing bus
            l0 = length(findall(bus_data.bus_id .== gen_data.bus_id[i]))
            if (l0 == 0)
                error("The STS unit ", string(i), " is situated at some non existing bus:", string(gen_data.bus_id[i]))
            end
            i0 = findall(bus_data.bus_id .== gen_data.bus_id[i])

            sts_data[i_sts, :Zone] = bus_data.bus_id[i0][1] #string.( bus_data.country[i0], "_", bus_data.bus_id[i0] )[1]
            sts_data[i_sts, :Name] = gen_data.unit_id[i] #gen_data.name[i]

            sts_data[i_sts, :MaxPower] = gen_data.capacity_mw[i]
            sts_data[i_sts, :MinPower] = -1.0 * gen_data.capacity_mw[i]

            tech_info = cSTS_l[gen_data.primary_fuel[i]]
            sts_data[i_sts, :MaxVolume] = tech_info[1] * gen_data.capacity_mw[i]
            sts_data[i_sts, :TurbineEfficiency] = tech_info[2]
            sts_data[i_sts, :PumpingEfficiency] = tech_info[3]
        end
    end

    #
    # The Seasonal Storage units
    #
    nb_ss = 0
    nb_ss_max_d = Dict([("FR", 0.0), ("ES", 0.0)])
    for i = 1:nT
        if (gen_data.primary_fuel[i] in SS_list)
            nb_ss += 1
            # first check if it actually exists at some existing bus
            l0 = length(findall(bus_data.bus_id .== gen_data.bus_id[i]))
            if (l0 == 0)
                error("The SS unit ", string(i), " is situated at some non existing bus:", string(gen_data.bus_id[i]))
            end
            i0 = findall(bus_data.bus_id .== gen_data.bus_id[i])

            nb_ss_max_d[bus_data.country[i0][1]] += gen_data.capacity_mw[i]
        end
    end

    ss_data = DataFrame(Name=Vector{String}(undef, nb_ss), Zone=Vector{String}(undef, nb_ss),
        HydroSystem=Vector{Int64}(undef, nb_ss), NumberUnits=Vector{Int64}(undef, nb_ss),
        MaxPower=Vector{Float64}(undef, nb_ss), MinPower=Vector{Float64}(undef, nb_ss),
        MaxVolume=Vector{Float64}(undef, nb_ss), MinVolume=Vector{Float64}(undef, nb_ss),
        Inflows=Vector{Float64}(undef, nb_ss), InflowsProfile=Vector{String}(undef, nb_ss),
        InitialVolume=Vector{Float64}(undef, nb_ss), TurbineEfficiency=Vector{String}(undef, nb_ss), PumpingEfficiency=Vector{String}(undef, nb_ss),
        AddPumpedStorage=Vector{Float64}(undef, nb_ss), WaterValues=Vector{String}(undef, nb_ss))

    ss_data[!, :NumberUnits] .= 1
    ss_data[!, :MinVolume] .= 0.0
    ss_data[!, :MinPower] .= 0.0
    ss_data[!, :TurbineEfficiency] .= 1.0
    ss_data[!, :PumpingEfficiency] .= 0.0
    ss_data[!, :AddPumpedStorage] .= 0.0
    ss_data[!, :WaterValues] .= "bellman_nutsx.csv"

    i_ss = 0
    for i = 1:nT
        if (gen_data.primary_fuel[i] in SS_list)
            i_ss += 1
            # Add this fellow - first check if it actually exists at some existing bus
            l0 = length(findall(bus_data.bus_id .== gen_data.bus_id[i]))
            if (l0 == 0)
                error("The SS unit ", string(i), " is situated at some non existing bus:", string(gen_data.bus_id[i]))
            end
            i0 = findall(bus_data.bus_id .== gen_data.bus_id[i])

            ss_data[i_ss, :Zone] = bus_data.bus_id[i0][1] #string.( bus_data.country[i0], "_", bus_data.bus_id[i0] )[1]
            ss_data[i_ss, :Name] = gen_data.unit_id[i] #gen_data.name[i]

            ss_data[i_ss, :MaxPower] = gen_data.capacity_mw[i]

            tech_info = cSS_l[gen_data.primary_fuel[i]]
            country_info = ss_d[bus_data.country[i0][1]]

            mx_tot = nb_ss_max_d[bus_data.country[i0][1]]

            # We will proportionally dispatch the stuff unto the units
            ss_data[i_ss, :HydroSystem] = country_info[2]
            ss_data[i_ss, :MaxVolume] = country_info[1] * (gen_data.capacity_mw[i] / mx_tot)

            ss_data[i_ss, :Inflows] = country_info[3] * (gen_data.capacity_mw[i] / mx_tot)
            ss_data[i_ss, :InflowsProfile] = string(tech_info[2], bus_data.country[i0][1], tech_info[3])
            ss_data[i_ss, :InitialVolume] = country_info[4] * ss_data[i_ss, :MaxVolume]
        end
    end

    # Check if the all units have been handled
    if (nb_ss + nb_thf + nb_res + nb_sts != nT)
        error("Some units have been lost...")
    end

    return (tu_thf_data, res_units_data, sts_data, ss_data)

end