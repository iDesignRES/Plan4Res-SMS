using DataFrames
using Dates

include("extract_ts_sequence.jl")

"""
    fname   : the name of the file to which to write
    ts_dir  : where to find the timeseries

    b_date  : the beginning Date
    e_date  : the end Date
    s_idx   : the scenario index if relevant
    st_idx  : Index of the relevant stochastic stage
    t_phi   : the 'tangent' phi value, with tan(phi)=Q/P => Q = tan(phi) P 

    z_data  : the table containing the node names
    zv_data : the table with load information
    ic_data : the table with interconnection data

    thf_data: the table with the thermal generators
    res_data: the table with the RES generators
    sts_data: the table with the pumped storage data
    ss_data : the seasonal storage data
"""
function write_smspp_file(fname, ts_dir, b_date, e_date, s_idx, st_idx, t_phi, z_data, zv_data, ic_data, thf_data, res_data, sts_data, ss_data)
    open(fname, "w") do fic
        # Write the Header part
        println(fic, "netcdf Block_0 {")
        println(fic, "")
        println(fic, "// global attributes:")
        println(fic, "        :SMS++_file_type = 1LL ;")
        println(fic, "")

        # Print some general informat
        nb_thf = size(thf_data)[1]
        nb_res = size(res_data)[1]
        nb_sts = size(sts_data)[1]
        nb_ss  = size(ss_data)[1]
        nb_ss_sys = length(unique(ss_data.HydroSystem))
        nbslack= size(z_data)[1] 
        println(string("Found a total of ", string(nb_thf), " Thermal plants"))
        println(string("Found a total of ", string(nb_res), " RES plants"))
        println(string("Found a total of ", string(nb_sts), " STS plants"))
        println(string("Found a total of ", string(nb_ss), " Seasonal Storage plants, in ", string(nb_ss_sys), " Hydro systems" ) )
        println(string("Found a total of ", string(nbslack), " Slack Units"))
        
        # Open the UC block part
        println(fic, "group: Block_0 {")
        println(fic, "  dimensions:")
        # Dimension stuff       
        println(fic, string("  	TimeHorizon = ",string(convert(Dates.Hour,(e_date - b_date)).value), " ; "))
        println(fic, string("  	NumberUnits = ",string(nb_thf+nb_res+nb_sts+nb_ss_sys + nbslack), " ; ") )
        println(fic, string("  	NumberElectricalGenerators = ",string(nb_thf+nb_res+nb_sts+nb_ss + nbslack), " ; ") )
        println(fic, string("  	NumberNodes = ",string(size(z_data)[1]), " ; "))
        println(fic, string("  	NumberLines = ",string(size(ic_data)[1]), " ; "))

        println(fic, "  variables:")
        # Variables
        println(fic, "  	uint GeneratorNode(NumberElectricalGenerators) ;")
        println(fic, "  	string NodeName(NumberNodes) ; ")
        println(fic, "  	double ActivePowerDemand(NumberNodes, TimeHorizon) ;")        
        # AC OPF related stuff
        if ( "NodeSusceptance" in names(z_data) )
            println(fic, "  	string NetworkBlockClassname ;")           
            println(fic, "  	double ReactivePowerDemand(NumberNodes, TimeHorizon) ;")
            println(fic, "  	double NodeConductance(NumberNodes) ;")
            println(fic, "  	double NodeSusceptance(NumberNodes) ;")
            println(fic, "  	double NodeMinVoltage(NumberNodes) ;")
            println(fic, "  	double NodeMaxVoltage(NumberNodes) ;")
        end
        println(fic, "  	uint StartLine(NumberLines) ;")
        println(fic, "  	uint EndLine(NumberLines) ;")
        println(fic, "  	double MinPowerFlow(NumberLines) ;")
        println(fic, "  	double MaxPowerFlow(NumberLines) ;")
        println(fic, "  	double LineSusceptance(NumberLines) ;")
        println(fic, "  	double NetworkCost(NumberLines) ;")
        println(fic, "  	string LineName(NumberLines) ;")
        # AC OPF related stuff
        if ( "LineReactance" in names(ic_data) )
            println(fic, "  	double LineResistance(NumberLines) ;")
            println(fic, "  	double LineReactance(NumberLines) ;")
            println(fic, "  	double LineMinAngle(NumberLines) ;")
            println(fic, "  	double LineMaxAngle(NumberLines) ;")
            println(fic, "  	double LineShiftAngle(NumberLines) ;") 
            println(fic, "  	double LineRATEA(NumberLines) ;") 
            println(fic, "  	double LineRatio(NumberLines) ;")            
        end

        #
        println(fic, "")
        println(fic, "  // group attributes:")
        println(fic, "  		:id = \"0\" ;")
        println(fic, "  		:type = \"UCBlock\" ;")
        if ( "NodeSusceptance" in names(z_data) )
            println(fic, "  		:baseMVA = \"1.0\" ;")
        end
        println(fic, "  data:")

        println(fic, "")

        if ( "NodeSusceptance" in names(z_data) )
            println(fic, "   NetworkBlockClassname = \"ACNetworkBlock\" ;" )
            println(fic, "")
        end

        # UC block main part
        #
        gnode=Vector{Int64}(undef,0)
        for ithf=1:nb_thf
            i0 = findall( z_data.Countries .== thf_data.Zone[ithf] )
            if ( isempty(i0) )
                error(string(" Thermal Generator ", thf_data.Name[ithf], " located at non existing node ", thf_data.Zone[ithf]))
            end
            append!(gnode,[i0[1]-1])
        end
        for ires=1:nb_res
            i0 = findall( z_data.Countries .== res_data.Zone[ires] )
            if ( isempty(i0) )
                error(string(" RES Generator ", res_data.Name[ires], " located at non existing node ", res_data.Zone[ires]))
            end
            append!(gnode,[i0[1]-1])
        end
        for ists=1:nb_sts
            i0 = findall( z_data.Countries .== sts_data.Zone[ists] )
            if ( isempty(i0) )
                error(string(" STS Generator ", sts_data.Name[ists], " located at non existing node ", sts_data.Zone[ists]))
            end
            append!(gnode,[i0[1]-1])
        end
        for iss=1:nb_ss
            i0 = findall( z_data.Countries .== ss_data.Zone[iss] )
            if ( isempty(i0) )
                error(string(" STS Generator ", ss_data.Name[iss], " located at non existing node ", ss_data.Zone[iss]))
            end
            append!(gnode,[i0[1]-1])
        end
        # Slacks can be done in one go
        append!( gnode, collect(1:nbslack).-1 )
        println(fic, string("   GeneratorNode = ", chop(string(gnode),head=1), "  ;" ) )
        println(fic, "")

        println(fic, string("   NodeName = ", strVtoList(z_data.Countries), " ;"))
        println(fic, "")
        #
        ReacPwrDem = Matrix{Float64}(undef, length(zp_data.Countries), convert(Dates.Hour,(e_date - b_date)).value)
        i_nd = 1
        # for each node
        println(fic, "   ActivePowerDemand = " )
        for nd in z_data.Countries
            i0 = findall( zv_data.Zone .== nd )
            if ( isempty(i0) )
                error(string("Node : ", nd, " not found!"))
            end
            ts_fname = string(ts_dir, "/", zv_data.Profile_Timeserie[i0][1])

            ts_vals  = extract_ts_sequence( ts_fname, b_date, e_date, s_idx )

            println( string("Node ", nd, " timeseries : ", ts_fname) )

            ActiveDemand = round.(ts_vals.*zv_data.value[i0], digits=5)
            ReacPwrDem[i_nd,:] = t_phi*ActiveDemand
            i_nd += 1

            if ( nd == zp_data.Countries[end] )
                println(fic, string("     ", chop(string(ActiveDemand),head=1), "  ;" ) )
            else
                println(fic, string("     ", chop(string(ActiveDemand),head=1), " ," ) )
            end
        end
        println(fic, "")
        # AC OPF related stuff
        if ( "NodeSusceptance" in names(z_data) )
            println(fic, "   ReactivePowerDemand = " )
            for ind=1:length(z_data.Countries)-1
                println(fic, string("     ", chop(string(ReacPwrDem[ind,:]),head=1), " ," ) )
            end
            println(fic, string("     ", chop(string(ReacPwrDem[end,:]),head=1), "  ;" ) )
            println(fic, "")

            println(fic, string("   NodeConductance = ", chop(string(z_data.NodeConductance),head=1), "  ;" ) )
            println(fic, "")
            println(fic, string("   NodeSusceptance = ", chop(string(z_data.NodeSusceptance),head=1), "  ;" ) )
            println(fic, "")
            println(fic, string("   NodeMinVoltage = ", chop(string(z_data.NodeMinVoltage),head=1), "  ;" ) )
            println(fic, "")
            println(fic, string("   NodeMaxVoltage = ", chop(string(z_data.NodeMaxVoltage),head=1), "  ;" ) )
            println(fic, "")
        end

        # 
        # Work out the lines
        sline=Vector{Int64}(undef,0)
        eline=Vector{Int64}(undef,0)
        mnflw=Vector{Float64}(undef,0)
        mxflw=Vector{Float64}(undef,0)
        suscp=Vector{Float64}(undef,0)
        netcs=Vector{Float64}(undef,0)
        # AC OPF line stuff
        resnc=Vector{Float64}(undef,0)
        reacc=Vector{Float64}(undef,0)
        mnang=Vector{Float64}(undef,0)
        mxang=Vector{Float64}(undef,0)
        lsang=Vector{Float64}(undef,0)
        lnratea=Vector{Float64}(undef,0)
        lnratio=Vector{Float64}(undef,0)
        for iln = 1:size(ic_data)[1]
            i0 = findall( z_data.Countries .== ic_data.StartLine[iln] )
            i1 = findall( z_data.Countries .== ic_data.EndLine[iln] )
            if ( isempty(i0) || isempty(i1) )
                error(" End or start bus not found")
            end
            append!(sline,[i0[1] - 1])
            append!(eline,[i1[1] - 1])
            append!(mnflw,ic_data.MinPowerFlow[iln][1])
            append!(mxflw,ic_data.MaxPowerFlow[iln][1])
            if ( "Susceptance" in names(ic_data) )
                append!(suscp,ic_data.Susceptance[iln][1])
            else
                append!(suscp, 0)
            end
            if ( "LineReactance" in names(ic_data) )
                append!(resnc,ic_data.LineResistance[iln][1])
                append!(reacc,ic_data.LineReactance[iln][1])
                append!(mnang,ic_data.LineMinAngle[iln][1])
                append!(mxang,ic_data.LineMaxAngle[iln][1])
                append!(lsang,ic_data.LineShiftAngle[iln][1])
                append!(lnratea,ic_data.LineRATEA[iln][1])
                append!(lnratio,ic_data.LineRatio[iln][1])                
            end
            if ( "NetworkCost" in names(ic_data) )
                append!(netcs,ic_data.NetworkCost[iln][1])
            else
                append!(netcs, 0)
            end
        end
        println(fic, string("   StartLine = ", chop(string(sline),head=1), "  ;" ) )
        println(fic, "")
        println(fic, string("   EndLine = ", chop(string(eline),head=1), "  ;" ) )
        println(fic, "")
        println(fic, string("   MinPowerFlow = ", chop(string(mnflw),head=1), "  ;" ) )
        println(fic, "")
        println(fic, string("   MaxPowerFlow = ", chop(string(mxflw),head=1), "  ;" ) )
        println(fic, "")
        println(fic, string("   LineSusceptance = ", chop(string(suscp),head=1), "  ;" ) )
        println(fic, "")
        if ( "LineReactance" in names(ic_data) )
            println(fic, string("   LineResistance = ", chop(string(resnc),head=1), "  ;" ) )
            println(fic, "")
            println(fic, string("   LineReactance = ", chop(string(reacc),head=1), "  ;" ) )
            println(fic, "")
            println(fic, string("   LineMinAngle = ", chop(string(mnang),head=1), "  ;" ) )
            println(fic, "")
            println(fic, string("   LineMaxAngle = ", chop(string(mxang),head=1), "  ;" ) )
            println(fic, "")
            println(fic, string("   LineShiftAngle = ", chop(string(lsang),head=1), "  ;" ) )
            println(fic, "")      
            println(fic, string("   LineRATEA = ", chop(string(lnratea),head=1), "  ;" ) )
            println(fic, "")
            println(fic, string("   LineRatio = ", chop(string(lnratio),head=1), "  ;" ) )
            println(fic, "")     
        end
        println(fic, string("   NetworkCost = ", chop(string(netcs),head=1), "  ;" ) )
        println(fic, "")
        println(fic, string("   LineName = ", strVtoList(ic_data.Name), " ;"))
        println(fic, "")

        # Hydro System Blocks
        isys = 0
        for hsys in unique(ss_data.HydroSystem)
            isys += 1
            # filter out all fellows belonging to this HSystem
            sub_tab = filter(row -> ( row.HydroSystem .== hsys ), ss_data )
        
            write_smspp_HBlocks(fic,ts_dir, isys - 1, s_idx, st_idx, t_phi, sub_tab, b_date, e_date)
        end        
        # Now write all Thermal unit blocks
        for ithf=1:nb_thf
            write_smspp_TBlocks(fic, ts_dir, nb_ss_sys + ithf - 1, t_phi, thf_data[ithf,:], b_date, e_date )
        end
        # Now write all the RES blocks
        for ires=1:nb_res
            write_smspp_RBlocks(fic, ts_dir, nb_ss_sys + nb_thf + ires - 1, s_idx, t_phi, res_data[ires,:], b_date, e_date )
        end
        # Now write all the Battery / STS blocks
        for ists=1:nb_sts
            write_smspp_BatBlocks(fic, nb_ss_sys + nb_thf + nb_res + ists - 1, t_phi, sts_data[ists,:], b_date, e_date )
        end
        # For each bus write a slack block
        for isl=1:nbslack
            naam = string("SlackUnit_", string(z_data.Countries[isl]))
            write_smspp_SlackBlocks(fic, nb_ss_sys + nb_thf + nb_res + nb_sts + isl - 1, naam, b_date, e_date )
        end

        #
        println(fic, "")
        println(fic, "  } // group Block_0")
        println(fic, "}")
    end
end

"""
    fic   : the file handler
    ts_dir: directory to the timeseries
    
    b_idx  : the block index
    s_idx  : the scenario index
    st_idx : the stage index
    t_phi  : the 'tangent' phi value, with tan(phi)=Q/P => Q = tan(phi) P 
    tab    : the data frame related to this HSystem

    b_date : begin Date
    e_date : end Date
"""
function write_smspp_HBlocks(fic,ts_dir, b_idx, s_idx, st_idx, t_phi, tab, b_date, e_date)
    println(fic,string("  group: UnitBlock_", string(b_idx)," { "))
    println(fic,"    dimensions:")
    println(fic,string("    	NumberHydroUnits = ", string(size(tab)[1]), " ; "))
    println(fic,"")
    println(fic,"    // group attributes:")
    println(fic,"    		:type = \"HydroSystemUnitBlock\" ;")
    println(fic,"")
    # Now write each of the units in this Block
    nb_tab = size(tab)[1]
    for itab=1:nb_tab
        write_smspp_HunitBlocks(fic, ts_dir, itab - 1, s_idx, t_phi, tab[itab,:], b_date, e_date)
    end
    wvalname=""
    if ( "WaterValues" in names(tab) )
        wvfile = unique(tab.WaterValues)
        if (length(wvfile)>1)
            error("Each Hydrosystem should have a unique Bellman file")
        end
        wvalname = string(ts_dir, "/", wvfile[1] )
    end
    # Write the PolyhedralFunctionBlock if needed ...
    write_smspp_PolyBlock(fic, wvalname, st_idx)
    
    println(fic,string("    } // group UnitBlock_",string(b_idx)))
    println(fic,"")
end


"""
    fic   : the file handler
    
    wvalname  : Path to watervalue / Bellman function file
    st_idx    : the stage index
"""
function write_smspp_PolyBlock(fic, wvalname, st_idx)
    println(fic,string("    group: PolyhedralFunctionBlock { "))
    println(fic,"      dimensions:")
    # Interrogate the Bellman fct file
    bvalst = DataFrame()
    # Read Bellman file
    if ( !isempty(wvalname) )
        bval = CSV.read(wvalname, DataFrame; delim=';')
        bvalst = filter( row -> (row.Timestep .== st_idx), bval )
    end
    (nCuts, nCols) = size(bvalst)
    nVars = nCols - 2
    println(fic,string("        PolyFunction_NumVar = ", string(max(nCols-2,0)), " ; "))       
    if ( nCuts > 0 )
        println(fic,string("        PolyFunction_NumRow = ", string(nCuts), " ; "))
        println(fic,"      variables:")
        println(fic,"        double PolyFunction_A(PolyFunction_NumRow, PolyFunction_NumVar) ;")
        println(fic,"        double PolyFunction_b(PolyFunction_NumRow) ;")
    end
    println(fic,"")    
    println(fic,"    // group attributes:")
    println(fic,"    		:type = \"PolyhedralFunctionBlock\" ;")
    if ( nCuts > 0 )
        println(fic,"      data:")
    end
    cutRow = Vector{Float64}(undef, nVars)
    println(fic, "       PolyFunction_A =")
    for ic=1:nCuts
        for iv=1:nVars
            cutRow[iv] = bvalst[ic,1+iv]
        end
        if ( ic == nCuts )
            println(fic, string("     ", chop(string(cutRow),head=1), "  ;" ) )
        else
            println(fic, string("     ", chop(string(cutRow),head=1), " ," ) )
        end
    end
    println(fic,"")
    if ( nCuts > 0 )
        println(fic, string("       PolyFunction_b = ", chop(string(bvalst.b),head=1), "  ;" ))
    end
    println(fic,"      } // group PolyhedralFunctionBlock")
end

"""
    fic   : the file handler
    ts_dir: directory to the timeseries
    
    b_idx  : the block index
    s_idx  : the scenario index
    t_phi  : the 'tangent' phi value, with tan(phi)=Q/P => Q = tan(phi) P 
    t_line : the data frame line

    b_date : begin Date
    e_date : end Date
"""
function write_smspp_HunitBlocks(fic, ts_dir, b_idx, s_idx, t_phi, t_line, b_date, e_date)
    println(fic,string("    group: HydroUnitBlock_", string(b_idx)," { "))
    println(fic,"      dimensions:")
    println(fic,string("      	NumberReservoirs = ", string(1), " ; "))
    println(fic,string("      	NumberArcs = ", string(1), " ; "))
    println(fic,string("      	NumberIntervals = ", string(convert(Dates.Hour,(e_date - b_date)).value), " ; "))
    println(fic,"      variables:")
    println(fic,"      	uint StartArc(NumberArcs) ;")
    println(fic,"      	uint EndArc(NumberArcs) ;")
    println(fic,"      	double MaxVolumetric(NumberReservoirs) ;")
    println(fic,"      	double MinVolumetric(NumberReservoirs) ;")
    println(fic,"      	double Inflows(NumberReservoirs, NumberIntervals) ;")
    println(fic,"      	double MaxFlow(NumberArcs) ;")
    println(fic,"      	double MaxPower(NumberArcs) ;")
    println(fic,"      	double MinFlow(NumberArcs) ;")
    println(fic,"      	double MinPower(NumberArcs) ;")
    println(fic,"      	double MaxReactivePower(NumberArcs) ;")
    println(fic,"      	double VoltageMagnitude(NumberArcs) ;")
    println(fic,"      	uint NumberPieces(NumberArcs) ;")
    println(fic,"      	double LinearTerm(NumberArcs) ;")
    println(fic,"      	double ConstantTerm(NumberArcs) ;")
    println(fic,"      	double InitialVolumetric(NumberReservoirs) ;")
    println(fic,"      	double InitialFlowRate(NumberArcs) ;")
    println(fic,"")
    println(fic,"      // group attributes:")
    println(fic,"    		:type = \"HydroUnitBlock\" ;")
    println(fic,string("    		:name = \"", string(t_line.Name),"\" ; ") )
    println(fic,"      data:")
    println(fic,"")
    println(fic,string("       StartArc = ",string(0), " ; "))
    println(fic,"")
    println(fic,string("       EndArc = ",string(1), " ; "))
    println(fic,"")
    println(fic,string("       MaxVolumetric = ",string(t_line.MaxVolume), " ; "))
    println(fic,"")
    println(fic,string("       MinVolumetric = ",string(t_line.MinVolume), " ; "))
    println(fic,"")
    # Inflows
    ts_fname = string(ts_dir, "/", t_line.InflowsProfile )
    println( string("Hunit @ Node ", t_line.Zone, " timeseries : ", ts_fname) )
    ts_vals  = extract_ts_sequence( ts_fname, b_date, e_date, s_idx )    
    infl = round.(ts_vals.*t_line.Inflows, digits=7)
    #
    println(fic,string("       Inflows = ", chop(string(infl),head=1), "  ;" ));
    println(fic,"")
    println(fic,string("       MaxFlow = ",string(t_line.MaxPower), " ; "))
    println(fic,"")
    println(fic,string("       MaxPower = ",string(t_line.MaxPower), " ; "))
    println(fic,"")
    println(fic,string("       MinFlow = ",string(0), " ; "))
    println(fic,"")
    println(fic,string("       MinPower = ",string(t_line.MinPower), " ; "))
    println(fic,"")
    println(fic,string("       MaxReactivePower = ",string(t_phi*t_line.MaxPower), " ; "))
    println(fic,"")
    println(fic,string("       VoltageMagnitude = ",string(1.0), " ; "))
    println(fic,"")
    println(fic,string("       NumberPieces = ",string(1), " ; "))
    println(fic,"")
    println(fic,string("       LinearTerm = ",string(1), " ; "))
    println(fic,"")
    println(fic,string("       ConstantTerm = ",string(0), " ; "))
    println(fic,"")
    println(fic,string("       InitialVolumetric = ",string(t_line.InitialVolume), " ; "))
    println(fic,"")
    println(fic,string("       InitialFlowRate = ",string(0), " ; "))
    println(fic,"")
    println(fic,string("      } // group HydroUnitBlock_",string(b_idx)))
    println(fic,"")
end

"""
    fic   : the file handler
    ts_dir: directory to the timeseries
    
    b_idx  : the block index
    t_phi  : the 'tangent' phi value, with tan(phi)=Q/P => Q = tan(phi) P
    t_line : the data frame line
    
    b_date : begin Date
    e_date : end Date
"""
function write_smspp_TBlocks(fic, ts_dir, b_idx, t_phi, t_line, b_date, e_date)
    println(fic,string("  group: UnitBlock_", string(b_idx)," { "))
    println(fic,"    dimensions:")
    println(fic,string("    	NumberIntervals = ", string(convert(Dates.Hour,(e_date - b_date)).value), " ; "))
    println(fic,"    variables:")
    if (isempty(t_line.MaxPowerProfile))
        println(fic,"    	double MaxPower ;")
        println(fic,"    	double MaxReactivePower ;")
    else
        println(fic,"    	double MaxPower(NumberIntervals) ;")
        println(fic,"    	double MaxReactivePower(NumberIntervals) ;")
    end
    println(fic,"    	double VoltageMagnitude ;")
    println(fic,"    	double LinearTerm ;")
    println(fic,"    	double ConstTerm ;")
    println(fic,"    	double InitialPower ;")
    println(fic,"    	uint InitUpDownTime ;")
    println(fic,"")
    println(fic,"    // group attributes:")
    println(fic,"    		:type = \"ThermalUnitBlock\" ;")
    println(fic,string("    		:name = \"", string(t_line.Name),"\" ; ") )
    println(fic,"    data:")
    if (isempty(t_line.MaxPowerProfile))
        println(fic,string("     MaxPower = ",string(t_line.MaxPower), " ; "))
        println(fic,string("     MaxReactivePower = ",string(t_phi*t_line.MaxPower), " ; "))
    else
        ts_fname = string(ts_dir, "/", t_line.MaxPowerProfile )
        ts_vals  = extract_ts_sequence( ts_fname, b_date, e_date, 0 )

        println( string("THF @ Node ", t_line.Zone, " timeseries : ", ts_fname) )
        mxP = round.(ts_vals.*t_line.MaxPower, digits=5)
        println(fic,string("     MaxPower = ", chop(string(mxP),head=1), "  ;" ));
        println(fic,string("     MaxReactivePower = ", chop(string(t_phi*mxP),head=1), "  ;" ));
    end
    println(fic,"")
    println(fic,string("     VoltageMagnitude = ",string(1.0), " ; "))
    println(fic,"")
    println(fic,string("     LinearTerm = ",string(t_line.VariableCost), " ; "))
    println(fic,"")
    println(fic,string("     ConstTerm = ",string(t_line.FixedCost), " ; "))
    println(fic,"")
    println(fic,string("     InitialPower = ",string(t_line.MaxPower), " ; "))
    println(fic,"")
    println(fic,string("     InitUpDownTime = ",string(1), " ; "))
    println(fic,"")
    println(fic,string("    } // group UnitBlock_",string(b_idx)))
    println(fic,"")
end

"""
    fic   : the file handler
    ts_dir: directory to the timeseries
    
    b_idx  : the block index
    s_idx  : the scenario index
    t_phi  : the 'tangent' phi value, with tan(phi)=Q/P => Q = tan(phi) P
    t_line : the data frame line

    b_date : begin Date
    e_date : end Date
"""
function write_smspp_RBlocks(fic, ts_dir, b_idx, s_idx, t_phi, t_line, b_date, e_date)
    println(fic,string("  group: UnitBlock_", string(b_idx)," { "))
    println(fic,"    dimensions:")
    println(fic,string("    	NumberIntervals = ", string(convert(Dates.Hour,(e_date - b_date)).value), " ; "))
    println(fic,"    variables:")
    if (isempty(t_line.MaxPowerProfile))
        println(fic,"    	double MaxPower ;")
        println(fic,"    	double MaxReactivePower ;")
    else
        println(fic,"    	double MaxPower(NumberIntervals) ;")
        println(fic,"    	double MaxReactivePower(NumberIntervals) ;")
    end
    println(fic,"    	double MinPower ;")
    println(fic,"    	double VoltageMagnitude ;")
    println(fic,"    	double Gamma ;")
    println(fic,"    	double Kappa ;")
    println(fic,"")
    println(fic,"    // group attributes:")
    println(fic,"    		:type = \"IntermittentUnitBlock\" ;")
    println(fic,string("    		:name = \"", string(t_line.Name),"\" ; ") )
    println(fic,"    data:")
    if (isempty(t_line.MaxPowerProfile))
        println(fic,string("     MaxPower = ",string(t_line.MaxPower), " ; "))
        println(fic,string("     MaxReactivePower = ",string(t_phi*t_line.MaxPower), " ; "))
    else        
        ts_fname = string(ts_dir, "/", t_line.MaxPowerProfile )

        println( string("RES @ Node ", t_line.Zone, " timeseries : ", ts_fname) )

        ts_vals  = extract_ts_sequence( ts_fname, b_date, e_date, 0 )
        
        mxP = round.(ts_vals.*t_line.MaxPower, digits=5)
        println(fic,string("     MaxPower = ", chop(string(mxP),head=1), "  ;" ));
        println(fic,string("     MaxReactivePower = ", chop(string(t_phi*mxP),head=1), "  ;" ));
    end
    println(fic,"")
    println(fic,string("     MinPower = ",string(t_line.MinPower), " ; "))
    println(fic,"")
    println(fic,string("     VoltageMagnitude = ",string(1.0), " ; "))
    println(fic,"")
    println(fic,string("     Gamma = ",string(0), " ; "))
    println(fic,"")
    println(fic,string("     Kappa = ",string(t_line.Kappa), " ; "))
    println(fic,"")
    println(fic,string("    } // group UnitBlock_",string(b_idx)))
    println(fic,"")
end

"""
    fic   : the file handler
    
    b_idx  : the block index
    t_phi  : the 'tangent' phi value, with tan(phi)=Q/P => Q = tan(phi) P
    t_line : the data frame line

    b_date : begin Date
    e_date : end Date
"""
function write_smspp_BatBlocks(fic, b_idx, t_phi, t_line, b_date, e_date)
    println(fic,string("  group: UnitBlock_", string(b_idx)," { "))
    println(fic,"    dimensions:")
    println(fic,string("    	NumberIntervals = ", string(convert(Dates.Hour,(e_date - b_date)).value), " ; "))
    println(fic,"    variables:")
    println(fic,"    	double MaxPower ;")
    println(fic,"    	double MinPower ;")
    println(fic,"    	double MaxReactivePower ;")
    println(fic,"    	double MinReactivePower ;")    
    println(fic,"    	double VoltageMagnitude ;")
    println(fic,"    	double MaxStorage ;")
    println(fic,"    	double MinStorage ;")
    println(fic,"    	double InitialPower ;")
    println(fic,"    	double InitialStorage ;")
    println(fic,"    	double StoringBatteryRho ;")
    println(fic,"    	double ExtractingBatteryRho ;")
    println(fic,"    	double Demand(NumberIntervals) ;")
    println(fic,"")
    println(fic,"    // group attributes:")
    println(fic,"    		:type = \"BatteryUnitBlock\" ;")
    println(fic,string("    		:name = \"", string(t_line.Name),"\" ; ") )
    println(fic,"    data:")
    println(fic,"")
    println(fic,string("     MaxPower = ",string(t_line.MaxPower), " ; "))
    println(fic,"")
    println(fic,string("     MinPower = ",string(t_line.MinPower), " ; "))
    println(fic,"")
    println(fic,string("     MaxReactivePower = ",string(t_phi*t_line.MaxPower), " ; "))
    println(fic,"")
    println(fic,string("     MinReactivePower = ",string(t_phi*t_line.MinPower), " ; "))
    println(fic,"")
    println(fic,string("     VoltageMagnitude = ",string(1.0), " ; "))
    println(fic,"")
    println(fic,string("     MaxStorage = ",string(t_line.MaxVolume), " ; "))
    println(fic,"")
    println(fic,string("     MinStorage = ",string(t_line.MinVolume), " ; "))
    println(fic,"")
    println(fic,string("     InitialPower = ",string(0), " ; "))
    println(fic,"")
    println(fic,string("     InitialStorage = ",string(t_line.InitialVolume), " ; "))
    println(fic,"")
    println(fic,string("     StoringBatteryRho = ",string(t_line.PumpingEfficiency), " ; "))
    println(fic,"")
    println(fic,string("     ExtractingBatteryRho = ",string(t_line.TurbineEfficiency), " ; "))
    println(fic,"")    
    dVec = Vector{Float64}(undef,convert(Dates.Hour,(e_date - b_date)).value)
    dVec .= 0.0
    println(fic,string("     Demand = ", chop(string(dVec),head=1), "  ;" ));
    println(fic,"")
    println(fic,string("    } // group UnitBlock_",string(b_idx)))
    println(fic,"")
end

"""
    fic   : the file handler
    
    b_idx  : the block index
    naam   : the name of this fellow

    b_date : begin Date
    e_date : end Date
"""
function write_smspp_SlackBlocks(fic, b_idx, naam, b_date, e_date)
    println(fic,string("  group: UnitBlock_", string(b_idx)," { "))
    println(fic,"    dimensions:")
    println(fic,string("    	NumberIntervals = ", string(convert(Dates.Hour,(e_date - b_date)).value), " ; "))
    println(fic,"    variables:")
    println(fic,"    	double MaxPower ;")
    println(fic,"    	double MaxReactivePower ;")
    println(fic,"    	double ActivePowerCost ;")
    println(fic,"    	double VoltageMagnitude ;")
    println(fic,"")
    println(fic,"    // group attributes:")
    println(fic,"    		:type = \"SlackUnitBlock\" ;")
    println(fic,string("    		:name = \"", naam,"\" ; ") )
    println(fic,"    data:")
    println(fic,"")
    println(fic,string("     MaxPower = ",string(1500000), " ; "))
    println(fic,"")
    println(fic,string("     MaxReactivePower = ",string(800000), " ; "))
    println(fic,"")    
    println(fic,string("     ActivePowerCost = ",string(10000), " ; "))
    println(fic,"")
    println(fic,string("     VoltageMagnitude = ",string(1.0), " ; "))
    println(fic,"")
    println(fic,string("    } // group UnitBlock_",string(b_idx)))
    println(fic,"")
end

function strVtoList(d::Vector{<:AbstractString})
    s = "";
    for el in d
        s = string(s,"\"",el,"\",")
    end
    s = chop(s)
    return s
end