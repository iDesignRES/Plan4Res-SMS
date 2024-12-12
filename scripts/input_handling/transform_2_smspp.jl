using CSV
using DataFrames
using Dates
using Statistics

include("make_smspp_uc.jl")

# GitHub Repository structure:
#
# -- Absolute or relative Path
github_local_d = string(@__DIR__, "/../..")
github_smsppin = "/smspp_in"
github_smsppout = "/smspp_out"
#
# File names as put together by NTNU
#
ntnu_bus_file  = string(github_local_d, "/input_data/","Bus_data.csv")
ntnu_line_file = string(github_local_d, "/input_data/","lines.csv") #"Line.csv"
ntnu_gen_file  = string(github_local_d, "/input_data/","generation.csv") #"Generators.csv"
ntnu_load_file = string(github_local_d, "/input_data/","load.csv")

#
# Other data
#
idx_scen      = 0 #Some scenario index
st_idx        = 0 #Stochastic stage belonging to this computation => to be deduced from ssv_stp
bgn_date      = DateTime("02/07/2030 00:00", "dd/mm/yyyy HH:MM")
#
uc_bgn_date   = DateTime("08/07/2030 00:00", "dd/mm/yyyy HH:MM")
uc_end_date   = DateTime("09/07/2030 00:00", "dd/mm/yyyy HH:MM")
#
generic_year  = 2050
# 
edf_impexp    = string(github_local_d, github_smsppout, "/nuts0/results_simul/Flows/Flows",string(idx_scen),".csv")
focus_country = "ES"
#
# Do we want to make a computation with DC opf stuff :
#
with_dcopf    = false #true
with_acopf    = false #supersedes with_dcopf flag
tan_phi       = 0.57 # assuming 30° phase angle
nb_clip       = 1 # 24 #Clip the timeperiod into multiple substeps in case of DC opf
nb_sstep      = div(convert(Dates.Hour, (uc_end_date - uc_bgn_date)).value, nb_clip)

# Bellman values
# 
ssv_step      = 168
bell_vals     = string(github_local_d, github_smsppout, "/nuts0/results_simul/BellmanValuesOUT.csv")
vol_file      = string(github_local_d, github_smsppout, "/nuts0/results_simul/Volume/Volume",string(idx_scen),".csv")

# Compute the stochastic stage index based on the temporal data
st_idx = div(convert(Dates.Hour, (uc_end_date - bgn_date)).value, ssv_step)

println(" -- Starting the transformation of the data -- ")

## Reading the underlying data DataFrames
#

# Load the data of the buses
#
bus_data  = CSV.read(ntnu_bus_file, DataFrame; delim=',')

# Load the data of the lines
#
lines_data = CSV.read(ntnu_line_file, DataFrame; delim=',')
f_to_km = 1.0 ; #Conversion factor of km to the unit of "length" ;

# Check if we have multiple lines with the same name
# 
if ( size(lines_data)[1] != length(unique(lines_data.line_id)) )    
    for linid in lines_data.line_id
        if (length( findall( lines_data.line_id .== linid ))>1 )
            println(string("Line with name ", linid, " has duplicates"))
        end
    end
    error("Duplicate lines found")
end

# Load the generator data
#
gen_data = CSV.read(ntnu_gen_file, DataFrame; delim=',')

# if The generator data has a unit_id field, then we can replace the names with that one
if ( "unit_id" in names(gen_data) )
    sbgd = filter( row -> ( ismissing(row.name) ), gen_data)
    for uid in sbgd.unit_id
        #println(uid)
        i0 = findall( gen_data.unit_id .== uid )
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
new_names = replace.(gen_data[:,:name],"Á"=> "A", "É"=>"E", "Í"=> "I", "Ñ"=>"N", "Ó"=>"O","Ú" => "U", "Ü"=>"U")
gen_data[:,:name] .= new_names

# Read the load data
#
load_data = CSV.read(ntnu_load_file, DataFrame; delim=',')

# Read the Import / Export data
# 
impexp_data = CSV.read( edf_impexp, DataFrame; delim=',')
coln = names(impexp_data)[2:end]
countries_to_cut = setdiff(unique(bus_data.country),[focus_country])

# Read the Bellman values
#
bell_data = CSV.read( bell_vals, DataFrame; delim=',')
vol_data  = CSV.read( vol_file, DataFrame; delim=',')

# Build the import data frame == thermal units generating the to be imported production
# Make the thermal unit file from this
nb_imp = 0
tu_imp_data = DataFrame(Zone=Vector{String}(undef, nb_imp), Name=Vector{String}(undef, nb_imp), NumberUnits=Vector{Int64}(undef, nb_imp),
                        MaxPower=Vector{Float64}(undef, nb_imp), MaxPowerProfile=Vector{String}(undef, nb_imp), VariableCost=Vector{Float64}(undef, nb_imp),	FixedCost=Vector{Float64}(undef, nb_imp),	
                        InvestmentCost=Vector{Float64}(undef, nb_imp),	Capacity=Vector{Float64}(undef, nb_imp),	Energy=Vector{Float64}(undef, nb_imp),	
                       	MaxAddedCapacity=Vector{Float64}(undef, nb_imp) )

# Some work may be needed in case the first column of the import export file is not strings but integers
tcol = Vector{String}(undef, 0)
if ( typeof(impexp_data[1,1]) == Int64 )
    dlist = bgn_date + Dates.Hour.(impexp_data[:,1])

    dlist2 = copy(dlist)
    bsyear = Dates.year(dlist[1])
    for yr in unique(Dates.year.(dlist))
        I = findall( Dates.year.(dlist) .== yr )
        dlist2[I] = replaceparts.( dlist[I], year = generic_year + (yr - bsyear) )
    end
    tcol  = Dates.format.(dlist, "dd/mm/yyyy HH:MM")
else
    tcol  = impexp_data[:,1]
end

nb_ts = size(impexp_data)[1]
imp_ts = DataFrame( Date=Vector{String}(undef, nb_ts), MaxPower=Vector{Float64}(undef, nb_ts) )
imp_ts[:,:Date] = tcol

exp_ts = DataFrame( Date=Vector{String}(undef, nb_ts), Exp=Vector{Float64}(undef, nb_ts) )
exp_ts[:,:Date] = tcol

# Make the Deterministic Timeseries file for the smspp scripts
#
det_ts = DataFrame( Date=Vector{String}(undef, nb_ts) )
det_ts[:,:Date] = Dates.format.(dlist2, "dd/mm/yyyy HH:MM")

for j=1:length(countries_to_cut)
    ln_name1 = string(focus_country,">",countries_to_cut[j])
    ln_name2 = string(countries_to_cut[j],">",focus_country)

    # Find all buses in target country
    I = findall( bus_data.country .== countries_to_cut[j])
    weight = Vector{Float64}(undef,length(I))
    weight[:] .= 0.0
    for l=1:length(I)
        # For each of these buses find all lines to and from this bus
        #JJ = findall( (lines_data.bus0 .== bus_data[I[l],:bus_id]) .| (lines_data.bus1 .== bus_data[I[l],:bus_id]) .& (lines_data.length .> 10) )
        JJ = findall( (lines_data.bus0 .== bus_data[I[l],:bus_id]) .| (lines_data.bus1 .== bus_data[I[l],:bus_id]) )
        # For all of those lines, find out some strength
        weight[l] = sum( (lines_data[JJ,:voltage].*lines_data[JJ,:Imax]) )
    end

    if ( ln_name1 in coln )
        kcol = findall( ln_name1 in coln )[1]
        # Negative values are importations to focus_country : these are productions in the target nodes ; positive values are demands        
        for l=1:length(I)
            wgh = weight[l]/sum(weight)

            #Import stuff 
            fic_name_i = string("imp_pmax_", bus_data.bus_id[[I[l]]][1], ".csv")
            imp_ts[:,:MaxPower] = -1.0*min.( wgh.*impexp_data[:,coln[kcol]] ,0.0)
            CSV.write(string(github_local_d, github_smsppin, "/ts/",fic_name_i), imp_ts; delim=',')

            #push!( tu_imp_data, [string(bus_data.country[I[l]],"_", bus_data.bus_id[I[l]]),"IMP",1,1.0,fic_name_i,
            #                     0.001, 0,0,0,0,0 ] )

            push!( tu_imp_data, [string(bus_data.bus_id[I[l]]),"IMP",1,1.0,fic_name_i,
                                 0.001, 0,0,0,0,0 ] )

            #Add little column to det_ts
            det_ts[:, string("imp_pmax_", bus_data.bus_id[[I[l]]][1])] = imp_ts[:,:MaxPower]

            #Export stuff
            fic_name_e = string(github_local_d, github_smsppin,"/ts/exp_load_", bus_data.bus_id[[I[l]]][1], ".csv")
            exp_ts[:,:Exp] = max.( wgh.*impexp_data[:,coln[kcol]] ,0.0)
            CSV.write(fic_name_e, exp_ts; delim=',')
        end
    end
    if ( ln_name2 in coln )
        kcol = findall( ln_name2 in coln )[1]
        # Negative values are exportations to focus_country : 
        for l=1:length(I)
            wgh = weight[l]/sum(weight)

            #Import stuff 
            fic_name_i = string("imp_pmax_", bus_data.bus_id[[I[l]]][1], ".csv")
            imp_ts[:,:MaxPower] = max.( wgh.*impexp_data[:,coln[kcol]] ,0.0)
            CSV.write(string(github_local_d, github_smsppin,"/ts/",fic_name_i), imp_ts; delim=',')

            #push!( tu_imp_data, [string(bus_data.country[I[l]],"_", bus_data.bus_id[I[l]]),"IMP",1, 1.0, fic_name_i,
            #                     0.001, 0,0,0,0,0 ] )
            
            push!( tu_imp_data, [string(bus_data.bus_id[I[l]]),"IMP",1, 1.0, fic_name_i,
                                 0.001, 0,0,0,0,0 ] )

            #Add little column to det_ts
            det_ts[:, string("imp_pmax_", bus_data.bus_id[[I[l]]][1])] = imp_ts[:,:MaxPower]

            #Export stuff
            fic_name_e = string(github_local_d, github_smsppin,"/ts/exp_load_", bus_data.bus_id[[I[l]]][1], ".csv")
            exp_ts[:,:Exp] = -1.0*min.( wgh.*impexp_data[:,coln[kcol]] ,0.0)
            CSV.write(fic_name_e, exp_ts; delim=',')
        end
    end
end
# Write the deterministic TS file
CSV.write(string(github_local_d, github_smsppin, "/ts/DeterministicTS.csv"), det_ts; delim=';')

# Make the ZP_Partition file
#
zp_data = DataFrame()
#both options below work
#zp_data[:, :Countries] = string.(bus_data[!, "country"], "_", bus_data[!, "bus_id"]) #
zp_data[:, :Countries] = bus_data.bus_id #string.(bus_data.country, "_", bus_data.bus_id)
zp_data[:, :Continent] .= "Continent"
if (with_acopf)
    insertcols!(zp_data, :NodeSusceptance=>0.0 )
    insertcols!(zp_data, :NodeConductance=>0.0 )
    insertcols!(zp_data, :NodeMinVoltage=>0.9 )
    insertcols!(zp_data, :NodeMaxVoltage=>1.1 )
    #
    zp_data[:, :NodeMinVoltage] .= bus_data.voltage*0.9
    zp_data[:, :NodeMaxVoltage] .= bus_data.voltage*1.1
end

CSV.write(string(github_local_d, github_smsppin, "/nutsx/ZP_ZonePartition.csv"), zp_data; delim=';')

# Make the ZV_ZoneValues file
# Total load comes from ENTSOE-E sources
# 
load_d   = Dict([("FR", (1410503232, "AggregatedTimeSerie_Total_BigFrance.csv")), ("ES", (228746533,"Profile-Iberia.csv")),
                 ("PT", (50567216,"Profile-Iberia.csv")) ])

nb_bus = length(bus_data.bus_id)
#
zv_zone_data = DataFrame(Type=Vector{String}(undef, nb_bus), Zone=Vector{String}(undef, nb_bus), 
                         value=Vector{Float64}(undef, nb_bus), Profile_Timeserie=Vector{String}(undef, nb_bus) )

zv_zone_data[:, :Type] .= "Total"
zv_zone_data[:, :Zone] = bus_data.bus_id #string.(bus_data.country, "_", bus_data.bus_id)

for i=1:nb_bus
    l_info = load_d[bus_data.country[i]]

    l0 = length(findall(load_data.busID .== bus_data.bus_id[i]))
    if ( l0 == 0)
        error("The bus ", bus_data.bus_id[i], " does not have load data ")
    end
    i0 = findall(load_data.busID .== bus_data.bus_id[i])

    if ( bus_data.country[i] in countries_to_cut )
        zv_zone_data[i,:value] = 1.0
        zv_zone_data[i,:Profile_Timeserie] = string("exp_load_", bus_data.bus_id[i], ".csv")
    else
        zv_zone_data[i,:value] = load_data.demand[i0][1]*l_info[1]
        zv_zone_data[i,:Profile_Timeserie] = l_info[2]
    end
end
CSV.write(string(github_local_d, github_smsppin, "/nutsx/ZV_ZoneValues.csv"), zv_zone_data; delim=';')

# Make the Aggregate ZV :
#
zv_agg_d = similar(zv_zone_data,0)
for cn in unique(bus_data.country)
    l_info = load_d[cn]
    push!( zv_agg_d, ["Total", cn, l_info[1], l_info[2]] )
end
CSV.write(string(github_local_d, github_smsppin,"/nuts0/ZV_ZoneValues.csv"), zv_agg_d; delim=';')

#
# Figure out the transformation buses
#
I_all = collect(1:nb_bus)
I_trns = []
nb_trans = 0
while ( length(I_all) > 0 )
    i_all = I_all[1]
    # Check if this bus has a transformer
    I = findall( (bus_data.x .== bus_data.x[i_all]) .& (bus_data.y .== bus_data.y[i_all] ) )
    # This of course contains at least i_all
    #
    global nb_trans += length(I) - 1
    if (length(I) > 1)
        global I_trns = union(I_trns, i_all)
    end
    global I_all = setdiff(I_all, I)
end
#nb_trans is not the lenght of I_trns if one bus can be split over multiple sub_buses with transofrmators

if ( !with_dcopf && !with_acopf )
    transf_data = DataFrame(Name=Vector{String}(undef, nb_trans), StartLine=Vector{String}(undef, nb_trans), EndLine=Vector{String}(undef, nb_trans),
                        MaxPowerFlow=Vector{Float64}(undef, nb_trans),	MinPowerFlow=Vector{Float64}(undef, nb_trans),
                        MaxAddedCapacity=Vector{Float64}(undef, nb_trans), MaxRetCapacity=Vector{Float64}(undef, nb_trans), InvestmentCost=Vector{Float64}(undef, nb_trans) )
else
    transf_data = DataFrame(Name=Vector{String}(undef, nb_trans), StartLine=Vector{String}(undef, nb_trans), EndLine=Vector{String}(undef, nb_trans),
                        MaxPowerFlow=Vector{Float64}(undef, nb_trans),	MinPowerFlow=Vector{Float64}(undef, nb_trans), Susceptance=Vector{Float64}(undef, nb_trans),
                        MaxAddedCapacity=Vector{Float64}(undef, nb_trans), MaxRetCapacity=Vector{Float64}(undef, nb_trans), InvestmentCost=Vector{Float64}(undef, nb_trans) )
end
if ( with_acopf )
    insertcols!(transf_data, :LineResistance=>Vector{Float64}(undef, nb_trans) )
    insertcols!(transf_data, :LineReactance=>Vector{Float64}(undef, nb_trans) )
    insertcols!(transf_data, :LineMinAngle=>Vector{Float64}(undef, nb_trans) )
    insertcols!(transf_data, :LineMaxAngle=>Vector{Float64}(undef, nb_trans) )
    insertcols!(transf_data, :LineShiftAngle=>Vector{Float64}(undef, nb_trans) )
    insertcols!(transf_data, :LineRATEA=>Vector{Float64}(undef, nb_trans) )
    insertcols!(transf_data, :LineRatio=>Vector{Float64}(undef, nb_trans) )
end

i_trns = 0
for i=1:length(I_trns)
    i_all = I_trns[i]
    I = findall( (bus_data.x .== bus_data.x[i_all]) .& (bus_data.y .== bus_data.y[i_all] ) )
    for j=2:length(I)
        global i_trns += 1            
        #transf_data[i_trns, :StartLine] = string.( bus_data.country[I[1]], "_", bus_data.bus_id[I[1]] )
        #transf_data[i_trns, :EndLine] = string.( bus_data.country[I[j]], "_", bus_data.bus_id[I[j]] )
        transf_data[i_trns, :StartLine] = string.( bus_data.bus_id[I[1]] )
        transf_data[i_trns, :EndLine] = string.( bus_data.bus_id[I[j]] )

        transf_data[i_trns, :Name] = string.("T_", bus_data.bus_id[I[1]], "_", bus_data.bus_id[I[j]] )
    end
end
transf_data[!, :MaxPowerFlow] .= 10000.0
transf_data[!, :MinPowerFlow] .= -10000.0
transf_data[!, :MaxAddedCapacity] .= 0.0
transf_data[!, :MaxRetCapacity] .= 0.0
transf_data[!, :InvestmentCost] .= 0.0
if ( with_dcopf || with_acopf )
    transf_data[!, :Susceptance] .= 0.0
end
if ( with_acopf )
    transf_data[!, :LineResistance] .= 0.0
    transf_data[!, :LineReactance] .= 0.0
    transf_data[!, :LineMinAngle] .= 0.0
    transf_data[!, :LineMaxAngle] .= 0.0
    transf_data[!, :LineShiftAngle] .= 0.0
    transf_data[!, :LineRATEA] .= round(sqrt(1+tan_phi^2)*10000 ; digits=3)
    transf_data[!, :LineRatio] .= 0.0
end

# Make the interconnection file
#
#	https://www.imse.iastate.edu/files/2021/03/EnergyProject_Capacity_of_Transmission_Lines.pdf
# 
#   The capacity is (kV)^2 * sin(30 °) / (length (in km) * 0.327) in MW
#
s_cst  = 0.5 # sin( 1/12*2*pi ) = 0.5
ps_cst = 0.327 # positive-sequence reactance per phase is 0.327 ohm/km

#Characteristics of the line resistance [Ohm/km], reactance [Ohm/km] and capacitance [nF/km] (R, X, C)
line_charac = Dict( [(132,(0.0949,0.38,9.2)), 
                     (220,(0.06,0.301,12.5)),
                     (225,(0.06,0.301,12.5)),
                     (250,(0.022,1e-6,200)),
                     (320,(0.04,0.265,13.2)),
                     (400,(0.03,0.246,13.8)) ]   )

# B = susceptance = -X / (R^2 + X^2)

b_shunt_cap  = 0.012 # mu S / km 
L_inductance = 1.6   # m H / km / phase 
C_capacitance= 10 # nF / km / phase
# The formula for the susceptance is
#
# B := b / Z_0 with Z_0 = sqrt( L / C ) : [ mu S / Omega] and Omega = sqrt( H / F )
#

#
# incon_data = DataFrame()
# Compute the true number of lines
(nR, nC) = size( lines_data )
Ilines = collect(1:nR)
I_real_lines = Vector{Int64}(undef,0)
nb_par_lines = Vector{Int64}(undef,0)
while ( length(Ilines) >= 1 )
    iR = Ilines[1]

    sub_line = filter( row -> ( ( (row.bus0 .== lines_data.bus0[iR]) && (row.bus1 .== lines_data.bus1[iR]) 
                                                                     && (row.voltage .== lines_data.voltage[iR]) )
                                            ||
                                ( (row.bus1 .== lines_data.bus0[iR]) && (row.bus0 .== lines_data.bus1[iR]) 
                                                                     && (row.voltage .== lines_data.voltage[iR]) ) ) , lines_data )

    for nm in sub_line.line_id
        i0 = findall( lines_data.line_id .== nm )
        setdiff!( Ilines, i0[1] )        
    end
    append!(I_real_lines, iR )
    append!(nb_par_lines, length(sub_line.line_id))
end

nb_icdata = length(I_real_lines)
if ( !with_dcopf && !with_acopf )
    incon_data = DataFrame(Name=Vector{String}(undef, nb_icdata), StartLine=Vector{String}(undef, nb_icdata), EndLine=Vector{String}(undef, nb_icdata),
                        MaxPowerFlow=Vector{Float64}(undef, nb_icdata),	MinPowerFlow=Vector{Float64}(undef, nb_icdata),
                        MaxAddedCapacity=Vector{Float64}(undef, nb_icdata), MaxRetCapacity=Vector{Float64}(undef, nb_icdata), InvestmentCost=Vector{Float64}(undef, nb_icdata) )
else
    incon_data = DataFrame(Name=Vector{String}(undef, nb_icdata), StartLine=Vector{String}(undef, nb_icdata), EndLine=Vector{String}(undef, nb_icdata),
                        MaxPowerFlow=Vector{Float64}(undef, nb_icdata),	MinPowerFlow=Vector{Float64}(undef, nb_icdata), Susceptance=Vector{Float64}(undef, nb_icdata),
                        MaxAddedCapacity=Vector{Float64}(undef, nb_icdata), MaxRetCapacity=Vector{Float64}(undef, nb_icdata), InvestmentCost=Vector{Float64}(undef, nb_icdata) )
end
if ( with_acopf )
    # Observe that LineRATEA, ShiftAngle, LineRatio are not used
    insertcols!(incon_data, :LineResistance=>Vector{Float64}(undef, nb_icdata) )
    insertcols!(incon_data, :LineReactance=>Vector{Float64}(undef, nb_icdata) )
    insertcols!(incon_data, :LineMinAngle=>Vector{Float64}(undef, nb_icdata) )
    insertcols!(incon_data, :LineMaxAngle=>Vector{Float64}(undef, nb_icdata) )
    insertcols!(incon_data, :LineShiftAngle=>Vector{Float64}(undef, nb_icdata) )
    insertcols!(incon_data, :LineRATEA=>Vector{Float64}(undef, nb_icdata) )
    insertcols!(incon_data, :LineRatio=>Vector{Float64}(undef, nb_icdata) )
end

incon_data[:, :StartLine] .= "NUL"
incon_data[:, :EndLine] .= "NUL"
for itl=1:nb_icdata
    i = I_real_lines[itl] 
    incon_data[itl, :Name] = string.( "L", lines_data.line_id[i] )    
    
    if ( ( length(findall(bus_data.bus_id .== lines_data.bus0[i])) == 0) || (length(findall(bus_data.bus_id .== lines_data.bus1[i])) == 0) )
        error("line : ", string(i), " has not existing start and/or end bus : ", lines_data.bus0[i], " or ", lines_data.bus1[i] )
    end
    i0 = findall(bus_data.bus_id .== lines_data.bus0[i])
    i1 = findall(bus_data.bus_id .== lines_data.bus1[i])
        
    incon_data[itl, :StartLine] = bus_data.bus_id[i0][1] #string.( bus_data.country[i0], "_", bus_data.bus_id[i0] )[1]
    incon_data[itl, :EndLine] = bus_data.bus_id[i1][1] #string.( bus_data.country[i1], "_", bus_data.bus_id[i1] )[1]

    # The following two are based on the above formula
    #incon_data[itl, :MaxPowerFlow] = nb_par_lines[itl]*round.(f_to_km*(s_cst*lines_data.voltage[i].^2)./(lines_data.length[i]*ps_cst);digits=3)
    #incon_data[itl, :MinPowerFlow] = -nb_par_lines[itl]*round.(f_to_km*(s_cst*lines_data.voltage[i].^2)./(lines_data.length[i]*ps_cst);digits=3)
    incon_data[itl, :MaxPowerFlow] = round(sqrt(3)/(sqrt(1 + tan_phi^2))*nb_par_lines[itl]*lines_data.Imax[i]*lines_data.voltage[i]; digits=3)
    incon_data[itl, :MinPowerFlow] = - incon_data[itl, :MaxPowerFlow]


    l_data = line_charac[lines_data.voltage[i]]
    if ( with_dcopf || with_acopf )
        incon_data[itl, :Susceptance] = 0.0
        if (lines_data.dc[i] == "f")
            # Use the formula
            # N.B. : (L * ℓ * 1e-3) / (C * ℓ * 1e-9 ) = 1e6 * (L/C)   )
            #incon_data[itl, :Susceptance] = 1e-3*( (b_shunt_cap*lines_data.length[i]*f_to_km ) / sqrt( L_inductance / C_capacitance ) )
            #
            # New formula
            incon_data[itl, :Susceptance] = nb_par_lines[itl]*( l_data[2] /( (l_data[1]^2 + l_data[2]^2)*lines_data.length[i]  ))
            if ( with_acopf )
                incon_data[itl, :LineResistance] = (l_data[1]*lines_data.length[i]) / nb_par_lines[itl] #Cf. Mise en // de lignes et mail de D. Croteau
                incon_data[itl, :LineReactance] = (l_data[2]*lines_data.length[i]) / nb_par_lines[itl]
                incon_data[itl, :LineMinAngle] = -30.0
                incon_data[itl, :LineMaxAngle] = 30.0
                incon_data[itl, :LineShiftAngle] = 0.0
                incon_data[itl, :LineRATEA] = round(sqrt(3)*nb_par_lines[itl]*lines_data.Imax[i]*lines_data.voltage[i]; digits=3) #round(sqrt(1 + tan_phi^2)*incon_data.MaxPowerFlow[itl][1]; digits=3)
                incon_data[itl, :LineRatio] = 0.0                
            end
        end
    end
end
# Investment data
incon_data[:, :MaxAddedCapacity] .= 500.0
incon_data[:, :MaxRetCapacity] .= 0.0
incon_data[:, :InvestmentCost] .= 100.0

# Append the transformers to the set
append!(incon_data, transf_data)

CSV.write(string(github_local_d, github_smsppin,"/nutsx/IN_Interconnections.csv"), incon_data; delim=';')

#
# Make the interconnection file at NUTS0 level as well
#
incon_agg_d = similar(incon_data, 0)
if ( with_dcopf || with_acopf)
    # Delete de Susceptance column not useful for NUTS0 computations
    select!(incon_agg_d, Not("Susceptance"))
end
if ( with_acopf )
    # Delete the other superfluous columns not useful for NUTS0 computations
    select!(incon_agg_d, Not("LineResistance"))
    select!(incon_agg_d, Not("LineReactance"))
    select!(incon_agg_d, Not("LineMinAngle"))
    select!(incon_agg_d, Not("LineMaxAngle"))
    select!(incon_agg_d, Not("LineShiftAngle"))
    select!(incon_agg_d, Not("LineRATEA"))
    select!(incon_agg_d, Not("LineRatio"))
end
# For each country we do stuff
c_list = unique(bus_data.country)
for i=1:length(c_list)
    cn = c_list[i]
    Ibus = findall( bus_data.country .== cn )
    for j=i+1:length(c_list)
        cn2 = c_list[j]
        Ibus2 = findall( bus_data.country .== cn2 )

        # filter out all lines with one bus in one country and the other in the other
        sub_t1 = filter(row -> ( ((row.bus0 in bus_data.bus_id[Ibus]) && (row.bus1 in bus_data.bus_id[Ibus2])) ||
                                 ((row.bus1 in bus_data.bus_id[Ibus]) && (row.bus0 in bus_data.bus_id[Ibus2])) ), lines_data)

        #I = findall( sub_t1.length .> 10 )
        mxFlow = round(sum(sqrt(3)*sub_t1.Imax.*sub_t1.voltage); digits=3) #round.( f_to_km*sum( (s_cst*sub_t1[I,:voltage].^2)./(sub_t1[I,:length]*ps_cst) );digits=3)

        push!( incon_agg_d, [string(cn,">",cn2), cn, cn2, mxFlow, -1.0*mxFlow, 500.0, 0.0, 100.0] )
    end
end
CSV.write(string(github_local_d, github_smsppin,"/nuts0/IN_Interconnections.csv"), incon_agg_d; delim=';')

#
# Make the Generation data from this stuff
#
# Type conversion on the column
gen_data[!, :primary_fuel] = convert.(String, gen_data[:, :primary_fuel])
# Some fix on the type of hydro
for ig = 1: size(gen_data)[1]
    if ( gen_data.primary_fuel[ig] == "Hydro" )
        #println( gen_data.unit_id[ig], "_", gen_data.technology[ig] )
        gen_data.primary_fuel[ig] = string(gen_data.primary_fuel[ig], "|",gen_data.technology[ig] )
    end
end
#
println("Found the following technologies: ", unique(gen_data.primary_fuel))

Thf_list = ["Coal", "Gas", "Biomass", "Waste", "Nuclear", "Oil", "Oil/Diesel", "Diesel", "Oil/Gas/Diesel", "Oil/Gas"];
Res_list = ["Wind", "Solar", "Hydro|run_of_river"]
STS_list = ["Hydro|pumped_storage"]
SS_list  = ["Hydro|reservoir"]

# List of technology related costs and information (prop cost, fixed cost, InvestmentCost)
#
cThf_l   = Dict([("Coal", (13.05, 42380, 1695200)), ("Gas", (31.24,	18011.5, 663122.3529)),
                 ("Biomass", (64.022, 55242.33,	2511015)), ("Waste", (62.00, 55242.33,	2511015)),
                 ("Nuclear", (12.42,	100122.75,	6357000)), ("Oil", (78.92,	6897.345, 688675)),
                 ("Oil/Diesel", (85.00,	6500.345, 688675)), ("Diesel", (90.00,	6500.345, 688675)),
                 ("Oil/Gas/Diesel", (83.00,	6500.345, 688675)), ("Oil/Gas", (82.00,	6500.345, 688675)) ])

# List of technology related information 
#
cRes_l   = Dict([("Wind",(1059500, "EDF__WindOnshore-LoadFactor-PresentClimate-","__13082019__13082019__v1.csv", 1.0)),
                 ("Solar",(476775,"EDF__PV-LoadFactor-PresentClimate-","__13082019__13082019__v1.csv", 1.0        )),
                 ("Hydro", (0 , "EDF__RunOfRiver-HourlyCoefficient-PresentClimate-","__18092019__18092019__v1.csv", 2500.0)), 
                 ("Hydro|run_of_river", (0 , "EDF__RunOfRiver-HourlyCoefficient-PresentClimate-","__18092019__18092019__v1.csv", 2500.0)) ])
                 
# List of technology related information maxVol, turEff, pumpEff
cSTS_l  = Dict([ ("Hydro|pumped_storage", (100, 1, 0.866) ),
                 ("Hydro (Pumped storage with natural inflow)", (100, 1, 0.866) ) ])

# List of technology related stuff
cSS_l = Dict([ ("Hydro|reservoir", (100, "EDF__Inflow-HourlyCoefficient-PresentClimate-","__18092019__18092019__v1.csv") ) ])
# We can rely on vol_data to find out the right percentage
# st_idx indicates the next stochastic stage, so we are in st_idx - 1
# println(string("vol init FR : ", string(vol_data[st_idx*ssv_step, 2]/16256280.72), " vol init ES : ", string(vol_data[st_idx*ssv_step, 3]/16409036.88) ))
fr_per = st_idx > 1 ? round(vol_data[(st_idx-1)*ssv_step, 2]/16256280.72, digits=3) : 0.3
es_per = st_idx > 1 ? round(vol_data[(st_idx-1)*ssv_step, 3]/16409036.88, digits=3) : 0.3
# maxVol, Hydrosystem name, inflows, %initial volume vs maxvol
ss_d   = Dict([("FR", (16256280.72,0,96118939.76,fr_per)), ("ES", (16409036.88,1,27075319.09, es_per)) ])

#
# The Thermal units
#
nb_thf = 0
(nT, nC) = size(gen_data)
for i=1:nT
    if ( gen_data.primary_fuel[i] in Thf_list )
        global nb_thf += 1
    end
end
# Make the thermal unit file from this
tu_thf_data = DataFrame(Zone=Vector{String}(undef, nb_thf), Name=Vector{String}(undef, nb_thf), NumberUnits=Vector{Int64}(undef, nb_thf),
                        MaxPower=Vector{Float64}(undef, nb_thf), MaxPowerProfile=Vector{String}(undef, nb_thf),	VariableCost=Vector{Float64}(undef, nb_thf),	FixedCost=Vector{Float64}(undef, nb_thf),	
                        InvestmentCost=Vector{Float64}(undef, nb_thf),	Capacity=Vector{Float64}(undef, nb_thf),	Energy=Vector{Float64}(undef, nb_thf),	
                       	MaxAddedCapacity=Vector{Float64}(undef, nb_thf) )

tu_thf_data[!,:NumberUnits] .= 1
tu_thf_data[!,:Energy] .= 0.0
tu_thf_data[!,:MaxPowerProfile] .= ""

i_thf = 0
for i=1:nT
    if ( gen_data.primary_fuel[i] in Thf_list )
        global i_thf += 1
        # Add this fellow - first check if it actually exists at some existing bus
        l0 = length(findall(bus_data.bus_id .== gen_data.bus_id[i]))
        if ( l0 == 0)
            error("The generator ", string(i)," is situated at some non existing bus:", string(gen_data.bus_id[i]))
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
tu_thf_data[!,:Capacity] .= tu_thf_data[!,:MaxPower]
tu_thf_data[!,:MaxAddedCapacity] .= ceil.( 0.1*tu_thf_data[!,:MaxPower] ) #round.( 0.1*tu_thf_data[!,:MaxPower]; digits=2 )

# Append the transformers to the set
append!(tu_thf_data, tu_imp_data)

#tu_thf_data[!,:Zone] .= "Nowhere"
CSV.write(string(github_local_d, github_smsppin,"/nutsx/TU_ThermalUnits.csv"), tu_thf_data; delim=';')

#
# Aggregate these into a NUTS0 version
#
tu_agg_d = similar(tu_thf_data,0)
# Delete de MaxPowerProfile column not useful for NUTS0 computations and add the MaxRetCapacity Column
select!(tu_agg_d, Not("MaxPowerProfile"))
insertcols!(tu_agg_d, :MaxRetCapacity=>Vector{Float64}(undef, 0) )

# For each country and each techno we add stuff
for cn in [focus_country] #unique(bus_data.country)
    Ibus = findall( bus_data.country .== cn )
    # Filter all generators in this country
    sub_t = filter(row -> row.bus_id in bus_data.bus_id[Ibus], gen_data)

    for i=1:length(Thf_list)
        Ig = findall( sub_t.primary_fuel .== Thf_list[i])

        cost_info = cThf_l[Thf_list[i]]
        Pmax = sum( sub_t.capacity_mw[Ig] )

        push!( tu_agg_d, [cn, Thf_list[i], 1, Pmax, cost_info[1], cost_info[2], cost_info[3], Pmax, 0.0, ceil.( 0.1*Pmax ), 0.0])
    end 
end
# Append the "base data" covering the other countries
tu_agg_elsewhere = CSV.read(string(github_local_d, github_smsppin,"/nuts0/TU_ThermalUnits_base.csv"), DataFrame; delim=';')
append!(tu_agg_d, tu_agg_elsewhere)
CSV.write(string(github_local_d, github_smsppin,"/nuts0/TU_ThermalUnits.csv"), tu_agg_d; delim=';')

#
# The RES units
#
nb_res = 0
for i=1:nT
    if ( gen_data.primary_fuel[i] in Res_list )
        global nb_res += 1
    end
end

# res DataFrame
res_units_data = DataFrame(Name=Vector{String}(undef, nb_res), Zone=Vector{String}(undef, nb_res), NumberUnits=Vector{Int64}(undef, nb_res),
                        MaxPower=Vector{Float64}(undef, nb_res),	MinPower=Vector{Float64}(undef, nb_res),	MaxPowerProfile=Vector{String}(undef, nb_res),	
                        Energy=Vector{Float64}(undef, nb_res),	Kappa=Vector{Float64}(undef, nb_res),	Capacity=Vector{Float64}(undef, nb_res),	
                       	MaxAddedCapacity=Vector{Float64}(undef, nb_res), MaxRetCapacity=Vector{Float64}(undef, nb_res), InvestmentCost=Vector{Float64}(undef, nb_res) )

res_units_data[!,:NumberUnits] .= 1
res_units_data[!,:Kappa] .= 1.0
res_units_data[!,:Energy] .= 0.0

i_res = 0
for i=1:nT
    if ( gen_data.primary_fuel[i] in Res_list )
        global i_res += 1
        # Add this fellow - first check if it actually exists at some existing bus
        l0 = length(findall(bus_data.bus_id .== gen_data.bus_id[i]))
        if ( l0 == 0)
            error("The generator ", string(i)," is situated at some non existing bus:", string(gen_data.bus_id[i]))
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
        res_units_data[i_res, :MaxPower] = gen_data.capacity_mw[i]*t_info[4]
    end
end
res_units_data[!,:MinPower] .= 0.0
res_units_data[!,:MaxRetCapacity] .= 0.0

#res_units_data[!,:Capacity] .= res_units_data[!,:MaxPower]
res_units_data[!,:MaxAddedCapacity] .= ceil.( 0.1*res_units_data[!,:MaxPower] )

CSV.write(string(github_local_d, github_smsppin,"/nutsx/RES_RenewableUnits.csv"), res_units_data; delim=';')

#
# Aggregate these into a NUTS0 version as well
#
res_agg_d = similar(res_units_data,0)
# For each country and each techno we add stuff
for cn in [focus_country] #unique(bus_data.country)
    Ibus = findall( bus_data.country .== cn )
    # Filter all generators in this country
    sub_t = filter(row -> row.bus_id in bus_data.bus_id[Ibus], gen_data)

    for i=1:length(Res_list)
        Ig = findall( sub_t.primary_fuel .== Res_list[i])

        t_info = cRes_l[Res_list[i]]
        Pmax = sum( sub_t.capacity_mw[Ig] )

        push!( res_agg_d, [Res_list[i], cn, 1, Pmax*t_info[4], 0.0, string(t_info[2], cn, t_info[3]), 0.0, 1.0, Pmax, ceil.( 0.1*Pmax ), 0.0, t_info[1] ])
    end 
end
# Append the "base data" covering the other countries
res_agg_elsewhere = CSV.read(string(github_local_d, github_smsppin,"/nuts0/RES_RenewableUnits_base.csv"), DataFrame; delim=';')
append!(res_agg_d, res_agg_elsewhere)
CSV.write(string(github_local_d, github_smsppin,"/nuts0/RES_RenewableUnits.csv"), res_agg_d; delim=';')

#
# Now the Pumped Storage file
#
nb_sts = 0
for i=1:nT
    if ( gen_data.primary_fuel[i] in STS_list )
        global nb_sts += 1
    end
end

sts_data = DataFrame(Name=Vector{String}(undef, nb_sts), Zone=Vector{String}(undef, nb_sts), NumberUnits=Vector{Int64}(undef, nb_sts),
                        MaxPower=Vector{Float64}(undef, nb_sts), MaxVolume=Vector{Float64}(undef, nb_sts), TurbineEfficiency=Vector{Float64}(undef, nb_sts),
                        PumpingEfficiency=Vector{Float64}(undef, nb_sts), MinPower=Vector{Float64}(undef, nb_sts), MinVolume=Vector{Float64}(undef, nb_sts),	
                        Energy=Vector{Float64}(undef, nb_sts), Inflows=Vector{Float64}(undef, nb_sts), InitialVolume=Vector{Float64}(undef, nb_sts),
                        AddPumpedStorage=Vector{Float64}(undef, nb_sts),	
                       	MaxAddedCapacity=Vector{Float64}(undef, nb_sts), MaxRetCapacity=Vector{Float64}(undef, nb_sts), InvestmentCost=Vector{Float64}(undef, nb_sts) )

sts_data[!,:NumberUnits] .= 1
sts_data[!,:MinVolume] .= 0.0
sts_data[!,:Energy] .= 0.0
sts_data[!,:Inflows] .= 0.0
sts_data[!,:InitialVolume] .= 0.0
sts_data[!,:AddPumpedStorage] .= 0.0
sts_data[!,:MaxAddedCapacity] .= 0.0
sts_data[!,:MaxRetCapacity] .= 0.0
sts_data[!,:InvestmentCost] .= 0.0

i_sts = 0
for i=1:nT
    if ( gen_data.primary_fuel[i] in STS_list )
        global i_sts += 1
        # Add this fellow - first check if it actually exists at some existing bus
        l0 = length(findall(bus_data.bus_id .== gen_data.bus_id[i]))
        if ( l0 == 0)
            error("The STS unit ", string(i)," is situated at some non existing bus:", string(gen_data.bus_id[i]))
        end
        i0 = findall(bus_data.bus_id .== gen_data.bus_id[i])

        sts_data[i_sts, :Zone] = bus_data.bus_id[i0][1] #string.( bus_data.country[i0], "_", bus_data.bus_id[i0] )[1]
        sts_data[i_sts, :Name] = gen_data.unit_id[i] #gen_data.name[i]

        sts_data[i_sts, :MaxPower] = gen_data.capacity_mw[i]
        sts_data[i_sts, :MinPower] = -1.0*gen_data.capacity_mw[i]

        tech_info = cSTS_l[gen_data.primary_fuel[i]]
        sts_data[i_sts, :MaxVolume] = tech_info[1]*gen_data.capacity_mw[i]
        sts_data[i_sts, :TurbineEfficiency] = tech_info[2]
        sts_data[i_sts, :PumpingEfficiency] = tech_info[3]
    end
end
CSV.write(string(github_local_d, github_smsppin,"/nutsx/STS_ShortTermStorage.csv"), sts_data; delim=';')

#
# Aggregate into a NUTS0 version as well
#
sts_agg_d = similar(sts_data,0)
# For each country and each techno we add stuff
for cn in [focus_country] #unique(bus_data.country)
    Ibus = findall( bus_data.country .== cn )
    # Filter all generators in this country
    sub_t = filter(row -> row.bus_id in bus_data.bus_id[Ibus], gen_data)

    for i=1:length(STS_list)
        Ig = findall( sub_t.primary_fuel .== STS_list[i])

        tech_info = cSTS_l[STS_list[i]]
        Pmax = sum( sub_t.capacity_mw[Ig] )

        push!( sts_agg_d, [STS_list[i], cn, 1, Pmax, tech_info[1]*Pmax, tech_info[2], tech_info[3], -1.0*Pmax, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ])
    end 
end
# Append the "base data" covering the other countries
sts_agg_elsewhere = CSV.read(string(github_local_d, github_smsppin,"/nuts0/STS_ShortTermStorage_base.csv"), DataFrame; delim=';')
append!(sts_agg_d, sts_agg_elsewhere)
CSV.write(string(github_local_d, github_smsppin,"/nuts0/STS_ShortTermStorage.csv"), sts_agg_d; delim=';')

#
# The Seasonal Storage units
#
nb_ss = 0
nb_ss_max_d = Dict([ ("FR",0.0), ("ES",0.0) ])
for i=1:nT
    if ( gen_data.primary_fuel[i] in SS_list )
        global nb_ss += 1
        # first check if it actually exists at some existing bus
        l0 = length(findall(bus_data.bus_id .== gen_data.bus_id[i]))
        if ( l0 == 0)
            error("The SS unit ", string(i)," is situated at some non existing bus:", string(gen_data.bus_id[i]))
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
                        AddPumpedStorage=Vector{Float64}(undef, nb_ss), WaterValues=Vector{String}(undef, nb_ss)  )

ss_data[!,:NumberUnits] .= 1
ss_data[!,:MinVolume] .= 0.0
ss_data[!,:MinPower] .= 0.0
ss_data[!,:TurbineEfficiency] .= 1.0
ss_data[!,:PumpingEfficiency] .= 0.0
ss_data[!,:AddPumpedStorage] .= 0.0
ss_data[!,:WaterValues] .= "bellman_nutsx.csv"

i_ss = 0
for i=1:nT
    if ( gen_data.primary_fuel[i] in SS_list )
        global i_ss += 1
        # Add this fellow - first check if it actually exists at some existing bus
        l0 = length(findall(bus_data.bus_id .== gen_data.bus_id[i]))
        if ( l0 == 0)
            error("The SS unit ", string(i)," is situated at some non existing bus:", string(gen_data.bus_id[i]))
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
        ss_data[i_ss, :MaxVolume] = country_info[1] *(gen_data.capacity_mw[i]/mx_tot)

        ss_data[i_ss, :Inflows] = country_info[3] *(gen_data.capacity_mw[i]/mx_tot)
        ss_data[i_ss, :InflowsProfile] = string(tech_info[2], bus_data.country[i0][1], tech_info[3])
        ss_data[i_ss, :InitialVolume] = country_info[4]*ss_data[i_ss, :MaxVolume]         
    end
end
CSV.write(string(github_local_d, github_smsppin,"/nutsx/SS_SeasonalStorage.csv"), ss_data; delim=';')

nb_cuts = size(bell_data)[1]
bell_nutsx_data = DataFrame(Timestep=Vector{Int64}(undef, nb_cuts) )
for i=1:nb_ss
    bell_nutsx_data[!, string("a_",string(i-1))] .= 0.0
end
bell_nutsx_data[!, string("b")] .= 0.0

bell_nutsx_data[!, :Timestep] .= bell_data[!,:Timestep]

## Filter out the target country => assumed to be some Hydrosystem, here the first one ;
stage = unique( bell_data.Timestep );
for st in stage
    t_step = (st+1)*ssv_step

    I = findall( bell_data.Timestep .== st )
    b_rhs = bell_data.b[I]

    # Loop over countries with Hydro but not focus country
    for cn in setdiff(keys(ss_d),[focus_country])
        country_info = ss_d[cn]
        h_idx = country_info[2]
        
        # a coefficients
        a_v = bell_data[I,2+h_idx]
        # volume at the end 
        vol = vol_data[t_step, 2 + h_idx]

        b_rhs .+= a_v*vol
    end
    h_idx_focus = ss_d[focus_country][2]
    bell_nutsx_data[I,2:1+nb_ss] .= bell_data[I, 2 + h_idx_focus]
    bell_nutsx_data[I,:b] .= b_rhs
end
CSV.write(string(github_local_d, github_smsppin,"/ts/bellman_nutsx.csv"), bell_nutsx_data; delim=';')

# Check if the all units have been handled
if ( nb_ss + nb_thf + nb_res + nb_sts != nT )
    error("Some units have been lost...")
end

# 
# Some entries for SettingsCreate
#
function strVtoTable(d::Vector{String})
    s = "[";
    for el in d
        s = string(s,"'",el,"',")
    end
    s = string(chop(s),"]")
    return s
end

println("    Countries: ", strVtoTable(zv_zone_data[:,:Zone]) )
#
println("    thermal: ", strVtoTable(unique(tu_thf_data[:,:Name])) )
#
println("    reservoir: ", strVtoTable(unique(ss_data[:,:Name])) )
#

# Type conversion on the column
gen_data[!, :unit_id] = convert.(String, gen_data[:, :unit_id])
sub_ror = filter( row -> row.primary_fuel .== "Hydro|run_of_river", gen_data )
sub_res = filter( row -> ( (row.primary_fuel .== "Wind") || (row.primary_fuel .== "Solar") ), gen_data )

println("    hydrostorage: ", strVtoTable(unique(sts_data[:,:Name])) )
println("    battery: ['Energy Storage System|Compressed Air','Energy Storage System|Lithium-Ion']")
println("    res: ", strVtoTable(unique(sub_res[:,:unit_id])) )
println("    runofriver: ", strVtoTable(unique(sub_ror[:,:unit_id])) )

# Using the direct UC script to write the data
#
#
#@time begin
smspp_bname = "spain_block"    
if ( with_dcopf && !with_acopf )
    smspp_bname = string(smspp_bname, "_dcopf")
end
if ( with_acopf )
    smspp_bname = string(smspp_bname, "_acopf")
end
smspp_bname = string(smspp_bname, "_", string(nb_sstep), "h")
uc_b_dt = uc_bgn_date
uc_e_dt = uc_b_dt + Dates.Hour(nb_sstep)
for iclip=1:nb_clip
    println(string("Generating data for time : ", string(uc_b_dt), " to ", string(uc_e_dt)))
    #
    smspp_bname_l = string(github_local_d, github_smsppin,"/nutsx/", smspp_bname, "_", string(iclip), ".txt")
    #write_smspp_file(smspp_bname, "C:/LocalDriveD/Tools/Spain/TimeSeries", uc_bgn_date, uc_end_date, idx_scen, st_idx, zp_data, zv_zone_data, incon_data, tu_thf_data, res_units_data, sts_data, ss_data )
    write_smspp_file(smspp_bname_l, string(github_local_d, github_smsppin,"/ts"), uc_b_dt, uc_e_dt, idx_scen, st_idx, tan_phi, zp_data, zv_zone_data, incon_data, tu_thf_data, res_units_data, sts_data, ss_data )
    #
    global uc_b_dt = uc_b_dt + Dates.Hour(nb_sstep)
    global uc_e_dt = uc_e_dt + Dates.Hour(nb_sstep)
end
#end

# Make a file for POMATO
nb_pomato = nb_thf + nb_res + nb_sts + nb_ss
pomato_data = DataFrame(index=Vector{String}(undef, nb_pomato), node=Vector{String}(undef, nb_pomato), 
                        mc=Vector{Float64}(undef, nb_pomato), mc_heat=Vector{Float64}(undef, nb_pomato),
                        g_max=Vector{Float64}(undef, nb_pomato), eta=Vector{Float64}(undef, nb_pomato),
                        plant_type=Vector{String}(undef, nb_pomato), storage_power=Vector{Float64}(undef, nb_pomato),
                        storage_capacity=Vector{Float64}(undef, nb_pomato) )
for i=1:nb_thf
    pomato_data[i,:index]= tu_thf_data.Name[i]
    pomato_data[i,:node]= tu_thf_data.Zone[i]
    pomato_data[i,:mc]  = tu_thf_data.VariableCost[i]
    pomato_data[i,:mc_heat]=0.0
    pomato_data[i,:g_max]=tu_thf_data.MaxPower[i]
    pomato_data[i,:eta]=1.0
    i0 = findall(gen_data.unit_id .== tu_thf_data[i,:Name])
    pomato_data[i,:plant_type]=gen_data.primary_fuel[i0[1]]
    pomato_data[i,:storage_power]=0
    pomato_data[i,:storage_capacity]=0
end
for i=1:nb_res
    ip = nb_thf + i
    pomato_data[ip,:index]= res_units_data.Name[i]
    pomato_data[ip,:node]= res_units_data.Zone[i]
    pomato_data[ip,:mc]  = 0.0
    pomato_data[ip,:mc_heat]=0.0
    pomato_data[ip,:g_max]=res_units_data.MaxPower[i]
    pomato_data[ip,:eta]=1.0
    i0 = findall(gen_data.unit_id .== res_units_data[i,:Name])
    pomato_data[ip,:plant_type]=gen_data.primary_fuel[i0[1]]
    pomato_data[ip,:storage_power]=0
    pomato_data[ip,:storage_capacity]=0
end
for i=1:nb_sts
    ip = nb_thf + nb_res + i
    pomato_data[ip,:index]= sts_data.Name[i]
    pomato_data[ip,:node]= sts_data.Zone[i]
    pomato_data[ip,:mc]  = 0.0
    pomato_data[ip,:mc_heat]=0.0
    pomato_data[ip,:g_max]=sts_data.MaxPower[i]
    pomato_data[ip,:eta]=sts_data.PumpingEfficiency[i]
    i0 = findall(gen_data.unit_id .== sts_data[i,:Name])
    pomato_data[ip,:plant_type]=gen_data.primary_fuel[i0[1]]
    pomato_data[ip,:storage_power]=sts_data.MaxPower[i]
    pomato_data[ip,:storage_capacity]=sts_data.MaxVolume[i]
end
for i=1:nb_ss
    ip = nb_thf + nb_res + nb_sts + i
    pomato_data[ip,:index]= ss_data.Name[i]
    pomato_data[ip,:node]= ss_data.Zone[i]

    subbell = filter(row->(row.Timestep == st_idx), bell_nutsx_data)

    pomato_data[ip,:mc]  = mean(abs.(subbell[:,"a_0"]))
    pomato_data[ip,:mc_heat]=0.0
    pomato_data[ip,:g_max]=ss_data.MaxPower[i]
    pomato_data[ip,:eta]=ss_data.PumpingEfficiency[i]
    i0 = findall(gen_data.unit_id .== ss_data[i,:Name])
    pomato_data[ip,:plant_type]=gen_data.primary_fuel[i0[1]]
    pomato_data[ip,:storage_power]=ss_data.MaxPower[i]
    pomato_data[ip,:storage_capacity]=ss_data.MaxVolume[i]
end
CSV.write(string(github_local_d, github_smsppin,"/nutsx/Pomato_input.csv"), pomato_data; delim=';')

#
#
#
println(" ---------------------------------------------------------------------------------------------------------- ")
println(" ------------------------------Observations on Usage------------------------------------------------------- ")
println(" ---------------------------------------------------------------------------------------------------------- ")
println(" > Run the NUTS0 SMS++ input files through the PLAN4RES python scripts")
println(" -- Within the properly installed plan4res environment -- ")
println(" To generate the netcdf files for SSV : ")
println(" > ./plan4res.sh FORMAT dataSetName simul ")
println(" To run SSV : ")
println(" > ./plan4res.sh SSV dataSetName simul ")
println(" To run the simulation and obtain the output in the smspp_out folder : ")
println(" > ./plan4res.sh SIM dataSetName")
println(" Then copy the output in the smspp_out/results_simul folder ")
println(" ---------------------------------------------------------------------------------------------------------- ")



