using CSV
using DataFrames
using Dates
using Colors

# For plots
using GMT

cur_dir  = pwd()
work_dir = "C:/LocalDriveD/Tools/SMSpp/julia"

# Import the needed data
#
#
# File names as put together by NTNU
#
ntnu_bus_file  = "Bus_Data.csv"
ntnu_line_file = "lines.csv"
ntnu_gen_file  = "generation.csv"
ntnu_load_file = "load.csv"

# Load the data of the buses
bus_data  = CSV.read(ntnu_bus_file, DataFrame; delim=',')

# Load the data of the lines
lines_data = CSV.read(ntnu_line_file, DataFrame; delim=',')

# Load the generator data
gen_data = CSV.read(ntnu_gen_file, DataFrame; delim=',')
# Do some cleaning on the names
# if The generator data has a unit_id field, then we can replace the names with that one
if ( "unit_id" in names(gen_data) )
    sbgd = filter( row -> ( ismissing(row.name) ), gen_data)
    for uid in sbgd.unit_id
        #println(uid)
        i0 = findall( gen_data.unit_id .== uid )
        gen_data.name[i0[1]] = uid
    end
end
#
#new_names = replace.(gen_data[:,:name],"Á"=> "A", "É"=>"E", "Í"=> "I", "Ñ"=>"N", "Ó"=>"O","Ú" => "U", "Ü"=>"U")
#gen_data[:,:name] .= new_names

incon_data = CSV.read("smspp_in/IN_Interconnections.csv", DataFrame; delim=';')

pbalNames = ["lightcyan", "paleturquoise1", "cyan", "darkturquoise", "turquoise4"]
for cnom in pbalNames
    if ( cnom ∉ keys(Colors.color_names) )
        error(string("undefined color : ", cnom))
    end
end

# Available colours
# https://juliagraphics.github.io/Colors.jl/stable/namedcolors/
# Colors.color_names:

# A color range for a heatmap kind of idea 
cRange = range(colorant"green", stop=colorant"red", length=11);
cNames = ["chartreuse4", "forestgreen", "green", "darkolivegreen", "darkgoldenrod4", "darkorange4", "orangered4", "firebrick", "red3", "red2", "red"]
for cnom in cNames
    if ( cnom ∉ keys(Colors.color_names) )
        error(string("undefined color : ", cnom))
    end
end
if (length(cRange) != length(cNames) )
    error("not allowed yet")
end

#
#
println(string("Found the following possible voltage levels", string(unique(lines_data.voltage)) ) )

v2c = Dict([(132, "papayawhip"), (220, "peachpuff1"), (225, "peachpuff2"), (250, "burlywood3"), (320, "tan3"), (400, "chocolate3") ])

#
# output directory
#
b_dir       = "//Atlas.edf.fr/in/retd-00ccr/D34552/Projects/iDesignRES/Spain_data/entsogridkit"
r_dir       = string(b_dir, "/", "smspp_in/nutsx/new_data_res_dcopf")
#r_dir       = string(b_dir, "/", "smspp_in/nutsx/results_dcopf")
outfile_ext = "OUT"

with_acopf = false
ts_prod = CSV.read(string(r_dir, "/ActivePower/ActivePower", outfile_ext, ".csv" ), DataFrame; delim=',')
ts_flow = CSV.read(string(r_dir, "/Flows/Flows", outfile_ext, ".csv" ), DataFrame; delim=',')
q_file  = string(r_dir, "/Flows/FlowsImag", outfile_ext, ".csv" )
if ( isfile(q_file) )
    with_acopf = true
    ts_flowQ= CSV.read(q_file, DataFrame; delim=',')
end
ts_dem  = CSV.read(string(r_dir, "/Demand/Demand", outfile_ext, ".csv" ), DataFrame; delim=',')

Im66    = Vector{Int64}(undef,0)
Im33    = Vector{Int64}(undef,0)
Inaught = Vector{Int64}(undef,0)
Ip33    = Vector{Int64}(undef,0)
Ip66    = Vector{Int64}(undef,0)
Idef    = Vector{Int64}(undef,0)
for nd in bus_data.bus_id
    i0 = findall( bus_data.bus_id .== nd )
    nodeNm = nd #string(bus_data.country[i0[1]], "_", nd)
    
    # Find the demand in this node
    dem = ts_dem[:, nodeNm]

    sgen = filter(row -> (row.bus_id .== nd), gen_data)
    totalP = Vector{Float64}(undef, size(ts_prod)[1] )
    totalP .= 0.0
    for g in sgen.unit_id
        totalP += ts_prod[:,g]
    end
    #
    ODEbal = dem - totalP

    imb_unit = string("SlackUnit_", nodeNm)
    def      = ts_prod[:,imb_unit]

    avg_d = mean(dem)
    avg_o = mean(ODEbal)
    avg_def = mean(def)

    #println(string(nd, " avg_d = ", string(avg_d), " avg_o = ", string(avg_o), " r = ", avg_o/avg_d))

    if ( ( avg_o / avg_d >= -0.33 ) && ( avg_o / avg_d <= 0.33 ) )
        append!(Inaught, i0)
    elseif ( ( avg_o / avg_d >= -0.66 ) && ( avg_o / avg_d <= -0.33 ) )
        append!(Im33, i0)
    elseif ( ( avg_o / avg_d >=  0.33 ) && ( avg_o / avg_d <= 0.66 ) )
        append!(Ip33, i0)
    elseif ( ( avg_o / avg_d <= -0.66 ) )
        append!(Im66, i0)
    elseif ( ( avg_o / avg_d >= 0.66 ) || ( abs(avg_d) < 1e-6 ) )
        append!(Ip66, i0)
    end

    if ( avg_def > 0.03*avg_d )
        append!(Idef, i0)
    end
end
println("The buses with imb are ", bus_data.bus_id[Idef] )

for i_d in Idef
    nd_name = bus_data.bus_id[i_d]
    imb_unit = string("SlackUnit_", nd_name)
    def      = ts_prod[:,imb_unit]
    println(string("The bus ", nd_name, " has max imbalance ", maximum(abs.(def))))
end

lost_bus=setdiff(collect(1:size(bus_data)[1]), union(Im66,Im33,Inaught,Ip33,Ip66))
if ( length(lost_bus) > 0)
    println(string("The following buses have been lost : ", string(lost_bus)))
end

#nodeNm = "ES00101"
#sgen = filter(row -> (row.bus_id .== nodeNm), gen_data)
#totalP = Vector{Float64}(undef, size(ts_prod)[1] )
#totalP .= 0.0
#for g in sgen.unit_id
#    global totalP += ts_prod[:,g]
#end


# Start plotting -> we need to change the working directory for some obscure reason

cd( work_dir )

gmtbegin("Spain", fmt=:png)

coast(
    region="-9/35/4/45+r", 
    proj=(name=:laea, center=[3,39]), 
    frame=:ag, 
    res=:full,
    area=500, 
    shore=:thin, 
    rivers=:thin,
    borders=1,
    
    # outline germany in red
    DCW=((country="ES", 
            pen=(1,:turquoise3))),
    figsize=100
)

latt = bus_data.x #bus_data.y
long = bus_data.y #bus_data.x

# Plots with export and import buses : colours : pbalNames
ndCols = [(Im66, pbalNames[1]), (Im33, pbalNames[2]), (Inaught, pbalNames[3]), (Ip33, pbalNames[4]), (Ip66, pbalNames[5])]

for ndc in ndCols
    (Ic, ndcolour) = ndc
    GMT.scatter!(
        latt[Ic], long[Ic], 
        fmt=:png, 
        marker=:circle,
        markeredgecolor=0, 
        size=0.25, 
        #markerfacecolor=:orchid,
        markerfacecolor=ndcolour, 
    )
end

# Defaillance
if ( length(Idef) > 0)
    GMT.scatter!(
        latt[Idef], long[Idef], 
        fmt=:png, 
        marker=:cross,
        markeredgecolor=0, 
        size=0.20, 
        markerfacecolor=:red, 
    )
end

for vt in unique(lines_data.voltage)
    
    Ivolt = findall( lines_data.voltage .== vt )
    xLines = Array{Float64}(undef, 0, 2)
    yLines = Array{Float64}(undef, 0, 2)
    cSc = Vector{Int64}(undef,0)
    
    for iln in Ivolt
        i0 = findall( bus_data.bus_id .== lines_data.bus0[iln])
        i1 = findall( bus_data.bus_id .== lines_data.bus1[iln])

        #if ( string("L", lines_data.line_id[iln]) ∉ names(ts_flow) )
        #    error(string("Line with line_id ", lines_data.line_id[iln], " not found in results "))
        #end
        # Some lines are doubled and have been handled as a single one
        if ( string("L", lines_data.line_id[iln]) in names(ts_flow) )
            i2 = findall( incon_data.Name .== string("L", lines_data.line_id[iln]) )
            if ( !with_acopf )
                mx_sat  = incon_data.MaxPowerFlow[i2[1]]
                avg_sat = mean(ts_flow[:, string("L", lines_data.line_id[iln])])
            else
                mx_sat  = incon_data.LineRATEA[i2[1]]
                avg_sat = mean(sqrt.(ts_flow[:, string("L", lines_data.line_id[iln])].^2 + ts_flowQ[:, string("L", lines_data.line_id[iln])].^2 ))
            end

            xLines = vcat(xLines, [ bus_data.x[i0[1]] bus_data.x[i1[1]] ] )
            yLines = vcat(yLines, [ bus_data.y[i0[1]] bus_data.y[i1[1]] ] )
            append!( cSc, [Int(round((length(cRange)-1)*(abs(avg_sat)/mx_sat)))+1] )

            if ( abs(avg_sat)/mx_sat > 1 )
                println(string(" L name ", lines_data.line_id[iln], " avg = ", string(avg_sat), " mx ", string(mx_sat)))
            end
        end
    end
    #println( string( cSc ))
    nb_lns = size(xLines)[1]

    for i=1:nb_lns
        GMT.lines!(
            xLines[i, :], yLines[i, :],
            fmt=:png,
            lw=1,
            #lc=:chocolate3,
            #lc=v2c[vt],
            lc=cNames[cSc[i]],
            #lc=(Int64(round(256*red(cRange[cSc[i]])))/255, Int64(round(256*green(cRange[cSc[i]])))/255, Int64(round(256*blue(cRange[cSc[i]])))/255 ),
            #lc=(min(Int64(round(256*red(cRange[cSc[i]]))),255), min(Int64(round(256*green(cRange[cSc[i]]))),255), min(Int64(round(256*blue(cRange[cSc[i]]))),255) ),
            #lc=(Int8(round(256*red(cRange[cSc[i]]))), Int8(round(256*green(cRange[cSc[i]]))), Int8(round(256*blue(cRange[cSc[i]]))) ),
            #lc= (238, 118, 33), #(red(cRange[cSc[i]]),green(cRange[cSc[i]]),blue(cRange[cSc[i]]) ),
        )
    end
end

gmtend(show=true)




cd(cur_dir)