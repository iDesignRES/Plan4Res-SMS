using CSV
using DataFrames
using Dates
using Colors
using Statistics

# For plots
using GLMakie
GLMakie.activate!()

# Import the needed data
#
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

# Current Results Name
res_type = "flows"
res_name = string("results_",res_type)

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

#
# output directory
#
r_dir       = string(github_local_d, github_smsppout, "/nutsx/", res_name)
outfile_ext = "OUT"

with_acopf = false
ts_prod = CSV.read(string(r_dir, "/ActivePower/ActivePower", outfile_ext, ".csv" ), DataFrame; delim=',')
mx_prod = CSV.read(string(r_dir, "/MaxPower/MaxPower", outfile_ext, ".csv" ), DataFrame; delim=',')
ts_flow = CSV.read(string(r_dir, "/Flows/Flows", outfile_ext, ".csv" ), DataFrame; delim=',')
q_file  = string(r_dir, "/Flows/FlowsImag", outfile_ext, ".csv" )
if ( isfile(q_file) )
    with_acopf = true
    ts_flowQ= CSV.read(q_file, DataFrame; delim=',')
end
# If the file exists but LineRATEA is not in the incond_data file, finally we are not dealing with ACOPF data
if ( "LineRATEA" ∉ names(incon_data) )
    with_acopf = false
end
ts_dem  = CSV.read(string(r_dir, "/Demand/Demand", outfile_ext, ".csv" ), DataFrame; delim=',')

nbgen = length(gen_data.unit_id)

fig = Figure()

xx = ts_prod.Timestep
ax = Axis(fig[1, 1], title="Empty")

sl_x = Slider(fig[2, 1], range = 1:1:nbgen, startvalue = 1)

gl = GridLayout(fig[3, 1], tellwidth = false)
toggle = Toggle(gl[1, 1], active = false)
Label(gl[1,2], "Emergency Clean")

tb = Textbox(fig[4,1], placeholder = "Enter a plant Name", tellwidth = false)
# Faire quelque chose avec le textbox
on(tb.stored_string) do s
    if (s in gen_data.unit_id)
        i0 = findall( gen_data.unit_id .== s )
        empty!(ax)
        lines!(ax, xx, ts_prod[!,s], color=:blue)
        lines!(ax, xx, mx_prod[!,s], color=:red)
        ax.title = string("Generation of Plant " ,s, " [", gen_data.primary_fuel[i0[1]], "]")
    end
end

# Code to activate the toggle
on(fig.scene.events.tick) do tick
    toggle.active[] || return
    #empty!(ax)
    #t[] += tick.delta_time
end

# Faire une action avec le Toggle
on(toggle.active) do val
    if val
        empty!(ax)
    end
end

# Fonction pour mettre à jour le plot en fonction de la colonne sélectionnée
function update_plot!(col)
    lines!(ax, xx, ts_prod[!,gen_data.unit_id[Int(col)]], color=:blue)
    lines!(ax, xx, mx_prod[!,gen_data.unit_id[Int(col)]], color=:red)
    ax.title = string("Generation of Plant " ,gen_data.unit_id[Int(col)], " [", gen_data.primary_fuel[Int(col)], "]")
    #Label(fig, "base")
    #title=string("Generation of Plant " ,gen_data.unit_id[Int(col)], " [", gen_data.primary_fuel[Int(col)], "]"))
    #title!(ax, string("Generation of Plant " ,gen_data.unit_id[Int(col)], " [", gen_data.primary_fuel[Int(col)], "]"))
end

# Connexion du slider à la fonction de mise à jour
on(sl_x.value) do col
    empty!(ax)
    update_plot!(Int(col))
end

# Affichage initial
update_plot!(1)

# Make another pair of axis for load
axd = Axis(fig[5, 1], title="Empty")
tbd = Textbox(fig[6,1], placeholder = "Enter a node Name", tellwidth = false)
# Faire quelque chose avec le textbox
on(tbd.stored_string) do s
    if (s in bus_data.bus_id)
        #i0 = findall( gen_data.unit_id .== s )
        empty!(axd)
        lines!(axd, xx, ts_dem[!,s], color=:blue)
        #lines!(ax, xx, mx_prod[!,s], color=:red)
        axd.title = string("Demand of Node " ,s)
    end
end

# Affichage de la figure
display(fig)

fig_name = string(r_dir, "/", "spain_some_bus.png")
#save(fig_name, fig ; px_per_unit=4)



#
#
function total_gen(fg, g_data, ts_prd, mx_prd)
"""

"""
ax = Axis(fg[7, 1], title="Generation by Tech")

xx = ts_prd.Timestep

tech_list = unique(g_data.primary_fuel)

for tech in tech_list
    tc_prod = Vector{Float64}(undef,length(xx))
    tmx_prod = Vector{Float64}(undef,length(xx))
    tc_prod .= 0.0
    tmx_prod .= 0.0
    s_g_data = filter( row -> (row.primary_fuel .== tech), g_data)
    # See if the technology field is filled in:
    sb_t_f = filter( row -> ( ismissing(row.technology) ), s_g_data)
    if ( length(sb_t_f.unit_id) > 0 )
        for uid in s_g_data.unit_id
            tc_prod += ts_prd[:,uid]        
        end
        lines!(ax, xx, tc_prod, label=[tech])
    else
        sub_tech = unique( s_g_data.technology )
        for stech in sub_tech
            tc_prod .= 0.0
            s_s_g_data = filter( row -> (row.technology .== stech), s_g_data )
            for uid in s_s_g_data.unit_id
                tc_prod += ts_prd[:,uid]        
            end
            lines!(ax, xx, tc_prod, label=[string(tech,"|",stech)])
        end
    end

    # Compute mx prod with availability factor
    for uid in s_g_data.unit_id
        tmx_prod += mx_prd[:,uid]        
    end

    println( string("Total installed capacity of this tech=", tech, " = ", sum(s_g_data.capacity_mw), " average max available ", mean(tmx_prod) ) )

end
axislegend(ax, position=:cb, orientation = :horizontal)

end

function total_res(fg, g_data, ts_prd, mx_prd)
    """
    
    """

    ax = Axis(fg[7, 1], title="Generation by Tech")
    
    xx = ts_prd.Timestep
    
    tech_list = ["Solar", "Wind"]
    
    for tech in tech_list
        tc_prod = Vector{Float64}(undef,length(xx))
        tmx_prod = Vector{Float64}(undef,length(xx))
        tc_prod .= 0.0
        tmx_prod .= 0.0
        s_g_data = filter( row -> (row.primary_fuel .== tech), g_data)
        for uid in s_g_data.unit_id
            tc_prod += ts_prd[:,uid]
            tmx_prod += mx_prd[:,uid]        
        end
        lines!(ax, xx, tc_prod, label=[tech])
        lines!(ax, xx, tmx_prod, label=[string(tech,"|max")])
    end
    axislegend(ax, position=:cb, orientation = :horizontal)
end