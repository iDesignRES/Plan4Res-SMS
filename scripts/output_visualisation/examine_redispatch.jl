#
#
# Check some plants and the Redispatch
using CSV
using DataFrames
using GLMakie
#using Makie.GeometryBasics

# Import the needed data
#
github_local_d = string(@__DIR__, "/../..")
github_smsppin = "/smspp_in"
github_smsppout = "/smspp_out"
#
# File names as put together by NTNU
#
ntnu_gen_file  = string(github_local_d, "/input_data/","generation.csv") #"Generators.csv"

# Current Results Name
res_type = "rdacopf"
res_name = string("results_",res_type)

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

# Import Pomatwo stuff
pomatwo_dir     = string(github_local_d, "/../POMATWO/results/dayahead")
pomatwo_res     = string(pomatwo_dir, "/pomatwo_DA_results_GEN.csv")
res_pomatwo_t = CSV.read(pomatwo_res, DataFrame; delim=',')
res_pomatwo   = unstack(res_pomatwo_t, :Time, :index,:GEN)

# Import output
r_dir       = string(github_local_d, github_smsppout, "/nutsx/", res_name)
outfile_ext = "OUT"
ts_prod = CSV.read(string(r_dir, "/ActivePower/ActivePower", outfile_ext, ".csv" ), DataFrame; delim=',')

nbT = size(ts_prod)[1]
dt = collect(1:nbT)
nbperLine = 7;

#Makie.inline!(true)
fig = Figure(size=(700, 700))
ga = fig[1, 1]

global line_number = 1
Tlist = unique(gen_data.technology);

for t in Tlist
    I = findall( gen_data.technology .== t)
    genIds = gen_data.unit_id[I];

    dev_prod = Vector{Float64}(undef,nbT)
    dev_prod .= 0.0

    for id in genIds
        dev_prod += abs.( ts_prod[:, id] - res_pomatwo[:,id]);
    end
    println( string("techno", t , maximum(dev_prod)))
    axis = Axis(ga[mod(line_number-1, nbperLine)+1,div(line_number-1,nbperLine)+1], ylabel=string(t),xlabel="Time Steps", xticks = 1:nbT)
    scatterlines!(axis, 1:nbT, dev_prod)
    global line_number += 1
end
display(fig)
fig_name = string(r_dir, "/", "redispatch_", res_type, ".png")
save(fig_name, fig)
