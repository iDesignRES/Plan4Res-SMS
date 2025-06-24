using CSV
using DataFrames

# GitHub Repository structure:
#
# -- Absolute or relative Path
github_local_d = string(@__DIR__, "/../..")
fcast_dir  = string(github_local_d, "/../NTNU/ES")

wind_base_file  = string(fcast_dir, "/Wind_Onshore/2_12_2024.csv")
solar_base_file = string(fcast_dir, "/Solar/2_12_2024.csv")
load_base_file  = string(fcast_dir, "/load/2_12_2024.csv")
wfct = CSV.read(wind_base_file, DataFrame; delim=',')
sfct = CSV.read(solar_base_file, DataFrame; delim=',')
lfct = CSV.read(load_base_file, DataFrame; delim=',')

wind_norm_factor = 30159
solr_norm_factor = 23867
load_norm_factor = 228746533 #remember total annual value

fct_frame = DataFrame(Zone=Vector{String}(undef, 24), Solar=Vector{Float64}(undef, 24), Wind=Vector{Float64}(undef, 24), Load=Vector{Float64}(undef, 24) )
fct_frame.Zone .= "ES"

# DA
fct_frame[:,:Wind]  .= wfct[:,2]./wind_norm_factor
fct_frame[:,:Solar] .= sfct[:,2]./solr_norm_factor
fct_frame[:,:Load]  .= lfct[:,2]./load_norm_factor
# Save
fct_fname = string("avail_ID_d2_da",".csv")
CSV.write(string(fcast_dir, "/forecast/",fct_fname), fct_frame; delim=',')

# g1
fct_frame[:,:Wind] .= wfct[:,3]./wind_norm_factor
fct_frame[:,:Solar] .= sfct[:,3]./solr_norm_factor
fct_frame[:,:Load]  .= lfct[:,3]./load_norm_factor
# Save
fct_fname = string("avail_ID_d2_g1",".csv")
CSV.write(string(fcast_dir, "/forecast/",fct_fname), fct_frame; delim=',')

# g2
fct_frame[:,:Wind]  .= wfct[:,4]./wind_norm_factor
fct_frame[:,:Solar] .= sfct[:,4]./solr_norm_factor
fct_frame[:,:Load]  .= lfct[:,4]./load_norm_factor
# Save
fct_fname = string("avail_ID_d2_g2",".csv")
CSV.write(string(fcast_dir, "/forecast/",fct_fname), fct_frame; delim=',')

#24 hours
for i=1:24
    fct_frame[:,:Wind]  .= wfct[:,4+i]./wind_norm_factor
    fct_frame[:,:Solar] .= sfct[:,4+i]./solr_norm_factor
    fct_frame[:,:Load]  .= lfct[:,4+i]./load_norm_factor
    # Save
    fct_fname = string("avail_ID_d2_c",lpad(i-1,2,"0"),".csv")
    CSV.write(string(fcast_dir, "/forecast/",fct_fname), fct_frame; delim=',')
end