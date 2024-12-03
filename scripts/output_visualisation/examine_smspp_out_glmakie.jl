using CSV
using DataFrames
using Dates
using Colors
using Statistics

# For plots
using GLMakie
GLMakie.activate!()

using Tyler
using Tyler.TileProviders
using Tyler.MapTiles
using Tyler.Extents

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

# Make some power plant icons
# The following strings are the svg code paths
batsymbol_string = "M96.84 141.998c-4.947-23.457-20.359-32.211-25.862-13.887-11.822-22.963-37.961-16.135-22.041 6.289-3.005-1.295-5.872-2.682-8.538-4.191-8.646-5.318-15.259-11.314-19.774-17.586-3.237-5.07-4.994-10.541-4.994-16.229 0-19.774 21.115-36.758 50.861-43.694.446-.078.909-.154 1.372-.231-22.657 30.039 9.386 50.985 15.258 24.645l2.528-24.367 5.086 6.52H103.205l5.07-6.52 2.543 24.367c5.842 26.278 37.746 5.502 15.414-24.429 29.777 6.951 50.891 23.936 50.891 43.709 0 15.136-12.406 28.651-31.609 37.267 14.842-21.822-10.867-28.266-22.549-5.549-5.502-18.325-21.147-9.341-26.125 13.886z"
windsymb_string  = "M48.7,35.8L34,26c-0.2-1.7-1.1-3.2-2.3-4.2c2.6-3.4,3.7-6.8,3.2-9.9c-0.7-4.4-4.2-7.2-6.2-8.5c-0.6-0.4-1.5-0.4-2.1-0.1  c-0.7,0.4-1.1,1.1-1.1,1.8v15.8c-2.3,0.9-4,3.2-4,5.9c0,0.3,0,0.5,0,0.7c-5.4-0.3-9.5,0.8-12.1,3.4c-3.3,3.3-3.3,8-3,10.5  c0.1,0.8,0.6,1.4,1.3,1.7c0.3,0.1,0.5,0.2,0.8,0.2c0.5,0,0.9-0.2,1.3-0.5l13.7-11.2V52c0,0.6,0.4,1,1,1H31c0.6,0,1-0.4,1-1V33.2  c1.4,3.2,3.4,5.5,5.9,6.7c1.5,0.7,3,1,4.4,1c2.5,0,4.8-0.8,6.2-1.4c0.7-0.3,1.2-1,1.2-1.7C49.7,36.9,49.4,36.2,48.7,35.8z M27.5,5  c0.1,0,0.1,0,0.1,0c1.8,1.1,4.8,3.6,5.3,7.2c0.4,2.7-0.6,5.6-3,8.6c-0.7-0.3-1.4-0.4-2.2-0.4c-0.1,0-0.2,0-0.3,0L27.5,5z M27.7,22.4  c2.4,0,4.3,1.9,4.3,4.3S30.1,31,27.7,31s-4.3-1.9-4.3-4.3S25.4,22.4,27.7,22.4z M8.6,41.3c0,0-0.1,0-0.2-0.1  C8.2,39,8.2,35,10.8,32.3c1.9-1.9,4.9-2.9,8.9-2.9c0.8,0,1.5,0,2.4,0.1c0.1,0.2,0.2,0.3,0.3,0.5L8.6,41.3z M25.5,51V32.6  c0.7,0.3,1.4,0.4,2.2,0.4c0.8,0,1.5-0.1,2.2-0.4V51H25.5z M47.6,37.6c-1.9,0.8-5.7,2-9,0.5c-2.5-1.2-4.5-3.8-5.7-7.8  c0.4-0.6,0.7-1.2,0.9-2l13.8,9.1C47.6,37.4,47.7,37.5,47.6,37.6z"
solarsymb_string = "M4.86261 5.07032L5.59077 4.5L18.409 4.5L19.1371 5.07032L21.728 15.5703L20.9999 16.5L12.7499 16.5L12.7499 18.75L15.7499 18.75L15.7499 20.25L8.24991 20.25L8.24991 18.75L11.2499 18.75L11.2499 16.5L2.99986 16.5L2.2717 15.5703L4.86261 5.07032ZM6.1782 6L5.25288 9.75L8.95789 9.75L9.88321 6L6.1782 6ZM3.95742 15L4.88275 11.25L8.58776 11.25L7.66243 15L3.95742 15ZM9.20742 15L10.1327 11.25L14.2079 11.25L15.1332 15L9.20742 15ZM12.9124 6L13.8378 9.75L10.5029 9.75L11.4282 6L12.9124 6ZM14.4574 6L15.3827 9.75L18.7469 9.75L17.8215 6L14.4574 6ZM16.6782 15L15.7529 11.25L19.117 11.25L20.0423 15L16.6782 15Z"
hydrosymb_string1 = "M10.323,414.019c7.78,0,11.155,1.701,15.823,4.056c5.537,2.792,12.428,6.267,25.119,6.267
				c12.691,0,19.582-3.475,25.119-6.267c4.669-2.354,8.043-4.056,15.825-4.056c7.781,0,11.156,1.701,15.826,4.056
				c5.537,2.792,12.428,6.267,25.12,6.267c12.692,0,19.583-3.475,25.12-6.266c4.67-2.355,8.044-4.057,15.827-4.057
				c7.781,0,11.155,1.701,15.825,4.056c5.537,2.792,12.429,6.267,25.12,6.267s19.581-3.475,25.118-6.267
				c4.669-2.354,8.043-4.056,15.822-4.056c7.781,0,11.155,1.701,15.826,4.056c5.537,2.792,12.428,6.267,25.12,6.267
				c12.693,0,19.584-3.475,25.121-6.267c4.67-2.354,8.044-4.056,15.826-4.056c7.782,0,11.157,1.701,15.827,4.057
				c5.537,2.791,12.429,6.266,25.121,6.266s19.584-3.475,25.121-6.266c4.67-2.355,8.044-4.057,15.827-4.057
				c7.783,0,11.159,1.701,15.829,4.057c5.538,2.791,12.43,6.266,25.123,6.266c12.694,0,19.586-3.475,25.124-6.266
				c4.67-2.355,8.045-4.057,15.829-4.057c5.7,0,10.323-4.622,10.323-10.323s-4.622-10.323-10.323-10.323
				c-0.5,0-0.974,0.014-1.456,0.025V55.045c0-5.7-4.622-10.323-10.323-10.323h-41.29c-5.7,0-10.323,4.622-10.323,10.323v22.71
				H76.996v-22.71c0-5.7-4.621-10.323-10.323-10.323h-41.29c-5.701,0-10.323,4.622-10.323,10.323v338.504
				c-1.486-0.113-3.059-0.177-4.738-0.177C4.621,393.372,0,397.995,0,403.695S4.621,414.019,10.323,414.019z M458.932,65.367h20.645
				v332.785c-1.062,0.502-2.065,1.005-3.022,1.487c-4.67,2.355-8.047,4.057-15.83,4.057c-0.629,0-1.219-0.016-1.793-0.037V65.367z
				 M185.503,362.024c2.358,4.352,3.915,7.227,3.915,14.138c0,6.911-1.558,9.785-3.914,14.137c-0.642,1.185-1.317,2.45-1.981,3.814
				c-2.752-0.457-5.851-0.74-9.422-0.74c-9.525,0-15.782,1.957-20.653,4.114c2.495-4.866,5.003-11.121,5.003-21.325
				c0.001-12.143-3.552-18.701-6.407-23.972c-2.358-4.352-3.914-7.227-3.914-14.137c0-6.909,1.558-9.784,3.914-14.134
				c2.855-5.269,6.408-11.827,6.408-23.97s-3.553-18.7-6.408-23.97c-2.358-4.351-3.914-7.226-3.914-14.135
				c0-6.909,1.558-9.785,3.914-14.135c1.916-3.536,4.146-7.652,5.406-13.646h30.006c-0.557,1.239-1.215,2.453-1.951,3.811
				c-2.855,5.269-6.408,11.827-6.408,23.97c0,12.143,3.553,18.7,6.408,23.97c2.358,4.351,3.914,7.226,3.914,14.135
				c0,6.909-1.558,9.785-3.914,14.135c-2.855,5.269-6.408,11.827-6.408,23.969C179.096,350.195,182.649,356.753,185.503,362.024z
				 M148.129,213.417c-5.701,0-10.323,4.622-10.323,10.323c0,6.909-1.558,9.784-3.914,14.134c-0.326,0.602-0.662,1.224-1,1.866
				h-7.472v-41.29h92.903v41.29H207.43c1.501-4.068,2.634-9.143,2.634-16c0-5.701-4.621-10.323-10.323-10.323H148.129z
				 M305.548,223.739c0,6.909-1.558,9.784-3.914,14.134c-0.326,0.602-0.662,1.224-1,1.866h-7.473v-41.29h92.903v41.29h-10.892
				c1.501-4.068,2.634-9.143,2.634-16c0-5.701-4.622-10.323-10.323-10.323h-51.613C310.171,213.417,305.548,218.039,305.548,223.739
				z M353.245,362.024c2.358,4.352,3.915,7.227,3.915,14.138c0,6.911-1.558,9.785-3.915,14.137c-0.758,1.399-1.564,2.894-2.34,4.557
				c-3.57-0.888-7.802-1.482-13.028-1.482c-6.635,0-11.676,0.954-15.776,2.259c2.16-4.603,4.092-10.516,4.092-19.47
				c0.001-12.145-3.551-18.702-6.406-23.973c-2.358-4.352-3.915-7.227-3.915-14.137c0-6.909,1.558-9.784,3.914-14.134
				c2.855-5.269,6.408-11.827,6.408-23.97s-3.553-18.7-6.408-23.97c-2.357-4.351-3.914-7.226-3.914-14.135
				c0-6.909,1.558-9.785,3.915-14.135c1.915-3.536,4.146-7.652,5.405-13.646h30.006c-0.556,1.239-1.215,2.453-1.951,3.811
				c-2.855,5.269-6.408,11.827-6.408,23.97c0,12.143,3.553,18.7,6.408,23.97c2.357,4.351,3.914,7.226,3.914,14.135
				c0,6.909-1.558,9.785-3.915,14.135c-2.854,5.269-6.407,11.827-6.407,23.969C346.838,350.195,350.39,356.753,353.245,362.024z
				 M76.997,141.159h180.645c5.7,0,10.323-4.622,10.323-10.323c0-5.7-4.622-10.323-10.323-10.323H76.997V98.399h361.289v22.114
				h-82.501c-5.7,0-10.323,4.622-10.323,10.323c0,5.7,4.622,10.323,10.323,10.323h82.501v255.441
				c-4.554-1.777-10.323-3.226-18.514-3.226c-12.692,0-19.584,3.475-25.121,6.266c-4.67,2.355-8.044,4.057-15.827,4.057
				c-3.717,0-6.425-0.391-8.753-1.054c0.406-0.806,0.849-1.627,1.325-2.508c2.855-5.27,6.407-11.828,6.407-23.972
				c0-12.144-3.552-18.702-6.407-23.973c-2.358-4.352-3.915-7.227-3.915-14.137c0-6.909,1.558-9.784,3.914-14.134
				c2.855-5.269,6.408-11.827,6.408-23.97s-3.553-18.7-6.408-23.97c-2.357-4.351-3.914-7.226-3.914-14.135
				c0-0.51,0.011-0.991,0.028-1.459h28.876c5.7,0,10.323-4.622,10.323-10.323v-61.935c0-5.7-4.622-10.323-10.323-10.323H282.839
				c-5.7,0-10.323,4.622-10.323,10.323v61.935c0,5.7,4.622,10.323,10.323,10.323h12.414c-0.011,0.482-0.027,0.958-0.027,1.459
				c0,12.143,3.553,18.7,6.408,23.97c2.357,4.351,3.914,7.226,3.914,14.135c0,6.909-1.558,9.785-3.915,14.135
				c-2.854,5.269-6.407,11.827-6.407,23.969c0,12.143,3.552,18.701,6.407,23.972c2.358,4.352,3.915,7.227,3.915,14.138
				c0,6.911-1.558,9.785-3.915,14.137c-1.887,3.483-4.072,7.54-5.341,13.388c-7.322-0.093-10.648-1.759-15.186-4.047
				c-5.537-2.791-12.428-6.266-25.119-6.266c-12.691,0-19.581,3.475-25.118,6.267c-4.669,2.354-8.043,4.056-15.822,4.056
				c-5.607,0-8.925-0.887-12.089-2.261c0.225-0.423,0.457-0.854,0.7-1.302c2.855-5.27,6.408-11.828,6.408-23.972
				c0-12.145-3.553-18.702-6.408-23.973c-2.358-4.352-3.914-7.227-3.914-14.137c0-6.909,1.558-9.784,3.914-14.134
				c2.855-5.269,6.408-11.827,6.408-23.97s-3.553-18.7-6.408-23.97c-2.358-4.351-3.914-7.226-3.914-14.135
				c0-0.51,0.011-0.991,0.028-1.459h28.873c5.701,0,10.323-4.622,10.323-10.323v-61.935c0-5.7-4.621-10.323-10.323-10.323H115.097
				c-5.701,0-10.323,4.622-10.323,10.323v61.935c0,5.7,4.621,10.323,10.323,10.323h12.414c-0.011,0.482-0.027,0.958-0.027,1.459
				c0,12.143,3.553,18.7,6.408,23.97c2.358,4.351,3.914,7.226,3.914,14.135s-1.558,9.785-3.914,14.135
				c-2.855,5.269-6.408,11.827-6.408,23.969c0,12.143,3.553,18.701,6.407,23.972c2.358,4.352,3.915,7.227,3.915,14.138
				c0,6.911-1.558,9.785-3.914,14.137c-1.861,3.436-4.014,7.432-5.29,13.156c-4.769-0.545-7.63-1.979-11.272-3.814
				c-5.537-2.792-12.428-6.267-25.12-6.267c-6.33,0-11.216,0.865-15.212,2.078V141.159z M35.706,65.367h20.645v338.027
				c-1.48,0.195-3.146,0.301-5.086,0.301c-7.631,0-11.027-1.639-15.559-3.923V65.367z"
nuclearsymb_string = "M44.1969,66c-5-6.1714-8-26.7429-5-36H16c3,9.2571,0,29.8286-5,36Z"
thfsymb_string = "M413.426,314.527l-5.764-266.46c1.726-0.629,2.965-2.27,2.965-4.214
		c0-2.011-1.328-3.694-3.149-4.271l-0.217-10.014c-0.053-2.446-2.052-4.402-4.499-4.402h-22.41c-2.447,0-4.446,1.956-4.499,4.402
		l-0.217,10.014c-1.82,0.576-3.149,2.259-3.149,4.271c0,1.944,1.239,3.584,2.965,4.214l-5.764,266.46h-24.117v19.15h-49.518
		c-8.948-41.729-13.6-84.278-13.907-127.021h1.962c2.485,0,4.5-2.015,4.5-4.5s-2.015-4.5-4.5-4.5h-2.002h-109.66h-1.927
		c-2.485,0-4.5,2.015-4.5,4.5s2.015,4.5,4.5,4.5h1.887c-0.201,27.741-2.241,55.696-6.116,83.171
		c-1.125,7.967-2.419,15.953-3.856,23.904c-7.216-38.43-10.989-77.456-11.257-116.646h2.019c2.485,0,4.5-2.015,4.5-4.5
		s-2.015-4.5-4.5-4.5H34.759c-2.485,0-4.5,2.015-4.5,4.5s2.015,4.5,4.5,4.5h2.019c-0.411,59.977-9.006,119.58-25.611,177.231
		L0,413.087h438.253v-98.56H413.426z M383.598,87.736h15.92l0.657,30.383h-17.234L383.598,87.736z M382.746,127.119h17.624
		l0.657,30.382h-18.938L382.746,127.119z M399.323,78.736h-15.53l0.657-30.383h14.216L399.323,78.736z M381.895,166.501h19.327
		l0.657,30.383h-20.642L381.895,166.501z M381.043,205.884h21.031l0.657,30.383h-22.346L381.043,205.884z M380.191,245.266h22.735
		l0.657,30.382h-24.05L380.191,245.266z M398.359,34.167l0.112,5.187h-13.826l0.112-5.187H398.359z M379.339,284.648h24.439
		l0.646,29.879h-25.732L379.339,284.648z M273.591,226.584h-92.624c0.06-1.543,0.101-3.085,0.149-4.628h92.323
		C273.487,223.499,273.532,225.042,273.591,226.584z M249.698,212.956v-6.299h8.919v6.299H249.698z M240.698,212.956h-8.919v-6.299
		h8.919V212.956z M222.779,212.956h-8.919v-6.299h8.919V212.956z M204.86,212.956h-8.919v-6.299h8.919V212.956z M273.223,212.956
		h-5.607v-6.299h5.508C273.14,208.758,273.188,210.856,273.223,212.956z M186.941,206.657v6.299h-5.606
		c0.036-2.1,0.083-4.201,0.098-6.299H186.941z M142.664,218.782h-97.37c0.075-1.896,0.126-3.794,0.185-5.691h97.001
		C142.538,214.988,142.589,216.886,142.664,218.782z M126.979,204.091h-9.75v-7.003h9.75V204.091z M135.979,197.088h6.176
		c0.016,2.335,0.067,4.669,0.107,7.003h-6.283V197.088z M79.729,204.091v-7.003h9.75v7.003H79.729z M98.479,197.087h9.75v7.003
		h-9.75V197.087z M70.729,197.087v7.004h-9.75v-7.004H70.729z M51.979,197.087v7.004h-6.284c0.04-2.334,0.091-4.668,0.107-7.004
		H51.979z M11.958,404.087l7.857-27.278c14.002-48.613,22.382-98.595,25.051-149.026h98.225
		c2.669,50.431,11.049,100.413,25.051,149.026L176,404.087H11.958z M185.366,404.087l-8.575-29.769
		c-3.613-12.545-6.825-25.187-9.68-37.901c3.224-14.922,5.952-30.165,8.093-45.332c2.589-18.356,4.363-36.924,5.343-55.502h93.461
		c2.539,48.278,10.552,96.138,23.956,142.688l7.436,25.814H185.366z M429.253,404.087H314.766l-8.153-28.305
		c-3.157-10.965-6.003-22.005-8.554-33.104h56.515v-19.15h74.68V404.087z" 

batsymbol = BezierPath(batsymbol_string, fit = true, flipy = true)
windsymb  = BezierPath(windsymb_string, fit = true, flipy = true)
solarsymb = BezierPath(solarsymb_string, fit = true, flipy = true)
hydrosymb1 = BezierPath(hydrosymb_string1, fit = true, flipy = true)
nucsymb    = BezierPath(nuclearsymb_string, fit = true, flipy = true)
thfsymb    = BezierPath(thfsymb_string, fit = true, flipy = true)

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

incon_data = CSV.read(string(github_local_d, github_smsppin,"/nutsx/IN_Interconnections.csv"), DataFrame; delim=';')

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
r_dir       = string(github_local_d, github_smsppout, "/nutsx/", res_name)
outfile_ext = "OUT"

with_acopf = false
ts_prod = CSV.read(string(r_dir, "/ActivePower/ActivePower", outfile_ext, ".csv" ), DataFrame; delim=',')
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

long = bus_data.y
latt = bus_data.x

function to_web_mercator(lat,lo)
    return Point2f(MapTiles.project((lat,lo), MapTiles.wgs84, MapTiles.web_mercator))
end

# Plots with export and import buses : colours : pbalNames
ndCols = [(Im66, pbalNames[1]), (Im33, pbalNames[2]), (Inaught, pbalNames[3]), (Ip33, pbalNames[4]), (Ip66, pbalNames[5])]

# Start plotting 

# Création de la figure
fig = Figure(; size = (1200, 1000))
#display(fig, px_per_unit = 2)

# Ajout de l'arrière-plan de la carte avec Tyler
ax = Axis(fig[1, 1], title = "Map of the Spanish Case Study")

prov = Esri(:WorldImagery)
extent = Rect2f(-9, 35, 13, 10)   #

# Charger les tuiles de la carte
tm = Tyler.Map(extent; provider=prov, figure=fig, axis=ax);

# Change the units of the axis:
#tm.axis.xticks = (0:10:180, string)
#tm.axis.yticks = (0:10:90, string)

# Ajout des nœuds
for ndc in ndCols
    (Ic, ndcolour) = ndc

    pts = to_web_mercator.(latt[Ic],long[Ic])
    objscatter = scatter!(tm.axis, pts; color = ndcolour, markersize = 5, marker = :circle)
    # The translate functions moves the cloud higher up in the z-axis
    translate!(objscatter, 0, 0, 2)
end

# Defaillance
if ( length(Idef) > 0)
    pts = to_web_mercator.(latt[Idef], long[Idef])

    objscatter = scatter!(tm.axis, pts; color = :red, markersize = 5, marker = :cross)
    # The translate functions moves the cloud higher up in the z-axis
    translate!(objscatter, 0, 0, 2)
end

# The power lines
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
        pt1 = to_web_mercator.(xLines[i, 1], yLines[i, 1])
        pt2 = to_web_mercator.(xLines[i, 2], yLines[i, 2])

        objL = lines!(tm.axis, [pt1, pt2], color=cNames[cSc[i]], linewidth=2)
        translate!(objL, 0, 0, 2)
    end
end

Tech_to_Symb = Dict([("Solar",solarsymb),("Wind",windsymb),("Hydro",hydrosymb1),
                     ("Coal",thfsymb),("Gas",thfsymb),("Biomass",thfsymb),("Oil/Diesel",thfsymb),("Diesel",thfsymb),("Oil/Gas/Diesel",thfsymb),("Oil",thfsymb),("Oil/Gas",thfsymb),
                     ("Nuclear",nucsymb) ])

#Plot all the generators on the map
p_scatters = []
for i=1:length(gen_data.unit_id)
    pt_plant = to_web_mercator(gen_data.x[i], gen_data.y[i] )
    o_s = scatter!(tm.axis, pt_plant, marker = Tech_to_Symb[gen_data.primary_fuel[i]], markersize = 10, color = :black)
    translate!(o_s, 0, 0, 1)

    append!(p_scatters, [o_s])

    t_s = text!(tm.axis, pt_plant, text = gen_data.unit_id[i]; fontsize = 10, color = :darkblue, align =(:center,:top))
    translate!(t_s, 0, 0, 1)
end

# Fonction pour mettre à jour la taille des marqueurs en fonction du zoom
function update_markersize!(zf)
    # Calcul de la nouvelle taille des marqueurs en fonction de l'échelle de l'axe
    for ob_s in p_scatters
        ob_s.markersize[] = 10*max(zf - 7,1)
    end
end

#Intercept the zoom events
# Création d'un observable pour suivre le zoom
zoom_observable = Observable(ax.scene.camera)
#on(ax.scene.events.mouseposition) do _
on(ax.finallimits) do _
    sf = tm.zoom.val
    update_markersize!(sf)
end

fig_name = string(r_dir, "/", "spain_", res_type, "_fig.png")
img_name = string(r_dir, "/", "spain_", res_type, "_tm.png")

#save( fig_name, fig; px_per_unit=2 )
save( img_name, tm )
# Affichage de la figure
fig

