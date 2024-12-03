using CSV
using DataFrames
using Dates

replaceparts(d::DateTime; year=Dates.year(d), month=Dates.month(d), day=Dates.day(d), hour=Dates.hour(d), min=Dates.minute(d)) = typeof(d)(year, month, day, hour, min)

"""
    fname   : the name of the file

    b_date  : the beginning Date
    e_date  : the end Date (excluded)
    s_idx   : the scenario idx

"""
function extract_ts_sequence(fname, b_date, e_date, s_idx)
    
    ts_data = CSV.read(fname, DataFrame; delim=',')

    dd = try
        DateTime.(ts_data[!,1], DateFormat("dd/mm/yyyy HH:MM"))
    catch
        DateTime.(ts_data[!,1], DateFormat("dd/mm/yyyy HH:MM:SS"))
    end
    if ( length(dd) != size(ts_data)[1] )
        error(string(" The length [",string(length(dd)),"] of the Dates in the times series has been converted longer than the info [",string(size(ts_data)[1]), "] - Error with file ", fname))
    end
    
    # Time step in hours
    deltat = Dates.Hour(dd[2]-dd[1]).value
    #println(string("deltat=",deltat))
    if ( deltat < 1 )
        error("Unsupported less than 1 hour tstep for now")
    end
    i0 = findall( dd .== b_date )
    i1 = findall( dd .== e_date )

    off_set = 0
    if ( isempty(i0) || isempty(i1) || (i1[1] < i0[1]) )
        # in this case it is more or less assumed that the TS has a generic year
        yr = Dates.Year( dd[1] ).value

        g_b_date = replaceparts( b_date, year=yr )
        g_e_date = replaceparts( e_date, year=yr )

        i0 = findall( dd .<= g_b_date )
        i1 = findall( dd .>= g_e_date )
        if ( isempty(i0) || isempty(i1) || (i1[1] < i0[end]) )
            error("Not a time series object - critical error")
        end
        off_set = convert(Dates.Hour,(dd[i0[end]] - g_b_date)).value 
        #println(string(" the indexes ", string(i0), " and ", string(i1) ) )
    else
        off_set = convert(Dates.Hour,(dd[i0[end]] - b_date)).value
    end
    nbH = convert(Dates.Hour,(e_date - b_date)).value
    tsv = Vector{Float64}(undef, nbH)
    if ( deltat == 1 )
        tsv = ts_data[collect(i0[end]:i1[1]-1), s_idx + 2]
    else
        i_tsv=1
        #println(string("the offset is : ", string(off_set)))
        for el in ts_data[collect(i0[end]:i1[1]-1), s_idx + 2]
            tsv[ collect(i_tsv:min(i_tsv + deltat + off_set - 1,nbH) ) ] .= el
            i_tsv += deltat
        end
    end

    return tsv
end

