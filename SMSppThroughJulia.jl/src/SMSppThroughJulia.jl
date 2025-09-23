module SMSppThroughJulia

using Libdl # This packet is used for handling dynamic libraries

"""
    This function assumes as input a block structured netcdf file in proper SMS++ format having notably the proper description of some Thermal unit
        See documentation :
            https://gitlab.com/smspp/smspp-project/-/tree/develop/doc/SMS++%20File%20Format%20Manual?ref_type=heads
    
        the string fname is the name of this file
    
    The second input is a (the negative of) price vector of appropriate dimension (the number of time steps)

    as output the user recovers the production output of the unit
"""
function value_nuclear_on_price( a::Vector{Float64}, fname::String)
    libpath = joinpath(pwd(), "cpp", "emx_smspp_lib")

    # Open Shared Library
    libc = Libdl.dlopen(libpath) # ; throw_error=false)
    if libc == C_NULL
        error("Impossible to load the library : $libpath")
    end
    
    try
        # Recover the pointer to the function
        emxThf = dlsym(libc, :emx_calling_thf_C)
        
        # vector a = - Market Price 
        #a = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
        n = length(a)
        result_size = Ref{Int32}(0)
        
        # The file with the data
        t_file = joinpath(pwd(), fname)

        # Calling the underlying C function
        values = ccall(emxThf, Ptr{Cdouble}, (Cstring, Ptr{Cdouble}, Csize_t, Ref{Cint}),
                       t_file, a, n, result_size)
        
        result_size[] <= 0 && error("Invalid size of result : $(result_size[])")
        
        # Création du tableau Julia à partir des données C (own=true means that Julia own the C ptr)
        result = unsafe_wrap(Array{Float64,1}, values, (result_size[],); own=true)
        
        return result
        
    catch e
        println("Some error occured when calling the C library : $(e)")
        throw(e)  # Répendre l'exception pour gérer l'erreur plus haut si nécessaire
    finally
        # close the lib
        Libdl.dlclose(libc)
    end
end


function test()

    rt = SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0] , "test/data/NBlock_ramp.nc4")

    rb = [855.0, 830.0, 805.0, 780.0, 755.0, 730.0, 705.0, 680.0, 655.0, 630.0, 605.0, 580.0, 555.0, 530.0, 505.0, 480.0, 455.0, 430.0, 405.0, 380.0, 355.0, 330.0, 305.0, 280.0]

    nbTest = 0
    if ( maximum(abs.(rt - rb)) < 1e-6 )
        nbTest += 1
    end

    rt = SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0] , "test/data/NBlock_ramp2.nc4")

    if ( maximum(abs.(rt - rb)) < 1e-6 )
        nbTest += 1
    end

    rt = SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0] , "test/data/NBlock_mindown.nc4")

    rb = [855.0, 830.0, 805.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if ( maximum(abs.(rt - rb)) < 1e-6 )
        nbTest += 1
    end

    rt = SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0] , "test/data/NBlock_mindown2.nc4")

    rb = [855.0, 830.0, 805.0, 800.0, 0.0, 0.0, 0.0, 0.0, 800.0, 830.0, 860.0, 890.0, 920.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0]

    if ( maximum(abs.(rt - rb)) < 1e-6 )
        nbTest += 1
    end

    rt = SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0] , "test/data/Nblock_mindown.nc4")

    rb = [855.0, 830.0, 805.0, 800.0, 800.0, 800.0, 830.0, 860.0, 890.0, 920.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0]

    if ( maximum(abs.(rt - rb)) < 1e-6 )
        nbTest += 1
    end

    rt = SMSppThroughJulia.value_nuclear_on_price( [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0] , "test/data/NBlock_mod.nc4")

    rb = [887.0, 886.0, 885.0, 884.0, 883.0, 882.0, 881.0, 856.0, 855.0, 854.0, 853.0, 828.0, 827.0, 826.0, 825.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if ( maximum(abs.(rt - rb)) < 1e-6 )
        nbTest += 1
    end

    rt = SMSppThroughJulia.value_nuclear_on_price( [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0] , "test/data/NBlock_mod2.nc4")

    rb = [905.0, 910.0, 909.0, 884.0, 883.0, 882.0, 881.0, 856.0, 855.0, 854.0, 853.0, 828.0, 827.0, 826.0, 825.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if ( maximum(abs.(rt - rb)) < 1e-6 )
        nbTest += 1
    end

    rt = SMSppThroughJulia.value_nuclear_on_price( [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0] , "test/data/Nblock_mod3.nc4")

    rb = [910.0, 911.0, 912.0, 913.0, 925.0, 925.0, 920.0, 910.0, 900.0, 890.0, 880.0, 855.0, 845.0, 835.0, 825.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if ( maximum(abs.(rt - rb)) < 1e-6 )
        nbTest += 1
    end

    println(string("Testing : ", nbTest, " / 8 - passed") )

end

end # module
