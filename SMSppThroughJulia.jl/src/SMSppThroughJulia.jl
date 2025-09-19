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
        a = [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
        n = length(a)
        result_size = Ref{Int32}(0)
        
        # The file with the data
        t_file = joinpath(pwd(), "NBlock_mod3.nc4")

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

end # module
