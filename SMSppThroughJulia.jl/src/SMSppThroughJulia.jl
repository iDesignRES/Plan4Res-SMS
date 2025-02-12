module SMSppThroughJulia

using Libdl # paquet utilisé pour la manipulation des librairies dynamiques

function test()
    libpath = joinpath(pwd(), "cpp", "emx_smspp_lib")

    # Ouvrir la bibliothèque partagée
    libc = Libdl.dlopen(libpath) # ; throw_error=false)
    if libc == C_NULL
        error("Impossible de charger la bibliothèque : $libpath")
    end
    
    try
        # Récupérer le pointeur de fonction
        emxThf = dlsym(libc, :emx_calling_thf_C)
        
        a = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
        n = length(a)
        result_size = Ref{Int32}(0)
        
        t_file = joinpath(pwd(), "NBlock_ramp.nc4")

        # Appel à la fonction C
        values = ccall(emxThf, Ptr{Cdouble}, (Cstring, Ptr{Cdouble}, Csize_t, Ref{Cint}),
                       t_file, a, n, result_size)
        
        result_size[] <= 0 && error("La taille du résultat est invalide : $(result_size[])")
        
        # Création du tableau Julia à partir des données C (own=true means that Julia own the C ptr)
        result = unsafe_wrap(Array{Float64,1}, values, (result_size[],); own=true)
        
        return result
        
    catch e
        println("Erreur lors de l'appel à la bibliothèque C : $(e)")
        throw(e)  # Répendre l'exception pour gérer l'erreur plus haut si nécessaire
    finally
        # Fermeture de la bibliothèque
        Libdl.dlclose(libc)
    end
end

end # module
