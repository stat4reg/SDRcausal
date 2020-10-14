.onAttach <- function(libname, pkgname)
   packageStartupMessage("SDRcausal v 0.3.0")

.onUnload <- function(libpath)
    library.dynam.unload("SDRcausal",  libpath)
