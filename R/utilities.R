.onAttach <- function(libname, pkgname)
   packageStartupMessage("SDRcausal v 0.2-0")

.onUnload <- function(libpath)
    library.dynam.unload("SDRcausal",  libpath)
