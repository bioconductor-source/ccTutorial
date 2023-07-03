
.onLoad <- function(libname, pkgname) {
    
  if(.Platform$OS.type == "windows" && interactive()
     && .Platform$GUI ==  "Rgui"){
    addVigs2WinMenu("ccTutorial")
    addVigs2WinMenu("ccTutorialSupplement")
  }

}
.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.19")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}

