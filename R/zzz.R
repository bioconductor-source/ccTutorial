
.onLoad <- function(libname, pkgname) {
    
  if(.Platform$OS.type == "windows" && interactive()
     && .Platform$GUI ==  "Rgui"){
    addVigs2WinMenu("ccTutorial")
    addVigs2WinMenu("ccTutorialSupplement")
  }

}
