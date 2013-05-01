phyml.path <- function(){
  libarch <- if (nzchar(.Platform$r_arch)) paste("libs", .Platform$r_arch, sep='') else "libs"
  if(.Platform$OS.type == "windows"){
    phymlexe <- system.file(libarch,"phyml.exe",package="phyexe")
  }
  else{
    phymlexe <- system.file(libarch,"phyml",package="phyexe")
  }
  phymlexe
}
