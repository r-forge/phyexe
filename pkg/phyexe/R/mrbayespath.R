mrbayes.path <- function(){
  libarch <- if (nzchar(.Platform$r_arch)) paste("libs", .Platform$r_arch, sep='') else "libs"
  if(.Platform$OS.type == "windows"){
    mrbayesexe <- system.file(libarch,"muscle.exe",package="phyexe")
  }
  else{
    mrbayesexe <- system.file(libarch,"muscle",package="phyexe")
  }
  mrbayesexe
}
