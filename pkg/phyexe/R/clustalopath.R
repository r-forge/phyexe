clustalo.path <- function(){
  libarch <- if (nzchar(.Platform$r_arch)) paste("libs", .Platform$r_arch, sep='') else "libs"
  if(.Platform$OS.type == "windows"){
    clustaloexe <- system.file(libarch,"clustalo.exe",package="phyexe")
  }
  else{
    clustaloexe <- system.file(libarch,"clustalo",package="phyexe")
  }
  clustaloexe
}
