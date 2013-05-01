muscle.path <- function(){
  libarch <- if (nzchar(.Platform$r_arch)) paste("libs", .Platform$r_arch, sep='') else "libs"
  if(.Platform$OS.type == "windows"){
    muscleexe <- system.file(libarch,"muscle.exe",package="phyexe")
  }
  else{
    muscleexe <- system.file(libarch,"muscle",package="phyexe")
  }
  muscleexe
}
