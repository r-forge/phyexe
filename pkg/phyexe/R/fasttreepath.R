fasttree.path <- function(){
  libarch <- if (nzchar(.Platform$r_arch)) paste("libs", .Platform$r_arch, sep='') else "libs"
  if(.Platform$OS.type == "windows"){
    fasttreeexe <- system.file(libarch,"FastTree.exe",package="phyexe")
  }
  else{
    fasttreeexe <- system.file(libarch,"FastTree",package="phyexe")
  }
  fasttreeexe
}
