raxml.path <- function(){
  libarch <- if (nzchar(.Platform$r_arch)) paste("libs", .Platform$r_arch, sep='') else "libs"
  if(.Platform$OS.type == "windows"){
    raxmlexe <- system.file(libarch,"raxmlHPC-PTHREADS-SSE3.exe",package="phyexe")
  }
  else{
    raxmlexe <- system.file(libarch,"raxmlHPC-PTHREADS-SSE3",package="phyexe")
  }
  raxmlexe
}
