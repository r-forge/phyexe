clustalo.cl <- function(input=NULL, output=NULL, MoreArgs="--outfmt=fasta"){
  libarch <- if (nzchar(.Platform$r_arch)) paste("libs", .Platform$r_arch, sep='') else "libs"
  if(.Platform$OS.type == "windows"){
    clustaloexe <- system.file(libarch,"clustalo.exe",package="phyexe")
  }
  else{
    clustaloexe <- system.file(libarch,"clustalo",package="phyexe")
  }
  if(is.null(input)){
    stop("Input file needed")
  }
  clustalocmd <- paste(clustaloexe,"-i",input,sep=" ")
  if(!is.null(output)){
    clustalocmd <- paste(clustalocmd,"-o",output,sep=" ")
  }
  clustalocmd <- paste(clustalocmd,MoreArgs,sep=" ")
  clustalocmd
}
