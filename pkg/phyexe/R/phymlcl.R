phyml.cl <- function(input=NULL,datatype="nt",pars=TRUE,bootstrap=0,model="GTR",frequencies="m",pinv=NULL,nclasses=NULL,alpha=NULL,searchstrategy="BEST",usertree=NULL,params="tlr",seed=NULL,quiet=TRUE){
  libarch <- if (nzchar(.Platform$r_arch)) paste("libs", .Platform$r_arch, sep='') else "libs"
  if(.Platform$OS.type == "windows"){
    phymlexe <- system.file(libarch,"phyml.exe",package="phyexe")
  }
  else{
    phymlexe <- system.file(libarch,"phyml",package="phyexe")
  }
  phymlcmd <- paste(phymlexe,"-i",input,"-d",datatype,"-b",bootstrap,"-m",model,"-f",frequencies,"-s",searchstrategy,"-o",params,sep=" ")
  if(pars){
    phymlcmd <- paste(phymlcmd,"-p",sep=" ")
  }
  if(!is.null(pinv)){
    phymlcmd <- paste(phymlcmd,"-v",pinv,sep=" ")
  }
  if(!is.null(nclasses)){
    phymlcmd <- paste(phymlcmd,"-c",nclasses,sep=" ")
  }
  if(!is.null(alpha)){
    phymlcmd <- paste(phymlcmd,"-a",alpha,sep=" ")
  }
  if(!is.null(usertree)){
    phymlcmd <- paste(phymlcmd,"-u",usertree,sep=" ")
  }
  if(!is.null(seed)){
    phymlcmd <- paste(phymlcmd,"--r_seed",seed,sep=" ")
  }
  if(quiet){
    phymlcmd <- paste(phymlcmd,"--quiet")
  }
  phymlcmd
}
