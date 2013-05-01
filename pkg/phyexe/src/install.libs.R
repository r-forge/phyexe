libarch <- if (nzchar(R_ARCH)) paste("libs", R_ARCH, sep='') else "libs"
dest <- file.path(R_PACKAGE_DIR, libarch)
files <- if (WINDOWS) c("phyml/src/phyml.exe","phytime/src/phytime.exe","raxml/raxmlHPC-PTHREADS-SSE3.exe","mrbayes/mb.exe","fasttree/FastTree.exe","lap/lap.exe","clustalo/src/clustalo.exe","muscle/muscle.exe") else c("phyml/src/phyml","phytime/src/phytime","raxml/raxmlHPC-PTHREADS-SSE3","mrbayes/mb","fasttree/FastTree","lap/lap","clustalo/src/clustalo","muscle/muscle")
# Only copy if built
files <- files[file.exists(files)]
if (length(files)) {
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(files, dest, overwrite = TRUE)
}
