dropbox.path <- paste(Sys.getenv('USERPROFILE'), '\\Dropbox', sep='')



.First <- function(){

}

# Run everything in lib/ when the project is opened

for (file in list.files("lib", pattern = "[.][Rr]$")) {
  source(paste("lib/", file, sep=''))
}
rm(file)
