#### Indigo Ag Maize Microbiome Project
#### Corey Schultz - UGA - Wallace Lab
#
### Script to parse apart fungal and bacteria sequences

library(dplyr)

getwd()

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")

list.dirs()

fastqs <- list.files("Indigo_RawData/g2f_indigo_combined")

crude <- read.csv("MetaData_Stuff/Crude_16s_vs_ITS.csv")

# assign samples bacteria or fungus

crude['assign'] <- NA

for(r in 1:nrow(crude)){
  bact <- crude[r,"crude_16s"]
  fung <- crude[r,"crude_its"]
  
  if(bact > fung){
    crude[r,"assign"] = "bact"
  } else {
    crude[r,"assign"] = "fung"
  }
}

# divide fastqs based on that assignment

bacts <- filter(crude, assign == "bact")
fungs <- filter(crude, assign == "fung")

# sort bacteria 
bfiles <- list()
for(r in 1:nrow(bacts)){
  tstring <- bacts[r,1]
  bfiles <- append(bfiles,fastqs[grep(tstring,fastqs)])
}

move_bact <- function(x){
  file.rename(from = file.path("\\\\wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome/Indigo_RawData/g2f_indigo_combined", x) ,
               to = file.path("\\\\wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome/Indigo_RawData/Bacteria_Sequences", x))
}

# sort fungus 
ffiles <- list()
for(r in 1:nrow(fungs)){
  tstring <- fungs[r,1]
  ffiles <- append(ffiles,fastqs[grep(tstring,fastqs)])
}

move_fung <- function(x){
  file.rename(from = file.path("\\\\wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome/Indigo_RawData/g2f_indigo_combined", x) ,
              to = file.path("\\\\wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome/Indigo_RawData/Fungal_Sequences", x))
}

lapply(ffiles, move_fung)



