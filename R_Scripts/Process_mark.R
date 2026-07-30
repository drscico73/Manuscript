#This script is used. It is correct. I remove all the rows that have atleat one NA

#Importing Libraries
library(dplyr)
library(stringr)
library(vcfR)

#Reading the Merged Markers vcf file
a<- read.vcfR("~/Desktop/Projects/HelOreSam/markers.vcf.gz")
sample_calc<- cbind(a@fix,a@gt)
sample_calc<- as.data.frame(sample_calc)
sample_calc<- distinct(sample_calc)
af_index<- 8

# Extracting only the AF information
af_values <- sample_calc
af_values[, 10:ncol(sample_calc)] <- apply(sample_calc[, 10:ncol(sample_calc)], c(1, 2), function(x) {
  fields <- strsplit(x, ":")[[1]]
  if (length(fields) >= af_index) {
    return(fields[af_index]) 
  } else {
    return(NA)  
  }
})
af_values<- af_values[-c(3,6,7,8,9)]
af_values$`Dmel_OregonR_non-inbred`<- as.numeric(af_values$`Dmel_OregonR_non-inbred`)
af_values$Dmel_Helsinki_inbred_F11<- as.numeric(af_values$Dmel_Helsinki_inbred_F11)
af_values$Dmel_Samarkand_inbred<- as.numeric(af_values$Dmel_Samarkand_inbred)
af_values<- af_values[rowSums(is.na(af_values)) == 0,]
af_values <- af_values[nchar(af_values$REF) == 1 & nchar(af_values$ALT) == 1, ]
af_values$POS<- as.numeric(af_values$POS)

#Generating unique SNP files

  #Alternate allele in the focal line
  sam_1<- af_values[af_values$Dmel_Samarkand_inbred==1&af_values$Dmel_Helsinki_inbred_F11==0&
                      af_values$`Dmel_OregonR_non-inbred`==0,]
  ore_1<- af_values[af_values$`Dmel_OregonR_non-inbred`==1&af_values$Dmel_Helsinki_inbred_F11==0&
                      af_values$Dmel_Samarkand_inbred==0,]
  hel_1<- af_values[af_values$Dmel_Helsinki_inbred_F11==1&af_values$`Dmel_OregonR_non-inbred`==0&
                      af_values$Dmel_Samarkand_inbred==0,]

  #Reference allele in the focal line
  sam_0<- af_values[af_values$Dmel_Samarkand_inbred==0&af_values$Dmel_Helsinki_inbred_F11==1&
                      af_values$`Dmel_OregonR_non-inbred`==1,]
  ore_0<- af_values[af_values$`Dmel_OregonR_non-inbred`==0&af_values$Dmel_Helsinki_inbred_F11==1&
                      af_values$Dmel_Samarkand_inbred==1,]
  hel_0<- af_values[af_values$Dmel_Helsinki_inbred_F11==0&af_values$`Dmel_OregonR_non-inbred`==1&
                      af_values$Dmel_Samarkand_inbred==1,]

  #Merging both
  hel<- rbind(hel_0, hel_1)
  hel$line<- paste("hel")
  sam<- rbind(sam_0, sam_1)
  sam$line<- paste("sam")
  ore<- rbind(ore_0, ore_1)
  ore$line<- paste("ore")
  ore_pos<- ore[c(1,2,2)]
  hel_pos<- hel[c(1,2,2)]
  sam_pos<- sam[c(1,2,2)]
   write.table(ore_pos, file = "~/Desktop/ore_pos_file.txt", sep="\t", 
           row.names = FALSE, col.names = FALSE, quote=FALSE)
   write.table(hel_pos, file = "~/Desktop/hel_pos_file.txt", sep="\t", 
               row.names = FALSE, col.names = FALSE, quote=FALSE)
   write.table(sam_pos, file = "~/Desktop/sam_pos_file.txt", sep="\t", 
               row.names = FALSE, col.names = FALSE, quote=FALSE)

   
   ore_marker<- ore[c(1,2,3,4)]
   hel_marker<- hel[c(1,2,3,4)]
   sam_marker<- sam[c(1,2,3,4)]
   write.table(ore_marker, file = "~/Desktop/ore_marker_file.txt", sep="\t", 
               row.names = FALSE, col.names = FALSE, quote=FALSE)
   write.table(hel_marker, file = "~/Desktop/hel_marker_file.txt", sep="\t", 
               row.names = FALSE, col.names = FALSE, quote=FALSE)
   write.table(sam_marker, file = "~/Desktop/sam_marker_file.txt", sep="\t", 
               row.names = FALSE, col.names = FALSE, quote=FALSE)
   
  marker<- rbind(sam,ore,hel)
  marker <- marker[order(marker$CHROM, marker$POS, decreasing = FALSE), ]
  
  pos<- marker[c(1,2,2)]
  # write.table(pos, file = "~/Desktop/marker_pos_file.txt", sep="\t", 
  #             row.names = FALSE, col.names = FALSE, quote=FALSE)
  mark<- marker[c(1,2,8)]
  mark$pos<- paste(mark$CHROM, mark$POS, sep="_")
mark<- mark[-c(1,2)]  
  write.table(mark, file = "~/Desktop/Projects/HelOreSam/Window_analysis/mark.txt", sep="\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)
