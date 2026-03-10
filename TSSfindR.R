#################################################################
##                                                             ##
##  ///////  //////   //////      ////              /   ////   ##
##    //     \\       \\         //                 /  \\   \  ##
##    \\     //       //         ////  //  /        /  //////  ##
##    //     \\\\\    \\\\\\     //        ////  ////  \\\     ##
##    \\         \\        \\    //    //  /  /  /  /  // \    ##
##    //         //        //    //    //  /  /  /  /  \\  \   ##
##    \\     \\\\\\    \\\\\\    //    //  /  /  ////  //   \  ##
##                                                             ##
##################################################### uillar ####

# Load packages & set wd
{
  library(tidyverse)
  library(dplyr)
  setwd(".../TSS/")
}

# Open data
#FANTOM5 <- read.delim("inputs/FANTOM5_filtered.bed", header = F)
#refTSS_v4.1 <- read.delim("inputs/refTSS_v4.1_annot.txt", header = T)

# Select columns and rename
#FANTOM5 <- FANTOM5[, c(1:6)]
#FANTOM5$DB <- c("FANTOM5")
#colnames(FANTOM5) <- c("chr", "start", "end", "ID", "diSco", "strand", "DB")
#refTSS_v4.1$DB <- c("refTSS_v4.1")
#colnames(refTSS_v4.1) <- c("chr", "start", "end", "ID", "diSco", "strand", "idents", "DB")

# Split databases in +/- strands
#FANTOM5 <- FANTOM5 %>% mutate(mean = round((start + end) / 2))
#FANTOM5P <- subset(FANTOM5, FANTOM5$strand == "+")
#FANTOM5M <- subset(FANTOM5, FANTOM5$strand == "-")
#rm(FANTOM5)
#refTSS_v4.1 <- refTSS_v4.1 %>% mutate(mean = round((start + end) / 2))
#refTSS_v4.1P <- subset(refTSS_v4.1, refTSS_v4.1$strand == "+")
#refTSS_v4.1M <- subset(refTSS_v4.1, refTSS_v4.1$strand == "-")
#rm(refTSS_v4.1)
#save.image("inputs/TSS.RData")

# Open previously saved data
load("inputs/TSS.RData")

# load candidate genes with lncRNAKB annotation with strand info
interest <- read.delim("inputs/TSS_input.txt", header = T, sep = "\t", row.names = 1)

# How much candidates per strand and biotype?
table(interest$strand)
table(interest$Biotype)
table(interest$strand, interest$Biotype)

# Separate interest list in strand +/-
interestP <- subset(interest, interest$strand == "+")
interestM <- subset(interest, interest$strand == "-")
interest.P <- subset(interest, interest$strand == ".")
interest.P$GeneID <- paste0(interest.P$GeneID, "_asP")
interest.M <- subset(interest, interest$strand == ".")
interest.M$GeneID <- paste0(interest.M$GeneID, "_asM")

##############
## + strand ##
##############

# chr column as character
interestP$chr <- as.character(interestP$chr)
refTSS_v4.1P$chr <- as.character(refTSS_v4.1P$chr)
FANTOM5P$chr <- as.character(FANTOM5P$chr)

# Empty dataframe to store the results
resultP <- data.frame()

# Iterate with each row in interestP
for (i in 1:nrow(interestP)) {
  # Actual row in interestP
  current_row <- interestP[i, ]
  
  # Check if interest gene is in refTSS annotated database by geneSymbol (idents)
  matching_rows_1 <- refTSS_v4.1P[refTSS_v4.1P$idents == current_row$idents, ]
  
  if (nrow(matching_rows_1) > 0) {
    # Encontrar el valor mínimo absoluto de diSco
    min_abs_diSco <- min(abs(matching_rows_1$diSco))
    
    # Filtrar las filas con el valor mínimo absoluto de diSco
    closest_rows <- matching_rows_1[abs(matching_rows_1$diSco) == min_abs_diSco, ]
    
    # Verificar si hay más de una fila con el mínimo valor absoluto
    if (nrow(closest_rows) > 1) {
      # Inicializar un contador para el mismo valor más cercano a 0
      count <- 0
    } else {
      # Si solo hay una fila, el contador se mantiene en 0
      count <- NA
    }
    
    # Iterar sobre las filas filtradas y guardar todas las instancias con el mismo valor más cercano a 0
    for (j in seq_len(nrow(closest_rows))) {
      current_row_closest <- closest_rows[j, ]
      
      # Modificar el nombre solo si hay más de una fila con el mínimo valor absoluto
      if (!is.na(count)) {
        # Modificar el nombre para reflejar el contador (A, B, ...)
        modified_name <- paste0(rownames(current_row), "_", letters[count + 1])
      } else {
        # Si solo hay una fila, mantener el nombre sin cambios
        modified_name <- rownames(current_row)
      }
      
      # Añadir GeneID y otras columnas necesarias
      current_row_closest$GeneID <- modified_name
      current_row_closest$mean_minus <- current_row_closest$mean - 200
      current_row_closest$mean_plus <- current_row_closest$mean + 1
      
      # Unir los resultados al dataframe resultP
      resultP <- rbind(resultP, current_row_closest)
      
      # Imprimir un mensaje indicando que el GeneID está anotado en refTSS
      cat("GeneID", modified_name, "is annotated in refTSS as", current_row$idents, "\n")
      
      # Incrementar el contador solo si hay más de una fila con el mínimo valor absoluto
      if (!is.na(count)) {
        count <- count + 1
      }
    }
    
    # Moverse al siguiente GeneID en interestP
    next
    
  } else {
    # Print a message indicating that the GeneID is not annotated in refTSS
    cat("GeneID", rownames(current_row), "is NOT annotated in refTSS\n")
  }
  
  # If no match in refTSS, filter FANTOM5P to get just same chr TSS entries
  matching_rows2 <- FANTOM5P[FANTOM5P$chr == current_row$chr, ]
  
  # Filter possible TSS entries based on distance
  filtered_rows <- matching_rows2[abs(matching_rows2$mean - current_row$start) < 200, ]
  if (nrow(filtered_rows) > 0) {
    # If there are matches, select the row with the highest score
    selected_row <- filtered_rows[which.max(filtered_rows$diSco), ]
    
    # Add GeneID & target region for FASTA sequence retrieval to the selected row
    selected_row$GeneID <- rownames(current_row)
    selected_row$mean_minus <- selected_row$mean - 200
    selected_row$mean_plus <- selected_row$mean + 1
    
    # Adapt "selected_row" for rbind
    selected_row$idents <- selected_row$GeneID
    selected_row<- selected_row[, c(1:6,12, 7:11)]
    
    # Bind individual gene results to the results dataframe
    resultP <- rbind(resultP, selected_row)
    
    # Print a message indicating that TSS peaks were found for the GeneID
    cat("TSS peaks found for GeneID", rownames(current_row), "\n")
  } else {
  
    # If no matches, add a row to resultP using the exact position
    new_row <- current_row
    new_row$GeneID <- rownames(current_row)
    new_row$mean_minus <- current_row$start - 200
    new_row$mean_plus <- current_row$start + 1
    new_row$diSco <- c("NA")
    new_row$DB <- c("lncRNAKB")
    new_row$mean <- current_row$start
    new_row <- new_row[, c(4,5,6,3,10,7,3,11,12,1,8,9)]
    colnames(new_row) <- colnames(resultP)

    # Bind the new row to the results dataframe
    resultP <- rbind(resultP, new_row)

    # Print a message indicating that no TSS peaks was found
    cat("No TSS peaks found within 200bp window for GeneID", rownames(current_row), "--> added the region by exact position\n")
    
    # Move to the next GeneID in interestP
    next
  }
}

write.table(resultP, "TSS_candidatesP.txt", sep = "\t")

resultP_bed <- resultP[, c(1,11,12,10,5,6)]
resultP_bed$chr <- sub("chr", "", resultP_bed$chr)

##############
## - strand ##
##############

# Same for - strand genes/TSS
# chr column as character
interestM$chr <- as.character(interestM$chr)
refTSS_v4.1M$chr <- as.character(refTSS_v4.1M$chr)
FANTOM5M$chr <- as.character(FANTOM5M$chr)

# Empty dataframe to store the results
resultM <- data.frame()

# Iterate with each row in interestM
for (i in 1:nrow(interestM)) {
  # Actual row in interestM
  current_row <- interestM[i, ]
  
  # Check if interest gene is in refTSS annotated database by geneSymbol (idents)
  matching_rows_1 <- refTSS_v4.1M[refTSS_v4.1M$idents == current_row$idents, ]
  
  if (nrow(matching_rows_1) > 0) {
    # Encontrar el valor mínimo absoluto de diSco
    min_abs_diSco <- min(abs(matching_rows_1$diSco))
    
    # Filtrar las filas con el valor mínimo absoluto de diSco
    closest_rows <- matching_rows_1[abs(matching_rows_1$diSco) == min_abs_diSco, ]
    
    # Verificar si hay más de una fila con el mínimo valor absoluto
    if (nrow(closest_rows) > 1) {
      # Inicializar un contador para el mismo valor más cercano a 0
      count <- 0
    } else {
      # Si solo hay una fila, el contador se mantiene en 0
      count <- NA
    }
    
    # Iterar sobre las filas filtradas y guardar todas las instancias con el mismo valor más cercano a 0
    for (j in seq_len(nrow(closest_rows))) {
      current_row_closest <- closest_rows[j, ]
      
      # Modificar el nombre solo si hay más de una fila con el mínimo valor absoluto
      if (!is.na(count)) {
        # Modificar el nombre para reflejar el contador (A, B, ...)
        modified_name <- paste0(rownames(current_row), "_", letters[count + 1])
      } else {
        # Si solo hay una fila, mantener el nombre sin cambios
        modified_name <- rownames(current_row)
      }
      
      # Añadir GeneID y otras columnas necesarias
      current_row_closest$GeneID <- modified_name
      current_row_closest$mean_minus <- current_row_closest$mean - 1
      current_row_closest$mean_plus <- current_row_closest$mean + 200
      
      # Unir los resultados al dataframe resultM
      resultM <- rbind(resultM, current_row_closest)
      
      # Imprimir un mensaje indicando que el GeneID está anotado en refTSS
      cat("GeneID", modified_name, "is annotated in refTSS as", current_row$idents, "\n")
      
      # Incrementar el contador solo si hay más de una fila con el mínimo valor absoluto
      if (!is.na(count)) {
        count <- count + 1
      }
    }
    
    # Moverse al siguiente GeneID en interestM
    next
    
  } else {
    # Print a message indicating that the GeneID is not annotated in refTSS
    cat("GeneID", rownames(current_row), "is NOT annotated in refTSS\n")
  }
  
  # If no match in refTSS, filter FANTOM5M to get just same chr TSS entries
  matching_rows2 <- FANTOM5M[FANTOM5M$chr == current_row$chr, ]
  
  # Filter possible TSS entries based on distance
  filtered_rows <- matching_rows2[abs(matching_rows2$mean - current_row$start) < 200, ]
  if (nrow(filtered_rows) > 0) {
    # If there are matches, select the row with the highest score
    selected_row <- filtered_rows[which.max(filtered_rows$diSco), ]
    
    # Add GeneID & target region for FASTA sequence retrieval to the selected row
    selected_row$GeneID <- rownames(current_row)
    selected_row$mean_minus <- selected_row$mean - 1
    selected_row$mean_plus <- selected_row$mean + 200
    
    # Adapt "selected_row" for rbind
    selected_row$idents <- selected_row$GeneID
    selected_row<- selected_row[, c(1:6,12, 7:11)]
    
    # Bind individual gene results to the results dataframe
    resultM <- rbind(resultM, selected_row)
    
    # Print a message indicating that TSS peaks were found for the GeneID
    cat("TSS peaks found for GeneID", rownames(current_row), "\n")
  } else {
  
    # If no matches, add a row to resultM using the exact position
    new_row <- current_row
    new_row$GeneID <- rownames(current_row)
    new_row$mean_minus <- current_row$start - 1
    new_row$mean_plus <- current_row$start + 200
    new_row$diSco <- c("NA")
    new_row$DB <- c("lncRNAKB")
    new_row$mean <- current_row$end
    new_row <- new_row[, c(4,5,6,3,10,7,3,11,12,1,8,9)]
    colnames(new_row) <- colnames(resultM)

    # Bind the new row to the results dataframe
    resultM <- rbind(resultM, new_row)

    # Print a message indicating that no TSS peaks was found
    cat("No TSS peaks found within 200bp window for GeneID", rownames(current_row), "--> added the region by exact position\n")
    
    # Move to the next GeneID in interestP
    next
  }
}

write.table(resultM, "TSS_candidatesM.txt", sep = "\t")

resultM_bed <- resultM[, c(1,11,12,10,5,6)]
resultM_bed$chr <- sub("chr", "", resultM_bed$chr)

##################
## . strand asP ##
##################

# For . strand genes we will try to find peaks trying both options --> here as + strand (asP)
# chr column as character
interest.P$chr <- as.character(interest.P$chr)
refTSS_v4.1P$chr <- as.character(refTSS_v4.1P$chr)
FANTOM5P$chr <- as.character(FANTOM5P$chr)

# Empty dataframe to store the results
result.P <- data.frame()

# Iterate with each row in interest.P
for (i in 1:nrow(interest.P)) {
  # Actual row in interest.P
  current_row <- interest.P[i, ]
  
  # Check if interest gene is in refTSS annotated database by geneSymbol (idents)
  matching_rows_1 <- refTSS_v4.1P[refTSS_v4.1P$idents == current_row$idents, ]
  
  if (nrow(matching_rows_1) > 0) {
    # Encontrar el valor mínimo absoluto de diSco
    min_abs_diSco <- min(abs(matching_rows_1$diSco))
    
    # Filtrar las filas con el valor mínimo absoluto de diSco
    closest_rows <- matching_rows_1[abs(matching_rows_1$diSco) == min_abs_diSco, ]
    
    # Verificar si hay más de una fila con el mínimo valor absoluto
    if (nrow(closest_rows) > 1) {
      # Inicializar un contador para el mismo valor más cercano a 0
      count <- 0
    } else {
      # Si solo hay una fila, el contador se mantiene en 0
      count <- NA
    }
    
    # Iterar sobre las filas filtradas y guardar todas las instancias con el mismo valor más cercano a 0
    for (j in seq_len(nrow(closest_rows))) {
      current_row_closest <- closest_rows[j, ]
      
      # Modificar el nombre solo si hay más de una fila con el mínimo valor absoluto
      if (!is.na(count)) {
        # Modificar el nombre para reflejar el contador (A, B, ...)
        modified_name <- paste0(rownames(current_row), "_", letters[count + 1])
      } else {
        # Si solo hay una fila, mantener el nombre sin cambios
        modified_name <- rownames(current_row)
      }
      
      # Añadir GeneID y otras columnas necesarias
      current_row_closest$GeneID <- modified_name
      current_row_closest$mean_minus <- current_row_closest$mean - 200
      current_row_closest$mean_plus <- current_row_closest$mean + 1
      
      # Unir los resultados al dataframe result.P
      result.P <- rbind(result.P, current_row_closest)
      
      # Imprimir un mensaje indicando que el GeneID está anotado en refTSS
      cat("GeneID", current_row$GeneID, "is annotated in refTSS as", current_row$idents, "\n")
      
      # Incrementar el contador solo si hay más de una fila con el mínimo valor absoluto
      if (!is.na(count)) {
        count <- count + 1
      }
    }
    
    # Moverse al siguiente GeneID en interest.P
    next
    
  } else {
    # Print a message indicating that the GeneID is not annotated in refTSS
    cat("GeneID", current_row$GeneID, "is NOT annotated in refTSS\n")
  }
  
  # If no match in refTSS, filter FANTOM5P to get just same chr TSS entries
  matching_rows2 <- FANTOM5P[FANTOM5P$chr == current_row$chr, ]
  
  # Filter possible TSS entries based on distance
  filtered_rows <- matching_rows2[abs(matching_rows2$mean - current_row$start) < 200, ]
  if (nrow(filtered_rows) > 0) {
    # If there are matches, select the row with the highest score
    selected_row <- filtered_rows[which.max(filtered_rows$diSco), ]
    
    # Add GeneID & target region for FASTA sequence retrieval to the selected row
    selected_row$GeneID <- rownames(current_row)
    selected_row$mean_minus <- selected_row$mean - 200
    selected_row$mean_plus <- selected_row$mean + 1
    
    # Adapt "selected_row" for rbind
    selected_row$idents <- selected_row$GeneID
    selected_row<- selected_row[, c(1:6,12, 7:11)]
    
    # Bind individual gene results to the results dataframe
    result.P <- rbind(result.P, selected_row)
    
    # Print a message indicating that TSS peaks were found for the GeneID
    cat("TSS peaks found for GeneID", current_row$GeneID, "\n")
  } else {
      
    # If no matches, add a row to result.P using the exact position
    new_row <- current_row
    new_row$GeneID <- rownames(current_row)
    new_row$mean_minus <- current_row$start - 200
    new_row$mean_plus <- current_row$start + 1
    new_row$diSco <- c("NA")
    new_row$DB <- c("lncRNAKB")
    new_row$mean <- current_row$end
    new_row <- new_row[, c(4,5,6,3,10,7,3,11,12,1,8,9)]
    colnames(new_row) <- colnames(result.P)

    # Bind the new row to the results dataframe
    result.P <- rbind(result.P, new_row)

    # Print a message indicating that no TSS peaks was found
    cat("No TSS peaks found within 200bp window for GeneID", rownames(current_row), "--> added the region by exact position\n")
    
    # Move to the next GeneID in interestP
    next
  }
}

write.table(result.P, "TSS_candidates.P.txt", sep = "\t")

result.P_bed <- result.P[, c(1,11,12,10,5,6)]
result.P_bed$chr <- sub("chr", "", result.P_bed$chr)

##################
## . strand asM ##
##################

# For . strand genes we will try to find peaks trying both options --> here as + strand (asM)
# chr column as character
interest.M$chr <- as.character(interest.M$chr)
refTSS_v4.1M$chr <- as.character(refTSS_v4.1M$chr)
FANTOM5M$chr <- as.character(FANTOM5M$chr)

# Empty dataframe to store the results
result.M <- data.frame()

# Iterate with each row in interest.M
for (i in 1:nrow(interest.M)) {
  # Actual row in interest.M
  current_row <- interest.M[i, ]
  
  # Check if interest gene is in refTSS annotated database by geneSymbol (idents)
  matching_rows_1 <- refTSS_v4.1M[refTSS_v4.1M$idents == current_row$idents, ]
  
  if (nrow(matching_rows_1) > 0) {
    # Encontrar el valor mínimo absoluto de diSco
    min_abs_diSco <- min(abs(matching_rows_1$diSco))
    
    # Filtrar las filas con el valor mínimo absoluto de diSco
    closest_rows <- matching_rows_1[abs(matching_rows_1$diSco) == min_abs_diSco, ]
    
    # Verificar si hay más de una fila con el mínimo valor absoluto
    if (nrow(closest_rows) > 1) {
      # Inicializar un contador para el mismo valor más cercano a 0
      count <- 0
    } else {
      # Si solo hay una fila, el contador se mantiene en 0
      count <- NA
    }
    
    # Iterar sobre las filas filtradas y guardar todas las instancias con el mismo valor más cercano a 0
    for (j in seq_len(nrow(closest_rows))) {
      current_row_closest <- closest_rows[j, ]
      
      # Modificar el nombre solo si hay más de una fila con el mínimo valor absoluto
      if (!is.na(count)) {
        # Modificar el nombre para reflejar el contador (A, B, ...)
        modified_name <- paste0(rownames(current_row), "_", letters[count + 1])
      } else {
        # Si solo hay una fila, mantener el nombre sin cambios
        modified_name <- rownames(current_row)
      }
      
      # Añadir GeneID y otras columnas necesarias
      current_row_closest$GeneID <- modified_name
      current_row_closest$mean_minus <- current_row_closest$mean - 1
      current_row_closest$mean_plus <- current_row_closest$mean + 200
      
      # Unir los resultados al dataframe result.M
      result.M <- rbind(result.M, current_row_closest)
      
      # Imprimir un mensaje indicando que el GeneID está anotado en refTSS
      cat("GeneID", current_row$GeneID, "is annotated in refTSS as", current_row$idents, "\n")
      
      # Incrementar el contador solo si hay más de una fila con el mínimo valor absoluto
      if (!is.na(count)) {
        count <- count + 1
      }
    }
    
    # Moverse al siguiente GeneID en interest.M
    next
    
  } else {
    # Print a message indicating that the GeneID is not annotated in refTSS
    cat("GeneID", current_row$GeneID, "is NOT annotated in refTSS\n")
  }
  
  # If no match in refTSS, filter FANTOM5M to get just same chr TSS entries
  matching_rows2 <- FANTOM5M[FANTOM5M$chr == current_row$chr, ]
  
  # Filter possible TSS entries based on distance
  filtered_rows <- matching_rows2[abs(matching_rows2$mean - current_row$start) < 200, ]
  if (nrow(filtered_rows) > 0) {
    # If there are matches, select the row with the highest score
    selected_row <- filtered_rows[which.max(filtered_rows$diSco), ]
    
    # Add GeneID & target region for FASTA sequence retrieval to the selected row
    selected_row$GeneID <- rownames(current_row)
    selected_row$mean_minus <- selected_row$mean - 1
    selected_row$mean_plus <- selected_row$mean + 200
    
    # Adapt "selected_row" for rbind
    selected_row$idents <- selected_row$GeneID
    selected_row<- selected_row[, c(1:6,12, 7:11)]
    
    # Bind individual gene results to the results dataframe
    result.M <- rbind(result.M, selected_row)
    
    # Print a message indicating that TSS peaks were found for the GeneID
    cat("TSS peaks found for GeneID", current_row$GeneID, "\n")
  } else {
      
    # If no matches, add a row to result.M using the exact position
    new_row <- current_row
    new_row$GeneID <- rownames(current_row)
    new_row$mean_minus <- current_row$start - 1
    new_row$mean_plus <- current_row$start + 200
    new_row$diSco <- c("NA")
    new_row$DB <- c("lncRNAKB")
    new_row$mean <- current_row$end
    new_row <- new_row[, c(4,5,6,3,10,7,3,11,12,1,8,9)]
    colnames(new_row) <- colnames(result.M)

    # Bind the new row to the results dataframe
    result.M <- rbind(result.M, new_row)

    # Print a message indicating that no TSS peaks was found
    cat("No TSS peaks found within 200bp window for GeneID", rownames(current_row), "--> added the region by exact position\n")
    
    # Move to the next GeneID in interestP
    next
  }
}

write.table(result.M, "TSS_candidates.M.txt", sep = "\t")

result.M_bed <- result.M[, c(1,11,12,10,5,6)]
result.M_bed$chr <- sub("chr", "", result.M_bed$chr)

##################
## FlashFry BED ##
##################

FlashFry_input <- rbind(resultP_bed, resultM_bed, result.P_bed, result.M_bed)

write.table(FlashFry_input, "FlashFry/FlashFry_input.bed", sep = "\t", col.names = F, row.names = F, quote = FALSE)
