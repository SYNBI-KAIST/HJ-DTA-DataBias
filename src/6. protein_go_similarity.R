library(parallel)
library(foreach)
library(doParallel)
library(data.table)
library(dplyr)
library(protr)
library(GOSemSim)
library(org.Hs.eg.db)


## Calculate average binding affinity for each Ikey in training data
dta <- fread('./Training_Set.csv', header = T, stringsAsFactors = F)
dta <- dta %>% distinct(Ikey, UNIPROT_AC, .keep_all = TRUE)
prn_info <- fread('./Input_Prn_Info.csv', header = T, stringsAsFactors = F) %>% filter(Organism == 'Homo sapiens (Human)') %>% dplyr::select(UNIPROT_AC, Sequence, GeneID) %>% na.omit()
goData <- godata('org.Hs.eg.db', ont = "BP")

library(parallel)
library(doParallel)
library(data.table)
library(protr)
library(foreach)

nodes <- paste0("node", 50:59)
cl <- makeCluster(nodes, type = "SOCK")

clusterEvalQ(cl, {
  library(doParallel)
  library(protr)
  library(foreach)
  library(parallel)
  library(GOSemSim)
})

low_cv_dta <- as.data.table(dta)
prn_info <- as.data.table(prn_info)
unique_mappings <- unique(dta, by = c("Ikey", "UNIPROT_AC"))
mapped_sequences <- merge(unique_mappings, prn_info, by = "UNIPROT_AC")

ikeys <- unique(mapped_sequences$Ikey)
ikeys_split <- split(ikeys, cut(seq_along(ikeys), length(cl)))

# calculateSimilarity
calculateSimilarity <- function(subset_keys, mapped_sequences) {
  results <- mclapply(subset_keys, function(ikey) {
    subset_data <- mapped_sequences[mapped_sequences$Ikey == ikey, ]
    temp_results <- data.frame(Ikey = character(), UNIPROT_AC1 = character(), UNIPROT_AC2 = character(), Similarity = numeric(), stringsAsFactors = FALSE)
    
    if (nrow(subset_data) > 1) {
      for (i in 1:(nrow(subset_data) - 1)) {
        for (j in (i + 1):nrow(subset_data)) {
          sim_score <- geneSim(as.character(subset_data$GeneID[i]), as.character(subset_data$GeneID[j]), semData = goData, measure = "Wang")
          
          if (length(sim_score) > 1){
            temp_results <- dplyr::bind_rows(temp_results, data.frame(Ikey = ikey, UNIPROT_AC1 = subset_data$UNIPROT_AC[i], UNIPROT_AC2 = subset_data$UNIPROT_AC[j], Similarity = sim_score$geneSim, stringsAsFactors = FALSE))
          }
        }
      }
    }
    return(temp_results)
  }, mc.cores = 24)  # Adjust mc.cores to utilize available cores effectively
  
  return(do.call(rbind, results))
}


results <- clusterApply(cl, ikeys_split, function(keys) calculateSimilarity(keys, mapped_sequences))

stopCluster(cl)

similarity_results <- do.call(rbind, results)
similarity_results <- as.data.table(similarity_results)