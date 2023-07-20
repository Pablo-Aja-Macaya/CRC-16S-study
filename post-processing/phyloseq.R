

# Libraries
library(purrr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(plyr)
library(vegan)
library("devtools")
library(data.table)
library(ggtree)
library(phyloseq)
library(plotly)
library(VennDiagram)
require(RColorBrewer)
library(egg)
library(randomcoloR)
library(glue)
library(tidyr)
library(ANCOMBC)
library("ggvenn")


# Options
options(dplyr.summarise.inform = FALSE)
set.seed(123)

# Input
feature_table_file <- "./data/feature-table.json"
contamination_file <- "./data/contaminants.tsv"
metadata_file <- "./data/non-ffpe_metadata.tsv"
taxmat_file <- "./data/separated_taxonomy.tsv"
tree_file <- "./data/tree.nwk"
output_folder <- "./output"
biosample_col <- "biosample_name"
variables_of_interest <- c("origin", "sample.type")

cat("\nStarting Phyloseq module!\n")

#######################
# ---- Load data ---- #
#######################

# ---- OTUs -----
cat("\nProccessing .biom...\n")

# Transform biome data to JSON.biome
# biom convert -i in.biom -o out.json.biom --table-type="OTU table" --to-json
otu_df <- import_biom(feature_table_file)

# Eliminar muestras donde toda la columna sea 0 (no tienen datos)
condition <- colSums(otu_df) != 0
otu_df_empty_cols <- colnames(otu_df)[!condition]
if (length(otu_df_empty_cols)>0){
  print(glue("WARNING: There are {length(otu_df_empty_cols)} empty columns in ASV dataframe (no data, this can happen in controls) and they will be removed: '{paste(otu_df_empty_cols, collapse='; ')}'"))
}
otu_df <- otu_df[, condition]

print(head(otu_df))

# ---- Metadata ----
cat("\nProccessing metadata...\n")

metadata <- read.delim(metadata_file, comment.char = "#")

# Check if there are columns not in otu_df besides the otu_df columns that were empty
condition <- metadata[,'sample.id'] %in% colnames(otu_df)
sample.ids.not.in.otus <- metadata[!condition,]$sample.id
if (length(setdiff(otu_df_empty_cols, sample.ids.not.in.otus))>0){
  stop("There are unexpected samples in metadata that do not exist in ASV dataframe")
}
metadata <- metadata[condition,] # eliminar muestras si no existen en otu_df
rownames(metadata) <- metadata[,'sample.id']

# Create column for grouping by variables_of_interest
if ("grouping_var" %in% names(metadata)){
  stop("Ups, metadata already has grouping_var as column")
} else {
  metadata$grouping_var <- do.call(paste, c(metadata[variables_of_interest], sep="__"))
}

# Make categorical columns
for (i in variables_of_interest){
  metadata[[i]] <- factor(metadata[[i]])
}

print(head(metadata))

# ---- Taxmat ----
cat("\nProccessing taxonomy...\n")

# Taxonomy must be split into columns
taxmat = as.matrix(read.delim(taxmat_file, sep='\t'))

# Propagate NAs
keep_unidentified <- function(tx_table){
  # Keep the ones that have one of c(NA, "", " ", "\t") in target_level column
  # and then find the last level that is known (can be Genus, Family...)
  # Assign to the target_level level the last known level + x__XXX NA
  # Family (last known)   Species
  # f__Lachnospiraceae    f__Lachnospiraceae NA
  elements_to_target <- colnames(tx_table)
  tx_table <- as.data.frame(tx_table)
  previous_level <- elements_to_target[1]
  for (level in elements_to_target[2:length(elements_to_target)]){
    empty_values <- tx_table[[level]] %in% c(NA, "", " ", "\t")
    tx_table[[level]][empty_values] <- paste(tx_table[[previous_level]][empty_values], "NA", sep="_")
    tx_table[[level]][empty_values] <- gsub("(_NA[> ]*)*", "\\1", tx_table[[level]][empty_values]) # Replace repetitions of NA by one NA

    previous_level <- level
  }

  return(tx_table)
}
taxmat <- as.matrix(keep_unidentified(taxmat))

# Set row names
rownames(taxmat) <- rownames(otu_df)

print(head(taxmat))

#####################
# ---- Control ---- #
#####################
# Calcular la suma de cada ASV en los controles
sum_rows <- function(input){
  # Si el input es dataframe hacer rowSums para obtener un vector
  # Si ya es un vector devolverlo
  if (is.null(ncol(input))){
    return (input)
  } else {
    return (rowSums(input))
  }
}

substract_control <- function(features, metadata, importance_multiplier=1){
  # Substract ASVs that appear in controls per sequencing-run
  # - features: ASV table
  # - metadata: dataframe with metadata
  # - importance_multiplier: increase this to a big number to fully remove ASVs that appear in controls (like 999999)
  cat("\nUsing controls for cleaning...\n")
  
  # Initialize clean dataframe with non-control ids
  clean.otu.df <- features[, subset(metadata, origin!="control")$sample.id]
  # clean.otu.df <- features
  
  # Split "control.applies.to.origin" column by ";"
  # The resulting elements will be searched in "origin" column 
  # and the control will be applied to corresponding samples
  metadata$control.applies.to_split <- strsplit(metadata$control.applies.to.origin, ";")
  
  # For each control ID
  for (control_id in subset(metadata, origin=="control")$sample.id ){
    # Get the metadata of the control ID
    r <- subset(metadata, sample.id==control_id)
    
    # Log
    control_id_seq_run = r$sequencing.run
    print(glue("\n==== Using sample.id '{control_id}' as control in sequencing.run='{control_id_seq_run}' ===="))
    
    # Check that if element "all" is present no other elements exist
    # and that an element is not repeated with unique
    applies_to_elements <- unique(r$control.applies.to_split[[1]]) # To access split elements: metadata$control.applies.to_split[65][[1]][1]
    if ("all" %in% applies_to_elements & length(applies_to_elements)>1){
      stop("Element 'all' was specified, but extra elements are present. Please choose either 'all' or the other elements")
    }
    
    # For each origin to which the control applies to
    for (applies_to in applies_to_elements){ 
      print(glue("\n~~~ Applying control to: '{applies_to}' ~~~"))
      
      tmp <- subset(metadata, origin!="control" & sequencing.run==control_id_seq_run)
      
      # If the control applies to all samples
      if (applies_to == "all"){
        tmp <- tmp
        # If it only applies to a certain origin type which IS in metadata$origin
      } else if (applies_to %in% metadata$origin){
        tmp <- subset(tmp, origin==applies_to)
        # If it is not "all" or the element doesnt exist in metadata$origin, raise error 
      # If it is not "all" or the element doesnt exist in metadata$origin, raise error 
        # If it is not "all" or the element doesnt exist in metadata$origin, raise error 
      } else {
        warning(glue("Origin '{applies_to}' is not in metadata$origin column, please fix"))
      }
      
      if (length(tmp$sample.id) >= 1){
        print(glue("Applying control.id={control_id} from sequencing.run={control_id_seq_run} to element={applies_to} in {length(tmp$sample.id)} samples"), "\n")
        clean.otu.df[,tmp$sample.id] <- clean.otu.df[,tmp$sample.id] - sum_rows(features[,c(control_id)])*importance_multiplier      
      } else {
        stop(glue("WARNING: Control control.id={control_id} from sequencing.run={control_id_seq_run} cant be applied to element={applies_to}, no samples match ({length(tmp$sample.id)} samples)\nThis in unexpected, every control should apply to at least one sample! Check metadata?"))
      }
    }
  }

  # Negative numbers to 0
  clean.otu.df[clean.otu.df<0] <- 0
  
  return(clean.otu.df)
}

clean.otu.df <- substract_control(otu_df, metadata)

########################
# ---- Physeq obj ---- #
########################
cat("Creating Physeq object...\n")

# ---- Phyloseq object ----
OTU = otu_table(clean.otu.df, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAM = sample_data(metadata)
TREE = read_tree(tree_file)

physeq = phyloseq(OTU, TAX, SAM, TREE)
physeq

# ---- Apply simple contamination filter using genus ----
pop_taxa = function(ps, bad_asvs){
  all_asvs = taxa_names(ps)
  keep_these <- all_asvs[!(all_asvs %in% bad_asvs)]
  return(prune_taxa(keep_these, ps))
}
# Only bacteria
print(glue("Filtering to only bacteria..."))

bad_asvs <- taxa_names(subset_taxa(tax_table(physeq), Domain!="d__Bacteria"))
physeq <- pop_taxa(physeq, bad_asvs)

# Remove typical contamination
print(glue("Removing taxa based on {contamination_file}..."))
typical_contamination <- read.csv(contamination_file, sep='\t')
typical_contamination$Genus <- paste("g__", typical_contamination$Genus, sep="")

print("Using genus in contamination file")
cont_genus <- c(typical_contamination$Genus, paste(typical_contamination$Genus, "NA", sep= "_"))
physeq <- subset_taxa(physeq, !Genus %in% cont_genus)

# Manual names
print("Manual names")
bad_asvs <- taxa_names(
  subset_taxa(
    tax_table(physeq), 
    Family %in% c("f__Dermacoccaceae",
                  "f__Mitochondria",
                  "f__Chloroplast"
                  ) |
    Genus %in% c(
      "g__Burkholderia-Caballeronia-Paraburkholderia",
      "g__Methylobacterium-Methylorubrum",
      "g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
      "g__Acinetobacter"
    )
  )
)
physeq <- pop_taxa(physeq, bad_asvs)

# Keep ASVs with more than one count
physeq <- prune_taxa(taxa_sums(physeq) > 1, physeq) 

# Remove 1 sample that is possibly wrong (studied through manual beta-diversity analysis, not shown)
physeq <- subset_samples(physeq, 
                         !sample.id %in% c(
                           "CCR-142-1-s" # In theory it is saliva, but looks like faeces
                         )
)


# ---- Simplify metadata ----
temp <- data.frame(sample_data(physeq))

# Simplify tnm
temp$has_affected_lymph_nodes <- as.integer(temp$tnm_n != "N0")
temp$has_big_tumor <- as.integer(temp$tnm_t %in% c("T3","T4"))
temp$has_metastasis <- as.integer(temp$tnm_m %in% c("M1","M2"))

# Simplify age
temp <- temp %>% mutate(age_group = case_when(
  age < 60 ~ 'under60',
  age >= 60  & age < 70 ~ '60to69',
  age >= 70 ~ 'over70',
)
) 

# Simplify tumor location
temp$simplified_location <- temp$localization
temp$simplified_location[temp$simplified_location %in% c("right-colon", "hepatic-flexure", "cecum")] <- "RIGHT"
temp$simplified_location[temp$simplified_location %in% c("left-colon", "sigmoid-colon", "splenic-flexure")] <- "LEFT"

# Add odontological consult metadata
# odontological_data <- read.delim("/home/usuario/Proyectos/Results/Kelly/KellyCCR/data/metadata/crc_consulta_dental_metadata.tsv")
# temp <- merge(temp, odontological_data[-2], by.x="subject", by.y="subject", all.x=TRUE)

# Return new columns to phyloseq object
sample_data(physeq)$has_affected_lymph_nodes <- temp$has_affected_lymph_nodes
sample_data(physeq)$has_big_tumor <- temp$has_big_tumor
sample_data(physeq)$has_metastasis <- temp$has_metastasis
sample_data(physeq)$simplified_location <- temp$simplified_location
sample_data(physeq)$age_group <- temp$age_group
# sample_data(physeq)$periodontal_disease_stage <- temp$periodontal_disease_stage
# sample_data(physeq)$silness_loe_index <- temp$silness_loe_index
# sample_data(physeq)$absent_teeth <- temp$absent_teeth
# sample_data(physeq)$caries_index <- temp$caries_index
# sample_data(physeq)$infectious_bone_pathologies <- temp$infectious_bone_pathologies

# ---- Set up factors ----
sample_data(physeq)$sample.type <- factor(sample_data(physeq)$sample.type,      # Reordering group factor levels
                                          levels = c("non-crc", "crc")
)
sample_data(physeq)$origin <- factor(sample_data(physeq)$origin,      # Reordering group factor levels
                                     levels = c("adenocarcinoma","faeces","normal-mucosa","saliva","subgingival-fluid")
)

more_factor_cols = c("simplified_location","localization","subject","age_group","sex",
                     "tnm_t","tnm_n","tnm_m",
                     "has_big_tumor","has_metastasis","has_affected_lymph_nodes"
                     )
sample_data(physeq)[,more_factor_cols] <- lapply(sample_data(physeq)[,more_factor_cols], factor)


# ---- Calculate rarified phyloseq object ----
rar_physeq <- rarefy_even_depth(physeq, sample.size = 10000,
                                  replace = FALSE, trimOTUs = TRUE, verbose = TRUE)



##########################################
# ---------- Sample stats -------------- #
##########################################
subject_metadata <- data.frame(sample_data(physeq)[,c("subject","age","sex","simplified_location","age_group","sample.type")])
rownames(subject_metadata) <- NULL
subject_metadata <- subject_metadata %>% dplyr::distinct()

# Age density
temp <- list(
  female_crc=length(subset(subject_metadata, sample.type=='crc' & sex=="female")$sex),
  male_crc=length(subset(subject_metadata, sample.type=='crc' & sex=="male")$sex),
  female_noncrc=length(subset(subject_metadata, sample.type=='non-crc' & sex=="female")$sex),
  male_noncrc=length(subset(subject_metadata, sample.type=='non-crc' & sex=="male")$sex)
)
p <- ggplot(subject_metadata, aes(x = age)) +
  geom_density(aes(color = sex, fill = sex), 
               position = "identity", alpha = 0.3) +
  facet_wrap(~sample.type) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  ggtitle(glue("Age distribution - (CRC: M={temp$male_crc}, F={temp$female_crc} | non-CRC: M={temp$male_noncrc}, F={temp$female_noncrc})")) +
  theme_bw()
p


###########################################################

save.physeq <- function(physeq_obj, output_file){
  # Merge taxonomy with ASVs
  merged.taxa.otus <- merge(otu_table(physeq), tax_table(physeq), by.x = 0, by.y = 0, all.x = TRUE, all.y = TRUE)
  names(merged.taxa.otus)[1] <- "ASV"
  
  # Merge tax into one column
  tax_cols <- colnames(tax_table(physeq))
  merged.taxa.otus$TAXONOMY <- do.call(paste, c(merged.taxa.otus[tax_cols], sep=";"))
  for (c in tax_cols) {merged.taxa.otus[c] <- NULL}
  
  # Make ASV into row index
  rownames(merged.taxa.otus) <- merged.taxa.otus$ASV
  
  # Count total
  merged.taxa.otus$TOTAL <- rowSums(merged.taxa.otus[,-which(colnames(merged.taxa.otus) %in% c("TAXONOMY","ASV"))])
  
  # Reorder columns
  merged.taxa.otus <- merged.taxa.otus[,c("ASV","TAXONOMY","TOTAL", 
                                         colnames(merged.taxa.otus)[!colnames(merged.taxa.otus) %in% c("TAXONOMY","TOTAL","ASV")]
                                       )
  ]
  merged.taxa.otus <- merged.taxa.otus[order(-merged.taxa.otus$TOTAL),]
  
  # Save
  write.table(merged.taxa.otus, file=output_file, sep='\t', row.names=F)
}
# save.physeq(physeq, "asv_table__clean__ccr_kelly_2022-11-08.tsv")


###########################################################

##########################################
# ---- Relative abundance functions ---- #
##########################################

check_rel_abund_sum <- function(data, group_cols, rel_abund_col){
  # Check relative abundances sum to 1
  sum_check <- data %>% 
    dplyr::rename("rel_abund_col" = rel_abund_col) %>%
    dplyr::group_by(across(all_of(c(group_cols)))) %>%
    dplyr::summarise(sum_relative_abundance = sum(rel_abund_col)) 
  # Check has to be done with "all.equal" because at these scales using "==" gives unexpected results
  check <- mapply(function(x) {isTRUE(all.equal(x, 1))}, sum_check$sum_relative_abundance)
  if (any(!check) == TRUE){
    print(sum_check)
    stop("ERROR: The sum of relative abundances for each group should be 1. At least one group fails this rule. This should never happen. Stopping...")
  }
  return(TRUE)
}

get.otus.to.keep.by.prevalence <- function(group, data, abundance_threshold, prevalence_threshold, sample_col="sample.id", abundance_col="Abundance"){
  starting.otus.count <- length(unique(data$OTU))
  
  print(glue("\nFiltering group composed of: '{paste(group, collapse='; ')}'"))
  
  # For each sample, filter to OTUs that have a relative abundance over N% 
  # and appear in at least M% of the other samples in its group
  abundance_over_list <- list()
  abundance_under_list <- list()
  for (i in group){
    # For each sample.id in group
    # Calculate relative abundance for the sample.id's data
    # Get OTUs over threshold into a list
    df <- data[ data[[sample_col]] == i , ]
    df$relative_abundance <- df[[abundance_col]]/sum(df[[abundance_col]])
    otus_over_thresh <- subset(df, relative_abundance>=abundance_threshold)$OTU
    otus_under_thresh <- subset(df, relative_abundance<abundance_threshold)$OTU
    abundance_over_list[[i]] <- otus_over_thresh
    abundance_under_list[[i]] <- otus_under_thresh
    print(glue("Sample {i} keeps {length(otus_over_thresh)}/{length(otus_under_thresh)+length(otus_over_thresh)} ASVs filtering by abundance={abundance_threshold}"))
  }
  
  # Transform list of vectors into a presence/absence dataframe by matching ASV names
  # leaving NAs if a sample does not have the ASV.
  # Each row is an ASV, each column is a sample
  # Example:
  #                             ccr11-1s                       CCR-59-1-s                       CCR-57-1-s
  # 1   101637aba474a6614786486ad42acb11 101637aba474a6614786486ad42acb11 101637aba474a6614786486ad42acb11
  # 2   00a2f7e2564becc62e71c46821779802                             <NA>                             <NA>
  # 3   f7cf33648c7d6bda19e991b97583cdab f7cf33648c7d6bda19e991b97583cdab f7cf33648c7d6bda19e991b97583cdab
  # 21  0a49a32752b611a144ddd506a2155134                             <NA> 0a49a32752b611a144ddd506a2155134
  # 22  5258515ff2816e94d24b2bd02f8ae803                             <NA> 5258515ff2816e94d24b2bd02f8ae803
  # 23  c773422b3dd3b7f057a6ec0d5fb93733 c773422b3dd3b7f057a6ec0d5fb93733                             <NA>
  otu_presence <- abundance_over_list %>%
    map(~ data.frame(col = ., ., stringsAsFactors = FALSE)) %>%
    reduce(full_join, by = "col") %>%
    select(-col) %>%
    setNames(names(abundance_over_list))
  
  # Set the OTUs as row names
  rownames(otu_presence) <- apply(otu_presence, 1, function(z) na.omit(z)[1])
  
  # Select which OTUs have to be kept by prevalence
  # The ASVs that appear in at least X% of the samples (row-wise)
  condition <- (rowSums(!is.na(otu_presence)) / dim(otu_presence)[2]) >= prevalence_threshold
  otus_to_keep <- rownames(otu_presence)[condition]
  otus_to_loose_by_prevalence <- rownames(otu_presence)[!condition]
  print(glue("Surviving ASVs to abundance filter kept by prevalence: {length(otus_to_keep)}/{length(otus_to_loose_by_prevalence)+length(otus_to_keep)}"))
  
  print(glue("Final: {length(otus_to_keep)}/{starting.otus.count} ASVs have passed the abundance+prevalence filter."))
  return(data %>% filter(OTU %in% otus_to_keep))
}

calculate_props <- function(data, variables_of_interest, target_level, threshold, abundance_col){
  # Convert to data.table
  data <- data.table(data)
  
  # Group by OTU, target_level, Type and Site
  grouped_df <- data %>%
    dplyr::rename("abundance" = abundance_col) %>%
    dplyr::group_by(across(all_of(c("OTU", target_level, variables_of_interest)))) %>%
    dplyr::summarise(mean_abundance = mean(abundance)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(relative_abundance = mean_abundance/sum(mean_abundance))
  
  # Check
  check_rel_abund_sum(grouped_df, group_cols=c(variables_of_interest), rel_abund_col="relative_abundance")
  
  # Reorder by abundance
  grouped_df <- grouped_df[order(grouped_df[,'relative_abundance'], decreasing=TRUE),]
  
  # Merge rows under a certain abundance threshold to "other"
  under_thresh_df <- grouped_df[grouped_df$relative_abundance < threshold,]
  over_thresh_df <- grouped_df[grouped_df$relative_abundance >= threshold,]
  
  # FIXME: error here if using timeline
  under_thresh_df <- data.frame(t(c('__OTHER__', '__Other__', variables_of_interest, sum(under_thresh_df$mean_abundance), sum(under_thresh_df$relative_abundance)))) 
  names(under_thresh_df) <- c('OTU', target_level, variables_of_interest, "mean_abundance", 'relative_abundance')
  under_thresh_df$relative_abundance <- as.numeric(under_thresh_df$relative_abundance)
  
  final_df <- rbind(over_thresh_df, under_thresh_df)
  final_df$relative_abundance <- as.numeric(final_df$relative_abundance)
  
  # Select columns (return value cannot include grouping variables in group_modify)
  return(data.frame(final_df)[,c("OTU", target_level, 'relative_abundance')])
}

physeq_proccessing <- function(physeq_obj, target_level, elements_to_target, filter_abundance_threshold, 
         prevalence_threshold, plot_abundance_threshold, metadata, variables_of_interest, biosample_col){
  # - Debug -
  # physeq_obj = phyloseq::rarefy_even_depth(physeq, sample.size=10000)
  # target_level = "Genus"
  # elements_to_target = c("Domain","Phylum", "Class", "Order", "Family")
  # filter_abundance_threshold = 0.0001
  # prevalence_threshold = 0.6
  # plot_abundance_threshold = 0.0001
  # metadata = metadata
  # biosample_col = "biosample_name"
  # variables_of_interest = c("water_type","community")
  
  # Check
  if (plot_abundance_threshold<filter_abundance_threshold){
    stop("Plot abundance threshold can not be lower than the abundance threshold. It could lead to bad interpretation of the results. Stopping...")
  }

  # Melt physeq
  print(glue("Melting phyloseq object to long format..."))
  data_phylo_long <- physeq_obj %>%
    tax_glom(taxrank = target_level, bad_empty = c()) %>%   # Aglomerate to taxnomic range
    psmelt()                                                # Melt into long format
  
  # ---- Process biosamples and their reps ----
  # Abundance + prevalence in replicates of each biosample
  biosample_data <- data_phylo_long %>%  
    dplyr::rename("biosample" = biosample_col) %>%
    dplyr::group_by(biosample) %>% 
    dplyr::mutate(sample.id.count = n_distinct(sample.id)) %>%
    dplyr::group_modify(~get.otus.to.keep.by.prevalence(unique(.x$sample.id), .x, filter_abundance_threshold, prevalence_threshold, "sample.id", "Abundance")) %>%
    dplyr::ungroup() 
  
  # Mean abundance in biosample (across replicates)
  biosample_data <- biosample_data %>%
    dplyr::group_by(across(all_of(c("OTU", "biosample", variables_of_interest, elements_to_target, target_level, "sample.id.count")))) %>%
    dplyr::summarise(biosample_mean_abundance = mean(Abundance)) 
  
  # Replace .*__ with an emtpy string in target column
  biosample_data[[target_level]] <- gsub(".*__", '', biosample_data[[target_level]])
  
  
  # ---- Grouping biosamples ----
  # Abundance + prevalence
  group_prevalence <- biosample_data %>%  
    dplyr::group_by(across(all_of(variables_of_interest))) %>%
    dplyr::group_modify(~get.otus.to.keep.by.prevalence(unique(.x$biosample), .x, filter_abundance_threshold, prevalence_threshold, "biosample", "biosample_mean_abundance")) %>%
    dplyr::ungroup() 
  
  # Calculate group relative abundance
  print(glue("Calculate relative abundance by group..."))
  proportions_per_group <- group_prevalence %>%  
    dplyr::group_by(across(all_of(variables_of_interest))) %>%
    dplyr::group_modify(~calculate_props(.x, variables_of_interest, target_level, plot_abundance_threshold, "biosample_mean_abundance"), .keep=TRUE)
  
  # Remove empty ASVs
  proportions_per_group <- subset(proportions_per_group, relative_abundance != 0)
  
  # Ensure sum of relative abundances is 1 in each group
  check_rel_abund_sum(proportions_per_group, group_cols=variables_of_interest, rel_abund_col="relative_abundance")
  
  # ---- Reform the phyloseq object (using biosamples) ----
  # ASVs
  new_otu_df <- biosample_data[,c("OTU", "biosample", "biosample_mean_abundance")] %>% 
    dplyr::mutate(across(c("biosample_mean_abundance"), round)) %>%
    tidyr::pivot_wider(names_from = biosample, values_from = biosample_mean_abundance, values_fill = 0)
  new_otu_df <- as.data.frame(new_otu_df)
  rownames(new_otu_df) <- new_otu_df[,"OTU"]
  new_otu_df[,"OTU"] <- NULL
  new_otu_df <- as.matrix(new_otu_df)
  
  # Taxmat
  new_taxmat <- biosample_data[,c("OTU", elements_to_target, target_level)] %>% dplyr::distinct(across(c("OTU", elements_to_target, target_level)), .keep_all = TRUE)
  new_taxmat <- as.data.frame(new_taxmat)
  rownames(new_taxmat) <- new_taxmat[,"OTU"]
  new_taxmat[,"OTU"] <- NULL
  new_taxmat <- as.matrix(new_taxmat)
  
  # Metadata
  new_metadata <- biosample_data[,c("biosample", variables_of_interest)] %>% dplyr::distinct(biosample, .keep_all = TRUE)
  new_metadata <- as.data.frame(new_metadata)
  rownames(new_metadata) <- new_metadata[,"biosample"]
  
  # Form biosample physeq object
  biosample_physeq = phyloseq(
    otu_table(new_otu_df, taxa_are_rows = TRUE), 
    tax_table(new_taxmat), 
    sample_data(new_metadata),
    phy_tree(physeq_obj)
  )

  return(list(
    proportions_per_group=proportions_per_group, 
    biosample_data=biosample_data,
    biosample_physeq=biosample_physeq
  ))
}

get_barplots <- function(
  proportions, 
  target_level, 
  x="label", y="relative_abundance", 
  form="relative", title="", 
  column_order=NULL, fwrap_formula=NULL, 
  modify_legend = FALSE
  ){
  
  # Colors
  color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] # Creando vector de colores contrastantes
  phylum_colors <- distinctColorPalette(length(levels(as.factor(proportions[[target_level]]))))
  
  # Change color of __Other__ to black
  others_pos <- levels(as.factor(proportions[[target_level]])) == "__Other__"
  phylum_colors[others_pos] <- "black"
  
  # Split names by '_'
  proportions[[target_level]] <- gsub("_", ' ', proportions[[target_level]])
  
  # Config plot
  geom_bar_pos <- NULL
  ylab <- NULL
  if (form=="relative"){
    geom_bar_pos <- "fill"
    ylab <- "Relative Abundance (RA) %"
  } else if (form=="as_is"){
    geom_bar_pos <- "stack"
    ylab <- "Abundance"
  } else {
    stop("Value for variable form not recognized (has to be 'relative' or 'as_is'")
  }
  
  # Main
  p <- ggplot(proportions, aes_string(x = x, y = y, fill = target_level)) + 
    geom_bar(stat = "identity", position = geom_bar_pos) + scale_fill_manual(values = phylum_colors)
  
  # If facet_wrap
  if (!is.null(fwrap_formula)){
    p <- p + facet_wrap(as.formula(fwrap_formula), scales="free_x")
  }

  # If relative put percentages in values
  if (form=="relative"){
    p <- p + scale_y_continuous(labels = scales::percent_format())
  }
  
  # If column_order order columns
  if (!is.null(column_order)){
    p <- p + scale_x_discrete(limits = column_order)
  }
  
  # Some formatting
  p <- p +
    ylab(ylab) + ggtitle(title) + xlab("") + theme_bw() + 
    theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) + 
    theme(axis.text.x=element_text(angle = 90, hjust=0.95, vjust=0.2), legend.key.size = unit(0.5,"line"))
  
  # This removes legend from ggplotly if active but is still there in ggplot
  if (modify_legend){
    p <- p + guides(fill=guide_legend(ncol=1, title=target_level)) # + theme(legend.position="none")
  }
  
  return(p)
  
}


#################################
# ---- Relative abundance  ---- #
#################################
cat("\nCreating barplots...\n")

# target_level <- "Species" # Ex: "Species"
# elements_to_target <- c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species") # Ex: c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species")
filter_abundance_threshold <- 0.0001
prevalence_threshold <- 0.3
plot_abundance_threshold <- 0.0001

# Get proportions
conf = list(
  list("Genus", c("Domain","Phylum", "Class", "Order", "Family")),
  list("Species", c("Domain","Phylum", "Class", "Order", "Family", "Genus")),
  list("Family", c("Domain","Phylum", "Class", "Order"))
)

result_list <- list()
for (i in conf){
  target_level <- i[[1]]
  elements_to_target <- i[[2]]
  out_folder = glue("{output_folder}/barplots/abund-{filter_abundance_threshold}__prev-{prevalence_threshold}__plot-{plot_abundance_threshold}/{target_level}")
  # dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
  prop <- NULL
  
  cat("\n")
  print(glue("====== {target_level} level ======"))
  
  # Get proportions
  prop <- physeq_proccessing(
    physeq_obj = rar_physeq, 
    target_level = target_level,
    elements_to_target = elements_to_target,
    filter_abundance_threshold = filter_abundance_threshold, 
    prevalence_threshold = prevalence_threshold,
    plot_abundance_threshold = plot_abundance_threshold, 
    metadata = metadata, 
    variables_of_interest = variables_of_interest, 
    biosample_col = biosample_col
  )
  
  proportions_per_group <- prop$proportions_per_group
  biosample_data <- prop$biosample_data

  # ---- Grouped barplot ----
  # Prepare labels
  proportions_per_group$label <- do.call(paste, c(proportions_per_group[variables_of_interest], sep="__"))
  article_labels <- list(
    c("saliva__non-crc", "S (non-CRC)"),
    c("saliva__crc", "S (CRC)"),
    c("subgingival-fluid__crc", "GCF (CRC)"),
    c("adenocarcinoma__crc", "Ac (CRC)"),
    c("normal-mucosa__crc", "NM (CRC)"),
    c("faeces__crc", "F (CRC)"),
    c("faeces__non-crc", "F (non-CRC)")
  )
  column_order <- c()
  for (i in article_labels){
    data_label <- i[[1]]
    article_label <- i[[2]]
    proportions_per_group$label <- gsub(data_label, article_label, proportions_per_group$label)
    column_order <- c(column_order, article_label)
  }
  
  # Plot
  barplot <- get_barplots(
    proportions=proportions_per_group,
    target_level=target_level,
    x="label",
    y="relative_abundance",
    form="relative",
    title=glue("{target_level} (Abundance:{filter_abundance_threshold}, Prevalence:{prevalence_threshold}, Plot:{plot_abundance_threshold})"),
    column_order=column_order,
  )
  static_barplot <- barplot + 
    scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) + 
    guides(fill=guide_legend(ncol=1, title=target_level)) +
    ggtitle("") +
    theme(
          axis.text.x=element_text(angle = 50, hjust=0.95, vjust=0.95, size=10),
          axis.text.y=element_text(size=10), 
          axis.title.y=element_text(size=10),
          legend.text=element_text(size=8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x=element_blank()
          # axis.line = element_line(colour = "black")
          ) 

  # # Save plot and tables
  ggsave(glue("{out_folder}/grouped_barplot.png"), plot=static_barplot, dpi="retina", units="px", width=2000, height=2500)

  for (i in unique(proportions_per_group$label)){
    temp <- subset(proportions_per_group, label==i)
    temp$relative_abundance <- temp$relative_abundance*100
    write.table(temp, gsub(" ", "_", glue("{out_folder}/proportions__{i}.tsv")), sep="\t", row.names=FALSE)
  }

  htmlwidgets::saveWidget(
    as_widget(ggplotly(barplot)),
    glue("{out_folder}/groups_barplot.html")
  )

  # Store data in list
  result_list[[target_level]] <- list(
    proportions_per_group=proportions_per_group, 
    biosample_data=biosample_data, 
    biosample_physeq=prop$biosample_physeq
  )
}


#############################################

# ---- Percentage of samples in faeces with each oral genus ----
selected_genus <- c("Actinomyces", "Aggregatibacter", "Alloprevotella", 
                    "Campylobacter", "Capnocytophaga", "Catonella", 
                    "Desulfobulbus", "Dialister", "Filifactor", 
                    "Fusobacterium", "Mogibacterium", "Parvimonas", 
                    "Peptococcus", "Peptostreptococcus", 
                    "Porphyromonas", "Prevotella", "Prevotella_7", "Prevotella_9",
                    "Pseudoramibacter", "Streptococcus", 
                    "Tannerella", "Treponema", "Veillonella"
)
# Filter to faeces and count number of biosamples in each sample.type
df <- subset(result_list$Genus$biosample_data, origin=="faeces")
sample_count <- list(
  "crc"=length(unique(subset(df, sample.type=="crc")$biosample)),
  "noncrc"=length(unique(subset(df, sample.type=="non-crc")$biosample))
)
# For each genus count the proportion of samples which have it in each category
pct_df <- data.frame()
for (g in selected_genus){
  crc_count <- length(subset(df, Genus==g & sample.type=="crc" & biosample_mean_abundance>0)$biosample)
  noncrc_count <- length(subset(df, Genus==g & sample.type=="non-crc" & biosample_mean_abundance>0)$biosample)
  pct_df <- rbind(pct_df, 
                data.frame(genus=g, 
                           crc_count=glue("{crc_count}/{sample_count$crc}"),
                           noncrc_count=glue("{noncrc_count}/{sample_count$noncrc}"),
                           crc_pct=round(100*crc_count/sample_count$crc, 2),
                           noncrc_pct=round(100*noncrc_count/sample_count$noncrc, 2)
                           )
                )
}
write.table(pct_df, glue("{output_folder}/prop_of_faeces_with_oral_genus.tsv"), sep="\t", row.names=FALSE)


# ---- Barplots for patients with Fusobacterium or Parvimonas in faeces ----
# Get patients that have Fuso/Parvi
df <- result_list$Genus$biosample_data
biosample.ids <- unique(subset(df, origin=="faeces" & Genus %in% c("Fusobacterium", "Parvimonas"))$biosample)
subject.ids <- sample_data(subset_samples(rar_physeq, biosample_name %in% biosample.ids))$subject

# Add patients ids to biosample dataframe
df <- merge(df, data.frame(sample_data(rar_physeq)[,c("biosample_name","subject")]), by.x="biosample", by.y="biosample_name")
df <- subset(df, subject %in% subject.ids)

# Make barplot
barplot <- get_barplots(
  proportions=df,
  target_level="Genus",
  x="biosample",
  y="biosample_mean_abundance",
  form="relative",
  title=glue("{target_level} (Abundance:{filter_abundance_threshold}, Prevalence:{prevalence_threshold}, Plot:{plot_abundance_threshold})"),
  column_order=NULL,
  fwrap_formula="~subject"
)

htmlwidgets::saveWidget(
  as_widget(ggplotly(barplot + 
                       theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
                       facet_wrap("~subject", scales="free_x", ncol=10), 
                     height=2000, width=1800
  )),
  glue("{out_folder}/patients_with_fuso_or_parvi.html")
)

# ---- Venn diagram for F,GCF and A in CRC ----
df <- result_list$Genus$proportions_per_group
df <- subset(df, sample.type=="crc")
l <- list(
  "Ac"=unique(subset(df, origin=="adenocarcinoma")$Genus),
  "GCF"=unique(subset(df, origin=="subgingival-fluid")$Genus),
  "Feces"=unique(subset(df, origin=="faeces")$Genus)
)

# Central
intersect(intersect(l$Ac, l$GCF), l$Feces) # shared in ac, gcf and feces

# Intermediate
setdiff(intersect(l$Ac, l$GCF), l$Feces) # only between ac and gcf
setdiff(intersect(l$Ac, l$Feces), l$GCF) # only between ac and feces
setdiff(intersect(l$GCF, l$Feces), l$Ac) # only between gcf and feces

# External
setdiff(l$Ac, c(l$GCF, l$Feces)) # only in ac
setdiff(l$GCF, c(l$Ac, l$Feces)) # only in gcf
setdiff(l$Feces, c(l$GCF, l$Ac)) # only in feces

# Plot
venn_p <- ggvenn(l, c("Ac","GCF", "Feces")) + 
  ggtitle("Genus Venn Diagram in CRC patients") + 
  theme(
    plot.title = element_text(hjust = 0.5, size=20)
  )
venn_p
ggsave(glue("{output_folder}/crc_venn.png"), dpi="retina", units="px")

#############################################

write.csv(result_list$Genus$proportions_per_group, "test.tsv")


###########################################
# ---- Alluvial plot for oral genus  ---- #
###########################################

library(ggalluvial)
require(RColorBrewer)
library(egg)
library(randomcoloR)
library(data.table)
proportions_per_group <- read.csv("test.tsv")
proportions_per_group$X <- NULL

# Subset commonly oral genus
selected_genus <- c("Actinomyces", "Aggregatibacter", "Alloprevotella", 
                    "Campylobacter", "Capnocytophaga", "Catonella", 
                    "Desulfobulbus", "Dialister", "Filifactor", 
                    "Fusobacterium", "Mogibacterium", "Parvimonas", 
                    "Peptococcus", "Peptostreptococcus", 
                    "Porphyromonas", "Prevotella", "Prevotella_7", "Prevotella_9",
                    "Pseudoramibacter", "Streptococcus", 
                    "Tannerella", "Treponema", "Veillonella"
)
# No hay: "Chloroflexi", "Cryptobacterium"
# Es intestinal: "Prevotella_9"
# Problemas: "Eubacterium"
oral_proportions <- subset(proportions_per_group, Genus %in% selected_genus)

# For each genus add a dummy entry in each group if it does not exist there, with RA=0
# That will make the ggalluvial plot be able to connect genus across columns
for (i in unique(oral_proportions$label)){
  for (j in selected_genus){
    temp <- subset(oral_proportions, label==i & Genus==j)
    if (nrow(temp)==0){
      oral_proportions <- rbind(oral_proportions, 
                                data.frame(
                                  origin="dummy",
                                  sample.type="dummy",
                                  OTU="randomid", 
                                  Genus=j,
                                  relative_abundance=0, 
                                  label=i
                                )
      )
    }
  }
}

# Colors
color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] # Creando vector de colores contrastantes
phylum_colors <- distinctColorPalette(length(levels(as.factor(oral_proportions[["Genus"]]))))
column_order <- c("S (non-CRC)", "S (CRC)","GCF (CRC)", "Ac (CRC)", "NM (CRC)", "F (CRC)", "F (non-CRC)")

p <- ggplot(oral_proportions, aes(x=label, y=relative_abundance, fill = Genus)) +
  geom_flow(aes(alluvium = Genus), alpha= 1, #color="gray40",
            curve_type = "xspline",
            width = .5) +
  geom_col(width = .5) +
  scale_x_discrete(limits = column_order, expand = c(0,0)) +
  scale_fill_manual(values = phylum_colors) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) +
  theme_classic() +
  ylab("Relative abundance (RA) %") + ggtitle("") + xlab("") + 
  guides(fill=guide_legend(ncol=1)) +
  theme(
    axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
    axis.text.x=element_text(angle = 50, hjust=0.95, vjust=0.95, size=10), legend.key.size = unit(0.5,"line"),
    axis.text.y=element_text(size=10), 
    axis.title.y=element_text(size=10),
    legend.text=element_text(size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    # axis.ticks.x=element_blank()
    # axis.line = element_line(colour = "black")
  ) 
p
# ggsave("alluvial_1.png", dpi="retina", units="px", width=2500, height=2000)




#############################################


##############################
# ---- Alpha diversity  ---- #
##############################
plot_alpha_diversity <- function(physeq_obj, measures, x, color_col){
  # ---- Plot alpha diversity ----
  alpha_div_plot <- plot_richness(physeq_obj, x=x, color=color_col, measures=measures) + 
    geom_boxplot() + theme_bw() +
    # facet_grid(variable~origin, scale="free") +
    # scale_x_discrete(limits = c("saliva__non-crc", "saliva__crc", "subgingival-fluid__crc", "faeces__non-crc", "faeces__crc", "normal-mucosa__crc", "adenocarcinoma__crc")) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
          axis.text = element_text(size = 12))

  # ---- Wilcoxon-test ----
  richness <- estimate_richness(physeq_obj, measures=measures)
  rownames(richness) <- gsub("\\.", '-', rownames(richness))
  richness <- merge(richness, data.frame(sample_data(physeq_obj)), by=0)
  richness$group <- as.factor(richness$group)
  
  # Run test for each combination of groups
  wilcox_list <- list()
  n = 1
  for (m in measures){
    for (i in sort(unique(richness$group))){
      for (j in sort(unique(richness$group))){
        res <- NULL
        # As long as the group is not the same (error if it is)
        if (i!=j){
          # Run the test between groups and print if significative
          res <- wilcox.test(as.formula(glue("{m} ~ group")), 
                      data=subset(richness, 
                                  group %in% c(i, j)
                      )
          )
          wilcox_list[[n]] <- c(i, j, m, res$p.value)
        } else {
          wilcox_list[[n]] <- c(i, j, m, NA)
        }
        n=n+1
      }
    }  
  }
  
  # Format Wilcox results
  wilcox_df <- as.data.frame(do.call(rbind, wilcox_list))
  names(wilcox_df) <- c("C1","C2","measure","pval")
  wilcox_df$pval <- as.numeric(wilcox_df$pval)
  
  # Plot each wilcoxon test
  significance_p_list <- list()
  for (m in measures){
    df <- subset(wilcox_df, measure==m)[c("C1","C2","pval")]
  
    # Turn to wide and block lower triangle
    df <- as.data.frame(df %>% tidyr::pivot_wider(names_from = C1, values_from = pval))
    df[-1][lower.tri(df[-1])] <- NA
    
    # Turn to long
    df <- df %>% pivot_longer(!C2, names_to = "C1", values_to = "pval")
    
    p <- ggplot(df, aes(C1, C2, fill=pval)) + 
      geom_tile() +
      scale_fill_gradient(low="darkgreen", high="white", na.value = "white", limits=c(0,0.1)) +
      # scale_y_discrete(position = "right") +
      theme_bw() + 
      ylab("") + ggtitle(glue("{m}")) + xlab("") + 
      theme(axis.text.x=element_text(angle = 60, hjust=0.95, vjust=0.95), legend.key.size = unit(0.5,"line"))
    significance_p_list[[m]] <- p
  }
  
  return(list(alpha_div_plot=alpha_div_plot, significance_p_list=significance_p_list))
  
}


measures <- c("Chao1","Shannon","Simpson")

# -- All samples --
ps <- subset_samples(physeq)
sample_data(ps)$group <- paste(
  sample_data(ps)$origin,
  sample_data(ps)$sample.type,
  sep='__'
)
res = plot_alpha_diversity(ps, measures, "group", "origin")
res$alpha_div_plot
res$significance_p_list$Chao1 + res$significance_p_list$Shannon + res$significance_p_list$Simpson + plot_layout(guides = "collect", ncol=3)

# -- Adenocarcinoma samples with location+tnms --
ps <- subset_samples(physeq, origin=="adenocarcinoma")
sample_data(ps)$group <- paste(
  sample_data(ps)$origin,
  paste("loc", sample_data(ps)$simplified_location, sep=":"),
  paste("ln", sample_data(ps)$has_affected_lymph_nodes, sep=":"),
  paste("bt", sample_data(ps)$has_big_tumor, sep=":"),
  paste("mt", sample_data(ps)$has_metastasis, sep=":"),
  sep='_'
)
res = plot_alpha_diversity(ps, measures, "group", "simplified_location")
res$alpha_div_plot
res$significance_p_list$Chao1 + res$significance_p_list$Shannon + res$significance_p_list$Simpson + plot_layout(guides = "collect", ncol=3)


# -- Normal mucosa samples with location --
ps <- subset_samples(physeq, origin %in% c("normal-mucosa"))
res = plot_alpha_diversity(ps, measures, "simplified_location", "simplified_location")
res$alpha_div_plot
res$significance_p_list$Chao1 + res$significance_p_list$Shannon + res$significance_p_list$Simpson + plot_layout(guides = "collect", ncol=3)




# ggsave(glue("{output_folder}/alphadiversity.png"), res$alpha_div_plot, dpi="retina", unit="px", height=2300, width=2800)



#############################
# ---- Beta diversity  ---- #
#############################
# Create beta-diversity plots for every combination of method+distance 
method_vector <- c("MDS","NMDS") #dPCoA
distance_vector <- c("jaccard","bray","jsd","wunifrac")
p_list <- list()
cbbPalette <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
for (i in asplit(crossing(method_vector, distance_vector), 1)){
  # Set up
  m <- i[["method_vector"]]
  d <- i[["distance_vector"]]
  dist <- NULL
  p <- NULL
  title = glue(glue("Method: {m}, Distance: {d}"))
  print(title)
  
  # Ordination
  if (d=="jaccard"){ # important for jaccard to add binary=TRUE
    dist <-ordinate(rar_physeq, m, d, binary = TRUE, formula=~origin+sample.type)
  } else {
    dist <-ordinate(rar_physeq, m, d, formula=~origin+sample.type)
  }
  
  # Plot
  p <- plot_ordination(rar_physeq, dist, 
                       type="sites", color="origin", shape="sample.type", 
                       title=title) + 
    stat_ellipse() + theme_classic() + scale_colour_manual(values=cbbPalette)
  
  p_list[[m]][[d]] <- p
}

# Plot
p_MDS <- p_list$MDS$bray + p_list$MDS$jaccard + p_list$MDS$jsd + p_list$MDS$wunifrac + plot_layout(guides = "collect")
ggsave(glue("{output_folder}/betadiversity_MDS.png"), p_MDS, dpi="retina", unit="px", height=3000, width=3000)

p_NDMS <- p_list$NMDS$bray + p_list$NMDS$jaccard + p_list$NMDS$jsd + p_list$NMDS$wunifrac + plot_layout(guides = "collect")
ggsave(glue("{output_folder}/betadiversity_NMDS.png"), p_NDMS, dpi="retina", unit="px", height=3000, width=3000)



p_mds_jsd <- p_list$MDS$jsd + scale_color_discrete() + theme_bw() + 
  theme(
    axis.text.x=element_text(angle = 0, hjust=0.95, vjust=0.95, size=10),
    axis.text.y=element_text(size=10), 
    axis.title.y=element_text(size=10),
    legend.text=element_text(size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # panel.background = element_blank(),
    # axis.ticks.x=element_blank()
    axis.line = element_line(colour = "black")
  )
ggsave(glue("{output_folder}/betadiversity_MDS_JSD.png"), p_mds_jsd, dpi="retina", unit="px", height=3000, width=3500)

# PcoA
# data_phylo_dist <-ordinate(rar_physeq, "PCoA", "jsd")
# pcoa_jsd = plot_ordination(rar_physeq, data_phylo_dist,
#                            type="sites", color="origin", shape="sample.type",
#                            title="pcoa_jsd") + theme_classic() + stat_ellipse()


#######################
# ---- PERMANOVA ---- #
#######################

# ----------------------------- #
# ------- Origin test --------- #
# ----------------------------- #
# Filter
rar_ps = subset_samples(rar_physeq)

# Matriz de distancias
rar_physeq_dist <- phyloseq::distance(rar_ps, "jsd")
sampledf <- data.frame(sample_data(rar_ps))

# Adonis test
permanova <- adonis2(
  rar_physeq_dist ~ origin + sample.type + age + sex, 
  data = sampledf, permutations = 9999,
  parallel = 60)

permanova

# adonis2(formula = rar_physeq_dist ~ origin + sample.type + age + sex, data = sampledf, permutations = 9999, parallel = 60)
#              Df SumOfSqs      R2       F Pr(>F)    
# origin        4   22.436 0.33965 46.5858 0.0001 ***
# sample.type   1    0.175 0.00266  1.4572 0.1208    
# age           1    0.318 0.00481  2.6370 0.0093 ** 
# sex           1    0.264 0.00400  2.1946 0.0213 *  
# Residual    356   42.862 0.64889                   
# Total       363   66.055 1.00000                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ----------------------------- #
# ---- Adenocarcinoma test ---- #
# ----------------------------- #

# Filter
rar_ps = subset_samples(rar_physeq, sample.type %in% c("crc") & origin %in% c("adenocarcinoma"))

# Matriz de distancias
rar_physeq_dist <- phyloseq::distance(rar_ps, "jsd")
sampledf <- data.frame(sample_data(rar_ps))

# Adonis test
permanova <- adonis2(
  rar_physeq_dist ~ simplified_location + age + sex + has_affected_lymph_nodes + has_big_tumor + has_metastasis, 
  data = sampledf, permutations = 9999,
  parallel = 60)

permanova

# adonis2(formula = rar_physeq_dist ~ simplified_location + age + sex + has_affected_lymph_nodes + has_big_tumor + has_metastasis, data = sampledf, permutations = 9999, parallel = 60)
#                          Df SumOfSqs      R2      F Pr(>F)
# simplified_location       3   0.6629 0.08559 1.6180 0.0081 **
# age                       1   0.3362 0.04340 2.4615 0.0021 **
# sex                       1   0.1993 0.02573 1.4593 0.0791 .
# has_affected_lymph_nodes  1   0.1715 0.02215 1.2559 0.1798
# has_big_tumor             1   0.1853 0.02393 1.3571 0.1267
# has_metastasis            1   0.1808 0.02335 1.3242 0.1343
# Residual                 44   6.0088 0.77585
# Total                    52   7.7447 1.00000
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# --------------------- #
# ---- Faeces test ---- #
# --------------------- #

# Filter
rar_ps = subset_samples(rar_physeq, 
                        sample.type %in% c("non-crc","crc") & 
                        origin %in% c("faeces")
                        )

# Matriz de distancias
rar_physeq_dist <- phyloseq::distance(rar_ps, "jsd")
sampledf <- data.frame(sample_data(rar_ps))

# Adonis test
permanova <- adonis2(
  rar_physeq_dist ~ sample.type + age + sex, 
  data = sampledf, permutations = 9999,
  parallel = 60)

permanova

# adonis2(formula = rar_physeq_dist ~ sample.type + age + sex, data = sampledf, permutations = 9999, parallel = 60)
#              Df SumOfSqs      R2      F Pr(>F)   
# sample.type   1   0.2207 0.01342 1.6170 0.0075 **
# age           1   0.2088 0.01270 1.5300 0.0149 * 
# sex           1   0.1806 0.01098 1.3232 0.0703 . 
# Residual    116  15.8304 0.96289                 
# Total       119  16.4404 1.00000                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# -------------------------------- #
# ---- Faeces test (only crc) ---- #
# -------------------------------- #
# Filter
rar_ps = subset_samples(rar_physeq, 
                        sample.type %in% c("crc") & 
                          origin %in% c("faeces")
)

# Matriz de distancias
rar_physeq_dist <- phyloseq::distance(rar_ps, "jsd")
sampledf <- data.frame(sample_data(rar_ps))

# Adonis test
permanova <- adonis2(
  rar_physeq_dist ~ age + sex + simplified_location + has_affected_lymph_nodes + has_big_tumor + has_metastasis, 
  data = sampledf, permutations = 9999,
  parallel = 60)

permanova

# adonis2(formula = rar_physeq_dist ~ age + sex + simplified_location + has_affected_lymph_nodes + has_big_tumor + has_metastasis, data = sampledf, permutations = 9999, parallel = 60)
#                          Df SumOfSqs      R2      F Pr(>F)   
# age                       1   0.1969 0.01549 1.4167 0.0343 * 
# sex                       1   0.1974 0.01553 1.4201 0.0310 * 
# simplified_location       4   0.7485 0.05888 1.3462 0.0031 **
# has_affected_lymph_nodes  1   0.1454 0.01144 1.0458 0.3715   
# has_big_tumor             1   0.1290 0.01015 0.9279 0.6127   
# has_metastasis            1   0.1740 0.01369 1.2517 0.1059   
# Residual                 80  11.1201 0.87483                 
# Total                    89  12.7112 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# -------------------------------- #
# ---- Subgingival-fluid test ---- #
# -------------------------------- #

# Filter
rar_ps = subset_samples(rar_physeq, 
                        sample.type %in% c("crc") & 
                        origin %in% c("subgingival-fluid") & 
                        infectious_bone_pathologies %in% c("Y","N")
                        )

# Matriz de distancias
rar_physeq_dist <- phyloseq::distance(rar_ps, "jsd")
sampledf <- data.frame(sample_data(rar_ps))

# Adonis test
permanova <- adonis2(
  rar_physeq_dist ~ age + sex + silness_loe_index + absent_teeth + caries_index + infectious_bone_pathologies, 
  data = sampledf, permutations = 9999,
  parallel = 60)

permanova

# adonis2(formula = rar_physeq_dist ~ age + sex + silness_loe_index + absent_teeth + caries_index + infectious_bone_pathologies, data = sampledf, permutations = 9999, parallel = 60)
#                             Df SumOfSqs      R2      F Pr(>F)
# age                          1  0.13061 0.11500 0.9395 0.5796
# sex                          1  0.14974 0.13184 1.0770 0.4309
# silness_loe_index            1  0.18541 0.16325 1.3336 0.2161
# absent_teeth                 1  0.14380 0.12661 1.0343 0.4582
# caries_index                 1  0.19594 0.17252 1.4094 0.1827
# infectious_bone_pathologies  1  0.19121 0.16835 1.3753 0.1893
# Residual                     1  0.13903 0.12241              
# Total                        7  1.13574 1.00000  



#########################
# ---- Saliva test ---- #
#########################
# Filter
rar_ps = subset_samples(rar_physeq, 
                        sample.type %in% c("crc","non-crc") & 
                          origin %in% c("saliva")
)

# Matriz de distancias
rar_physeq_dist <- phyloseq::distance(rar_ps, "jsd")
sampledf <- data.frame(sample_data(rar_ps))

# Adonis test
permanova <- adonis2(
  rar_physeq_dist ~ sample.type + age + sex, 
  data = sampledf, permutations = 9999,
  parallel = 60)

permanova

# adonis2(formula = rar_physeq_dist ~ sample.type + age + sex, data = sampledf, permutations = 9999, parallel = 60)
#              Df SumOfSqs      R2      F Pr(>F)
# sample.type   1   0.1014 0.01064 1.2626 0.1655
# age           1   0.1059 0.01110 1.3181 0.1377
# sex           1   0.0899 0.00943 1.1194 0.2744
# Residual    115   9.2403 0.96883              
# Total       118   9.5376 1.00000 



####################################
# ---- Differential abundance ---- #
####################################

run_ancombc_1vs1 <- function(physeq_obj, target_level, target_variables, prv_cut=0.3, qval_min=0.1){
  # IMPORTANT: ANCOM-BC v2.0.1 was used. Its output format varies depending on version
  if (target_level=="ASV"){
    ancom_target <- NULL
  } else {
    ancom_target <- target_level
  }
  
  # Run ANCOM-BC
  out <- ancombc(data=physeq_obj, tax_level=ancom_target, prv_cut=prv_cut, 
                            formula = paste(target_variables, collapse="+"), n_cl=60,
                            p_adj_method="holm", alpha = 0.05) 
  # Get results
  res = out$res
  
  # To long
  lfc_long <- res$lfc[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "lfc")
  diff_abn_long <- res$diff_abn[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "diff_abn")
  se_long <- res$se[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "se")
  qval_long <- res$q_val[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "qval")
  
  # Merge long
  merged <- Reduce(function(x, y) merge(x, y, by=c("taxon","var")), 
         list(lfc_long, diff_abn_long, se_long, qval_long))
  
  # Get taxons where any of the tests achived a certain qval
  selected_taxons <- res$q_val[rowSums(res$q_val[,c(-1,-2),drop=F]<qval_min)!=0,]$taxon
  merged <- subset(merged, taxon %in% selected_taxons) 
  merged
  
  # Change column name
  names(merged)[1] <- "taxon_id"

  # Add names if ASV was selected as target_level
  if (target_level=="ASV"){
    merged$taxon_id <- data.frame(tax_table(physeq_obj))[merged$taxon_id,][,"Species"]
  }
  
  # Format taxon_id to remove things like "g__"
  # merged$taxon_id <- gsub(".*__", '', merged$taxon_id)
  
  # Format results for plotting
  df_fig <- merged %>%
    dplyr::arrange(desc(lfc)) %>%
    dplyr::mutate(direction = ifelse(lfc > 0, "Positive LFC", "Negative LFC"))
  df_fig$direction <- factor(df_fig$direction, levels = c("Positive LFC", "Negative LFC"))
  
  # Significance label
  # - Corrected p-value (qval) under 0.001 will have "***"
  # - Corrected p-value (qval) under 0.01 will have "**"
  # - Corrected p-value (qval) under 0.05 will have "*"
  df_fig$qval.txt <- cut(df_fig$qval, c(-Inf,0.001,0.01,0.05,Inf), labels=c('***','**','*','-'))
  
  # Orientation for significance level
  df_fig$orientation <- 0
  df_fig$orientation[df_fig$lfc > 0] <- 1
  df_fig$orientation[df_fig$lfc < 0] <- -1  
  
  # Add target level
  df_fig$target_level <- target_level
  
  return(df_fig)
}


# - Use original physeq object without prevalence/abundance filters
#   (but controls applied and contamination removed)
# - Select faeces samples for comparison
# ps <- subset_samples(physeq, 
#                      sample.type %in% c("crc","non-crc") & 
#                        origin %in% c("faeces") 
#                      
# )
ps <- subset_samples(physeq, origin %in% c("faeces"))
target_variables <- c("sample.type")

# sample_data(ps)$test <- mapply(
#   paste0,
#   as.character(get_variable(ps, "sample.type")),
#   as.character(get_variable(ps, "sex")),
#   as.character(get_variable(ps, "age_group")),
#   collapse = "_"
# )

# Add previous level
temp <- as.data.frame(tax_table(ps))
temp$Species <- paste(temp[,"Genus"], temp[,"Species"], sep=';')
temp$Genus <- paste(temp[,"Family"], temp[,"Genus"], sep=';')
tax_table(ps) <- tax_table(as.matrix(temp))

# Group by sample.type, origin and subject
sample_data(ps)$group <- mapply(
  paste0,
  as.character(get_variable(ps, "sample.type")),
  as.character(get_variable(ps, "origin")),
  as.character(get_variable(ps, "subject")),
  as.character(get_variable(ps, target_variables)),
  collapse = "_"
)
derep_ps <- phyloseq::merge_samples(ps, "group", fun=median)

# The function merge_samples transforms string factors to numeric
# so we need to replace the metadata of derep_ps
# based on: https://github.com/joey711/phyloseq/issues/608
df <- sample_data(ps) #Extract sample data dataframe 
rmreps <- subset(df, !duplicated(group)) # Remove rows with duplicate values by the factor "SampleType" (it will retain 1 of each SampleType value -if associated variables are unique it will pick one)
sorted_rmreps <- rmreps[order(rmreps$group),] #Sort dataframe in order of "SampleType" to match merged phyloseq object
rownames(sorted_rmreps) <- sorted_rmreps$group #Rename rows of dataframe to match the sample_names in the merged phyloseq object
sample_data(derep_ps) <- sorted_rmreps # Replace metadata in merged phyloseq object with populated metadata file

# table(sample_data(derep_ps)$simplified_location)
# table(sample_data(derep_ps)$has_big_tumor)
# table(sample_data(derep_ps)$has_metastasis)
# table(sample_data(derep_ps)$has_affected_lymph_nodes)

# Run on each level
proc_family_ancombc <- NULL
proc_genus_ancombc <- NULL
proc_species_ancombc <- NULL
proc_asv_ancombc <- NULL

proc_family_ancombc <- run_ancombc_1vs1(derep_ps, "Family", target_variables)
proc_genus_ancombc <- run_ancombc_1vs1(derep_ps, "Genus", target_variables)
proc_species_ancombc <- run_ancombc_1vs1(derep_ps, "Species", target_variables)
proc_asv_ancombc <- run_ancombc_1vs1(derep_ps, "ASV", target_variables)

# Concatenate results
df_fig <- rbind(proc_family_ancombc, proc_genus_ancombc, proc_species_ancombc, proc_asv_ancombc)
df_fig$target_level <- factor(df_fig$target_level,      # Reordering group factor levels
                               levels = c("Family","Genus","Species","ASV")
                              )



### Modo multi ###
get_scale_alpha_values <- function(vec){
  l <- length(unique(vec))
  alpha_vec <- rep(1, l)
  if ("-" %in% vec){
    alpha_vec[length(alpha_vec)] <- 0.4
  }
  return(alpha_vec)
}
p <- NULL
p <- ggplot(data = subset(df_fig), aes(x = taxon_id, y = lfc, fill=var, alpha=qval.txt)) + 
  scale_alpha_manual(values=get_scale_alpha_values(df_fig$qval.txt)) +
  geom_bar(stat = "identity",
           position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width=1,
                position = position_dodge2(width = 0.9, preserve = "single"), 
                color = "gray36") + 
  geom_text(aes(label = qval.txt, y=lfc+orientation*(se+0.22)), vjust = 0.7, color = "gray36",
            position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_vline(xintercept = (1:length(df_fig$taxon_id)-1)+0.5, alpha=0.5, linetype="dotted") +
  geom_hline(yintercept = 0, alpha=0.5) +
  facet_wrap(~target_level, scale="free_y", ncol=1) +
  labs(x = NULL, y = "Log fold change", title = "") + 
  # scale_fill_discrete(name = NULL) +
  # scale_color_discrete(name = NULL) +
  # scale_x_discrete(limits = unique(df_fig$taxon_id)) +
  coord_flip() +
  theme_bw() + 
  theme(
    axis.text.x=element_text(angle = 0, hjust=0.95, vjust=0.95, size=10),
    axis.text.y=element_text(size=10), 
    axis.title.y=element_text(size=10, face = "italic"),
    legend.text=element_text(size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks.x=element_blank()
    # axis.line = element_line(colour = "black")
  )
p
ggsave(glue("ancombc.png"), 
       dpi="retina", unit="px", height=4000, width=5000)



# df_fig <- subset(df_fig, !taxon_id %in% c("Incertae_Sedis","uncultured_bacterium_209","Incertae_Sedis_NA"))

# # Vertical
# p <- ggplot(data = df_fig, aes(x = taxon_id, y = lfc, fill=direction)) + 
#   geom_bar(stat = "identity", width = 0.7, 
#            position = position_dodge(width = 0.4), color="gray36") +
#   geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), 
#                 width = 0.2, position = position_dodge(0.05), color = "gray36") + 
#   geom_text(aes(label = qval.txt, y=lfc+orientation*(se+0.15)), vjust = 0.7, color = "gray36") +
#   coord_flip() +
#   facet_wrap(~target_level, scale="free_y", nrow=4) + 
#   labs(x = NULL, y = "Log fold change", title = "") + 
#   scale_fill_discrete(name = NULL) +
#   scale_color_discrete(name = NULL) +
#   # scale_x_discrete(limits = df_fig$taxon_id) +
#   theme_bw() + 
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.minor.y = element_blank(),
#         axis.title = element_text(color = "black", size = 12),
#         axis.text.x = element_text(angle = 0, color = "black"),
#         axis.text.y = element_text(angle = 0, color = "black", face = "italic"),
#         axis.text = element_text(size = 12))
# p
# ggsave(glue("{output_folder}/ancombc.png"), dpi="retina", unit="px", height=2000, width=2500)

















###################################
# Exclusion zone
###################################

# # ---- Biosamples barplot ----
# # Reorder facets
# biosample_data$group <- do.call(paste, c(biosample_data[variables_of_interest], sep="__"))
# biosample_data$group <- factor(biosample_data$group,      # Reordering group factor levels
#                                levels = c("saliva__non-crc", "saliva__crc", "subgingival-fluid__crc", "adenocarcinoma__crc", "normal-mucosa__crc", "faeces__crc", "faeces__non-crc")
#                                )
# 
# barplot <- get_barplots(
#   proportions=subset(biosample_data, origin %in% c("faeces")),
#   target_level=target_level,
#   x="biosample",
#   y="biosample_mean_abundance",
#   form="relative",
#   title=glue("{target_level} (Abundance:{filter_abundance_threshold}, Prevalence:{prevalence_threshold}, Plot:{plot_abundance_threshold})"),
#   column_order=NULL,
#   fwrap_formula="~group"
# )
# ggplotly(barplot + theme(axis.text.x=element_blank(),
#                          axis.ticks.x=element_blank()))
# 
#   htmlwidgets::saveWidget(
#     as_widget(ggplotly(barplot)),
#     glue("{out_folder}/groups_barplot.html")
#   )

####################################33
# # Separate dataset
# crc_faeces_physeq <- subset_samples(physeq, origin=="faeces" & sample.type=="crc")
# noncrc_faeces_physeq <- subset_samples(physeq, origin=="faeces" & sample.type=="non-crc")
# 
# # Select
# n.selected.ids <- length(sample_data(noncrc_faeces_physeq)$sample.id)
# selected.ids <- sample(sample_data(crc_faeces_physeq)$sample.id, n.selected.ids)
# 
# # Merge
# ps <- merge_phyloseq(
#   subset_samples(crc_faeces_physeq, sample.id %in% selected.ids),
#   noncrc_faeces_physeq
# )
# 
# table(data.frame(sample_data(ps))[,c("sample.type","sex")])
# table(data.frame(sample_data(ps))[,c("sample.type","age_group")])
# 

######################################
# # Rename columns (WARNING: only appropiate for one to one comparisons)
# names(res$lfc) <- c("taxon_id","lfc_intercept", "lfc")
# names(res$diff_abn) <- c("taxon_id", "diff_abn_intercept", "diff_abn")
# names(res$se) <- c("taxon_id","se_intercept", "se")
# names(res$q_val) <- c("taxon_id","qval_intercept", "qval")
# 
# # Merge results into one dataframe by taxon_id
# merged <- Reduce(function(x, y) merge(x, y, by.x="taxon_id", by.y="taxon_id"),
#                  list(res$lfc, res$diff_abn, res$se, res$q_val))
# merged <- merged[,c("taxon_id","lfc","diff_abn","se","qval")]
# 
# # Select the significant ones
# merged <- subset(merged, qval<0.1) #diff_abn==TRUE

# df_fig <- merged %>% 
#   dplyr::arrange(desc(lfc)) %>%
# dplyr::mutate(direct = ifelse(lfc > 0, "Positive LFC", "Negative LFC"))
# df_fig$taxon_id <- factor(df_fig$taxon_id, levels = df_fig$taxon_id)
# df_fig$direct <- factor(df_fig$direct, levels = c("Positive LFC", "Negative LFC"))

#####################################
# # Horizontal
# p <- ggplot(data = df_fig, aes(x = taxon_id, y = lfc, fill=direct)) + 
#   geom_bar(stat = "identity", width = 0.7, 
#            position = position_dodge(width = 0.4), color="gray36") +
#   geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), 
#                 width = 0.2, position = position_dodge(0.05), color = "gray36") + 
#   geom_text(aes(label = qval.txt, y=lfc+orientation*(se+0.1)), color = "gray36") +
#   # coord_flip() +
#   facet_wrap(~target_level, scale="free_x", nrow=1) + 
#   labs(x = NULL, y = "Log fold change", title = "") + 
#   scale_fill_discrete(name = NULL) +
#   scale_color_discrete(name = NULL) +
#   # scale_x_discrete(limits = df_fig$taxon_id) +
#   theme_bw() + 
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.minor.y = element_blank(),
#         axis.title = element_text(color = "black", size = 12),
#         axis.text.y = element_text(angle = 0, color = "black"),
#         axis.text.x = element_text(angle = 90, color = "black", face = "italic", hjust=0.95, vjust=0.2),
#         axis.text = element_text(size = 12))
# p
# 















library(mia)

ps <- subset_samples(rar_physeq, origin %in% c("normal-mucosa"))
tse <- makeTreeSEFromPhyloseq(ps)
tse <- agglomerateByRank(tse, rank = "Genus", agglomerateTree=TRUE)


########################################################################

tse_dmn <- mia::runDMN(tse, name = "DMN", k = 1:7)

library(miaViz)
plotDMNFit(tse_dmn, type = "laplace")
getBestDMNFit(tse_dmn, type = "laplace")



DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn))
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn)))
head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn)))


prob <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn))
# Add column names
colnames(prob) <- c("comp1", "comp2")

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec <- colnames(prob)[max.col(prob,ties.method = "first")]


# Does clr transformation. Pseudocount is added, because data contains zeros.
tse <- transformCounts(tse, method = "relabundance", pseudocount = 1)
tse <- transformCounts(tse, "relabundance", method = "clr")

library(scater)

# Does principal coordinate analysis
df <- calculateMDS(tse, exprs_values = "clr", method = "euclidean")

# Creates a data frame from principal coordinates
euclidean_pcoa_df <- data.frame(pcoa1 = df[,1], 
                                pcoa2 = df[,2])



# Creates a data frame that contains principal coordinates and DMM information
euclidean_dmm_pcoa_df <- cbind(euclidean_pcoa_df,
                               dmm_component = vec)

euclidean_dmm_pcoa_df <- merge(euclidean_dmm_pcoa_df, metadata, by=0)

# Creates a plot
euclidean_dmm_plot <- ggplot(data = euclidean_dmm_pcoa_df, 
                             aes(x=pcoa1, y=pcoa2,
                                 color=dmm_component, shape=origin)) +
  geom_point() +
  facet_wrap(~origin)+
  labs(x = "Coordinate 1",
       y = "Coordinate 2",
       title = "PCoA with Aitchison distances") +  
  theme(title = element_text(size = 12)) # makes titles smaller

euclidean_dmm_plot


########################################################################333


library(bluster)
library(patchwork) # For arranging several plots as a grid
library(scater)

tse <- transformCounts(tse, method = "rclr")

# Performing and storing UMAP
tse <- runUMAP(tse, name="UMAP", exprs_values="rclr")

k <- c(2,3,5,11)
ClustAndPlot <- function(x) {
  # Creating the graph and running the short random walks algorithm  
  graph_clusters <- clusterRows(t(assays(tse)$rclr), NNGraphParam(k=x))
  
  # Results of the clustering as a color for each sample
  plotUMAP(tse, colour_by = I(graph_clusters), shape_by="origin") +
    labs(title = paste0("k = ", x))
}

# Applying the function for different k values
plots <- lapply(k,ClustAndPlot)

# Displaying plots in a grid
(plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]])

ggplotly(plots[[1]])
