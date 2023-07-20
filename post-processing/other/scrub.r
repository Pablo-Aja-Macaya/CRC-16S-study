# test in scrub_env: mamba install 'r-base>=3.6.3' bioconductor-phyloseq r-torch r-glmnet r-tidyverse r-devtools + devtools::install_github("shenhav-and-korem-labs/SCRuB")


library(SCRuB)
library(glue)
library(phyloseq)

feature_table_file <- as.character("/home/usuario/Proyectos/CRC-16S-study/data/feature-table.json")
taxmat_file <- as.character("/home/usuario/Proyectos/CRC-16S-study/data/taxonomy.tsv")
metadata_file <- as.character("/home/usuario/Proyectos/CRC-16S-study/data/non-ffpe_metadata.tsv")


data <- import_biom(feature_table_file)
data <- as.matrix(as.data.frame(t(data)))

metadata <- read.csv(metadata_file, sep="\t")
rownames(metadata) <- metadata$sample.id

metadata$is_control <- metadata$sample.type=="control"
metadata$sample_type <- metadata$sample.type
metadata <- metadata[,c("is_control","sample_type")] # ERROR if there are more columns than these

# Remove empty ASVs
# data <- data[,colSums(data)!=0]


# Check if there are columns not in otu_df
condition <- rownames(metadata) %in% rownames(data)
sample.ids.not.in.otus <- metadata[!condition,]$sample.id
metadata <- metadata[condition,]

# Check in the other way
condition <- rownames(data) %in% rownames(metadata)
sample.ids.not.in.md <- rownames(data[!condition,])
data <- data[condition,]

# order
metadata <- metadata[ order(match(rownames(metadata), rownames(data))), ]

# run scrub
scr_out <- SCRuB(data, metadata)

# see contaminants
threshold <- 0.01
contaminants<-colnames( scr_out$decontaminated_samples)[ 
                  which( scr_out$inner_iterations$`control`$gamma > threshold ) ]

taxonomy <- read.csv(taxmat_file, sep="\t")
res<-subset(taxonomy, Feature.ID %in% contaminants)$Taxon
res

# plot relative abundance of contaminants
temp <- data.frame(scr_out$inner_iterations$`control`$gamma)
temp <- merge(temp, taxonomy, by.x=0, by.y="Feature.ID")
temp <- subset(temp, scr_out.inner_iterations.control.gamma>0)
plot(temp$scr_out.inner_iterations.control.gamma)
