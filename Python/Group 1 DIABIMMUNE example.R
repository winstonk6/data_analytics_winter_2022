# Set working directory
setwd("C:/OneDrive/OneDrive - Hunter - CUNY/Folder Name")

# Load libraries
p1 <- c("tidyverse", "vegan", "BiocManager")
p2 <- c("phyloseq", "ANCOMBC", "DESeq2", "ComplexHeatmap")
load_package <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    ifelse(p %in% p1, 
           install.packages(p, repos = "http://cran.us.r-project.org/"), 
           BiocManager::install(p))
  }
  library(p, character.only = TRUE, quietly = TRUE)
}
invisible(lapply(c(p1,p2), load_package))
rm(p1, p2, load_package)

# Load 16S counts
data_16s <- read.table("https://diabimmune.broadinstitute.org/diabimmune/uploads/attachments/69/diabimmune_karelia_16s_otu_table.txt",
                       sep = "\t", header = T, comment.char = "", check.names = FALSE)

# Filter out rows that aggregate other rows
data_16s_all_taxa <- data_16s$sample # Get taxa names
otu <- c() # New vector
for (tax in data_16s_all_taxa) {
  if (length(str_split(tax, pattern="\\|", simplify = TRUE)) == 8) {
    otu <- append(otu, tax)
  }
}
data_16s_otu <- data.frame(otu)
data_16s <- left_join(data_16s_otu, data_16s, by=c("otu" = "sample"))
rownames(data_16s) <- otu
data_16s <- data_16s[-1]
rm(data_16s_otu, tax, data_16s_all_taxa)


# Load metadata
load(url("http://diverge.hunter.cuny.edu/~weigang/datasets-bio47120/DIABIMMUNE_Karelia_metadata.RData"))

# Remove rows from metadata where the SampleID is not in the OTU table
SampleID <- colnames(data_16s) # The otu table's columns are the SampleIDs
samples <- data.frame(SampleID)
metadata <- left_join(samples, metadata, by=c("SampleID" = "SampleID")) # Left join to only get rows where the SampleID is in both tables
rownames(metadata) <- metadata$SampleID  # Rename the row indices to be the SampleID
metadata <- metadata[-grep('SampleID', colnames(metadata))] # Remove SampleID column
rm(SampleID, samples)


# Make taxonomy table
kingdom <- c()
phylum <- c()
class <- c()
order <- c()
family <- c()
genus <- c()
species <- c()
otu_num <- c()

for (o in otu) {
  otu_split <- str_split(o, pattern="\\|", simplify = TRUE)
  kingdom <- append(kingdom, str_replace(otu_split[1], "k__", ""))
  phylum <- append(phylum, str_replace(otu_split[2], "p__", ""))
  class <- append(class, str_replace(otu_split[3], "c__", ""))
  order <- append(order, str_replace(otu_split[4], "o__", ""))
  family <- append(family, str_replace(otu_split[5], "f__", ""))
  genus <- append(genus, str_replace(otu_split[6], "g__", ""))
  species <- append(species, str_replace(otu_split[7], "s__", ""))
  otu_num <- append(otu_num, otu_split[8])
}
taxonomy <- data.frame(
  Kingdom = kingdom,
  Phylum = phylum,
  Class = class,
  Order = order,
  Family = family,
  Genus = genus,
  Species = species,
  OTU = otu_num,
  row.names = otu
)
rm(kingdom, phylum, class, order, family, genus, species, otu_num, o, otu, otu_split)


# Build phyloseq object
OTU = otu_table(as.matrix(data_16s), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxonomy))
SAMPLE <- sample_data(metadata)
ps <- phyloseq(OTU, TAX, SAMPLE)


# Rarefying: normalizing the library size (the counts of sequence reads) for samples in order to compare them.
# This minimizes the bias from varying sequencing depth.
set.seed(111)
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=200, replace=F)


# Alpha diversity metrics assess the species diversity within the ecosystems, telling you how diverse a community is.
plot_richness(ps.rarefied, x="allergy_peanut", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  geom_jitter(shape = 1, color="gray") +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))




# Beta diversity
# Assesses dissimilarity between ecosystem, telling you to what extent individuals are different.

# PCoA (Principal coordinates analysis) reduces the dimensionality of the data.
# Each point represents 1 sample, and distance between points indicates similarity in microbiome composition.
dist = phyloseq::distance(ps.rarefied, method="bray")
ordination = ordinate(ps.rarefied, method="PCoA", distance=dist)
plot_ordination(ps.rarefied, ordination, color="country") + 
  theme_classic() +
  facet_wrap(~gender) +
  theme(strip.background = element_blank())


# Does consumption of hydrosylated formula reduce allergies (including milk) in children of varying ages and backgrounds?

allergies = metadata %>% 
  mutate(Allergies = allergy_milk | allergy_egg | allergy_peanut | allergy_dustmite | allergy_cat | allergy_dog | allergy_birch | allergy_timothy) %>% 
  select(Allergies, hydrosylated_formula) %>% 
  filter(!is.na(Allergies))

allergies.table <- table(allergies) # Contingency table for observed allergies

# Chi squared test
chi.allergies <- chisq.test(allergies.table)
chi.allergies$expected

chi.allergies

# Tidy data frame for actual counts
totals.obs <- as.data.frame(allergies.table) %>% 
  mutate(ob.exp = "observed") # Add categorical column

# Tidy data frame for expected counts
totals.exp <- as.data.frame(chi.allergies$expected) %>% 
  mutate(Allergies = rownames(chi.allergies$expected)) %>% 
  pivot_longer(1:2, names_to = "hydrosylated_formula", values_to = "Freq") %>% 
  mutate(ob.exp = "expected")

totals.all <- bind_rows(totals.obs, totals.exp)

# Plot
totals.all %>% ggplot(aes(x = ob.exp, y = Freq, fill = hydrosylated_formula)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Allergies) +
  theme_bw()
