#general functions

parse_gff3 <- function(gff_file) {
  # Read the GFF3 file
  gff_data <- readLines(gff_file)
  # Extract locus_tag and GO terms
  parsed_data <- lapply(gff_data, function(line) {
    if (!startsWith(line, "#")) {
      cols <- unlist(strsplit(line, "\t"))
      if (length(cols) >= 9) {
        attributes <- cols[9]
        locus_tag <- ifelse(grepl("locus_tag=", attributes),
                            sub(".*locus_tag=([^;]+);.*", "\\1", attributes),
                            NA)
        go_terms <- ifelse(grepl("Ontology_term=", attributes),
                           sub(".*Ontology_term=([^;]+);.*", "\\1", attributes),
                           NA)
        return(data.frame(locus_tag = locus_tag, go_terms = go_terms))
      }
    }
    return(NULL)
  })
  # Combine and clean the data
  gff_df <- do.call(rbind, parsed_data) %>%
    filter(!is.na(go_terms)) %>%
    separate_rows(go_terms, sep = ",")  # Split GO terms into separate rows
  return(gff_df)
}

perform_go_enrichment <- function(gff_df, gene_list) {
  # Background GO term frequencies
  all_go_counts <- gff_df %>%
    group_by(go_terms) %>%
    summarise(background_count = n(), .groups = "drop")
  # Selected gene GO term frequencies
  selected_go_counts <- gff_df %>%
    filter(locus_tag %in% gene_list) %>%
    group_by(go_terms) %>%
    summarise(selected_count = n(), .groups = "drop")
  # Merge counts
  enrichment_df <- left_join(all_go_counts, selected_go_counts, by = "go_terms") %>%
    replace_na(list(selected_count = 0)) %>%
    mutate(
      not_selected = sum(selected_count) - selected_count,
      not_in_background = sum(background_count) - background_count
    )
  # Ensure all values are nonnegative
  enrichment_df <- enrichment_df %>%
    filter(
      selected_count >= 0 & not_selected >= 0 & 
      background_count >= 0 & not_in_background >= 0
    )
  # Perform Fisher's exact test row by row
  enrichment_df <- enrichment_df %>%
    rowwise() %>%
    mutate(
      p_value = tryCatch(
        fisher.test(
          matrix(
            c(selected_count, not_selected, background_count, not_in_background), 
            nrow = 2
          ), 
          alternative = "greater"
        )$p.value,
        error = function(e) NA  # Handle errors
      )
    ) %>%
    ungroup()
  # Adjust p-values for multiple testing
  enrichment_df <- enrichment_df %>%
    mutate(adjusted_p_value = p.adjust(p_value, method = "fdr")) %>%
    arrange(adjusted_p_value)
  return(enrichment_df)
}




#####################hflu


dev.new()

gff_file <- "hflu.gff3"  # Replace with your GFF file

gff_df <- parse_gff3(gff_file)

go_genes <- gff_df %>% 
  filter(grepl("GO:0003735", go_terms)) %>%
  pull(locus_tag)


# Read and preprocess data
tb.hflu <- read.delim("hflu_df.tsv", header = TRUE, sep = '\t')

print(length(which(tb.hflu[,3] >= 10)))
length(which(tb.hflu[,3] >= 10)) / nrow(tb.hflu)

goList = match(go_genes,tb.hflu[,1])

tb.hflu <- data.frame(Gene = tb.hflu[, 2], Expression = tb.hflu[, 3])  # Use column 2 for gene names
tb.hflu$LogExpression <- log(tb.hflu$Expression + 1, 10)  # Log-transform the expression data

# Specify the genes to highlight
highlight_genes <- c("ompA","pal","hfq","luxS","sodA","znuA","hxuB")
tb.hflu$Highlight <- ifelse(tb.hflu$Gene %in% highlight_genes, "Highlighted", "Normal")
tb.hflu$Highlight[goList] = "Ribosomal"

# Compute quartiles
quartiles <- quantile(tb.hflu$LogExpression, probs = seq(0, 1, by = 0.25), na.rm = TRUE)

# Create the jitter plot
jitter_plot <- ggplot(tb.hflu, aes(x = 1, y = LogExpression)) +
  # Add normal points (gray)
  geom_jitter(data = subset(tb.hflu, Highlight == "Normal"), color = "gray", width = 0.2, size = 3) +
  # Add ribosomal genes (blue)
  geom_jitter(data = subset(tb.hflu, Highlight == "Ribosomal"), color = "blue", width = 0.2, size = 3) +
  # Add highlighted genes (red)
  geom_jitter(data = subset(tb.hflu, Highlight == "Highlighted"), color = "red", width = 0.2, size = 3) +
  # Add quartile division lines
  geom_hline(yintercept = quartiles[-c(1, 5)], linetype = "dashed", color = "blue") +
  # Add labels for highlighted genes
  geom_text(
    data = subset(tb.hflu, Highlight == "Highlighted"),
    aes(label = Gene),
    hjust = -0.2, vjust = 0, size = 4, color = "red"
  ) +
  theme_minimal() +
  labs(x = "Samples", y = "Log10(Expression)", title = "hflu with Highlighted Genes and Labels") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank())

# Add marginal histogram
ggExtra::ggMarginal(
  jitter_plot,
  type = "histogram",
  margins = "y",          # Add histogram to the y-axis
  binwidth = 0.1,         # Adjust bin width
  fill = "gray",          # Neutral fill for histogram
  color = "black"         # Border color for histogram
)





### Go enrichment


tb.hflu = read.delim("hflu_df.tsv",header=T,sep='\t')

gene_list = tb.hflu[1:200,1] ## define gene list of interest

gff_file <- "hflu.gff3"  # Replace with your GFF file

gff_df <- parse_gff3(gff_file)


enrichment_results <- perform_go_enrichment(gff_df, gene_list)
write.csv(enrichment_results, "hflu-go_enrichment_results.csv", row.names = FALSE)





################### SPN

dev.new()

gff_file <- "spn.gff3"  # Replace with your GFF file

gff_df <- parse_gff3(gff_file)

go_genes <- gff_df %>% 
  filter(grepl("GO:0003735", go_terms)) %>%
  pull(locus_tag)


# Read and preprocess data
tb.spn <- read.delim("spn_df.tsv", header = TRUE, sep = '\t')

print(length(which(tb.spn[,3] >= 10)))
length(which(tb.spn[,3] >= 10)) / nrow(tb.spn)


goList = match(go_genes,tb.spn[,1])

tb.spn <- data.frame(Gene = tb.spn[, 2], Expression = tb.spn[, 3])  # Use column 2 for gene names
tb.spn$LogExpression <- log(tb.spn$Expression + 1, 10)  # Log-transform the expression data

# Specify the genes to highlight
highlight_genes <- c("psaA", "ply", "pspA", "cbpB", "lytA", "spxB", "pavB", "bgaA")
tb.spn$Highlight <- ifelse(tb.spn$Gene %in% highlight_genes, "Highlighted", "Normal")
tb.spn$Highlight[goList] = "Ribosomal"

# Compute quartiles
quartiles <- quantile(tb.spn$LogExpression, probs = seq(0, 1, by = 0.25), na.rm = TRUE)

# Create the jitter plot
jitter_plot <- ggplot(tb.spn, aes(x = 1, y = LogExpression)) +
  # Add normal points (gray)
  geom_jitter(data = subset(tb.spn, Highlight == "Normal"), color = "gray", width = 0.2, size = 3) +
  # Add ribosomal genes (blue)
  geom_jitter(data = subset(tb.spn, Highlight == "Ribosomal"), color = "blue", width = 0.2, size = 3) +
  # Add highlighted genes (red)
  geom_jitter(data = subset(tb.spn, Highlight == "Highlighted"), color = "red", width = 0.2, size = 3) +
  # Add quartile division lines
  geom_hline(yintercept = quartiles[-c(1, 5)], linetype = "dashed", color = "blue") +
  # Add labels for highlighted genes
  geom_text(
    data = subset(tb.spn, Highlight == "Highlighted"),
    aes(label = Gene),
    hjust = -0.2, vjust = 0, size = 4, color = "red"
  ) +
  theme_minimal() +
  labs(x = "Samples", y = "Log10(Expression)", title = "SPN with Highlighted Genes and Labels") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank())

# Add marginal histogram
ggExtra::ggMarginal(
  jitter_plot,
  type = "histogram",
  margins = "y",          # Add histogram to the y-axis
  binwidth = 0.1,         # Adjust bin width
  fill = "gray",          # Neutral fill for histogram
  color = "black"         # Border color for histogram
)




### Go enrichment

tb.spn = read.delim("spn_df.tsv",header=T,sep='\t')

gene_list = tb.spn[1:200,1] ## define gene list of interest

gff_file <- "spn.gff3"  # Replace with your GFF file

gff_df <- parse_gff3(gff_file)
enrichment_results <- perform_go_enrichment(gff_df, gene_list)
write.csv(enrichment_results, "spn-go_enrichment_results.csv", row.names = FALSE)



#################### mcat

dev.new()

gff_file <- "mcat.gff3"  # Replace with your GFF file

gff_df <- parse_gff3(gff_file)

go_genes <- gff_df %>% 
  filter(grepl("GO:0003735", go_terms)) %>%
  pull(locus_tag)


# Read and preprocess data
tb.mcat <- read.delim("mcat_df.tsv", header = TRUE, sep = '\t')

print(length(which(tb.mcat[,3] >= 10)))
length(which(tb.mcat[,3] >= 10)) / nrow(tb.mcat)


goList = match(go_genes,tb.mcat[,1])

tb.mcat <- data.frame(Gene = tb.mcat[, 2], Expression = tb.mcat[, 3])  # Use column 2 for gene names
tb.mcat$LogExpression <- log(tb.mcat$Expression + 1, 10)  # Log-transform the expression data

# Specify the genes to highlight
highlight_genes <- c("ompR", "ompA","pal","hfq","sodA")
tb.mcat$Highlight <- ifelse(tb.mcat$Gene %in% highlight_genes, "Highlighted", "Normal")
tb.mcat$Highlight[goList] = "Ribosomal"

# Compute quartiles
quartiles <- quantile(tb.mcat$LogExpression, probs = seq(0, 1, by = 0.25), na.rm = TRUE)

# Create the jitter plot
jitter_plot <- ggplot(tb.mcat, aes(x = 1, y = LogExpression)) +
  # Add normal points (gray)
  geom_jitter(data = subset(tb.mcat, Highlight == "Normal"), color = "gray", width = 0.2, size = 3) +
  # Add ribosomal genes (blue)
  geom_jitter(data = subset(tb.mcat, Highlight == "Ribosomal"), color = "blue", width = 0.2, size = 3) +
  # Add highlighted genes (red)
  geom_jitter(data = subset(tb.mcat, Highlight == "Highlighted"), color = "red", width = 0.2, size = 3) +
  # Add quartile division lines
  geom_hline(yintercept = quartiles[-c(1, 5)], linetype = "dashed", color = "blue") +
  # Add labels for highlighted genes
  geom_text(
    data = subset(tb.mcat, Highlight == "Highlighted"),
    aes(label = Gene),
    hjust = -0.2, vjust = 0, size = 4, color = "red"
  ) +
  theme_minimal() +
  labs(x = "Samples", y = "Log10(Expression)", title = "mcat with Highlighted Genes and Labels") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank())

# Add marginal histogram
ggExtra::ggMarginal(
  jitter_plot,
  type = "histogram",
  margins = "y",          # Add histogram to the y-axis
  binwidth = 0.1,         # Adjust bin width
  fill = "gray",          # Neutral fill for histogram
  color = "black"         # Border color for histogram
)



#### Go enrichment

tb.mcat = read.delim("mcat_df.tsv",header=T,sep='\t')

gene_list = tb.mcat[1:200,1] ## define gene list of interest

gff_file <- "mcat.gff3"  # Replace with your GFF file

gff_df <- parse_gff3(gff_file)
enrichment_results <- perform_go_enrichment(gff_df, gene_list)
write.csv(enrichment_results, "mcat-go_enrichment_results.csv", row.names = FALSE)
