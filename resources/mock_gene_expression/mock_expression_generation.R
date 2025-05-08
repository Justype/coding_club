library(tidyverse)

generate_mock_expression_matrix <- function(num_genes = 200, 
                                            degs = 20,
                                            samples_per_group = 5, 
                                            dispersion = 0.01, 
                                            fold_change = 3,
                                            random_seed = 42) {
  # Check if degs is divisible by 2
  if (degs %% 2 != 0) {
    stop("Number of DEGs must be even")
  }
  
  # Check if num_genes is divisible by degs
  if (num_genes %% degs != 0) {
    stop("num_genes must be divisible by DEGs")
  }
  
  set.seed(random_seed)
  
  total_samples <- samples_per_group * 2
  
  # Baseline mean expression for all genes
  base_means <- runif(num_genes, min = 5, max = 500)
  
  # Identify DEGs
  gene_base <- num_genes / degs
  upregulated_genes <- seq(gene_base, num_genes/2, by = gene_base)
  downregulated_genes <- seq(num_genes/2 + gene_base, num_genes, by = gene_base)
  
  # Define group-specific means
  mu_control <- base_means
  mu_treatment <- base_means
  
  # Apply fold changes
  mu_treatment[upregulated_genes] <- mu_control[upregulated_genes] * fold_change
  mu_treatment[downregulated_genes] <- mu_control[downregulated_genes] / fold_change
  
  # Initialize matrix
  expression_matrix <- matrix(nrow = num_genes, ncol = total_samples)
  
  # Simulate counts
  for (i in 1:num_genes) {
    control_counts <- rnbinom(samples_per_group, size = 1/dispersion, mu = mu_control[i])
    treatment_counts <- rnbinom(samples_per_group, size = 1/dispersion, mu = mu_treatment[i])
    expression_matrix[i, ] <- c(control_counts, treatment_counts)
  }
  
  # Label rows/columns
  rownames(expression_matrix) <- paste0("Gene", 1:num_genes)
  colnames(expression_matrix) <- c(paste0("Control", 1:samples_per_group),
                                   paste0("Treatment", 1:samples_per_group))
  
  
  return(expression_matrix)
}

gene_wise_wilcox_test <- function(expression_matrix) {
  require(dplyr)
  expression_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
    mutate(Group = ifelse(grepl("Control", Sample), "Control", "Treatment")) %>%
    group_by(Gene) %>%
    summarise(
      log2FC = log2(mean(Expression[Group == "Treatment"]) / mean(Expression[Group == "Control"])),
      p.val = wilcox.test(Expression ~ Group)$p.value
    ) %>%
    mutate(p.adj = p.adjust(p.val, method = "fdr")) %>% 
    return
}

gene_wise_log_t_test <- function(expression_matrix) {
  require(dplyr)
  expression_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
    mutate(Group = ifelse(grepl("Control", Sample), "Control", "Treatment")) %>%
    group_by(Gene) %>%
    summarise(
      log2FC = log2(mean(Expression[Group == "Treatment"]) / mean(Expression[Group == "Control"])),
      p.val = t.test(log2(Expression) ~ Group)$p.value
    ) %>%
    mutate(p.adj = p.adjust(p.val, method = "fdr")) %>% 
    return
}


experiment1 <- generate_mock_expression_matrix(num_genes = 200, degs = 10, random_seed = 42)
experiment1_wilcox_results <- gene_wise_wilcox_test(experiment1)
experiment1_log_t_results <- gene_wise_log_t_test(experiment1)
write_csv(experiment1_wilcox_results, "experiment1_wilcox_results.csv")
write_csv(experiment1_log_t_results, "experiment1_log_t_results.csv")


experiment2 <- generate_mock_expression_matrix(num_genes = 200, degs = 20, random_seed = 34)
experiment2_wilcox_results <- gene_wise_wilcox_test(experiment2)
experiment2_log_t_results <- gene_wise_log_t_test(experiment2)
write_csv(experiment2_wilcox_results, "experiment2_wilcox_results.csv")
write_csv(experiment2_log_t_results, "experiment2_log_t_results.csv")


