# Author: Kinga Golebiewska
# This is a function for calculation of normalized expression data from
# qPCR results. The output is similar to the easyqPCR's one - a list with expression
# values and sd.
      #1. Import CT values from qPCR in xlsx file. The first column of the table 
          #should be named "sample" with names of samples and technical replicates.
          #The second column of the table should be named "group" with sample type. 
          #The rest cols are the CT values for each analyzed gene
          #(check the file "sample_data.xlsx", sheet = Ct).
              #e.g. read_excel("C:/R/sample_data.xlsx", sheet = "Ct")
      #2. Import LinReg values with calculated primer amplification efficiency in xlsx
          #file (check the file "sample_data.xlsx", sheet = LinReg).
              #e.g. read_excel("C:/R/sample_data.xlsx", sheet = "LinReg")
      #3. Type name of the control group (for normalization):
              #e.g. normalize <- "N0_0"
      #4. Create a vector with names of the reference genes (at least 2):
              #e.g. reference <- c("ACT2", "UPL7")
      #5. Run the qPCR.expression function:


library(readxl)
my_reads <- read_excel("sample_data.xlsx", sheet = "Ct")
linreg_effic <- read_excel("sample_data.xlsx", sheet = "LinReg")
normalize <- "N0_0.4"
reference <- c("RRN16S", "RRN23S")


qPCR.expression <- function(my_reads, linreg_effic, normalize, reference) {
  require(dplyr); require(psych); require(data.table)
  my_reads$sample <- factor(my_reads$sample)
  my_reads$group <- factor(my_reads$group)
  genes <- colnames(my_reads[,3:length(my_reads)])
  for (i in 3:ncol(my_reads)) {
    y <- my_reads[,c(1,2,i)]
    x <- linreg_effic %>% filter(gene == genes[i-2])
    y <- cbind(y, effic=x$value)
    g_eff <- y[,4]^y[,3]
    y[ , ncol(y) + 1] <- g_eff
    colnames(y)[ncol(y)] <- "g_eff"
    colnames(y)[3] <- "gene"
    assign(genes[i-2], y)
  }
  
  m1 <- data.table(get(reference[1])[1])
  for (i in 1:length(reference)) {
    x <- get(reference[i])[5]
    names(x)[1] <- reference[i]
    m1 <- cbind(m1, x)
  }
  m2 <- as.matrix(m1[,2:ncol(m1)])
  rownames(m2) <- m1$sample
  geo_ref <- exp(rowMeans(log(m2)))
  geo_ref2 <- cbind(m1[,1], geo_ref)
  
  genes_list <- list()
  for (i in 1:length(genes)) {
    x <- get(genes[i])
    genes_list[[genes[i]]] <- x
  }
  
  names <- c("sample", names(genes_list))
  exp_summary <- data.frame(levels(my_reads$group))
  exp_sd <- data.frame(levels(my_reads$group))
  for (g in genes_list) {
    g <- g %>% left_join(geo_ref2, by = "sample") %>% group_by(sample) %>% 
      mutate(georef_eg = geo_ref/g_eff)
    n <- g %>% group_by(group) %>% 
      summarise(geo_reps = geometric.mean(georef_eg))
    g <- merge(x = g, y = n, by = "group", all.x = T)
    normal <- g %>% filter(group == normalize)
    normal <- normal[1,c(1,8)]
    g <- g %>% mutate(expression = georef_eg/normal$geo_reps)
    e <- g %>% group_by(group) %>% 
      summarise(express_mean = geometric.mean(expression), sd = sd(expression))
    exp_summary <- cbind(exp_summary, e$express_mean)
    exp_sd <- cbind(exp_sd, e$sd)
  }
  colnames(exp_summary) <- names
  colnames(exp_sd) <- names
  results <- list(expression = exp_summary, sd = exp_sd)
  return(results)
}


qPCR.expression(my_reads = my_reads, linreg_effic = linreg_effic, normalize = normalize, reference = reference)
