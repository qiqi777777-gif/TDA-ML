# ================================================================
# Load libraries
# ================================================================
library(TDA)
library(fda)
library(ggplot2)
library(ggrepel)
library(patchwork)
set.seed(123)

# ================================================================
# Path to files
# ================================================================
ph_folder <- "output/PH" 
files <- list.files(ph_folder, pattern = "_PH.csv$", full.names = TRUE)

#extract the numeric ID from filename
get_number <- function(fname){
 as.numeric(sub(".*adjmatrix([0-9]+)_PH.*", "\\1", fname))
}

#order files numerically
file_numbers <- sapply(files, get_number)
files <- files[order(file_numbers)]

# ================================================================
# Compute persistence landscapes
# ================================================================
DataSymp <- c()
DataAsymp <- c()
tseq <- seq(0, 1, length.out = 500)

read_landscape <- function(file_path, dim = 1, tseq = tseq) {
  PH <- read.csv(file_path)
  colnames(PH) <- c("Birth", "Death", "dimension")
  PH <- PH[c("dimension", "Birth", "Death")]
  PH <- as.matrix(PH)
  x <- landscape(PH, dimension = dim, KK = 1, tseq)
  return(x)
}

for(f in files){
  x1 <- read_landscape(f, dim = 1, tseq = tseq)
  
  # detect group from filename
  fname <- basename(f)
  if(grepl("^Asymptomatic", fname)){
    DataAsymp <- cbind(x1, DataAsymp)
  } else if(grepl("^Symptomatic", fname)){
    DataSymp <- cbind(x1, DataSymp)
  } else {
    warning(paste("File not assigned to group:", fname))
  }
}


DataSymp <- t(DataSymp)
DataAsymp <- t(DataAsymp)
DataBoth <- rbind(DataAsymp, DataSymp)

# ================================================================
# Create FD objects
# ================================================================
create_fd <- function(Data){
  x <- seq(0, 1, length.out = ncol(Data))
  basis <- create.bspline.basis(c(0,1), nbasis = ncol(Data)-1, norder = 2)
  fdobj <- smooth.basis(x, t(Data), basis)
  return(fdobj$fd)
}

fdSymp <- create_fd(DataSymp)
fdAsymp <- create_fd(DataAsymp)
fdBoth <- create_fd(DataBoth)

plot(mean.fd(fdSymp), main="Mean Symptomatic Landscape")
plot(mean.fd(fdAsymp), main="Mean Asymptomatic Landscape")

# ================================================================
# Functional PCA
# ================================================================
pca <- pca.fd(fdBoth, nharm=2)
m <- c("A1","A2","A3","A4","A5","A6","A7","A8",
       "S1","S2","S3","S4","S5","S6","S7","S8","S9")
col.group <- c(rep("black", 8), rep("blue", 9))
group <- c(rep("Asymptomatic", 8), rep("Symptomatic", 9))

pca_df <- data.frame(
  PC1 = pca$scores[, 1],
  PC2 = pca$scores[, 2],
  ID = m,
  Group = group
)

  ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = ID)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 4, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 4, show.legend = FALSE, max.overlaps = Inf) +
    scale_color_manual(values = c("Asymptomatic" = "black", "Symptomatic" = "blue")) +
    labs(
      x = "PC1 Score",
      y = "PC2 Score",
      title = ""
    ) +
    coord_fixed() +  
    theme_classic(base_size = 16) +
    theme(
      legend.position = ifelse(T, "right", "none"),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )


# ================================================================
# Permutation hypothesis test
# ================================================================
#Test of asymptomatic and symptomatic
ptest1 <- tperm.fd(fdAsymp, fdSymp, nperm = 10000, q = 0.05, plotres = TRUE)
ptest1$pval

#Test of symptomatic and symptomatic
splits <- combn(1:9, 4)
n_splits <- ncol(splits)
pvals <- numeric(n_splits)

for(i in 1:n_splits){
  
  group1_idx <- splits[, i]
  group2_idx <- setdiff(1:9, group1_idx)
  
  fdSymp1 <- create_fd(DataSymp[group1_idx, ])
  fdSymp2 <- create_fd(DataSymp[group2_idx, ])
  
  ptest <- tperm.fd(fdSymp1, fdSymp2,
                    nperm = 10000,
                    q = 0.05,
                    plotres = FALSE)
  
  pvals[i] <- ptest$pval
  print(i)
}

mean_p <- mean(pvals)
sd_p   <- sd(pvals)

mean_p
sd_p

#Test of asymptomatic and asymptomatic
splits <- combn(1:8, 4)
n_splits <- ncol(splits)
pvals <- numeric(n_splits)

for(i in 1:n_splits){
  
  group1_idx <- splits[, i]
  group2_idx <- setdiff(1:8, group1_idx)
  
  fdSymp1 <- create_fd(DataAsymp[group1_idx, ])
  fdSymp2 <- create_fd(DataAsymp[group2_idx, ])
  
  ptest <- tperm.fd(fdSymp1, fdSymp2,
                    nperm = 10000,
                    q = 0.05,
                    plotres = FALSE)
  
  pvals[i] <- ptest$pval
  print(i)
}

mean_p <- mean(pvals)
sd_p   <- sd(pvals)

mean_p
sd_p
