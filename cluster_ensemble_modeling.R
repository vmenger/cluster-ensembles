library(readr)
library(beepr)
library(fmsb)
library(ggplot2)
library(NbClust)
library(diceR)
library(sigclust)

# Define data path
data_path <- ""

### 1. Read file ###

ysr <- read_delim(paste(data_path, "ysr.csv", sep=""), ";", escape_double = FALSE, trim_ws = TRUE)

### 2. Compute PCA coordinates for YSR data ###

pca.ysr <- prcomp(ysr, center = TRUE,scale. = TRUE)

### 3. Determine number of clusters ### 

# We use the NbClust package to compute several validity indices for cluster sizes 2-8, and then 
# fix the number of clusters based on a majority vote. Since we will be applying various algorithms
# to the dataset in the next step, in this step we include both k-means clustering, and hierarchical 
# clustering with complete-linkage and with wards method. Based on the output of the three commands
# combined, we select three as the number of clusters in this dataset.

nbclust.kmeans <- NbClust(data=ysr, 
                          distance="euclidean", 
                          min.nc=2, 
                          max.nc=7, 
                          method=c("kmeans"), 
                          index="all")

nbclust.kmeans <- NbClust(data=ysr, 
                          distance="euclidean", 
                          min.nc=2, 
                          max.nc=7, 
                          method=c("ward.D"), 
                          index="all")

nbclust.kmeans <- NbClust(data=ysr, 
                          distance="euclidean", 
                          min.nc=2, 
                          max.nc=7, 
                          method=c("complete"), 
                          index="all")

### 4. Determine clustering algorithms ###

# We apply all algorithms to the dataset, and then pick a diverse subset
all_algs <- c("hc", "nmf", "diana", "km", "pam", "ap", "sc", "gmm", "block", "som", "cmeans", "hdbscan")
cc_objs <- c()

# Iterate over all algorithms
for (algorithm_name in all_algs) {
  
  # Use the consensus cluster function to apply the single algorithm to the dataset
  # Number of clusters is set to three, as described above  
  cc <- diceR::consensus_cluster(ysr, 
                                 nk=3, 
                                 p.item=1, 
                                 reps=1, 
                                 nmf.method=c("brunet"),
                                 algorithms=c(algorithm_name))
  
  # Append the clustering to a list so we can access it later
  cc_objs[[algorithm_name]] <- cc

  # Plot title
  plot.title = sprintf("%s algorithm", algorithm_name)
  
  # Create pca plot and save it
  png(sprintf(paste(data_path, "pca_%s.png", sep=""), algorithm_name)) 
  plot(pca.ysr$x[,'PC1'], pca.ysr$x[,'PC2'], col=cc[,,,])
  title(plot.title)
  dev.off()  

}

### 5. Apply cluster ensemble ###

# We select the k-means, gaussian mixture models, and affinity propagation algorithms 
# based on the pca plots. Furthermore, we will trim the cluster portfolio using Rank Aggregation
# based on several internal validity indices, and use CSPA to obtain one single clustering.
ysr.ensemble <- dice(ysr, 
                     nk=3, 
                     algorithms=c('ap', 'km', 'gmm'), 
                     reps=5, 
                     cons.funs=c("kmodes", "majority", "CSPA", "LCE"),
                     trim=TRUE,
                     reweigh=FALSE,
                     seed=1)

### 6. Write results so analysis can be continued in Python ###

# Gather results, together with PCA coordinates and single algorithm clusters in a dataframe
cluster_result_df <- cbind(pca.ysr$x[,'PC1'], 
                           pca.ysr$x[,'PC2'],
                           cc_objs[['km']][,,,],
                           cc_objs[['gmm']][,,,],
                           cc_objs[['ap']][,,,], 
                           col=ysr.ensemble$clusters[,'CSPA']
                          )

# Convert to dataframe and add appropriate column names
cluster_result_df <- as.data.frame(cluster_result_df)

colnames(cluster_result_df) <- c("PCA-1", 
                                 "PCA-2", 
                                 "partition_km", 
                                 "partition_gmm", 
                                 "partition_ap", 
                                 "partition_ensemble")

# Write output to file
write.table(cluster_result_df,
            file=paste(data_path, "clusters_result.csv", sep=""),
            sep=";"
            )

### 6. Assess statistical significance using siglust package ###

sig.pval <- sigclust(ysr, k=3, nsim=100, labflag=0, label=ysr.ensemble$clusters[,'CSPA'])@pval