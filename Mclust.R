setwd("~/Desktop/YNC/Otter:Bat Stuff")
library(mclust)
spectral_output <- read.csv("Manual Cluster Clean Data.csv", header = T)
spectral_output$index <- 1:nrow(spectral_output) #reference tag to extract calls
mdata <- spectral_output[spectral_output$Manual.Cluster == 1, -c(1:2)]
#selecting just the first homology that ZY identified, consistent with what we've been doing
mdata_cut <- data.frame(cbind(mdata$meandom, mdata$entropy, mdata$index))
colnames(mdata_cut) = c("meandom", "entropy", "index")

#intial clustering
bayes <- mclustBIC(mdata_cut[, -3])
plot(bayes)
summary(bayes)

mmod <- Mclust(mdata_cut[, -3], x = bayes)
summary(mmod, parameters = TRUE)
plot(mmod, what = "classification")

mmod$classification #who's in each cluster
#finds 4 clusters, which matches our result from silhouette


#checking if those clusters contain sub clusters
clus1 <- mdata[mmod$classification == 1, ]
bayes1 <- mclustBIC(clus1[, -30])
plot(bayes1)
summary(bayes1)
mmod1 <- Mclust(clus1[, -30], x = bayes1)
summary(mmod1, parameters = TRUE)
mmod1$classification
#no further cluster

clus2 <- mdata[mmod$classification == 2, ]
bayes2 <- mclustBIC(clus2[, -30])
plot(bayes2)
summary(bayes2)
mmod2 <- Mclust(clus2[, -30], x = bayes2)
summary(mmod2, parameters = TRUE)
mmod2$classification
#no further cluster

clus3 <- mdata[mmod$classification == 3, ]
bayes3 <- mclustBIC(clus3[, -30])
plot(bayes3)
summary(bayes3)
mmod3 <- Mclust(clus3[, -30], x = bayes3)
summary(mmod3, parameters = TRUE)
mmod3$classification
#no further cluster

clus4 <- mdata[mmod$classification == 4, ]
bayes4 <- mclustBIC(clus4[, -30])
plot(bayes4)
summary(bayes4)
mmod4 <- Mclust(clus4[, -30], x = bayes4)
summary(mmod4, parameters = TRUE)
mmod4$classification
#no further cluster




#now onto extracting calls and creating spectrograms
#calculate euclidean distance to determine proximity to central coords
summary(mmod, parameters = T) #gives center coords

#cluster 1
clus1$cen <- sqrt((4.373880-clus1$meandom)^2+(0.744392-clus1$entropy)^2)
centro1 <- clus1[clus1$cen == min(clus1$cen), ]
centro1 <- spectral_output[spectral_output$index %in% centro1$index, ]

#cluster 2
clus2$cen <- sqrt((0.1183985-clus2$meandom)^2+(0.7389623-clus2$entropy)^2)
centro2 <- clus2[clus2$cen == min(clus2$cen), ]
centro2 <- spectral_output[spectral_output$index %in% centro2$index, ]

#cluster 3
clus3$cen <- sqrt((0.4249154-clus3$meandom)^2+(0.7569739-clus3$entropy)^2)
centro3 <- clus3[clus3$cen == min(clus3$cen), ]
centro3 <- spectral_output[spectral_output$index %in% centro3$index, ]

#cluster 4
clus4$cen <- sqrt((1.9035132-clus4$meandom)^2+(0.7720732-clus4$entropy)^2)
centro4 <- clus4[clus4$cen == min(clus4$cen), ]
centro4 <- spectral_output[spectral_output$index %in% centro4$index, ]



setwd("~/Desktop/YNC/Otter:Bat Stuff/Raw Audio")

spectro_gen <- function (x) 
{
  for (i in 1:length(x)) 
  {
    k <- readWave(as.character(x[i]))
    spectro(k, k@samp.rate, 512, ovlp = 87.5, osc = T, main = as.character(i), flim = c(0, 10))
  }
}

#create spectrograms
spectro_gen(centro1$sound.files)
spectro_gen(centro2$sound.files)
spectro_gen(centro3$sound.files)
spectro_gen(centro4$sound.files)












#extracting calls [OLD]
ex1 <- spectral_output[spectral_output$index %in% clus1$index, ]
ex2 <- spectral_output[spectral_output$index %in% clus2$index, ]
ex3 <- spectral_output[spectral_output$index %in% clus3$index, ]
ex4 <- spectral_output[spectral_output$index %in% clus4$index, ]

set.seed(123)
ex1_sample <- ex1[sample(nrow(ex1), 10, replace = F), ]
ex2_sample <- ex2[sample(nrow(ex2), 10, replace = F), ]
ex3_sample <- ex3[sample(nrow(ex3), 10, replace = F), ]
ex4_sample <- ex4[sample(nrow(ex4), 10, replace = F), ]


spectro_gen(ex1_sample$sound.files)

spectro_gen(ex2_sample$sound.files)
spectro_gen(ex3_sample$sound.files)
spectro_gen(ex4_sample$sound.files)




