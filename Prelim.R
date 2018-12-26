library(caret)
library(warbleR)
library(knitr)
library(factoextra)
library(cluster)
library(clValid)
library(glmnet)
library(seewave)
library(tuneR)
library(MASS)
library(ggplot2)
library(mclust)

set.seed(123)
# setwd("~/Desktop/YNC/Otter:Bat Stuff/Raw Audio")
# sound.files <- list.files(pattern = "wav$")
# sound.file.sample <- sound.files[sample(1:469, 20, replace = FALSE)]
setwd("~/Desktop/YNC/Otter:Bat Stuff/Unknowns")
sound.files <- list.files(pattern = "wav$")
sound.file.sample <- sound.files
autodetec.output <- autodetec(flist = sound.file.sample, ssmooth = 350, power = 1,
                              bp = c(3, 13), ls = TRUE, wl = 1024, 
                              flim = c(0, 18), mindur = 0.05, maxdur = 0.5) # 34 min runtime :(
call.features <- specan(autodetec.output)
call.features <- na.omit(call.features)
call.features <- call.features[sample(1:nrow(call.features), 20, replace = F), ]
call.features$ID <- seq_len(nrow(call.features))

call.features$accuracy <- as.factor(c(0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0))

indices <- createDataPartition(call.features$ID, p = 0.3, list = FALSE)
train <- call.features[indices, ]
test <- call.features[-indices, ]

glm.model

coefs <- cor(call.features[ ,3:31], method = "spearman")

coef.df <- as.data.frame(coefs, header = TRUE)
coef.df$name <- row.names(coefs)
#WORK IN PROGRESS 
desc <- coef.df[order(coef.df$duration, decreasing = TRUE), ]
View(desc)


control <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
glm.model <- train(train[, 3:31], 
                  # The second function is the accuracy column
                  train[, 33], 
                  method = "glmnet", trControl = control)
glm.predictions <- predict(glm.model, newdata = test[, 3:31])
1-sum(rf.predictions == test$accuracy)/nrow(test)
print(glm.model)
#adaboost and rf test error rate is 16%. On a very small sample. 2 non-calls classified as calls
#xgb was awful. glment had test error rate of 8.3% (1 wrong). I recommend this.

#Moving on to the unsupervised learning part
setwd("~/Desktop/YNC/Otter:Bat Stuff")
#this is similar to the output of specan()
unsup <- read.csv("Manual Cluster Clean Data.csv", header = T)
str(unsup)
manual.cluster.1.subset <- unsup[unsup$Manual.Cluster == "1", ]
kable(head(manual.cluster.1.subset[1:5]))

plot(manual.cluster.1.subset$meandom, manual.cluster.1.subset$entropy)



#selecting subgroups with clustering rather than visual analysis of spectrograms
initial <- data.frame(cbind(manual.cluster.1.subset$meandom, 
          manual.cluster.1.subset$entropy, 1:nrow(manual.cluster.1.subset)))
initial.std <- scale(initial)
clusters.init <- fviz_nbclust(initial.std, pam, method = "silhouette")
print(clusters.init)
mod.init <- pam(initial.std, k = 2)
#plot(mod.init)          

#trying to initially cluster based on just two variables
init.data <- cbind(cluster = mod.init$clustering, initial.std)
init.pam <- lda(formula = cluster ~ ., data = data.frame(init.data))
coord.pam <- initial.std %*% init.pam$scaling
coord.pam <- data.frame(mod.init$clustering, coord.pam)
ggplot(coord.pam, aes(LD1, LD2)) + 
  geom_point(aes(colour = mod.init.clustering)) + 
  scale_color_gradientn(colours = rainbow(2), name = "Clusters") + theme_bw() +
  labs(title = "LDA Cluster Plot of Initial Data")
plot(mod.init)

Dist.init <- dist(initial.std, method = "euclidean")
tree.init <- hclust(Dist.init, method = "ward.D")
plot(tree.init, main = "Cluster Dendrogram of Initial Data")
rect.hclust(tree.init, k = 4, border = 1:4)

#then cluster each of the three groups using the other methods (LDA and dendro)
mod.init$clustering

g1 <- manual.cluster.1.subset[mod.init$clustering == 1, ]
g2 <- manual.cluster.1.subset[mod.init$clustering == 2, ]
g3 <- manual.cluster.1.subset[mod.init$clustering == 3, ]
g4 <- manual.cluster.1.subset[mod.init$clustering == 4, ]


#Group 1
standard.g1 <- scale(g1[, 4:31])
clustersdb.g1 <- fviz_nbclust(standard.g1, pam, method = "silhouette")
print(clustersdb.g1)
#silhouette suggests k=8, but k = 3 gives a good increase without making too many clusters given data size
mod.g1 <- pam(x = standard.g1, k = 3, metric = "euclidean", stand = FALSE)
mod.g1.data <- cbind(cluster = mod.g1$clustering, standard.g1)

LDA.mod.g1 <- lda(formula = cluster ~ ., data = data.frame(mod.g1.data))
coord.mod.g1 <- standard.g1 %*% LDA.mod.g1$scaling
coord.mod.g1 <- data.frame(mod.g1$clustering, coord.mod.g1)
ggplot(coord.mod.g1, aes(LD1, LD2)) + 
  geom_point(aes(colour = mod.g1.clustering)) + 
  scale_color_gradientn(colours = rainbow(3), name = "Clusters") + theme_bw() +
  labs(title = "LDA Cluster Plot of First Group")
#dendrogram
Dist.g1 <- dist(standard.g1, method = "euclidean")
tree.g1 <- hclust(Dist.g1, method = "ward.D")
plot(tree.g1, main = "Cluster Dendrogram of First Group", xlab = "")
rect.hclust(tree.g1, k = 3, border = 1:3)




#Group 2
standard.g2 <- scale(g2[, 4:31])
clustersdb.g2 <- fviz_nbclust(standard.g2, pam, method = "silhouette")
print(clustersdb.g2)
#silhouette suggests k=9, but k=4 gives good increase in silhouette without over clustering
mod.g2 <- pam(x = standard.g2, k = 4, metric = "euclidean", stand = FALSE)
mod.g2.data <- cbind(cluster = mod.g2$clustering, standard.g2)

LDA.mod.g2 <- lda(formula = cluster ~ ., data = data.frame(mod.g2.data))
coord.mod.g2 <- standard.g2 %*% LDA.mod.g2$scaling
coord.mod.g2 <- data.frame(mod.g2$clustering, coord.mod.g2)
ggplot(coord.mod.g2, aes(LD1, LD2)) + 
  geom_point(aes(colour = mod.g2.clustering)) + 
  scale_color_gradientn(colours = rainbow(4), name = "Clusters") + theme_bw() +
  labs(title = "LDA Cluster Plot of Second Group")
#dendrogram
Dist.g2 <- dist(standard.g2, method = "euclidean")
tree.g2 <- hclust(Dist.g2, method = "ward.D")
plot(tree.g2, main = "Cluster Dendrogram of Second Group", xlab = "")
rect.hclust(tree.g2, k = 4, border = 1:4)




#Group 3
standard.g3 <- scale(g3[, 4:31])
standard.g3 <- standard.g3[-62, ] #outlier observed at rownum 273, row 62
clustersdb.g3 <- fviz_nbclust(standard.g3, pam, method = "silhouette")
print(clustersdb.g3)

#silhouette suggests k=2, but it cannot be visualized by LDA due to collinearities
mod.g3 <- pam(x = standard.g3, k = 2, metric = "euclidean", stand = FALSE)
mod.g3.data <- cbind(cluster = mod.g3$clustering, standard.g3)
#will have to use components to visualize instead
clusplot(mod.g3, main = "Component Cluster Plot of Third Group")


# LDA.mod.g3 <- lda(formula = cluster ~ ., data = data.frame(mod.g3.data))
# coord.mod.g3 <- standard.g3 %*% LDA.mod.g3$scaling
# coord.mod.g3 <- data.frame(mod.g3$clustering, coord.mod.g3)
# ggplot(coord.mod.g3, aes(LD1, LD2)) + 
#   geom_point(aes(colour = mod.g3.clustering)) + 
#   scale_color_gradientn(colours = rainbow(2), name = "Clusters") + theme_bw()

#dendrogram
Dist.g3 <- dist(standard.g3, method = "euclidean")
tree.g3 <- hclust(Dist.g3, method = "ward.D")
plot(tree.g3, main = "Cluster Dendrogram of Third Group")
rect.hclust(tree.g3, k = 2, border = 1:2)


#Group 4
standard.g4 <- scale(g4[, 4:31])
standard.g4 <- data.frame(standard.g4)
standard.g4 <- standard.g4[-13, ] #outlier observed at rownum 178, row 13
standard.g4 <- standard.g4[, -22]
standard.g4 <- as.matrix(standard.g4)
clustersdb.g4 <- fviz_nbclust(standard.g4, pam, method = "silhouette")
print(clustersdb.g4)
#silhouette suggests k=10, but k=4 gives good increase in silhouette without over clustering
mod.g4 <- pam(x = standard.g4, k = 4, metric = "euclidean", stand = FALSE)
mod.g4.data <- cbind(cluster = mod.g4$clustering, standard.g4)

LDA.mod.g4 <- lda(formula = cluster ~ ., data = data.frame(mod.g4.data))
coord.mod.g4 <- standard.g4 %*% LDA.mod.g4$scaling
coord.mod.g4 <- data.frame(mod.g4$clustering, coord.mod.g4)
ggplot(coord.mod.g4, aes(LD1, LD2)) + 
  geom_point(aes(colour = mod.g4.clustering)) + 
  scale_color_gradientn(colours = rainbow(4), name = "Clusters") + theme_bw() +
  labs(title = "LDA Cluster Plot of Fourth Group")
#dendrogram
Dist.g4 <- dist(standard.g4, method = "euclidean")
tree.g4 <- hclust(Dist.g4, method = "ward.D")
plot(tree.g4, main = "Cluster Dendrogram of Fourth Group")
rect.hclust(tree.g4, k = 4, border = 1:4)







row.names(manual.cluster.1.subset) <- 1:nrow(manual.cluster.1.subset) #rownames were all wrong

standard.1 <- scale(manual.cluster.1.subset[, 4:31])
clustersdb.1 <- fviz_nbclust(standard.1, pam, method = "silhouette")
print(clustersdb.1)

pam.1 <- pam(x = standard.1, k = 4, metric = "euclidean", stand = FALSE)
print(pam.1)
plot(pam.1)
pam.data.1 <- cbind(cluster = pam.1$clustering, standard.1)

LDA.pam <- lda(formula = cluster ~ ., data = data.frame(pam.data.1))
coord.pam <- standard.1 %*% LDA.pam$scaling
coord.pam <- data.frame(pam.1$clustering, coord.pam)
ggplot(coord.pam, aes(LD1, LD2)) + 
  geom_point(aes(colour = pam.1.clustering)) + 
  scale_color_gradientn(colours = rainbow(4), name = "Clusters") + theme_bw()

coord.pam$LD1
coord.pam$LD2
coord.pam$index <- 1:nrow(coord.pam)



#Model Validation
pam.data.1 <- as.data.frame(pam.data.1)
sil <- silhouette(pam.1, full = T)
#Connectivity corresponds to what extent items are placed in the same cluster 
#as their nearest neighbors in the data space. 
#The connectivity has a value between 0 and infinity. Minimize.
connectivity(Data = pam.data.1, method = "euclidean", clusters = pam.data.1$cluster) # 88.95437

#Dunn index is the ratio of the smallest distance between observations 
#not in the same cluster to the largest intra-cluster distance. Maximize.
dunn(Data = pam.data.1, method = "euclidean", clusters = pam.data.1$cluster) # 0.08753131





#a possible way to comibine the validation and training methods
stability.valid <- clValid(standard.1, 2:10, clMethods = c("pam","fanny","clara",
                          "agnes","kmeans","hierarchical"), 
                           validation = c("internal"))

slotNames(stability.valid)
summary(stability.valid) #this suggests a few different k values for pam, usually 2, 4, or the max
clusters(stability.valid)


#stability validation for all columns, removing one at a time, means given for all k
stab <- matrix(0,nrow=ncol(standard.1),ncol=4)
colnames(stab) <- c("APN","AD","ADM","FOM")
output <- data.frame(matrix(0,nrow=10,ncol=4))
colnames(output) <- c("APN","AD","ADM","FOM")
k_val <- 1:10 

for (i in k_val) {
  Dist <- dist(standard.1, method = "euclidean")
  clusterObj <- pam(Dist, k = i)
  #nc <- 4 ## number of clusters
  cluster <- clusterObj$clustering
  for (j in 1:ncol(standard.1)) 
  {
    matDel <- standard.1[,-j]
    DistDel <- dist(matDel,method = "euclidean")
    clusterObjDel <- pam(DistDel, k = i)
    clusterDel <- clusterObjDel$clustering
    stab[j,] <- stability(standard.1, Dist, j, cluster, clusterDel)
  }
  output[i, ] <- colMeans(stab)
}


#validation of pam using stability measures suggests either 2 or 10 clusters. Not good.
for (i in 1:ncol(standard.1)) 
{
  matDel <- standard.1[,-i]
  DistDel <- dist(matDel,method = "euclidean")
  clusterObjDel <- pam(DistDel, k = 2)
  clusterDel <- clusterObjDel$clustering
  stab[i,] <- stability(standard.1, Dist, i, cluster, clusterDel)
}
colMeans(stab)





#trying out a few different kinds of method for clustering
#pam.2 <- sota(data = standard.1, maxCycles = 3, metric = "euclidean")
#pam.2 <- pam(x = standard.1, k = 4, metric = "manhattan", stand = F)
#This should come up with 4 relatively equal clusters at some point to mesh with LDA
Dist.1 <- dist(standard.1, method = "euclidean")
pam.2 <- hclust(Dist.1, method = "average")
plot(pam.2)
rect.hclust(pam.2, h = 35, border = 1:4)

pam.2$data
pam.2$height
#this analysis using a heirarchical algo suggests an outlier, point 273
#Also have to settle on method for hclust. Euclidean is typical for dist.
View(standard.1)
standard.2 <- standard.1
standard.2 <- standard.2[-273, ]
#note that rowname 273 is gone as well, not been reindexed

#hclust w/o outlier point
Dist.2 <- dist(standard.2, method = "euclidean")
pam.3 <- hclust(Dist.2, method = "ward.D")
plot(pam.3, labels = F, frame.plot = T)
rect.hclust(pam.3, k = 4, border = 1:4)
#compare the number found in these clusters to the numbers in LDA groups
#also the correspondance is unlikely to be 1:1, i.e. not same calls in same grps



#a different thing, that I don't remember now...
pam.data.2 <- data.frame(cbind(cluster = pam.2$clust, standard.1))

LDA.pam.2 <- lda(formula = cluster ~ ., data = data.frame(pam.data.2))
coord.pam.2 <- standard.1 %*% LDA.pam.2$scaling
coord.pam.2 <- data.frame(pam.2$clustering, coord.pam.2)

ggplot(coord.pam.2, aes(LD1, LD2)) + 
  geom_point(aes(colour = pam.2.clustering)) + 
  scale_color_gradientn(colours = rainbow(10), name = "Clusters")


temp <- eclust(x = standard.1, FUNcluster = "hclust", hc_metric = "euclidean")
temp$nbclust
fviz_cluster(temp, geom = "point",
            ggtheme = theme_minimal())








