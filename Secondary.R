#Raunak Vijayakar
#Credit to Yeo Zhi Yi, Philip Johns, and all who have worked on this project.
#24/12/18
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
library(dbscan)


#Here I am just working on selecting columns from spectral 
#analysis results to reduce collinearity. 
#This can be added before the call identification classifier. 
#I will also work on clustering the clean data Zhi Yi prepared.
#I will use an unsupervised learner to create initial homologies.
#I will then validate the clustering using either discriminant function analysis, 
#as in Mumm 2014, or using techniques in clValid.
#I will then investigate if those broad homologies have similar
#representative signals, or if they may benefit from further clustering.
setwd("~/Desktop/YNC/Otter:Bat Stuff")
spectral_output <- read.csv("Manual Cluster Clean Data.csv", header = T)
length(names(spectral_output)) #we have a maximum of 31 variables as input
length(spectral_output$selec) #and 335 observations
#We do not see p >> N, so I'm not too worried about dimensionality. 
#But variance is a concern for overfit. Now, that is mostly dealt with 
#by using k-medoids, which reduces variance.
#Removal is unnecessary for the first classifier in Prelim, 
#as we have used a lasso method, whereby only useful variables are kept (Hastie et al. 2017).
spearman_matrix <- cor(spectral_output[, -c(1:3)], method = "spearman")
#We see that temporal and frequency variables are positively correlated to each other.
#Skew and kurtosis as well. However, few are greater than r = 0.9.
#If we perform principal component analysis, we can summarize the data, and reduce redundancies.

#Though this is an unsupervised learning problem, we can actually use supervised clustering in 
#finding initial homologies, as Zhi Yi has manually checked this subset of data.
#This way we could see if PCA for dimension reduction improves accuracy, and should be applied for the unsupervised problem.
#But, it may be the case that other homologies exist in other data, so we do need an unsupervised method as well.

setwd("~/Desktop/YNC/Otter:Bat Stuff/Raw Audio")
otter_1 <- readWave("Crystal 2 June PR12 - 2.wav")
oscillo(otter_1, 44100, colwave = "cadetblue")
spectro(otter_1, 44100, 512, ovlp = 50, osc = T, main = "Spectrogram of Homology")
#I am sampling 10 calls for each homology ZY found, and making spectrograms of them to have an idea of what we're working with.
homo_1 <- spectral_output[spectral_output$Manual.Cluster == 1, ]
homo_2 <- spectral_output[spectral_output$Manual.Cluster == 2, ]
set.seed(123)
homo_1_sample <- homo_1[sample(nrow(homo_1), 10, replace = F), ]
homo_2_sample <- homo_2[sample(nrow(homo_2), 10, replace = F), ]

#input to spectrogram generation function
homo_1_sample$sound.files

spectro_gen <- function (x) 
{
  for (i in 1:length(x)) 
  {
    k <- readWave(as.character(x[i]))
    spectro(k, k@samp.rate, 512, ovlp = 87.5, osc = T, main = as.character(i), flim = c(0, 10))
  }
}
#frequency resolution of 44100/512 = 86.1 Hz. 
#Mumm used a resolution of 46.9 Hz. There is an inherent tradeoff between frequency and time resolution.
#I have chosen to prioritise time resolution for easier visual interpretation.
spectro_gen(homo_1_sample$sound.files)
spectro_gen(homo_2_sample$sound.files)


#Going to do a rework for the unsupervised clustering we were doing before
#First thing is dimension reduction using PCA prior to learning.
#The second is the use of HDBSCAN as our call clustering method.
set.seed(123)
scramble <- spectral_output[sample(nrow(spectral_output), nrow(spectral_output), replace = F), ]
pca_trainset <- scramble[1:201, -(1:3)]
pca_testset <- scramble[202:335, -(1:3)]
pca <- prcomp(pca_trainset, scale. = T)
pca_var <- (pca$sdev)^2
prop_var <- pca_var/sum(pca_var)
plot(cumsum(prop_var), xlab = "Principal Component", ylab = "Cumulative Proportion of Variance Explained", type = "b")
View(cumsum(prop_var))
#using the first 15 components provides 98% of the variance, with roughly half the variables
View(pca$x)
pcs <- data.frame(pca$x)
#now, we can use this as a prior to the initial call identification model.

set.seed(123)
setwd("~/Desktop/YNC/Otter:Bat Stuff/Unknowns")
sound.files <- list.files(pattern = "wav$")
sound.file.sample <- sound.files
autodetec.output <- autodetec(flist = sound.file.sample, ssmooth = 350, power = 1,
                              bp = c(3, 13), ls = TRUE, wl = 1024, 
                              flim = c(0, 18), mindur = 0.05, maxdur = 0.5) # 34 min runtime :(
call.features <- specan(autodetec.output)
call.features <- na.omit(call.features)
pca_call <- prcomp(call.features[, -(1:2)], scale. = T)
pca_call_var <- (pca_call$sdev)^2
prop_call_var <- pca_call_var/sum(pca_call_var)
plot(cumsum(prop_call_var), xlab = "Principal Component", ylab = "Cumulative Proportion of Variance Explained", type = "b")
View(cumsum(prop_call_var))
#Taking the first 14 components gives 97.8% of the variance.

#the section below is just so we can test the model
set.seed(123)
call.features.man <- call.features
call.features.man <- call.features.man[sample(1:nrow(call.features.man), 20, replace = F), ]
call.features.man$ID <- seq_len(nrow(call.features.man))
call.features.man$accuracy <- as.factor(c(0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0))

#This is where I will select another 20-30 calls to manually check through, and add to the previous set





set.seed(123)
indices <- createDataPartition(call.features.man$ID, p = 0.4, list = FALSE)
train <- call.features.man[indices, ]
test <- call.features.man[-indices, ]

#selecting new data with pca
pca_out <- pca_call$x
train_new <- data.frame(train$accuracy, pca_out[rownames(pca_out) %in% rownames(train), ])
test_new <- data.frame(test$accuracy, pca_out[rownames(pca_out) %in% rownames(test), ])
#LEAVING NOTE: when you come back, pick the columns for training and testing glm
train_new <- train_new[, 1:15]
test_new <- test_new[, 2:15]


control <- trainControl(method = "LOOCV", number = 10, repeats = 5)
glm.model <- train(train_new[, -1], 
                   # The second function is the accuracy column
                   train_new[, 1], 
                   method = "glmnet", tuneLength = 12, trControl = control)

glm.predictions <- predict(glm.model, newdata = test_new)

1-sum(glm.predictions == test$accuracy)/nrow(test) #test error rate

#Consider that the accuracy here does not need to be great if hdbscan can filter noise calls
#ideally though, i would have maybe 50-100 identified calls and not calls to properly train the model


#HDBSCAN examples from dbscan vignette
data("moons")
plot(moons, pch = 20)
cl <- hdbscan(moons, minPts = 5)
cl
plot(moons, col=cl$cluster+1, pch = 20)
plot(cl$hc, main = "Dendrogram")
plot(cl, show_flat = T)

data("DS3")
plot(DS3, pch = 20)
cl2 <- hdbscan(DS3, minPts = 25)
plot(DS3, col=cl2$cluster+1, 
     pch=ifelse(cl2$cluster == 0, 8, 1), # Mark noise as star
     cex=ifelse(cl2$cluster == 0, 0.5, 0.75), # Decrease size of noise
     xlab=NA, ylab=NA)
colors <- sapply(1:length(cl2$cluster), 
                 function(i) adjustcolor(palette()[(cl2$cluster+1)[i]], alpha.f = cl2$membership_prob[i]))
points(DS3, col=colors, pch=20)

plot(cl2, show_flat = T)

#HDBSCAN applied to ZY's data
hdb_data <- na.omit(spectral_output[spectral_output$Manual.Cluster == 1, -(1:3)]) #selecting clean calls

hdb <- hdbscan(hdb_data, minPts = 15)
hdb
plot(hdb, show_flat = T)
#concurs with ZY's manual selection of 2 clusters, but finds a lot of noise
plot(hdb$hc)
vis_data <- cbind(cluster = hdb$cluster, hdb_data)
lda.fit <- lda(cluster ~ ., data = vis_data)
coord.hdb <- as.matrix(hdb_data) %*% lda.fit$scaling
coord.hdb <- data.frame(hdb$cluster, coord.hdb)

ggplot(coord.hdb, aes(LD1, LD2)) + 
  geom_point(aes(colour = hdb$cluster)) + 
  scale_color_gradientn(colours = rainbow(3), name = "Clusters")

plot(coord.hdb$LD1, coord.hdb$LD2, col=hdb$cluster+1, 
     pch=ifelse(hdb$cluster == 0, 8, 1), # Mark noise as star
     cex=ifelse(hdb$cluster == 0, 0.5, 0.75), # Decrease size of noise
     xlab="LD1", ylab="LD2")

pairs(hdb_data, col = hdb$cluster+1) #scatterplot matrix
plot(hdb_data$skew, hdb_data$meandom, col = hdb$cluster+1) #close to a clear clustering in spectral vars
#it's ok to have a selective algorithm, at least from data quantity perspective

#HDBSCAN applied directly to the autodetec output, w/o the selection model (since it's not yet been trained on enough data)
hdb.full_data <- call.features[, -(1:2)]
hdb.full <- hdbscan(hdb.full_data, minPts = 5)
hdb.full
#turns out that it's way too noisy without any prior selection. Thinks it's all noise

#LEAVING NOTE: Consider tsne and pca as ways to visualize the clustering results. 
#Determine how hdbscan defines noise. 
#Perhaps ask it to recluster after removing everything it thinks it noise









library(bclust)

#we will attempt the bayesian cluster on ZY's cleaned data since we know what it tends to look like


head(spectral_output)
b_data <- as.matrix(spectral_output[spectral_output$Manual.Cluster == 1, -c(1:2)])
#cannot assume replicated data, so no x.id, no labels either
meansumsq <- meancss(b_data)
optimfunc <- function(phi) { -loglikelihood(x.mean = meansumsq$mean, x.css = meansumsq$css,
                    repno = meansumsq$repno, transformed.par = phi, var.select = FALSE)
}
xinit.tpar <- optim(rep(0, 5), optimfunc, method = "BFGS")$par
optimfunc <- function(phi) { -loglikelihood(x.mean = meansumsq$mean, x.css = meansumsq$css,
                     repno = meansumsq$repno, transformed.par = c(xinit.tpar[1:4], phi))
}
x.tpar <- c(xinit.tpar[1:4], optim(rep(0, 2), optimfunc,
                                   method = "BFGS")$par)

x.labels <- c(names(spectral_output[ , -c(1,2)]))

bclust.obj <- bclust(b_data, transformed.par = x.tpar)

dptplot(bclust.obj, scale = 10, horizbar.plot = TRUE,
        varimp = imp(bclust.obj)$var, horizbar.distance = 5, dendrogram.lwd = 2)

bclust.obj$var.select

par(mfrow=c(2,1))

plot(bclust.obj$clust.number, bclust.obj$logposterior, type = "b")
abline(h=max(bclust.obj$logposterior))
abline(v=obj$clust.number)

plot(spectral_output$meandom, spectral_output$entropy)

#apparently the cut occured at 2428.627. Seems unnecessary
#so we have the classic situation of a function (log post) that we are maximizing.
bclust.obj$optim.clustno
#apparently the optimal number of clusters is 55...
#The max value does occur at 55, but 
#we need a function that penalizes increasing clusters in relation to gain in loglik
#I'm not sure how to do that. It's a more objective way to prevent this kind of madness
#Rather than me saying that the cluster number should be 5 or whatever.
#Maybe some kind of percent increase in loglik per added cluster?

ditplot(bclust.obj, xlab = colnames(bclust.obj$data),
        ylab = bclust.obj$labels, dendrogram.lwd = 1, dendrogram.size = 2,
        xlab.mar = 3, ylab.mar = 3, image.col = rainbow(20),
        horizbar.plot = FALSE,
        horizbar.col = rev(c(heat.colors(5)[-4], "white")),
        horizbar.distance = 4,varimp = rep(0, ncol(bclust.obj$data)),
        horizbar.size = 0.5, vertbar = NULL,
        vertbar.col = rainbow(max(vertbar)),
        teeth.size = 0.25, plot.width = 10)

bclust.obj$merge

loglikelihood(means$mean, means$css, transformed.par = )
bayes <- bclust(b_data, transformed.par = )

#run the bayesian clustering on the 4 homologies from


library(fpc)
#here's an idea. Use a hierarchical clusterer, then bootstrap(?) to determine stability of results.
#I already did a method like that (stability validation, which just cfm pam is best)
useful <- as.matrix(spectral_output[spectral_output$Manual.Cluster == 1, -c(1:2)])
matrix <- scale(useful)
pcenter <- attr(matrix, "scaled:center")  
pscale <- attr(matrix, "scaled:scale")
#NOTE: i'm just doing this as proof of concept rn, but I would do the bclust here instead
d <- dist(matrix, "euclidean")
fit <- hclust(d, "ward.D")
plot(fit)
#let's just say that this shows 4 clusters, bclust would actually tell us, 
#but if we don't use bclust we can try 
#thorndike method to choose. For now though, just 4
rect.hclust(fit, k=4)
kbest <- 4
cboot.hclust <- clusterboot(matrix, clustermethod = hclustCBI, method = "ward.D",
                             k=kbest)

#the clusterboot func has a lot of cluster options, but bclust is not one of them, so if we want it,
#we'll need to write a custom interface
groups <- cboot.hclust$result$partition  
cboot.hclust$bootmean #want this closer to 1
cboot.hclust$bootbrd
#for interpretation, see details of clusterboot documentation page




#I think we've finally got it with the following
library(mclust)
spectral_output <- read.csv("Manual Cluster Clean Data.csv", header = T)
mdata <- spectral_output[spectral_output$Manual.Cluster == 1, -c(1:2)]
#selecting just the first homology that ZY identified, consistent with what we've been doing
mdata_cut <- data.frame(cbind(mdata$meandom, mdata$entropy))
colnames(mdata_cut) = c("meandom", "entropy")

#intial clustering
bayes <- mclustBIC(mdata_cut)
plot(bayes)
summary(bayes)

mmod <- Mclust(mdata_cut, x = bayes)
summary(mmod, parameters = TRUE)
plot(mmod, what = "classification")
mmod$classification #who's in each cluster
#finds 4 clusters, which matches our result from silhouette

#checking if thos clusters contain sub clusters
clus1 <- mdata[mmod$classification == 1, ]
bayes1 <- mclustBIC(clus1)
plot(bayes1)
summary(bayes1)
mmod1 <- Mclust(clus1, x = bayes1)
summary(mmod1, parameters = TRUE)
mmod1$classification
#no further cluster


clus2 <- mdata[mmod$classification == 2, ]
bayes2 <- mclustBIC(clus2)
plot(bayes2)
summary(bayes2)
mmod2 <- Mclust(clus2, x = bayes2)
summary(mmod2, parameters = TRUE)
mmod2$classification
#no further cluster

clus3 <- mdata[mmod$classification == 3, ]
bayes3 <- mclustBIC(clus3)
plot(bayes3)
summary(bayes3)
mmod3 <- Mclust(clus3, x = bayes3)
summary(mmod3, parameters = TRUE)
mmod3$classification
#no further cluster

clus4 <- mdata[mmod$classification == 4, ]
bayes4 <- mclustBIC(clus4)
plot(bayes4)
summary(bayes4)
mmod4 <- Mclust(clus4, x = bayes4)
summary(mmod4, parameters = TRUE)
mmod4$classification
#no further cluster




