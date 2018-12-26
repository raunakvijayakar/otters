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
View(spearman_matrix)
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


















