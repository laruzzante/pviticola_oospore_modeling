require(lubridate)
require(lme4)
require(broom.mixed)
require(effects)
require(dplyr)
require(FactoMineR)
require(factoextra)
require(missMDA)
require(ggbiplot)

## DATASET PROCESSING
setwd("~/mnt/Data-Work-CH/22_Plant_Production-CH/222.6_Mycologie_protected/Projets de recherche/38_SMALA/Livio - oospore modeling/pviticola_oospore_modeling")
df_all <- read.table("Oosp_not_all_2003-2024_v9.csv", sep = ";", header = T)
df_all$BBCH <- as.numeric(df_all$BBCH)
df_all$date <- as_datetime(df_all$date, format = "%d.%m.%Y")
df_all$MTG <- as.numeric(df_all$MTG)
df_all$nb_germ_oosp_1d <- as.numeric(df_all$nb_germ_oosp_1d)
df_all$cumul_precipit_1Jan <-  as.numeric(df_all$cumul_precipit_1Jan)
df_all$nb_days_rainfall_30d <-  as.numeric(df_all$nb_days_rainfall_30d) 
df_all$solar_radiation_1Jan <- as.numeric(df_all$solar_radiation_1Jan)
df_all$VPD <- as.numeric(df_all$VPD)
df_all$RH <- as.numeric(df_all$RH)
df_all$temp <- as.numeric(df_all$temp)
df_all$TDD <- as.numeric(df_all$TDD)

# solar_radiation variables were ultimately not included in the model variable selection
# because they were strongly correlated with TDD, thus biasing the predictions.
# Also, they included a lot of missing values, thus making TDD a better variable choice.

### PCA FUNCTION
pca <- function(df){
  dataPCA <- cbind(df$cumul_precipit_1Jan, df$nb_days_rainfall_30d, df$VPD,
                   df$RH, df$temp, df$TDD)
  dataPCA <- matrix(as.numeric(unlist(dataPCA)),nrow = nrow(dataPCA))
  colnames(dataPCA) <- (colnames(subset(df, select = c(cumul_precipit_1Jan, nb_days_rainfall_30d, VPD,
                                                       RH, temp, TDD))))
  pca <- prcomp(dataPCA, scale. = T)
  print(summary(pca))
  print(pca$rotation)
  ## PLOTS
  # specifying MTG categories for PCA groups
  MTG_cat <- df$MTG
  for (i in 1:length(MTG_cat)) {
    if (MTG_cat[i] < 3) {
      MTG_cat[i] <- "1-2"
    }
    if (MTG_cat[i] > 2) {
      MTG_cat[i] <- "3-10"
    }
  }
  p <- ggbiplot(pca, groups = MTG_cat, choices = c(1,2), ellipse = T, ellipse.prob = 0.4) + theme_bw()
  print(p)
}


### MODEL FUNCTIONS, NEEDS DATASET AS INPUT.
## the two functions creates distinct models: one with MGT as response variable,
## the other with Nspores as response variable
## they then plot the model partial plots, the QQ-residuals, the table statistics

### Average oospore maturation day
model_MGT <- function(df){
  MGT_model <- glm(data = df, formula =  MTG ~ cumul_precipit_1Jan + nb_days_rainfall_30d + 
                     + VPD + RH + temp + TDD, family = "poisson")
  
  # SHOWING DISTRIBUTION OF MAIN RESPONSE VARIABLES OF INTEREST
  hist(df$MTG)
  
  # MODEL INFO AND PARTIAL EFFECTS PLOTS
  plot(MGT_model)
  plot(allEffects(MGT_model))
  
  # MODEL STATISTICS TABLES
  tidy(MGT_model)
  # glance(MGT_model)
}

### Number of spores 1 day after first germination
model_Nspores1d <- function(df){
  Nspores_model <- glm(data = df, formula =  nb_germ_oosp_1d ~ cumul_precipit_1Jan + nb_days_rainfall_30d + 
                         VPD + RH + temp + TDD, family = "poisson") 

  # SHOWING DISTRIBUTION OF MAIN RESPONSE VARIABLES OF INTEREST
  hist(df$nb_germ_oosp_1d)
  
  # MODEL INFO AND PARTIAL EFFECTS PLOTS
  plot(Nspores_model)
  plot(allEffects(Nspores_model))
  
  # MODEL STATISTICS TABLES
  tidy(Nspores_model)
  # glance(Nspores_model)
}

### Number of spores 10 days after first germination
model_Nspores10d <- function(df){
  Nspores_model <- glm(data = df, formula =  nb_germ_oosp_10d ~ cumul_precipit_1Jan + nb_days_rainfall_30d + 
                         VPD + RH + temp + TDD, family = "poisson") 
  
  # SHOWING DISTRIBUTION OF MAIN RESPONSE VARIABLES OF INTEREST
  hist(df$nb_germ_oosp_10d)
  
  # MODEL INFO AND PARTIAL EFFECTS PLOTS
  plot(Nspores_model)
  plot(allEffects(Nspores_model))
  
  # MODEL STATISTICS TABLES
  tidy(Nspores_model)
  # glance(Nspores_model)
}


### Random Forest and Partition Trees (decision trees)
random_forest <- function(df){
  
  require(randomForest)
  require(caret)
  
  ### RANDOM FOREST COMPUTAION ON ALL EXPLANATORY VARIABLES
  
  # dataRF <- subset(df, select = -c(date, solar_radiation_30d, solar_radiation_1Jan, MTG, maturity, nb_germ_oosp_1d, nb_germ_oosp_10d))
  dataRF <- as.data.frame(cbind(df$cumul_precipit_1Jan, df$nb_days_rainfall_30d, df$VPD,
                                df$RH, df$temp, df$TDD))
  # specifying MTG categories for PCA groups
  MTG_cat <- df$MTG
  for (i in 1:length(MTG_cat)) {
    if (MTG_cat[i] < 3) {
      MTG_cat[i] <- "1-2"
    }
    if (MTG_cat[i] > 2) {
      MTG_cat[i] <- "3-10"
    }
  }
  
  colnames(dataRF) <- colnames(subset(df, select = c(cumul_precipit_1Jan, nb_days_rainfall_30d, VPD,
                                                      RH, temp, TDD)))
  
  dataRF$MTG_cat <- MTG_cat
  
  
  dataRF$MTG_cat <- as.factor(dataRF$MTG_cat)
  set.seed(111)
  ind <- sample(2, nrow(dataRF), replace = TRUE, prob = c(0.7, 0.3))
  train <- dataRF[ind == 1,]
  train$MTG_cat <- factor(train$MTG_cat)
  test <- dataRF[ind == 2,]
  test$MTG <- factor(test$MTG_cat)
  
  rf <- randomForest(MTG_cat~., data = train, proximity = TRUE, mtry = 3)
  print(rf)
  plot(rf)
  
  p1 <- predict(rf, train)
  confusionMatrix(p1, train$MTG_cat)
  
  p2 <- predict(rf, test)
  confusionMatrix(p2, test$MTG_cat)
  
  
  t <- tuneRF(subset(train, select = -c(MTG_cat)), train[,"MTG_cat"],
              stepFactor = 0.5,
              plot = TRUE,
              ntreeTry = 150,
              trace = TRUE,
              improve = 0.05)
  
  hist(treesize(rf),
       main = "No. of Nodes for the Trees",
       col = "green")
  # Variable Importance
  varImpPlot(rf,
             sort = T,
             n.var = 6,
             main = "Ranked Variable Importance")
  print(importance(rf))
  #MeanDecreaseGini
  
  partialPlot(rf, train, TDD, as.factor("1-2"))
  partialPlot(rf, train, TDD, "3-10")
  partialPlot(rf, train, nb_days_rainfall_30d, "1-2")
  partialPlot(rf, train, nb_days_rainfall_30d, "3-10")
  partialPlot(rf, train, TDD, "1-2")
  partialPlot(rf, train, TDD, "3-10")
  
  # devtools::install_github("araastat/reprtree")
  
  require(devtools)
  install_github("araastat/reprtree")
  library(reprtree)
  
  ## Repartition tree / Decision tree
  reprtree:::plot.getTree(rf)
  
  ## MDS plot
  MDSplot(rf, train$MTG_cat)
  
  require("rpart")
  require("rpart.plot")
  
  ### PARTITION TREES COMPUTAION ON SELECTED VARIABLES ONLY
  
  dataRpart <- subset(df, select = c(cumul_precipit_1Jan, nb_days_rainfall_30d, VPD,
                                     RH, temp, TDD))
  dataRpart <- matrix(as.numeric(unlist(dataRpart)),nrow = nrow(dataRpart))
  colnames(dataRpart) <- colnames(subset(df, select = c(cumul_precipit_1Jan, nb_days_rainfall_30d, VPD,
                                                        RH, temp, TDD)))
  dataRpart <- as.data.frame(dataRpart)
  
  ## Decision tree for MTG
  rf1 <- rpart(df$MTG ~., data = dataRpart, method = "poisson")
  rpart.plot(rf1)
  
  ## Decision tree for Nspores1d
  rf2 <- rpart(df$nb_germ_oosp_1d ~., data = dataRpart, method = "poisson")
  rpart.plot(rf2)
  
  ## Decision tree for Nspores10d
  rf3 <- rpart(df$nb_germ_oosp_10d ~., data = dataRpart, method = "poisson")
  rpart.plot(rf3)
  
}


## ALL BBCH DATASET
model_MGT(df_all)
model_Nspores1d(df_all)
model_Nspores10d(df_all)
pca(df_all)
random_forest(df_all)

## DATASET BBCH 0:12
df <- df_all %>% filter(df_all$BBCH < 13)
model_MGT(df)
model_Nspores1d(df)
model_Nspores10d(df)
pca(df)
random_forest(df)

## DATASET BBCH 13:+
df <- df_all %>% filter(df_all$BBCH >= 13)
model_MGT(df)
model_Nspores1d(df)
model_Nspores10d(df)
pca(df)
random_forest(df)
