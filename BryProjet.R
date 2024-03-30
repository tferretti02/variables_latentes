library(readr)
library(pscl)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(cluster)
library(glmnet)
library(dplyr)
library(tidyr)
data = read.csv("C:/Users/thoma/Downloads/Datagenus.csv",sep = "\t")
data = subset(data , select = -c(code))
convert_to_numeric <- function(column) {
  column <- gsub(",", ".", column)  # Remplacer les virgules par des points
  as.numeric(column)               # Convertir en nombre
}
# Appliquer cette fonction à toutes les colonnes de votre dataframe
data[] <- lapply(data, convert_to_numeric)
gen_columns <- grep("^gen", names(data), value = TRUE)
# Vérifier les types des colonnes pour s'assurer qu'ils sont maintenant numériques
str(data)
View(data)



X = data[names(data)[!names(data) %in% gen_columns]]

for(i in 1:27) {
  species_col_name <- paste("gen", i, sep="")
  new_col_name <- paste("Y", i, sep="")
  data[[new_col_name]] <- data[[species_col_name]] / data$surface
}



Y = data[grep("^Y", names(data), value = TRUE)]
logY = log1p(Y)


#####Tout ensemble
# Créer une formule qui combine toutes les variables Y
response_formula <- paste(names(Y), collapse = " + ")
full_formula <- as.formula(paste(response_formula, "~ ."))

# Ajuster le modèle de régression linéaire multivariée
multivariate_model <- lm(full_formula, data = cbind(Y, X))

# Voir le résumé du modèle
summary(multivariate_model)

# Créer une formule qui combine toutes les variables Y
response_formula <- paste(names(logY), collapse = " + ")
full_formula <- as.formula(paste(response_formula, "~ ."))

# Ajuster le modèle de régression linéaire multivariée
multivariate_model <- lm(full_formula, data = cbind(logY, X))

# Voir le résumé du modèle
summary(multivariate_model)

residuals <- resid(multivariate_model)

# Histogramme des résidus
par(mfrow = c(1, 2))
hist(residuals, breaks = 30, main = "Histogramme des Résidus", xlab = "Résidus")
# Diagramme Q-Q des résidus
qqnorm(residuals)
qqline(residuals, col = "red")
shapiro.test(residuals)

######Meilleur R² avec le logY en multivarié + residus suivent loi normal
par(mfrow = c(1, 1))

# Fonction pour créer un histogramme pour chaque gen avec le pourcentage de zéros
create_histograms <- function(data) {
  # Obtenir les noms des colonnes qui commencent par "gen"
  gen_columns <- grep("^gen", names(data), value = TRUE)
  
  # Parcourir chaque colonne gen et créer un histogramme
  for (gen_col in gen_columns) {
    # Calculer le pourcentage de zéros
    zero_percentage <- sum(data[[gen_col]] == 0) / nrow(data) * 100
    
    # Créer l'histogramme avec ggplot2
    p <- ggplot(data, aes_string(x = gen_col)) + 
      geom_histogram(binwidth = 1, fill = "blue", color = "black") +
      ggtitle(paste("Histogramme de", gen_col)) +
      xlab(gen_col) + ylab("Fréquence") +
      annotate("text", x = Inf, y = Inf, label = paste("Pourcentage de zéros:", round(zero_percentage, 2), "%"), 
               vjust = 1.5, hjust = 1.1, size = 5)
    
    # Afficher l'histogramme
    print(p)
  }
}

create_histograms(data)
density_columns = grep("^Y", names(data), value = TRUE)
other_columns <- names(data)[!names(data) %in% gen_columns] # Autres colonnes




library(corrplot)
correlations <- cor(logY[density_columns])
corrplot(correlations, method = "circle")
# Calcul de la matrice de distance
distance_matrix <- as.dist(1 - abs(correlations))

# Clustering hiérarchique
hc <- hclust(distance_matrix)

# Découpe du dendrogramme à un certain niveau pour identifier les clusters
groups <- cutree(hc, k = 2)

# Visualiser le dendrogramme
plot(hc)



########Loi 0infleted ?

library(pscl)
library(VGAM)
library(MASS)

# Exemple avec la variable Y1
# Modèle de Poisson
df = read.csv("C:/Users/thoma/Downloads/Datagenus.csv",sep = "\t")
df = subset(df , select = -c(code))
# Appliquer cette fonction à toutes les colonnes de votre dataframe
df[] <- lapply(df, convert_to_numeric)
Ydf = df[gen_columns]
Xdf = df[names(df)[!names(df) %in% gen_columns]]


poisson_model <- glm(Ydf$gen27 ~ ., data = X, family = poisson())

# Modèle binomial négatif
negbin_model <- glm.nb(Ydf$gen27 ~ ., data = X)

# Ajuster un modèle Zero-Inflated Poisson
zip_model <- zeroinfl(Ydf$gen27 ~ . | 1, data = Xdf, dist = "poisson")

# Ou un modèle Zero-Inflated Binomial Négatif
zinb_model <- zeroinfl(Ydf$gen27 ~ . | 1, data = Xdf, dist = "negbin")

# Comparer les modèles
# Vous pouvez utiliser des critères comme l'AIC pour comparer les modèles
AIC(poisson_model, negbin_model, zip_model,zinb_model)

###Que sur les gen et ca marche bien pour la 27 ou y a beaucoup de 0

#####Stats desciptive
hist(data$gen1, main = "Histogramme de gen1", xlab = "Valeurs", breaks = 30)
data$geology <- as.factor(data$geology)
barplot(table(data$geology), main = "Diagramme en Bâton de Géologie", xlab = "Catégories", ylab = "Fréquence")

data$forest <- as.factor(data$forest)
barplot(table(data$forest), main = "Diagramme en Bâton de forest", xlab = "Catégories", ylab = "Fréquence")





########ACP par theme
# Convertir 'geology' en variables indicatrices
data$geology = as.factor(data$geology)
geology_dummies <- model.matrix(~ geology - 1, data = data)

# Vous pouvez les ajouter au DataFrame original ou les garder séparément
data <- cbind(data, geology_dummies, forest_dummies)
data$forest = NULL
data$geology = NULL



library(FactoMineR)
library(factoextra)
gen_columns <- grep("^gen", names(data), value = TRUE)
env_columns = grep("^evi",names(data),value = TRUE)
geo_columns = c('lon','lat','surface')

data_acp <- data

# Obtenir les noms des colonnes 'gen'
gen_columns <- grep("gen", names(data), value = TRUE)

# Renommer les colonnes 'gen' dans le nouveau DataFrame
new_names <- as.character(1:length(gen_columns))
names(data_acp)[names(data_acp) %in% gen_columns] <- new_names

# Réaliser l'ACP sur les colonnes renommées du nouveau DataFrame
acp_gen <- PCA(data_acp[new_names], graph = FALSE)

# Visualisation
fviz_pca_var(acp_gen, 
             col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             title = "ACP - Variables Génétiques",
             repel = TRUE)

set.seed(123) # Pour la reproductibilité
# Extraction des scores de composantes principales pour les variables
scores_gen <- get_pca_var(acp_gen)$coord

# Clustering hiérarchique
hc_gen <- hclust(dist(scores_gen), method = "ward.D2")

# Coupe dendrogramme pour obtenir un nombre spécifique de clusters, disons 4
clusters_gen <- cutree(hc_gen, k = 4)

# Visualisation du dendrogramme avec les clusters
plot(hc_gen)
rect.hclust(hc_gen, k = 4, border = "red")

km_gen <- kmeans(scores_gen, centers = 4)
library(ClustOfVar)
result_clustofvar <- hclustvar(X.quanti = scores_gen)
plot(result_clustofvar)








# Identifier les colonnes des variables environnementales
evi_columns <- grep("evi", names(data), value = TRUE)


# Réaliser l'ACP sur les colonnes sélectionnées
acp_env <- PCA(data[,evi_columns], graph = FALSE)

# Visualisation sans étiquettes pour éviter l'encombrement
fviz_pca_var(acp_env, 
             col.var = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             title = "ACP - Variables Environnementales",
             repel = TRUE)
# Extraction des scores de composantes principales pour les variables
scores_gen <- get_pca_var(acp_env)$coord


# Effectuez un clustering hiérarchique
# Utilisez la méthode que vous préférez, par exemple 'ward.D2' pour la méthode Ward
hc <- hclust(dist(scores_gen), method = "ward.D2")
plot(hc)
rect.hclust(hc, k = 3, border = "red")
fviz_eig(acp_env, addlabels = TRUE, ylim = c(0, 100))


pcevi <- data.frame(matrix(ncol = 3, nrow = nrow(data)))
names(pcevi) <- c("PCEVI_1", "PCEVI_2", "PCEVI_3")


#Premier cluster
evi_columns = (hc$order[c(1,2,3,4,5,6,7,8,9,10,11,12)])
acp_cluster <- PCA(data[paste0("evi_", evi_columns)], graph = FALSE)
# Extraction des valeurs propres
valeurs_propres <- acp_cluster$eig
# Appliquer la règle de Kaiser
nombre_composantes <- sum(valeurs_propres[, 1] > 1)
print(nombre_composantes)#1 seule
scores_acp <- acp_cluster$ind$coord
pcevi <- scores_acp[, 1:1]
pcevi <- data.frame(pcevi)
names(pcevi) <- c("pc_evi_1")
head(pcevi)


##Deuxieme cluster
evi_columns = (hc$order[c(13,14,15,16,17)])
acp_cluster <- PCA(data[paste0("evi_", evi_columns)], graph = TRUE)
# Extraction des valeurs propres
valeurs_propres <- acp_cluster$eig
# Appliquer la règle de Kaiser
nombre_composantes <- sum(valeurs_propres[, 1] > 1)
print(nombre_composantes)
scores_acp <- acp_cluster$ind$coord
pcevi2 <- scores_acp[, 1:1]
pcevi2 <- data.frame(pcevi2)
pcevi["pc_evi_2"] = scores_acp[, 1:1] 

##Troisieme cluster
evi_columns = (hc$order[c(18,19,20,21,22,23)])
acp_cluster <- PCA(data[paste0("evi_", evi_columns)], graph = TRUE)
# Extraction des valeurs propres
valeurs_propres <- acp_cluster$eig
# Appliquer la règle de Kaiser
nombre_composantes <- sum(valeurs_propres[, 1] > 1)
print(nombre_composantes)
scores_acp <- acp_cluster$ind$coord
pcevi2 <- scores_acp[, 1:1]
pcevi2 <- data.frame(pcevi2)
pcevi["pc_evi_3"] = scores_acp[, 1:1] 


#####Pluvio
# ACP

pluvio_columns <- grep("pluvio_", names(data), value = TRUE)
acp_geo <- PCA(data[,pluvio_columns], graph = FALSE)
# Visualisation
fviz_pca_var(acp_geo, 
             col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             title = "ACP - Coordonnées Pluviométrique",
             repel = TRUE
)

#Premier cluster
acp_cluster <- PCA(data[,c("pluvio_1","pluvio_2","pluvio_12","pluvio_11")], graph = TRUE)
# Extraction des valeurs propres
valeurs_propres <- acp_cluster$eig
# Appliquer la règle de Kaiser
nombre_composantes <- sum(valeurs_propres[, 1] > 1)
print(nombre_composantes)#1 seule
scores_acp <- acp_cluster$ind$coord
pluvio <- scores_acp[, 1:1]
pluvio <- data.frame(pluvio)
names(pluvio) <- c("pluvioDF")
head(pluvio)
##Deuxieme cluster
acp_cluster <- PCA(data[,c("pluvio_6","pluvio_8","pluvio_7")], graph = TRUE)
# Extraction des valeurs propres
valeurs_propres <- acp_cluster$eig
# Appliquer la règle de Kaiser
nombre_composantes <- sum(valeurs_propres[, 1] > 1)
print(nombre_composantes)
scores_acp <- acp_cluster$ind$coord
pluvio2 <- scores_acp[, 1:1]
pluvio2 <- data.frame(pluvio2)
pluvio["pluvioJA"] = scores_acp[, 1:1] 
head(pluvio)






#####On rajoute les cluster et on eleve les evi de X + formation du X

X = data[names(data)[!names(data) %in% gen_columns]]
evi_columns <- grep("evi", names(data), value = TRUE)
X[evi_columns] = NULL
X[gen_columns] = NULL
X[names(pcevi)] = pcevi
X[names(pluvio)] = pluvio
X$geology = as.factor(X$geology)
geology_indicatrices <- model.matrix(~ geology - 1, data = X)
X <- cbind(X, geology_indicatrices)
X[c("pluvio_6","pluvio_8","pluvio_7","pluvio_1","pluvio_2","pluvio_12","pluvio_11")] = NULL
X$geology = NULL
X$geology6 = NULL #Pour ne pas l'avoir dans les interaction car plus forte modalité
X$forest = NULL

numerics <- sapply(X, is.numeric)
X_squares <- X[, numerics]**2
names(X_squares) <- paste0(names(X_squares), "^2")

# Ajouter des termes d'interaction
# Note : Cela peut créer un très grand nombre de nouvelles colonnes si vous avez beaucoup de variables
combinations <- combn(names(X[, numerics]), 2, simplify = FALSE)
interactions <- data.frame(
  lapply(combinations, function(cols) X[[cols[1]]] * X[[cols[2]]])
)
names(interactions) <- sapply(combinations, function(cols) paste(cols, collapse = ":"))

# Combinez le data frame original 'X', les carrés 'X_squares', et les interactions 'interactions'
X_enriched <- cbind(X, X_squares, interactions)
cols_to_remove <- colSums(X_enriched) == 0
X_enriched <- X_enriched[, !cols_to_remove]
X_enriched <- X_enriched[, !grepl("^2$", names(X_enriched)) & !grepl("geology.*\\^2$", names(X_enriched))]


######Pour Y
library(FactoMineR)
library(cluster)
data_acp <- data

# Obtenir les noms des colonnes 'gen'
gen_columns <- grep("gen", names(data), value = TRUE)

# Renommer les colonnes 'gen' dans le nouveau DataFrame
new_names <- as.character(1:length(gen_columns))
names(data_acp)[names(data_acp) %in% gen_columns] <- new_names
# Réaliser l'ACP
acp_gen <- PCA(data_acp[new_names], graph = FALSE)

# Clustering des variables basé sur les scores de l'ACP
# Choisissez le nombre de clusters (ici, 4 ou 5)
k <- 5 
clusters <- kmeans(acp_gen$var$coord, centers = k)

# Ajouter les clusters aux variables
variables_clustered <- data.frame(variable = rownames(acp_gen$var$coord), 
                                  cluster = clusters$cluster)

# Trouver la variable la plus corrélée à CP1 dans chaque cluster
most_correlated_vars <- sapply(split(variables_clustered, variables_clustered$cluster), function(cluster) {
  vars <- cluster$variable
  correlations <- apply(acp_gen$var$coord[vars, ], 1, function(x) x[1]) # Corrélation avec CP1
  vars[which.max(correlations)]
})

# Afficher les variables les plus corrélées à CP1 dans chaque cluster
print(most_correlated_vars)



Y = data_acp[c("12","23","11","24","4")]/data_acp$surface
Y = log1p(Y)####On a mis au log car meilleur représentation 




######REGRESSION
library(glmnet)

# Préparation des données
x <- as.matrix(X_enriched)
y <- as.matrix(Y)

# Appliquer la régression Ridge
set.seed(123) # Pour la reproductibilité
cv_ridge <- cv.glmnet(x, y, alpha = 0, family = "mgaussian")
best_lambda_ridge <- cv_ridge$lambda.min
ridge_model <- glmnet(x, y, alpha = 0, lambda = best_lambda_ridge,family = "mgaussian")
ridge_coeffs <- coef(ridge_model)

cv_lasso <- cv.glmnet(x, y, alpha = 1,family = "mgaussian") # alpha = 1 pour LASSO
best_lambda_lasso <- cv_lasso$lambda.min
lasso_model <- glmnet(x, y, alpha = 1, lambda = best_lambda_lasso,family = "mgaussian")
lasso_coeffs <- coef(lasso_model)
# Comparaison des erreurs de validation croisée
print(cv_ridge$cvm[cv_ridge$lambda == best_lambda_ridge])
print(cv_lasso$cvm[cv_lasso$lambda == best_lambda_lasso])

# Comparaison des coefficients
print(ridge_coeffs)
print(lasso_coeffs)

# Comparaison du nombre de coefficients non nuls
print(sum(ridge_coeffs != 0))
print(sum(lasso_coeffs != 0))

# Diviser les données en ensembles d'entraînement et de test
set.seed(123) # Pour la reproductibilité
train_indices <- sample(1:nrow(x), size = floor(0.8 * nrow(x)))
x_train <- x[train_indices, ]
y_train <- y[train_indices, ]
x_test <- x[-train_indices, ]
y_test <- y[-train_indices, ]

# Entraîner les modèles sur l'ensemble d'entraînement
cv_ridge <- cv.glmnet(x_train, y_train, alpha = 0,family = "mgaussian")
cv_lasso <- cv.glmnet(x_train, y_train, alpha = 1,family = "mgaussian")

# Prédire sur l'ensemble de test
predictions_ridge <- predict(cv_ridge, s = "lambda.min", newx = x_test)
predictions_ridge_adj <- matrix(predictions_ridge, nrow = nrow(y_test), ncol = ncol(y_test))
predictions_lasso <- predict(cv_lasso, s = "lambda.min", newx = x_test)
predictions_lasso_adj <- matrix(predictions_lasso, nrow = nrow(y_test), ncol = ncol(y_test))


# Calculer R^2
r2_ridge <- 1 - sum((y_test - predictions_ridge_adj)^2) / sum((y_test - mean(y_train))^2)
r2_lasso <- 1 - sum((y_test - predictions_lasso_adj)^2) / sum((y_test - mean(y_train))^2)

# Afficher les résultats
print(r2_ridge)
print(r2_lasso)

colnames(Y) <- paste0("Y_", colnames(Y))
response_formula <- paste(names(Y), collapse = " + ")
full_formula <- as.formula(paste(response_formula, "~ ."))

# Ajuster le modèle de régression linéaire multivariée
multivariate_model <- lm(full_formula, data = cbind(Y, X_enriched))

# Voir le résumé du modèle
summary(multivariate_model)




library(pls)
y
pls_model <- plsr(y ~ ., data = X_enriched,scale = TRUE, validation = "CV")

# Résumé du modèle
summary(pls_model)

# Choix du nombre de composantes
validationplot(pls_model, val.type = "MSEP")
pls_model <- plsr(y ~ ., data = X_enriched,scale = TRUE, ncomp = 50)
r2 <- R2(pls_model)
print(r2)


# Modèle de Poisson
df = read.csv("C:/Users/thoma/Downloads/Datagenus.csv",sep = "\t")
Ydf = df[c('gen12','gen23','gen11','gen24','gen4')]


poisson_model <- glm(Ydf$gen12 ~ ., data = X_enriched, family = poisson())

# Modèle binomial négatif
negbin_model <- glm.nb(Ydf$gen12 ~ ., data = X_enriched)

# Ajuster un modèle Zero-Inflated Poisson
zip_model <- zeroinfl(Ydf$gen12 ~ . | 1, data = X_enriched, dist = "poisson")

# Ou un modèle Zero-Inflated Binomial Négatif
zinb_model <- zeroinfl(Ydf$gen12 ~ . | 1, data = X_enriched, dist = "negbin")

# Comparer les modèles
# Vous pouvez utiliser des critères comme l'AIC pour comparer les modèles
AIC(poisson_model, negbin_model, zip_model,zinb_model)

# Prédire les valeurs
predictions_zinb <- predict(zinb_model, type = "response")

# Calculer le pseudo R^2
r2_zinb <- pseudoR2(Ydf$gen12, predictions_zinb)

# Créer le graphique
plot(Ydf$gen12, predictions_zinb, xlab = "Valeurs Observées", ylab = "Valeurs Prédites",
     main = "Prédictions du Modèle ZINB")
abline(0, 1, col = "red")  # Ligne y=x pour référence

# Vérifier les résumés des modèles
summary(zinb_model)

# Vérifier la distribution des prédictions et des observations
hist(predictions_zinb, main = "Distribution des Prédictions")
hist(Ydf$gen12, main = "Distribution des Observations")















################







#################

library(randomForest)


Y
X_enriched


library(caret)
set.seed(123) # Pour la reproductibilité
trainIndex <- createDataPartition(data$gen1, p = .8, 
                                  list = FALSE, 
                                  times = 1)
dataTrain <- data[trainIndex, ]
dataTest <- data[-trainIndex, ]
Y_Train <- dataTrain[,gen_columns] # Remplacer ... par les autres noms de colonnes gen
X_Train <- dataTrain[, !(names(dataTrain) %in% gen_columns)]
Y_Test <- dataTest[,gen_columns] # Remplacer ... par les autres noms de colonnes gen
X_Test <- dataTest[, !(names(dataTest) %in% gen_columns)]

models <- lapply(names(Y_Train), function(gen) {
  randomForest(X_Train, Y_Train[[gen]], ntree = 100)
})

predictions <- lapply(models, function(model) {
  predict(model, newdata = X_Test)
})

errors <- lapply(1:length(predictions), function(i) {
  mean((predictions[[i]] - Y_Test[[i]])^2)  # MSE
  # Ou pour MAE :
  # mean(abs(predictions[[i]] - Y_test[[i]]))
})



library(ggplot2)
plot_data <- data.frame(Real = unlist(Y_Test), Predicted = unlist(predictions))
ggplot(plot_data, aes(x = Real, y = Predicted)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Ligne y = x
  labs(x = "Valeurs Réelles", y = "Valeurs Prédites") +
  ggtitle("Comparaison des Valeurs Réelles et Prédites")



# Assurez-vous que X et Y sont correctement définis
Y <- as.matrix(data[, grep("gen", names(data))])
X <- as.matrix(data[, !names(data) %in% names(Y)])
X <- scale(X) # Standardisation des prédicteurs

lasso_models = cv.glmnet(X, data[,'gen2'],alpha =1) 
plot(lasso_models)
coef(lasso_models)


###################################


## data
data("bioChemists", package = "pscl")
View(bioChemists)

## without inflation
## ("art ~ ." is "art ~ fem + mar + kid5 + phd + ment")
fm_pois <- glm(art ~ ., data = bioChemists, family = poisson)
fm_qpois <- glm(art ~ ., data = bioChemists, family = quasipoisson)
fm_nb <- MASS::glm.nb(art ~ ., data = bioChemists)

## with simple inflation (no regressors for zero component)
fm_zip <- zeroinfl(art ~ . | 1, data = bioChemists)
fm_zinb <- zeroinfl(art ~ . | 1, data = bioChemists, dist = "negbin")

## inflation with regressors
## ("art ~ . | ." is "art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment")
fm_zip2 <- zeroinfl(art ~ . | ., data = bioChemists)
fm_zinb2 <- zeroinfl(art ~ . | ., data = bioChemists, dist = "negbin")

