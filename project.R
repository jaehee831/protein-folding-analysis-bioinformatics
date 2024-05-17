# Load Necessary Libraries
library(glmnet)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(readr)



# Load Datasets
load("~/Downloads/microbiome/trainMicrobiomeData.RData") # Load the dataset
load("~/Downloads/microbiome/trainMicrobiomeMetadata.RData") # Load the meta_datasets



# Data Transposing
trainMicrobiomeData_t <- as.data.frame(t(trainMicrobiomeData))
colnames(trainMicrobiomeData_t) <- rownames(trainMicrobiomeData)
rownames(trainMicrobiomeData_t) <- NULL
trainMicrobiomeData_t$sample.ID <- rownames(t(trainMicrobiomeData))
merged_data <- merge(trainMicrobiomeData_t, trainMicrobiomeMetadata, by = "sample.ID") # Merge the data
merged_data$Total_Microbiome <- rowSums(merged_data[,2:223], na.rm = TRUE) # Total Sum Scaling - proportional normalization



# TSS - Proportional Normalization
merged_data_ratio <- merged_data
for (i in 2:223) {merged_data_ratio[, i] <- merged_data[, i] / merged_data$Total_Microbiome}
merged_data_ratio$Total_Microbiome <- NULL # Remove the 'Total_Microbiome' column as it's no longer needed

sampled_rows <- sample(1:nrow(merged_data_ratio), size = floor(4 * nrow(merged_data_ratio) / 5)) # Random Sampling 4 out of 5
merged_ratio_sampled <- merged_data_ratio[sampled_rows, ] # Create new data frames based on the sampled rows



## Parameters for Elastic Model Fitting
afeatures <- c("Barnesiella", "Odoribacter", "Streptococcus", "Anaerostipes", "Intestinimonas","Ruminococcus",
              "Allisonella", "Evtepia", "Faecalimonas","Anaerobutyricum","Negativibacillus","Escherichia")


X_train <- as.matrix(merged_ratio_sampled[, features, drop = FALSE]) # bring feature matrix to model
y_train <- merged_ratio_sampled$Age # bring target variable(age)
X_question <- as.matrix(merged_data_ratio[, features, drop = FALSE]) # Prepare the data for prediction 

lasso_model <- glmnet(X_train, y_train, alpha = 0.75, nlambda = 100) # Fit the Lasso model
age_estimations <- predict(lasso_model, s = 0.001, newx = X_question) # Make age estimations. s=0.01 is convention

merged_data_ratio$Estimated_Age <- as.numeric(age_estimations) # Convert the matrix to a numeric vector



## Lasso Bootstrapping
n_bootstrap <- 300 # Initialize a matrix to s
age_estimations_lasso_bootstrap <- matrix(0, nrow = nrow(X_question), ncol = n_bootstrap) # tore bootstrap estimates

for (i in 1:n_bootstrap) {
  bootstrap_indices <- sample(1:nrow(X_train), replace = TRUE) # Resample the indices
  
  X_bootstrap <- X_train[bootstrap_indices, ] # feature matrix_ bring each Iterations
  y_bootstrap <- y_train[bootstrap_indices] # bring target variable each iterations
  
  lasso_model_bootstrap <- glmnet(X_bootstrap, y_bootstrap, alpha = 0.75, nlambda = 100)  # Fit the Lasso model
  age_estimations_lasso_bootstrap[, i] <- as.numeric(predict(lasso_model_bootstrap, s = 0.001, newx = X_question)) # Age estimations Iterations
}

# Average across all bootstrap estimates
merged_data_ratio$Estimated_Age_Lasso_Bootstrap <- rowMeans(age_estimations_lasso_bootstrap)

# Extracts the coefficients of the Lasso model that was trained on the original dataset
lasso_coefficients <- coef(lasso_model, s = 0.001)
print("Lasso Coefficients:")
print(lasso_coefficients)

# Calculates the average of the coefficients from the Bootstrap Lasso model
average_lasso_bootstrap_coefficients <- rowMeans(as.matrix(lasso_model_bootstrap$beta))
print("Average Lasso Bootstrap Coefficients:")
print(average_lasso_bootstrap_coefficients)



## Verification Process
# Establish the data to long format
df_long <- merged_data_ratio %>%
  select(Age, Estimated_Age, Estimated_Age_Lasso_Bootstrap) %>%
  gather(key = "Estimation_Type", value = "Estimated_Age_Value", -Age)

# Calculate R-squared values
r_squared_values <- df_long %>%
  group_by(Estimation_Type) %>%
  summarise(R_squared = summary(lm(Estimated_Age_Value ~ Age))$r.squared)

# Calculate R-squared and RMSE values
metrics_values <- df_long %>%
  group_by(Estimation_Type) %>%
  summarise(
    R_squared = summary(lm(Estimated_Age_Value ~ Age))$r.squared,
    RMSE = sqrt(mean((Estimated_Age_Value - Age)^2))
  )

# Create the scatter plot with facet grid
p <- ggplot(df_long, aes(x = Age, y = Estimated_Age_Value)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") + # y = x linear line for age reference line
  facet_wrap(~ Estimation_Type, scales = "free") +
  ggtitle("Age vs Estimated Age (Lasso and Lasso Bootstrap)") +
  guides(color = FALSE)+ # Remove color legend
  geom_label(fill = "lightblue", color = "black", alpha = 0.8, hjust = 0, vjust = 0, data = metrics_values, 
             aes(label = paste("R-squared = ", round(R_squared, 2), "\nRMSE = ", round(RMSE, 2)), x = 30, y = 60), size = 4) +
  theme(plot.title = element_text(hjust = 0.5),  # Center the plot title
        axis.title.x = element_text(margin = margin(t = 10)),  # Adjust x-axis title margin
        axis.title.y = element_text(margin = margin(r = 10)),  # Adjust y-axis title margin
        strip.text = element_text(size = 12, face = "bold"),  # Customize facet labels
        strip.background = element_blank(),  # Remove facet labels background
        legend.position = "none")  # Remove legend

# Add R-squared and RMSE values
metrics_values$Estimation_Type <- factor(metrics_values$Estimation_Type, levels = unique(df_long$Estimation_Type)) 
print(p) # Display the plot