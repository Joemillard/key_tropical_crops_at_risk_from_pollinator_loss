##Function for standardizing and centering predictor variables on very different scales, e.g. percNH
StdCenterPredictor <- function(x) {
  variable <- x
  sd <- sd(na.omit(variable))
  mean <- mean(na.omit(variable))
  variable.s <- (variable - mean)/sd
  return(variable.s)
}

# build base map for fertiliser/climate plot
get_basemap <- function(){
  
  # download full basemap
  base_map <- getMap(resolution = "high")
  
  # convert to correction projection
  proj4string(base_map) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  
  # return basemap
  return(base_map)
}

# function for predicting continuous variable
predict_continuous <- function(model,
                               model_data, 
                               response_variable,
                               categorical_variable,
                               continuous_variable,
                               continuous_transformation,
                               random_variable,
                               colour_palette){
  
  # set up the prediction dataframe
  prediction_data <- model_data[, c(response_variable, 
                                    random_variable[1], 
                                    random_variable[2], 
                                    random_variable[3], 
                                    categorical_variable[1], 
                                    continuous_variable[1])]
  
  # remove any incomplete rows (NAs) from the prediction data
  prediction_data <- prediction_data[complete.cases(prediction_data),]
  
  # predict the values for the model
  y_value <- c(StatisticalModels::PredictGLMER(model, data = prediction_data, se.fit = TRUE, seMultiplier = 1.96))[[1]]
  y_value_plus <- c(StatisticalModels::PredictGLMER(model, data = prediction_data, se.fit = TRUE, seMultiplier = 1.96))[[2]]
  y_value_minus <- c(StatisticalModels::PredictGLMER(model, data = prediction_data, se.fit = TRUE, seMultiplier = 1.96))[[3]]
  
  # bind the predicted values to the prediction data
  bound_values <- data.frame(cbind(prediction_data,
                                   y_value, 
                                   y_value_plus, 
                                   y_value_minus, 
                                   metric = response_variable,
                                   prediction_data[, continuous_variable[1]]))
  
  # rename last column after transformation
  colnames(bound_values)[ncol(bound_values)] <- paste(continuous_variable[1], "transform", sep = "_")
  
  # rename the response variable column "response_variable" for later plot function
  bound_values <- bound_values %>%
    rename("response_variable" = all_of(response_variable))
  
  # print the final dataframe of predicted values
  print(bound_values)
  
}

# function for plotting the output from a continuous variable glmer - predict_continuous()
# need to amend to read in categorical variable rather than zone/order
plot_fert_response <- function(data_set, categorical_variable){
  plot_obj <- ggplot(data_set) +
    geom_line(aes_string(x = "standard_anom_transform", y = "y_value", colour = categorical_variable), size = 1.5) +
    geom_ribbon(aes_string(x = "standard_anom_transform", ymin = "y_value_minus", ymax = "y_value_plus", fill = categorical_variable), alpha = 0.4) +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 2.39794, 2.69897, 3, 3.39794), labels = c(0.1, 1, 10, 100, 250, 500, 1000, 2500)) +
    theme_bw() +
    theme(panel.grid = element_blank())
  return(plot_obj)
}