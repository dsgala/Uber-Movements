###### Set-up
libraries <- c("tidyverse", "lubridate", "ggplot2", "ggfortify", "fpp",
               "stringr", "forecast", "scales", "zoo", "imputeTS", "vars",
               "MuMIn", "keras", "TSLSTM")
load <- lapply(libraries, require, character.only = TRUE)


###### 1. Load Data
get_file_names <- function(data_path){
  file_names <- list.files(path = data_path, pattern = "*.csv")
  file_times <- file.mtime(file.path(data_path, file_names))
  file_names_sorted <- file_names[order(file_times, decreasing = FALSE)]
}

load_data <- function(data_path){
  destination_dirs <- list.files(path = data_path)
  df_full <- data.frame()
  for (destination in destination_dirs){
    dest_path = paste0(data_path,"/",destination)
    file_list <- get_file_names(dest_path)
    for (file in file_list){
      tryCatch({
        temp <- read.csv(paste0(dest_path,"/",file))
        re <- "(?<=\\.\\s)\\D+(?=\\s*\\()|(?<=\\.\\s)\\D+$"
        temp$dest_name <- rep(trimws(str_extract_all(destination, re)),nrow(temp))
        df_full <- rbind(df_full, temp)
      }, error = function(e) {
        cat(paste("Warning:", "No data present in this file, continuing..."), "\n")
      })
    }
  }
  df <- df_full[,c(1,2,4,6,ncol(df_full))]
  names(df) <- c("date", "origin", "destination", "mean_time", "neighborhood")
  df$date <- as.Date(df$date, format = "%m/%d/%Y")
  df_list <- split(df, df$neighborhood)
  df_list_sorted <- lapply(df_list, function(x) x[order(x$date),])
  return(df_list_sorted)
}

###### 2. Preliminary Data Cleaning
# Apply weekly averaging
get_weekly_avg <- function(route){
  route$week <- floor_date(route$date, "week")
  weekly_df <- route %>%
    group_by(week) %>%
    summarize(mean_travel_time = mean(mean_time, na.rm = TRUE))
  return(weekly_df)
}

# Filter routes that have a high number of missing values
missing_breach <- function(weekly_df, threshold = 0.15){
  prop_missing = sum(is.na(weekly_df$mean_travel_time))/nrow(weekly_df)
  cat(paste0("Missing Values (%): ", round(prop_missing,4)*100, "%", "\n"))
  if (prop_missing > threshold){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Apply some basic cleaning
clean_df <- function(weekly_df){
  # 1. Subset series to end before March 2020 to avoid covid impact
  end <- "2020-02-23"
  weekly_df <- weekly_df[weekly_df$week <= end, ]
  
  # 2. Convert to data to minutes from seconds
  weekly_df$mean_travel_time <- weekly_df$mean_travel_time/60
  
  # # 3. Impute missing values (if necessary) using last value carried forward
  # cat(paste0("Missing Values (%): ", round(sum(is.na(weekly_df$mean_travel_time))/nrow(weekly_df),4)*100,"%","\n"))
  data_missing <- sum(is.na(weekly_df$mean_travel_time)) > 0
  weekly_df$mean_travel_time <- if(data_missing){na_locf(weekly_df$mean_travel_time)} else {weekly_df$mean_travel_time}
  cat("Imputing missing values using LOCF, if necessary...\n")
  cat(paste("Confirm No Missing Values:", sum(is.na(weekly_df$mean_travel_time))==0), "\n")
  return(weekly_df)
}

# Plot series
plot_series <- function(weekly_df){
  p <- ggplot(weekly_df) + 
    geom_line(aes(x = week, y = mean_travel_time)) +
    scale_x_date(labels = function(x) format(x, "%b-%Y")) + 
    ggtitle(paste0("Weekly Average Travel Times from LAX (734) to ",
                   unique(route))) +
    xlab("") + 
    ylab("Travel Time (Minutes)") +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 9))
  print(p)
}

###### 3. Diagnosing Stationarity
# Check for stationarity
check_stationarity <- function(weekly_df, adf_thresh = 0.05){
  
  is_stationary <- (adf.test(weekly_df$mean_travel_time)$p.value < adf_thresh)
  return(is_stationary)
}

###### 4. Modeling
model <- function(df, features, stationary){
  ts_data <- ts(df$mean_travel_time,
                start = c(year(df$week[1]), as.numeric(format(as.Date(df$week[1]), "%U"))),
                frequency = 52)
  
  # ARIMA
  cat("Running ARIMA...\n")
  arima_fit <- auto.arima(ts_data)
  arima_order <- arimaorder(arima_fit)
  arima_spec <- paste0("ARIMA(",arima_order[1], ",", arima_order[2], ",", arima_order[3],")",
                 "(", arima_order[4], ",", arima_order[5], ",", arima_order[6], ")",
                 "[", arima_order[7], "]")
  arima_df <- data.frame(`Model Type` = "ARIMA",
                         "Specification" = arima_spec,
                         AICc = arima_fit$aicc,
                         AIC = arima_fit$aic, 
                         BIC = arima_fit$bic,
                         RMSE = accuracy(arima_fit)[,'RMSE'],
                         correlated_residuals = as.numeric(checkresiduals(arima_fit, plot = FALSE)$p.value < 0.05))
  
  # Exponential Smoothing
  cat("Running Exponential Smoothing...\n")
  ets_fit <- ets(df$mean_travel_time)
  ets_spec <- ets_fit$method
  ets_df <- data.frame(`Model Type` = "Exponential Smoothing",
                       "Specification" = ets_spec,
                       AICc = ets_fit$aicc,
                       AIC = ets_fit$aic,
                       BIC = ets_fit$bic,
                       RMSE = accuracy(ets_fit)[,'RMSE'],
                       correlated_residuals = as.numeric(checkresiduals(ets_fit, plot = FALSE)$p.value < 0.05))

  # Linear Regression
  cat("Running Linear Regression...\n")
  df <- if(stationary){df} else {data.frame("week" = df$week[2:nrow(df)], "mean_travel_time" = diff(df$mean_travel_time))}
  regression_df <- left_join(df, features, by = 'week')
  full_lm <- lm(mean_travel_time~.-week, data = regression_df)
  best_lm <- step(full_lm, direction = 'both', trace = FALSE)
  vars_to_remove <- c("week", "mean_travel_time", "regression_df")
  chosen_vars <- all.vars(best_lm$call)[!all.vars(best_lm$call) %in% vars_to_remove]
  chosen_vars_read <- paste("X =", paste(chosen_vars, collapse = ", "))
  lm_df <- data.frame(`Model Type` = "Linear Regression",
                      "Specification" = chosen_vars_read,
                      AICc = AICc(best_lm),
                      AIC = AIC(best_lm),
                      BIC = BIC(best_lm),
                      RMSE = accuracy(best_lm)[,'RMSE'],
                      correlated_residuals = as.numeric(checkresiduals(best_lm, plot = FALSE)$p.value < 0.05))

  # Linear Regression with Arima Errors
  regression_ts <- ts(regression_df,
                      start = c(year(regression_df$week[1]),
                                as.numeric(format(as.Date(regression_df$week[1]), "%U"))),
                      frequency = 52)
  cat("Running Linear Regression with Arima Errors...\n")
  arimax_fit <- auto.arima(regression_ts[,'mean_travel_time'],
                           xreg = regression_ts[,chosen_vars])
  arimax_order <- arimaorder(arimax_fit)
  arimax_spec <- paste0("ARIMA(",arimax_order[1], ",", arimax_order[2], ",", arimax_order[3],")",
                        "(", arimax_order[4], ",", arimax_order[5], ",", arimax_order[6], ")",
                        "[", arimax_order[7], "]")
  arimax_df <- data.frame(`Model Type` = "Linear Regression with Arima Errors",
                          "Specification" = paste0(chosen_vars_read, " (Errors: ",arimax_spec, ")"),
                          AICc = AICc(arimax_fit),
                          AIC = AIC(arimax_fit),
                          BIC = BIC(arimax_fit),
                          RMSE = accuracy(arimax_fit)[,'RMSE'],
                          correlated_residuals = as.numeric(checkresiduals(arimax_fit, plot = FALSE)$p.value < 0.05))

  # Spectral Analysis
  tbats_fit <- tbats(ts_data)
  tbats_spec <- paste0("TBATS(",
                       ifelse(is.null(tbats_fit$lambda), 1, tbats_fit$lambda),
                       ",{", tbats_fit$parameters$control$p, ",",
                       tbats_fit$parameters$control$q, "},",
                       ifelse(is.null(tbats_fit$damping.parameter), "-", round(tbats_fit$damping.parameter,2)), ",",
                       ifelse(is.null(tbats_fit$seasonal.periods), "NULL", tbats_fit$seasonal.periods), ")")
  tbats_df <- data.frame(`Model Type` = "TBATS",
                         "Specification" = tbats_spec,
                         AICc = NA,
                         AIC = tbats_fit$AIC,
                         BIC = NA,
                         RMSE = accuracy(tbats_fit)[,'RMSE'],
                         correlated_residuals = as.numeric(checkresiduals(tbats_fit, plot = FALSE)$p.value < 0.05))

  model_summary_df <- rbind(arima_df, ets_df, lm_df, arimax_df, tbats_df)
  return(model_summary_df)
}

###### 5. Execute Pipeline
data_path <- paste0(getwd(),"/Data")
df_list_all <- load_data(data_path)

features <- read.csv("features.csv")
features$week <- as.Date(features$week, format = "%m/%d/%y")

df_list <- df_list_all
neighborhoods <- names(df_list)
models_df <- rbind()
for (route in neighborhoods){
  cat("***************************************************\n")
  cat(paste0("Building models for LAX to ", route, "..."), "\n")
  weekly_df <- get_weekly_avg(df_list[[route]])
  if(missing_breach(weekly_df)){
    cat(paste0("Removing ", route, " from analysis, as too many missing values\n"))
    next
  }
  else {
    weekly_df <- clean_df(weekly_df)
    # plot_series(weekly_df)
    is_stationary <- check_stationarity(weekly_df)
    cat(paste("Series is stationary:",is_stationary), "\n")
    temp <- cbind(Destination = route,
                  Stationary = is_stationary,
                  model(weekly_df, features = features, stationary = is_stationary))
    models_df <- rbind(models_df, temp)
  }
}
# write.csv(models_df_2, "models_df_2.csv")

###### 6. Calculating Test Error on Candidate Models (Cross-validation)

#### Arima Candidate Models
arima_1 <- function(ts){Arima(ts, order = c(0,1,1), seasonal = c(1,1,1))}
arima_2 <- function(ts){Arima(ts, order = c(0,1,1), seasonal = c(1,1,0))}
arima_3 <- function(ts){Arima(ts, order = c(0,1,2), seasonal = c(1,0,0))}
candidate_models_arima <- c("ARIMA(0,1,1)(1,1,1)[52]" = arima_1,
                            "ARIMA(0,1,1)(1,1,0)[52]" = arima_2)

arima_cv <- function(df, model, h = 12){
  ts_data <- ts(df$mean_travel_time,
                start = c(year(df$week[1]), as.numeric(format(as.Date(df$week[1]), "%U"))),
                frequency = 52)
  fit <- model(ts_data)
  ljung_box <- checkresiduals(fit)$p.value < 0.05
  arima_cv <- tsCV(ts_data, function(x, h){forecast(fit, h=h)}, h = h)
  arima_rmse <- sqrt(colMeans(arima_cv^2, na.rm = TRUE))
  return(c('train_rmse' = accuracy(fit)[,'RMSE'],
           'correlated_residuals' = ljung_box,
           'aic' = fit$aic,
           arima_rmse))
}

#### Arimax Candidate Models
arimax_1 <- function(ts){Arima(ts[,'mean_travel_time'], order = c(0,0,1), seasonal = c(1,0,0),
                               xreg = ts[,c("holiday")])}
arimax_2 <- function(ts){Arima(ts[,'mean_travel_time'], order = c(1,0,1), seasonal = c(1,0,0),
                               xreg = ts[,c("holiday")])}
arimax_3 <- function(ts){Arima(ts[,'mean_travel_time'], order = c(0,0,1), seasonal = c(1,1,1),
                               xreg = ts[,c("holiday")])}
candidate_models_arimax <- c("X = holiday (Errors: ARIMA(0,0,1)(1,0,0)[52])" = arimax_1,
                             "X = holiday (Errors: ARIMA(1,0,1)(1,0,0)[52])" = arimax_2)

arimax_cv <- function(df, features, model, h = 12){
  regression_df <- left_join(df, features, by = 'week')
  regression_ts <- ts(regression_df,
                      start = c(year(regression_df$week[1]),
                                as.numeric(format(as.Date(regression_df$week[1]), "%U"))),
                      frequency = 52)
  fit <- model(regression_ts)
  check_resid <- checkresiduals(fit, plot = FALSE)
  ljung_box <- check_resid$p.value < 0.05
  arimax_cv <- tsCV(regression_ts[,'mean_travel_time'],
                    function(x, h){forecast(fit, h=h, xreg = regression_ts[,'holiday'])}, h = h)
  arimax_rmse <- sqrt(colMeans(arimax_cv^2, na.rm = TRUE))
  return(c('train_rmse' = accuracy(fit)[,'RMSE'],
           'correlated_residuals' = ljung_box,
           'aic' = fit$aic,
           arimax_rmse))
}
candidate_models <- c(candidate_models_arima, candidate_models_arimax)
model_labels = c("ARIMA(0,1,1)(1,1,1)[52]",
                 "X = holiday (Errors: ARIMA(0,0,1)(1,0,0)[52])",
                 "ARIMA(0,1,1)(1,1,0)[52]", 
                 "X = holiday (Errors: ARIMA(1,0,1)(1,0,0)[52])")

#### Run Cross-validation
neighborhoods <- unique(read.csv("models_df.csv")$Destination)
cv_res_df <- c()
for (route in neighborhoods){
  cat("***************************************************\n")
  cat(paste0("Testing candidate models for LAX to ", route, "..."), "\n")
  weekly_df <- get_weekly_avg(df_list[[route]])
  weekly_df <- clean_df(weekly_df)
  cv_res <- c()
  for (i in 1:2){
    cat(paste0("Running cross validation for: ", names(candidate_models_arima[i])), "\n")
    arima_cv_rmse <- arima_cv(df = weekly_df, model = candidate_models_arima[[i]], h = 12)
    cat(paste0("Running cross validation for: ", names(candidate_models_arimax[i])), "\n")
    arimax_cv_rmse <- arimax_cv(df = weekly_df, features = features, model = candidate_models_arimax[[i]], h = 12)
    
    cv_res <- rbind(cv_res, arima_cv_rmse, arimax_cv_rmse)
  }
  res <- cbind(Model = model_labels,
               Destination = rep(route, 4),
               cv_res)
  cv_res_df <- rbind(cv_res_df, res)
}
cv_res_df <- data.frame(cv_res_df)
row.names(cv_res_df) <- NULL
# write.csv(cv_res_df, "cv_arima.csv")

###### 7. Appendix

# Plot of all routes
library(tidyr)
full_df_weekly <- data.frame('week' = as.Date(weekly_df$week))
for (route in neighborhoods){
  weekly_df <- get_weekly_avg(df_list[[route]])
  weekly_df <- clean_df(weekly_df)
  names(weekly_df)[2] <- route
  full_df_weekly <- left_join(full_df_weekly, weekly_df, by = 'week')
}

existing_colors <- c("#5d0000", "#800000", "#767676", "#D6D6CE")
additional_colors <- colorRampPalette(c("#5d0000", "#800000", "#767676", "#D6D6CE"))(46)
color_palette <- c(existing_colors, additional_colors)
df_long <- full_df_weekly %>%
  gather(key = "route", value = "time", -week)
p <- ggplot(data = df_long, aes(x = week, y = time, color = route)) +
  geom_line(alpha = 0.5) +
  labs(title = "Average Weekly Travel Time between LAX and 50 different locations in LA (Jan 2016 - Feb 2020)",
       x = "",
       y = "Travel Time (minutes)") +
  scale_color_manual(values = color_palette) +
  theme_minimal() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 12, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = 10, face = 'bold'))
p
# Save the plot
# ggsave(filename = "time_series_plot.png", plot = p, dpi = 600, width = 8, height = 5.5)

## Get Average Travel Times Across the 50 Routes
neighborhoods <- names(df_list)
all_times <- data.frame('week' = as.Date(weekly_df$week))
for (route in neighborhoods){
  cat("***************************************************\n")
  weekly_df <- get_weekly_avg(df_list[[route]])
  if(missing_breach(weekly_df)){
    cat(paste0("Removing ", route, " from analysis, as too many missing values\n"))
    next
  }
  else {
    weekly_df <- clean_df(weekly_df)
    names(weekly_df)[2] <- route
    all_times <- left_join(all_times, weekly_df, by = 'week')
  }
}
avgs <- rowMeans(all_times[,-1], na.rm = TRUE)
# write.csv(data.frame(all_times$week, avgs), "avgs.csv")


# Seasonal Naive Cross-Validation for Comparison
snaive_cv <- function(df, h = 12){
  ts_data <- ts(df$mean_travel_time,
                start = c(year(df$week[1]), as.numeric(format(as.Date(df$week[1]), "%U"))),
                frequency = 52)
  snaive_cv <- tsCV(ts_data, snaive, h = h)
  snaive_rmse <- sqrt(colMeans(snaive_cv^2, na.rm = TRUE))
  return(c('train_rmse'= accuracy(snaive(ts_data, h = h))[,'RMSE'],
           'test_rmse' = snaive_rmse))
}


cv_res_df_snaive <- c()
for (route in neighborhoods){
  cat("***************************************************\n")
  cat(paste0("Testing candidate models for LAX to ", route, "..."), "\n")
  weekly_df <- get_weekly_avg(df_list[[route]])
  weekly_df <- clean_df(weekly_df)
  snaive_cv_rmse <- snaive_cv(weekly_df, h = 12)
  cv_res_df_snaive <- rbind(cv_res_df_snaive, snaive_cv_rmse)
}
cv_res_df_snaive <- data.frame(cv_res_df_snaive)
row.names(cv_res_df_snaive) <- NULL
write.csv(cv_res_df_snaive, "snaive_cv.csv")
plot(colMeans(cv_res_df_snaive))

mod = snaive(weekly_df$mean_travel_time, h = 12)
accuracy(mod)
