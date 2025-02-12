#Testing data processing

#load libraries
require(chillR)
require(dormancyR)
require(pracma)
require(lubridate)
require(tidyr)
require(dplyr)
require(readr)
require(zoo)
require(fruclimadapt)
require(weathermetrics)
require(ggplot2)

trait_data = read_csv("1_Submission_Template_2024.csv")
meta_data = read_csv("2_Testing_Meta_Data_2024.csv")
soil_data = read_csv('3_Testing_Soil_Data_2024.csv')
weather_data = read_csv('4_Testing_Weather_Data_2024_full_year.csv')
geno_data = read_csv('5_Genotype_Data_All_2014_2025_Hybrids_numerical.csv')
ec_data = read_csv('6_Testing_EC_Data_2024.csv')

#meta_data processing
meta_data = meta_data %>%
  rename(lon = `Weather_Station_Longitude (in decimal numbers NOT DMS)`,
         lat = `Weather_Station_Latitude (in decimal numbers NOT DMS)`)

meta_data <- meta_data %>%
  # First, fill missing values within Experiment Code
  group_by(`Experiment_Code`) %>%
  mutate(
    lat = if_else(is.na(lat), first(lat[!is.na(lat)]), lat),
    lon = if_else(is.na(lon), first(lon[!is.na(lon)]), lon)
  ) %>%
  ungroup() %>%
  # Next, fill missing values within City
  group_by(City) %>%
  mutate(
    lat = if_else(is.na(lat), first(lat[!is.na(lat)]), lat),
    lon = if_else(is.na(lon), first(lon[!is.na(lon)]), lon)
  ) %>%
  ungroup()

###Date of harvest prediction
weather_data$Date = ymd(weather_data$Date)
weather_data$lon = meta_data$lon[match(weather_data$Env, meta_data$Env)]
weather_data$lat = meta_data$lat[match(weather_data$Env, meta_data$Env)]
weather_data = filter(weather_data, Env %in% unique(trait_data$Env))
weather_data = filter(weather_data,!is.na(T2M_MIN),!is.na(T2M_MAX))

columns_for_EWMA_REWMA = colnames(weather_data)[-c(1,2,19,20)]
EWMA_REWMA_windows = c(2,3,4,5,6,7,10,14,21,30,45,60,90)
columns_for_cumsum = c('ALLSKY_SFC_SW_DWN','PRECTOTCORR','ALLSKY_SFC_SW_DNI')


##EWMA_REWMA computation

{
  weather_data_EWMA_REWMA = weather_data[,-c(19,20)]
  
  #
  inverse_EWA = function(datalist,window=10){
    datalist_rev = rev(datalist)
    datalist_rev_EWA = movavg(datalist_rev,n = window,type = 'e')
    datalist_EWA = rev(datalist_rev_EWA)
    datalist_EWA = c(rep(NA,window), datalist_EWA)
    datalist_EWA = datalist_EWA[1:length(datalist)]
    return(datalist_EWA)
  }
  
  compute_ewma_and_rewma <- function(df, columns, window_sizes) {
    # Sort by Date within each group
    df <- df %>% arrange(Date)
    
    for (col in columns) {
      for (window in window_sizes) {
        # Compute EWMA
        ewma_col_name <- paste0(col, "_EWMA_", window)
        df[[ewma_col_name]] <- movavg(df[[col]], n = window, type = "e")
        
        # Compute REWMA using inverse_EWA
        rewma_col_name <- paste0(col, "_REWMA_", window)
        df[[rewma_col_name]] <- inverse_EWA(df[[col]], window = window)
      }
    }
    return(df)
  }
  
  weather_data_EWMA_REWMA_list <- split(weather_data_EWMA_REWMA, weather_data_EWMA_REWMA$Env)
  
  # Apply the combined EWMA and REWMA function to each dataframe
  processed_list <- lapply(weather_data_EWMA_REWMA_list, compute_ewma_and_rewma, 
                           columns = columns_for_EWMA_REWMA, 
                           window_sizes = EWMA_REWMA_windows)
  
  # Combine the processed dataframes back into one
  weather_data_EWMA_REWMA <- bind_rows(processed_list)
}

##cumsum computation

{
  
  weather_data_cumsum = weather_data[,c(1,2,which(colnames(weather_data) %in% columns_for_cumsum))]
  weather_data_cumsum <- weather_data_cumsum %>%
    arrange(Env, Date) %>%  # Ensure data is ordered by Env and Date
    group_by(Env) %>%
    mutate(across(
      .cols = all_of(columns_for_cumsum),
      .fns = ~ cumsum(.),
      .names = "{.col}_cumsum"  # Rename columns dynamically
    )) %>%
    ungroup()
  
}


#chilling and heat computation

{
  weather_data_chill_heat = weather_data[,c(1,2,which(colnames(weather_data) %in% c('T2M_MAX','T2M_MIN')),19,20)]
  weather_data_chill_heat = weather_data_chill_heat %>%
    rename(Tmax = T2M_MAX, Tmin = T2M_MIN)

  temperature_feature_generation = function(df){
    
    df <- df %>% 
      arrange(Date)
    df$DOY = yday(df$Date)
    df$Year = year(df$Date)
    df$Month = month(df$Date)
    df$Day = day(df$Date)
    df$Hour = hour(df$Date)
    
    # colnames(df)[which(names(df) == "hourly_temperature_2m")] <- "Temp"
    df = stack_hourly_temps(df,latitude = df$lat[1])[[1]]
    
    #chilling_GDH_GDD_computation
    CU <- chilling_units(df$Temp, summ = F)
    Utah<- modified_utah_model(df$Temp, summ = F)
    NC<-north_carolina_model(df$Temp, summ = F)
    DP <- Dynamic_Model(df$Temp, summ = F)
    GDH_10 <- GDH_linear(df[,-which(colnames(df) %in% c("datetime","Date",'Env'))], Tb = 10, Topt = 25, Tcrit = 36)
    GDH_7 <- GDH_linear(df[,-which(colnames(df) %in% c("datetime","Date",'Env'))], Tb = 7, Topt = 25, Tcrit = 36)
    GDH_4 <- GDH_linear(df[,-which(colnames(df) %in% c("datetime","Date",'Env'))], Tb = 4, Topt = 25, Tcrit = 36)
    GDH_0 <- GDH_linear(df[,-which(colnames(df) %in% c("datetime","Date",'Env'))], Tb = 0, Topt = 25, Tcrit = 36)
    GDD_0 <- GDD(df$Temp,summ = F,Tbase = 0)
    GDD_4 <- GDD(df$Temp,summ = F,Tbase = 4)
    GDD_7 <- GDD(df$Temp,summ = F,Tbase = 7)
    GDD_10 <- GDD(df$Temp,summ = F,Tbase = 10)
    
    CU = if_else(CU < 0, 0, CU)
    Utah = if_else(Utah < 0, 0 ,Utah)
    NC = if_else(NC < 0,0, NC)
    
    All_chilling_data <- data.frame(Date = df$Date,
                                    Month = format(df$Date,format = "%b"),
                                    Year = format(df$Date,format = "%Y"),
                                    CU = CU,
                                    Utah = Utah,
                                    NC = NC,
                                    DP = DP,
                                    GDD_0 = GDD_0,
                                    GDD_4 = GDD_4,
                                    GDD_7 = GDD_7,
                                    GDD_10 = GDD_10)
    
    All_chilling_data <- All_chilling_data[order(All_chilling_data$Date),]
    
    
    
    GDHs = data.frame(Date = as.Date(GDH_10$Date),
                      GDH10 = GDH_10$GDH,
                      GDH_7 = GDH_7$GDH,
                      GDH_4 = GDH_4$GDH,
                      GDH_0 = GDH_0$GDH)
    #delete the last row
    GDHs <- GDHs[-nrow(GDHs),]
    
    Chilling_data_summary_daily <-  All_chilling_data %>%
      group_by(Date) %>%
      summarise(
        CU = sum(CU),
        NC = sum(NC),
        Utah = sum(Utah),
        DP = sum(DP),
        GDD_0 = sum(GDD_0),
        GDD_4 = sum(GDD_4),
        GDD_7 = sum(GDD_7),
        GDD_10 = sum(GDD_10)
      ) %>%
      arrange(Date) %>%
      #delete the last two rows
      slice(1:(n()-2))
    
    Chilling_data_summary_daily <- left_join(Chilling_data_summary_daily,GDHs, by = "Date")
    
    apply_rollsum <- function(data, column_names, window_lengths) {
      for (column in column_names) {
        for (window in window_lengths) {
          new_column_name <- paste0(column, "_", window, "days")
          data <- data %>%
            mutate(!!new_column_name := rollsum(get(column), window, fill = NA, align = "right"))
        }
      }
      return(data)
    }
    
    # Define a function to compute rolling sums for the dormant season (Sept to May)
    apply_rollsum_dormant <- function(data, column_names) {
      data <- data %>%
        filter(format(Date, "%m") %in% c("09", "10", "11", "12"))  # Select Sept to May months
      for (column in column_names) {
        new_column_name <- paste0(column, "_dormant_sum")
        data <- data %>%
          group_by(season = ifelse(month(Date) >= 9, year(Date), year(Date) - 1)) %>%
          mutate(!!new_column_name := cumsum(get(column))) %>%
          ungroup()
      }
      return(data)
    }
    
    # Define a function to compute cumulative sum from the start of each year
    apply_cumsum_by_year <- function(data, column_names) {
      data <- data %>%
        mutate(year = year(Date)) %>%  # Extract year from the date column
        group_by(year)  # Group by year
      for (column in column_names) {
        new_column_name <- paste0(column, "_cumsum_year")
        data <- data %>%
          mutate(!!new_column_name := cumsum(get(column)))
      }
      data <- data %>% ungroup()  # Ungroup after operation
      return(data)
    }
    
    # Define the columns and window lengths
    columns_for_rollsum <- c("CU", "NC", "Utah", "DP", "GDD_0", "GDD_4", "GDD_7", "GDD_10", "GDH10", "GDH_7", "GDH_4", "GDH_0")
    dormant_columns <- c("CU", "NC", "Utah", "DP")
    columns_for_cumsum_year <- c("GDD_0", "GDD_4", "GDD_7", "GDD_10", "GDH10", "GDH_7", "GDH_4", "GDH_0")
    window_lengths <- c(3, 7, 14, 30, 60, 90)
    
    # Apply the rolling sum for all the columns
    Chilling_data_summary_daily <- apply_rollsum(Chilling_data_summary_daily, columns_for_rollsum, window_lengths)
    
    # Apply the cumulative sum from the start of each year for GDD and GDH columns
    Chilling_data_summary_daily <- apply_cumsum_by_year(Chilling_data_summary_daily, columns_for_cumsum_year)
    
    # Apply the rolling sum for the dormant season (Sept to May)
    Chilling_data_summary_daily_dormant <- apply_rollsum_dormant(Chilling_data_summary_daily, dormant_columns)
    Chilling_data_summary_daily_dormant = Chilling_data_summary_daily_dormant %>% 
      dplyr::select(Date, CU_dormant_sum, NC_dormant_sum, Utah_dormant_sum, DP_dormant_sum)
    Chilling_data_summary_daily = left_join(Chilling_data_summary_daily, Chilling_data_summary_daily_dormant, by = "Date")
    
    #Output

    return(Chilling_data_summary_daily)
    
  }
  
  weather_data_chill_heat_split <- split(weather_data_chill_heat, weather_data_chill_heat$Env)
  
  # Apply the temperature_feature_generation function to each dataframe
  weather_data_chill_heat_split <- lapply(names(weather_data_chill_heat_split), function(env) {
    df <- weather_data_chill_heat_split[[env]]
    df <- temperature_feature_generation(df)  # Apply your function
    df$Env <- env  # Add the `Env` column back to the dataframe
    return(df)
  })
  
  # Combine the processed dataframes back into a single dataframe
  final_weather_data_chill_heat <- bind_rows(weather_data_chill_heat_split)
  
}

weather_data_assembled = left_join(weather_data,
                                   weather_data_EWMA_REWMA %>%
                                     dplyr::select(-all_of(columns_for_EWMA_REWMA)), 
                                   by = c("Date","Env"))

weather_data_assembled = left_join(weather_data_assembled,
                                   weather_data_cumsum %>%
                                     dplyr::select(-all_of(columns_for_cumsum)), 
                                   by = c("Date","Env"))

weather_data_assembled = left_join(weather_data_assembled,
                                   final_weather_data_chill_heat, 
                                   by = c("Date","Env"))


weather_data_assembled$plant_date = mdy(meta_data$Date_Planted[match(weather_data_assembled$Env,meta_data$Env)])
weather_data_assembled = filter(weather_data_assembled, !is.na(plant_date))
weather_data_assembled = weather_data_assembled %>%
  group_by(Env) %>%
  mutate(season_GDD_0_cumsum = GDD_0_cumsum_year - GDD_0_cumsum_year[which(Date == plant_date)],
         season_GDD_4_cumsum = GDD_4_cumsum_year - GDD_4_cumsum_year[which(Date == plant_date)],
         season_GDD_7_cumsum = GDD_7_cumsum_year - GDD_7_cumsum_year[which(Date == plant_date)],
         season_GDD_10_cumsum = GDD_10_cumsum_year - GDD_10_cumsum_year[which(Date == plant_date)],
         season_GDH_0_cumsum = GDH_0_cumsum_year - GDH_0_cumsum_year[which(Date == plant_date)],
         season_GDH_4_cumsum = GDH_4_cumsum_year - GDH_4_cumsum_year[which(Date == plant_date)],
         season_GDH_7_cumsum = GDH_7_cumsum_year - GDH_7_cumsum_year[which(Date == plant_date)],
         season_GDH_10_cumsum = GDH10_cumsum_year - GDH10_cumsum_year[which(Date == plant_date)]
  )

weather_data_assembled_harvest_testing <-  weather_data_assembled %>%
  group_by(Env) %>% 
  filter(Date >= plant_date[1]) %>%
  ungroup() %>%
  dplyr::select(-c('year','plant_date'))


write_csv(weather_data_assembled_harvest_testing, file = 'harvest_data_test.csv')