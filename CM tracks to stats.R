rm(list=ls())
library(geosphere)
library(dplyr)
library(tidyr)
library(circular)

directory <- "path/CM_trimmed/" #Cleaned data directory
meta_data <- read.csv("path/CM info.csv") #Metadata
save_path <- "path/"


siteSH_coords <- c(-1.42375, 51.7527778) # Release site coords
siteLH_coords <- c(-1.3915, 51.8303056) # Release site coords
home_coords <- c(-1.3173753, 51.7828602)

csv_files <- paste(directory, list.files(directory), sep="")
track_files <- csv_files[grep("-", csv_files)]


full_meta_data <- meta_data %>%
  mutate(bird=Bird, partner=Partner, site=Additional_training_site, add_train=T) %>%
  bind_rows(meta_data %>% mutate(bird=Bird, partner=Partner, site=ifelse(Additional_training_site == "LH", "SH", "LH"),
                                 add_train=F))%>%
  bind_rows(meta_data %>% mutate(bird=Partner, partner=Bird, site=Additional_training_site,
                                 add_train=T))%>%
  bind_rows(meta_data %>% mutate(bird=Partner, partner=Bird, site=ifelse(Additional_training_site == "LH", "SH", "LH"),
                                 add_train=F))%>%
  select(-Bird, -Partner, -Additional_training_site)%>%
  arrange(Pair_ID)
full_meta_data

#Loading data and extracting info
gps_list <- lapply(track_files, function(file, full_meta_data) {
  print(file)
  # Read each CSV file
  data <- read.csv(file)
  # Ensure that time is in POSIXct format
  data$DT <- as.POSIXct(data$DT, tz="UTC")
  data$Course <- data$Course - 180
  ind_50_first <- which(data$DistRS>50)[1] #Leaving 50m of release site for first time
  ind_500h_first <- which(data$DistHome<500)[1] #Entering 500m of home site for first time
  ind_500h_first <- ifelse(is.na(ind_500h_first), nrow(data), ind_500h_first)
  data <- data[c(ind_50_first
                 :ind_500h_first),]
  data <- data[data$Speed.km.h.> 30,] # Put in speed filter of 30
  split <- strsplit(file, "[/_ .]")[[1]]
  data$bird <- split[14]
  data$treatment <- split[19] #Treatment here just refers to training or memory testing
  data$release_cond <- split[17]
  data$round <- NULL
  data$date <- split[20]
  data$time <- split[21]
  return(data)
})

gps_list2 <- lapply(gps_list, function(df) {
  merge(df, full_meta_data, by = c("bird", "site"), all.x=T, all.y=F)
})
gps_list <- NULL

results_df  <- do.call(rbind, lapply(gps_list2, function(df) {
  df[1, c("bird", "Pair_ID", "partner", "site", "treatment", "release_cond", 
          "date", "time", "add_train"),
     drop = FALSE]  # Extract the first row with specified columns
}))
results_df$dt1 <- as.POSIXct(paste(results_df$date, results_df$time), format="%Y-%m-%d %H:%M:%OS", tz="UTC") #Datetime
#Pair* treatment relates to a small number of instances in which the birds were meant to be released solo, but the two individuals in the pair joined up
# -> reclassified into the pair treatment
results_df$release_cond <- ifelse(results_df$release_cond=="pair*", "pair", results_df$release_cond)

results_df <- results_df %>%
  # Save the original row order
  mutate(original_order = row_number()) %>%
  # Group by 'bird', 'site', 'treatment', and arrange by 'dt1' (date-time of release)
  #Treatment here just refers to training vs memory testing
  group_by(bird, site, treatment) %>%
  arrange(bird, site, treatment, dt1) %>% 
  # Mutate to add the order and total_n columns
  mutate(
    order = row_number(),        # Sequential order within group
    total_n = n(),               # Total rows in each group
    order2 = ifelse(treatment != "CM", row_number() - n() - 1, NA) #Order2 is reverse order of training releases
  ) %>%
  # Ungroup and reorder based on the original row order
  ungroup() %>%
  arrange(original_order) %>%   # Restore original order
  select(-original_order)      # Remove the `original_order` column


results_df <- results_df %>%
  mutate(original_order = row_number()) %>%
  group_by(bird, site, treatment, release_cond) %>% #Release_cond is pair vs solo release
  arrange(bird, site, treatment, release_cond, dt1) %>%
  mutate(
    order3 = ifelse(treatment != "CM", row_number() - n() - 1, NA) #order3 is reverse order of releases by release_cond
  ) %>%
  ungroup() %>%
  arrange(original_order) %>%   # Restore original order
  select(-original_order)      # Optionally remove the `original_order` column


results_df$list_n <- c(1:nrow(results_df))
results_df$date <- as.POSIXct(results_df$date, tz="UTC")
results_df$treatment2 <- as.factor(ifelse(results_df$date > "2024-07-01" & 
                                            results_df$date < "2024-07-11",
                                          "ATR", as.character(results_df$treatment)))
#ATR -> Additional training / extra training treatment

#Function to find indices of tracks of partner for each track
assign_partner_list <- function(df) {
  # Create an empty column
  df$partner_list_n <- NA
  
  # Loop through each row
  for (i in 1:nrow(df)) {
    print(i)
    # Find the matching indices
    if (df$release_cond[i]!="solo"){
      match_index_partner <- which(df$Pair_ID == df$Pair_ID[i] & 
                                     df$site == df$site[i] & 
                                     df$bird == df$partner[i] &
                                     df$date == df$date[i] &
                                     df$release_cond!="solo")
      partner_matches <- results_df[match_index_partner,]
      partner_matches_order <- partner_matches[order(partner_matches$dt1),]
      
      
      match_index_bird <- which(df$Pair_ID == df$Pair_ID[i] & 
                                  df$site == df$site[i] & 
                                  df$bird == df$bird[i] &
                                  df$date == df$date[i] &
                                  df$release_cond!="solo")
      bird_matches <- results_df[match_index_bird,]
      bird_matches_order <- bird_matches[order(bird_matches$dt1),]
      bird_match_match <- which(bird_matches_order$dt1==results_df$dt1[i])
      partner_match_match <- partner_matches_order[bird_match_match,]
      
      if (length(match_index_bird)==length(match_index_partner)){
        df$partner_list_n[i] <- partner_match_match$list_n
      }else{
        df$partner_list_n[i] <- "CHECK" #Partner not found
      }
      
    }
  }
  
  return(df)  # Return the modified data frame
}
results_df <- assign_partner_list(results_df)

results_df$partner_list_n

results_df[which(results_df$partner_list_n=="CHECK"),]
#1 -> 20-C02 data excluded as it joined with another bird
#2-3 -> missing data from 20-C05 from 2024-07-04 (DEVICE LOST)
#4 - lost bird
results_df$partner_list_n[which(results_df$partner_list_n=="CHECK")] <- NA
results_df$partner_list_n <- as.numeric(results_df$partner_list_n)

#Assign order2 and order3 values to the memory testing
#In the analysis, only the first memory testing release was used (so where order2 and order3 =1)
results_df$order2[results_df$treatment == "CM" &
                    results_df$site=="SH" & results_df$dt1<as.POSIXct("2024-08-17 12:00:00", tz="UTC")] <- 1
results_df$order2[results_df$treatment == "CM" &
                    results_df$site=="SH" & results_df$dt1>as.POSIXct("2024-08-17 12:00:00", tz="UTC")] <- 2
results_df$order2[results_df$treatment == "CM" &
                    results_df$site=="LH" & results_df$dt1<as.POSIXct("2024-08-16 12:00:00", tz="UTC")] <- 1
results_df$order2[results_df$treatment == "CM" &
                    results_df$site=="LH" & results_df$dt1>as.POSIXct("2024-08-16 12:00:00", tz="UTC")] <- 2

results_df$order3[results_df$treatment == "CM"] <- results_df$order2[results_df$treatment == "CM"]

#Function gets distances between tracks of a pair
get_gps_distances_between_tracks <- function(gps_list, ind1, ind2, N, split_threshold_in) {
  print(ind1)
  if (is.na(ind2)) {
    updated_track1 <- gps_list[[ind1]]
    updated_track1$Dist_partner <- NA
    updated_track1$split <- NA
    updated_track1$Dist_ahead <- NA
    updated_track1$Dist_sideways <- NA
    updated_track1$mean_pair_course <- NA
    updated_track1$course_diff <- NA
    updated_track1$partner_list_ind <- NA
    updated_track1$track_ind <- ind1
    return(updated_track1)
  }
  
  # Extract the two tracks
  track1 <- gps_list[[ind1]]
  track2 <- gps_list[[ind2]]
  
  if (track1$release_cond[1] == "solo" | track1$release_cond[1] == "group"){
    updated_track1 <- gps_list[[ind1]]
    updated_track1$Dist_partner <- NA
    updated_track1$split <- NA
    updated_track1$Dist_ahead <- NA
    updated_track1$Dist_sideways <- NA
    updated_track1$mean_pair_course <- NA
    updated_track1$course_diff <- NA
    updated_track1$partner_list_ind <- NA
    updated_track1$track_ind <- ind1
    return(updated_track1)
  }
  
  DTs <- c(track1$DT, track2$DT)
  
  all_times <- do.call(c, lapply(gps_list, function(gps) {
    seq(from = min(DTs, na.rm = TRUE), to = max(DTs, na.rm = TRUE), by = N)
  }))
  
  all_times <- data.frame(time = unique(all_times))  # Ensure unique timestamps
  
  # Join track1 and track2 based on the time
  combined <- merge(all_times, track1, by.x = "time", by.y = "DT", all = FALSE)
  combined <- merge(combined, track2, by.x = "time", by.y = "DT", all = FALSE, 
                    suffixes = c("_track1", "_track2"))
  
  if (nrow(combined)==0){
    updated_track1 <- gps_list[[ind1]]
    updated_track1$Dist_partner <- NA
    updated_track1$split <- NA
    updated_track1$Dist_ahead <- NA
    updated_track1$Dist_sideways <- NA
    updated_track1$mean_pair_course <- NA
    updated_track1$course_diff <- NA
    updated_track1$partner_list_ind <- NA
    updated_track1$track_ind <- ind1
    return(updated_track1)
  }
  
  # Initialize Dist_partner columns
  combined$Dist_partner <- NA
  
  
  # Calculate distances only if there are matching rows
  if (nrow(combined) > 0) {
    combined$Dist_partner <- distHaversine(
      cbind(combined$Longitude_track1, combined$Latitude_track1),
      cbind(combined$Longitude_track2, combined$Latitude_track2)
    )
  }
  
  combined$Dist_ahead <- NA
  
  inds_together <- which(combined$Dist_partner<split_threshold_in & 
                           combined$Speed.km.h._track1 >30 & 
                           combined$Speed.km.h._track2 >30)
  last_ind_together <- inds_together[length(inds_together)]
  
  inds_apart <- which(combined$Dist_partner>=split_threshold_in)
  inds_split <- inds_apart[which(inds_apart>last_ind_together)]
  
  combined$split <- NA
  combined$split[1:last_ind_together] <- F
  if (length(inds_split)>0){
    combined$split[(last_ind_together+1):nrow(combined)] <- T
  }
  
  combined$Course_track1 <- circular(combined$Course_track1, units='degrees')
  combined$Course_track2 <- circular(combined$Course_track2, units='degrees')
  
  combined$mean_pair_course <- apply(data.frame(combined$Course_track1, combined$Course_track2), 
                                     1, function(row) {
                                       # Calculate the circular mean for the row (averaging the angles)
                                       row <- circular(row, units='degrees')
                                       return(mean.circular(row))
                                     })
  course_diff <- combined$Course_track1 - combined$Course_track2
  course_diff <- ifelse(course_diff > 180, -360 + course_diff, course_diff)
  combined$course_diff <- ifelse(course_diff < -180, 360 + course_diff, course_diff)
  
  
  bearing_between_birds <- deg2rad(bearing(cbind(combined$Longitude_track1, combined$Latitude_track1),
                                           cbind(combined$Longitude_track2, combined$Latitude_track2)))
  mean_pair_course_rad <- deg2rad(combined$mean_pair_course)
  
  
  delta_bearing <- abs(mean_pair_course_rad - bearing_between_birds)
  
  # Ensure the angle difference is within 0 to pi (0 to 180 degrees)
  delta_bearing <- ifelse(delta_bearing > pi, 2 * pi - delta_bearing, delta_bearing)
  
  combined$Dist_ahead <- combined$Dist_partner*cos(delta_bearing)
  
  
  delta_bearing2 <- abs(mean_pair_course_rad + pi/2 - bearing_between_birds)
  delta_bearing2 <- ifelse(delta_bearing2 > pi, 2 * pi - delta_bearing2, delta_bearing2)
  combined$Dist_sideways <- combined$Dist_partner*cos(delta_bearing2)
  
  first_common_time <- min(combined$time)
  
  # Create updated data frames for each individual
  updated_track1 <- merge(track1, combined[, c("time","Dist_partner", "split", "Dist_ahead", "Dist_sideways",
                                               "mean_pair_course", "course_diff")], 
                          by.x = "DT", by.y = "time", all.x = T)
  updated_track2 <- merge(track2, combined[, c("time","Dist_partner", "split", "Dist_ahead", "Dist_sideways",
                                               "mean_pair_course", "course_diff")], 
                          by.x = "DT", by.y = "time", all.x = T)
  
  updated_track1 <- updated_track1[updated_track1$DT >= first_common_time, ]  
  updated_track2 <- updated_track2[updated_track2$DT >= first_common_time, ]
  
  updated_track1$partner_list_ind <- ind2
  updated_track1$track_ind <- ind1
  
  # Return a list containing the two updated data frames
  return(updated_track1)
}

deg2rad <- function(degrees) {
  return(degrees * pi / 180)
}

#Threshold for birds splitting in metres (if the leave this threshold and never come back to within this distance over the course of the tracks)
split_threshold <- 150

gps_list3 <- lapply(c(1:nrow(results_df)), function(i) 
{get_gps_distances_between_tracks(gps_list2, results_df$list_n[i], 
                                  results_df$partner_list_n[i],
                                  1, split_threshold)})


gps_list2 <- NULL

results_df$partner <- ifelse(results_df$release_cond=="solo", NA, results_df$partner)


##################
results <- lapply(c(1:length(gps_list3)), function(i){
  df0 <- gps_list3[[i]]
  
  left_RS_wholetrack <- which(df0$DistRS<2000)[length(which(df0$DistRS<2000))] #When did it leave vicinity of release site for last time
  gps_list3[[i]]$LeftRS <<- T
  gps_list3[[i]]$LeftRS[1:left_RS_wholetrack] <<- F
  
  #Did the birds in pairs split - analyse paired tracks before split
  if (!is.na(results_df$partner_list_n[i])){
    split_any <- any(df0$split, na.rm=T)
    df <- df0[which(df0$split==F),]
  }else{
    split_any <- NA
    df <- df0
  }
  if (!nrow(df)){
    return(c(df0$bird[1], df0$site[1], df0$release[1], results_df$list_n[i], results_df$partner_list_n[i], NA, split_any))
  }
  leftRS <- which(df$DistRS<2000)[length(which(df$DistRS<2000))]
  df_leaving <- df[1:leftRS,]
  
  if (leftRS<nrow(df)){
    df_left <- df[(leftRS+1):nrow(df),]
    home_approach <- distHaversine(c(df_left$Longitude[1], df_left$Latitude[1]), home_coords) - 
      distHaversine(c(df_left$Longitude[nrow(df_left)], df_left$Latitude[nrow(df_left)]), home_coords)
    hei <- home_approach/sum(df_left$Distance.m., na.rm=T) #Homing efficiency index (for after leaving vicinity of RS)
  } else{
    hei <- NA
  }
  
  return(c(df$bird[1], df$site[1], results_df$release_cond[i], 
           results_df$list_n[i], results_df$partner_list_n[i], hei, split_any
  ))
})

results0 <- do.call(rbind, lapply(results, function(x) as.data.frame(t(x))))
results0
colnames(results0) <- c("bird", "site", "release_cond", "list_n", "partner_list_n",
                        "HEI", 'split_any')

results_df1 <- merge(results_df, results0, by = c("bird", "site", "list_n", "release_cond", "partner_list_n"), all.x=T, all.y=F)
results_df1$HEI <- as.numeric(results_df1$HEI)

results_df1$site <- as.factor(results_df1$site)
results_df1$Pair_ID <- as.factor(results_df1$Pair_ID)
results_df1$treatment <- as.factor(results_df1$treatment)
results_df1$partner <- ifelse(results_df1$release_cond=="solo",
                              ",", results_df1$partner)
results_df1$split_any <- as.logical(unlist(results_df1$split_any))

results_df1$order1 <- as.factor(ifelse(results_df1$treatment!="CM", 0, results_df1$order))
results_df1$order2 <- as.factor(results_df1$order2)
results_df1$order3 <- as.factor(results_df1$order3)


find_nearest_distance <- function(point, track) {
  distances <- distHaversine(c(point[2], point[1]), 
                             cbind(track$Longitude, track$Latitude))
  return(min(distances))  # Return the minimum distance
}
#Find distances to across previous baseline tracks (for route fidelity measure)
track_distances_average_across_tracks <- function(ind_in, left_RS=T, remove_splitting=T){
  list_n_solo <- results_df1$list_n[ind_in]
  
  previous_track <- results_df1$list_n[which(results_df1$bird == results_df1$bird[ind_in] & 
                                               results_df1$site==results_df1$site[ind_in] &
                                               results_df1$release_cond=="pair" &
                                               results_df1$treatment=="TR" &
                                               results_df1$date < 
                                               results_df1$date[ind_in])]
  
  if (length(previous_track)==0){
    return(NA)
  }
  print(previous_track)
  df <- gps_list3[[list_n_solo]]
  
  if (left_RS){
    df <- df[df$LeftRS==T,]
  }else{
    df <- df[df$LeftRS==F,]
  }
  
  if (remove_splitting){
    if (df$release_cond[1]=="pair"){
      df <- df[df$split==F & !is.na(df$split),]
    }
  }
  df <- df[df$Speed.km.h.>30,]
  
  all_mean_dists <- c()
  for (tr in 1:length(previous_track)){
    df_prev <- gps_list3[[previous_track[tr]]]
    if (left_RS){
      df_prev <- df_prev[df_prev$LeftRS==T,]
    }else{
      df_prev <- df_prev[df_prev$LeftRS==F,]
    }
    df_prev <- df_prev[df_prev$Speed.km.h.>30,]
    df_prev <- df_prev[!is.na(df_prev$split) & df_prev$split==F,]
    
    if (nrow(df)==0 | nrow(df_prev)==0){
      all_mean_dists[tr] <- NA
    }else{
      # Calculate the nearest neighbor distance for each point in track1
      nearest_distances <- apply(data.frame(df$Latitude, df$Longitude), 1,
                                 find_nearest_distance, track = df_prev)
      all_mean_dists[tr] <- mean(nearest_distances)
    }
  }
  # Calculate mean distance
  overall_mean_distance <- mean(all_mean_dists, na.rm=T)
  return(overall_mean_distance)
}
#Find average distance to closest previous baseline tracks (for route fidelity measure)
track_distances_closest_track <- function(ind_in, left_RS=T, remove_splitting=T){
  list_n_solo <- results_df1$list_n[ind_in]
  
  previous_track <- results_df1$list_n[which(results_df1$bird == results_df1$bird[ind_in] & 
                                               results_df1$site==results_df1$site[ind_in] &
                                               results_df1$release_cond=="pair" &
                                               results_df1$treatment=="TR" &
                                               results_df1$date < 
                                               results_df1$date[ind_in])]
  
  if (length(previous_track)==0){
    return(NA)
  }
  print(previous_track)
  df <- gps_list3[[list_n_solo]]
  
  df_prev <- do.call(rbind, gps_list3[previous_track])
  
  if (left_RS){
    df <- df[df$LeftRS==T,]
    df_prev <- df_prev[df_prev$LeftRS==T,]
  }else{
    df <- df[df$LeftRS==F,]
    df_prev <- df_prev[df_prev$LeftRS==F,]
  }
  
  if (remove_splitting){
    if (df$release_cond[1]=="pair"){
      df <- df[df$split==F & !is.na(df$split),]
    }
  }
  
  df <- df[df$Speed.km.h.>30,]
  df_prev <- df_prev[df_prev$Speed.km.h.>30,]
  df_prev <- df_prev[!is.na(df_prev$split) & df_prev$split==F,]
  
  if (nrow(df)==0 | nrow(df_prev)==0){
    return(NA)
  }
  # Calculate the nearest neighbor distance for each point in track1
  nearest_distances <- apply(data.frame(df$Latitude, df$Longitude), 1,
                             find_nearest_distance, track = df_prev)
  # Calculate mean distance
  mean_distance <- mean(nearest_distances)
  return(mean_distance)
}


results_df1$distance_training_tracks_average <- NA
results_df1$distance_training_tracks_average <- unlist(
  lapply(1:nrow(results_df1), function(i) {
    print(i)
    track_distances_average_across_tracks(i)}))
results_df1$distance_training_tracks_average
length(which(!is.na(results_df1$distance_training_tracks_average)))

results_df1$distance_training_tracks_closest <- NA
results_df1$distance_training_tracks_closest <- unlist(
  lapply(1:nrow(results_df1), function(i) {
    print(i)
    track_distances_closest_track(i)}))
results_df1$distance_training_tracks_closest
length(which(!is.na(results_df1$distance_training_tracks_closest)))

#Subset so that there is only one row per pair (and solos retained)
pair_results_df1 <- results_df1[which(results_df1$partner_list_n>results_df1$list_n |
                                        results_df1$release_cond=="solo"),]
pair_results_df1$HEI_analysis <- NA
pair_results_df1$Distance_average_analysis <- NA
pair_results_df1$Distance_closest_analysis <- NA

for (i in 1:nrow(pair_results_df1)){
  if (pair_results_df1$release_cond[i]=="pair"){
    partner_ind <- pair_results_df1$partner_list_n[i]
    #Take average across pair
    pair_results_df1$HEI_analysis[i] <- mean(c(pair_results_df1$HEI[i], results_df1$HEI[partner_ind])) 
    pair_results_df1$Distance_average_analysis[i] <- mean(c(pair_results_df1$distance_training_tracks_average[i], 
                                                            results_df1$distance_training_tracks_average[partner_ind]))
    pair_results_df1$Distance_closest_analysis[i] <- mean(c(pair_results_df1$distance_training_tracks_closest[i], 
                                                            results_df1$distance_training_tracks_closest[partner_ind]))  
  }else{
    pair_results_df1$HEI_analysis[i] <- pair_results_df1$HEI[i]
    pair_results_df1$Distance_average_analysis[i] <- pair_results_df1$distance_training_tracks_average[i]
    pair_results_df1$Distance_closest_analysis[i] <- pair_results_df1$distance_training_tracks_closest[i]
  }
}

#In the analysis, only the first memory testing release is used (so order2 and order3 =1)
analysis_data <- pair_results_df1[pair_results_df1$order3!=2,]

#Variable distinguishes between memory testing (CM) with additional/extra training (CM_A), and without (CM_N)
analysis_data$treatment3 <- as.factor(ifelse(analysis_data$treatment2=="CM", 
                                             ifelse(analysis_data$add_train==T, "CM_A", "CM_N"),
                                             as.character(analysis_data$treatment2)))


all_processed_data <- do.call(rbind, gps_list3)
all_processed_data$experiment <- "Collective_memory"

all_processed_data$release_cond[all_processed_data$release_cond == "pair*"] <- "pair"

all_processed_data$order3 <- NA
for (i in 1:nrow(analysis_data)){
  list_i <- analysis_data$list_n[i]
  order3_i <- as.numeric(as.character(analysis_data$order3[i]))
  all_processed_data$order3[all_processed_data$track_ind == list_i] <- order3_i
}

#All data for analysis
all_processed_data <- all_processed_data[which(all_processed_data$order3 != 2),]
all_processed_data$order3 <- NULL

head(all_processed_data)
#write.csv(all_processed_data, paste(save_path, "CM_all_processed_data.csv", sep=""))

#Clean data - removes unnecessary additional variables
analysis_data$order <- NULL
analysis_data$order1 <- NULL
analysis_data$order2 <- NULL
analysis_data$order3 <- NULL
analysis_data$total_n <- NULL
analysis_data$treatment2 <- NULL
analysis_data$treatment_addtrain <- analysis_data$treatment3
analysis_data$treatment3 <- NULL
head(analysis_data)

#write.csv(analysis_data, paste(save_path, "CM_analysis_df.csv", sep=""))


