library(data.table)
library(tidyverse)
library(zoo)

setwd("/rawdata")

# Read in pop_table files
pop_files <- list.files(pattern = "pop_table1")

# rbind pop files to a single data frame
pop_table <- rbindlist(lapply(pop_files, fread), fill = TRUE)

# Standardize County column, bring in county names from cor_county
pop_table <- pop_table %>% 
  mutate(County = ifelse(is.na(cor_county), County, cor_county))
# Check for NA's
sum(is.na(pop_table$County))#Nil

# Check unique county names
unique(pop_table$County)#few inconsistencies need to be reconciled

# Using case_when to standardize county names
pop_table <- pop_table %>% 
  mutate(County = case_when(County == "revard" ~ "Brevard", 
                            County == "Orange/Reedy" ~ "Orange_Reedy", 
                            County == "Palm Beach" ~ "Palm_Beach", 
                            County == "St. Johns" ~ "St_Johns", 
                            County == "St. Lucie" ~ "St_Lucie", 
                            County == "Indian River" ~ "Indian_River", 
                            County == "Santa Rosa" ~ "Santa_Rosa", 
                            County %in% c("Walton, North", "Walton, North ") ~ 
                              "Walton_North", 
                            County %in% c("South_Walton", "Walton, South", 
                                          "Walton, South ") ~ "Walton_South",
                            County %in% c("au", "sau") ~ "Nassau", 
                            County %in% c("ion", "rion") ~ "Marion", 
                            County %in% c("artin", "irtin", "vartin", "Mortin") ~ 
                              "Martin", 
                            County == "ay" ~ "Bay", 
                            TRUE ~ pop_table$County))

# Select useful columns
pop_table1 <- pop_table %>% 
  select(County, rep_date, n_birds_Flavi_sus, birds_flavi_susc, n_birds_Alpha_sus, 
         birds_alpha_susc, birds_susc)

# Creating a column that merges data from n_birds_Flavi_sus and birds_flavi_susc
pop_table1 <- pop_table1 %>% 
  mutate(birds_susc_flavi = ifelse(is.na(n_birds_Flavi_sus) & is.na(birds_flavi_susc), NA, 
                                   ifelse(n_birds_Flavi_sus > 0 & is.na(birds_flavi_susc), 
                                          n_birds_Flavi_sus, birds_flavi_susc)), 
         birds_susc_alpha = ifelse(is.na(n_birds_Alpha_sus) & is.na(birds_alpha_susc), NA, 
                                   ifelse(n_birds_Alpha_sus > 0 & is.na(birds_alpha_susc), 
                                          n_birds_Alpha_sus, birds_alpha_susc)))

# birds_susc col is not filled for earlier years, imputing that data from the newly
# created columns as the max value of birds_susc_flavi and birds_susc_alpha. This is 
# how this column has been filled for later years.
# Convert to numeric
pop_table1$birds_susc <- as.numeric(pop_table1$birds_susc)

# Rewrite column values
pop_table1 <- pop_table1 %>% 
  mutate(birds_susc = ifelse(is.na(birds_susc) & birds_susc_flavi >= 0 & birds_susc_alpha >= 0,
                             pmax(birds_susc_flavi, birds_susc_alpha, na.rm = TRUE), 
                             birds_susc))

# Selecting useful columns
pop_table1 <- pop_table1 %>% 
  select(County, rep_date, birds_susc, birds_susc_flavi, birds_susc_alpha)

# Check NA's in the date column before moving forward,
which(is.na(pop_table1$rep_date))
# Seems like for one full week between 6/11/2004 - 6/25/2004 the dates have not 
# been entered, same for between 11/6/2009 - 11/20/2009 these dates can be imputed

# Filling dates for the 6/11/2004 - 6/25/2004 gap with 6/18/2004
pop_table1[7333:7375, 2] <- "6/18/2004"

# Filling dates for the 11/6/2009 - 11/20/2009 gap with 11/13/2009
pop_table1[19538:19584, 2] <- "11/13/2009"

# Check NA's
sum(is.na(pop_table1$rep_date))#Nil

# Chanage NA's to 0 in the rest of the data frame
pop_table1[is.na(pop_table1)] <- 0

# Format date
pop_table1$rep_date <- parse_date_time(pop_table1$rep_date, orders = c('mdy', 'ymd'))

# Check NA's
which(is.na(pop_table1$rep_date))#Nil

#Set date
pop_table1$rep_date <- as.POSIXct(pop_table1$rep_date)

# Creating day, week, month, year column
pop_table1 <- pop_table1 %>% 
  mutate(day = day(pop_table1$rep_date), 
         week = isoweek(pop_table1$rep_date), 
         month = month(pop_table1$rep_date), 
         year = year(pop_table1$rep_date))

# Read in positive chicken data for wnv and eeev
wnv <- read.csv("WNV2001-2020_FDOH (3).csv")
eeev <- read.csv("EEEV_2001_2020_FDOH.csv")
str(wnv)
str(eeev)

# Create list of sites and counties
site_list_wnv <- unique(wnv$Site)#525 sites
site_list_eeev <- unique(eeev$Site)#277 sites

county_list_wnv <- unique(wnv$County)#43 counties
county_list_eeev <- unique(eeev$County)#38 counties

# Format column names
colnames(wnv) <- c("County", "CollDate_WNV", "Site", "LabNo_WNV", "Status_WNV", "Lon", "Lat")
colnames(eeev) <- c("County", "CollDate_EEEV", "Site", "LabNo_EEEV", "Status_EEEV", "Lon", "Lat", "year")

# Format date
wnv$CollDate_WNV <- parse_date_time(wnv$CollDate_WNV, orders = c('mdy', 'ymd'))
eeev$CollDate_EEEV <- parse_date_time(eeev$CollDate_EEEV, orders = c('mdy', 'ymd'))

# Check NA's
which(is.na(wnv$CollDate_WNV))#nil
which(is.na(eeev$CollDate_EEEV))#nil

# Set date
wnv$CollDate_WNV <- as.POSIXct(wnv$CollDate_WNV)
eeev$CollDate_EEEV <- as.POSIXct(eeev$CollDate_EEEV)

# Creating week and year column
wnv <- wnv %>% 
  mutate(week = isoweek(wnv$CollDate_WNV), 
         year = year(wnv$CollDate_WNV))

eeev <- eeev %>% 
  mutate(week = isoweek(eeev$CollDate_EEEV), 
         year = year(eeev$CollDate_EEEV))

# Check site names between the two df's in case there are different spellings
setdiff(eeev$Site, wnv$Site)#14 sites in eeev are not present in wnv

# Manually search for these in wnv to see if there are spelling errors

# Change site names in eeev N. Berry -> Newberry, MIKE -> Mike, 3-4 bear creek -> 3-4 Bear Creek
eeev <- eeev %>% 
  mutate(Site = case_when(Site == "N. Berry " ~ "Newberry", 
                          Site == "MIKE" ~ "Mike", 
                          Site %in% c("3-4 bear creek", "3-4 Bear creek") ~ 
                            "3-4 Bear Creek", 
                          TRUE ~ eeev$Site))

# Change site name in wnv 3-4 Bear creek -> 3-4 Bear Creek
wnv <- wnv %>% 
  mutate(Site = case_when(Site == "3-4 Bear creek" ~ "3-4 Bear Creek", 
                          TRUE ~ wnv$Site))

# Since there are still major inconsistencies between site names within each of
# these sheets writing these files to do clustering and cleaning of site names on 
# OpenRefine.
#write wnv
#write.csv(wnv, "C:/Users/Akshay/UFL Dropbox/Akshay Vinod Anand/chickens/data/chickens/raw_data/wnv_filtered.csv", row.names = FALSE)
#write eeev
#write.csv(eeev, "C:/Users/Akshay/UFL Dropbox/Akshay Vinod Anand/chickens/data/chickens/raw_data/eeev_filtered.csv", row.names = FALSE)

# Read in filtered wnv and eeev files
wnv_clean <- read.csv("wnv-filtered-OR.csv")
eeev_clean <- read.csv("eeev-filtered-OR.csv")
str(wnv_clean)
str(eeev_clean)

# Checking for inconsistencies
setdiff(eeev_clean$Site, wnv_clean$Site)#much better only unique sites are showing

# Calculate summarized number of positive chickens per site week combination
# for each virus to join with individual df's

# WNV
positive_wnv <- wnv_clean %>% 
  group_by(County, Site, week, year, Lon, Lat) %>% 
  summarise(n_positive_wnv = n())

# EEEV
positive_eeev <- eeev_clean %>% 
  group_by(County, Site, week, year, Lon, Lat) %>% 
  summarise(n_positive_eeev = n())

# Join to form positive chickens data
positive_chickens <- merge(positive_wnv, positive_eeev, all = TRUE)

# Check county names between pop_table1 and positive chickens
unique(pop_table1$County)
unique(positive_chickens$County)

# Some county names need to be changed in positive_chickens to merge data
positive_chickens <- positive_chickens %>% 
  mutate(County = case_when(County == "Indian River" ~ "Indian_River", 
                            County == "Orange/Reedy" ~ "Orange_Reedy", 
                            County == "Palm Beach" ~ "Palm_Beach", 
                            County == "St. Johns" ~ "St_Johns", 
                            County == "St. Lucie" ~ "St_Lucie", 
                            County %in% c("Walton (North)", "Walton (north)")
                            ~ "Walton_North", 
                            County %in% c("Walton (South)", "Walton (south)")
                            ~ "Walton_South", 
                            TRUE ~ positive_chickens$County))

# Check for differences
setdiff(positive_chickens$County, pop_table1$County)#Nil

# Creating a county + site list for all available sites
county_site <- positive_chickens %>% 
  group_by(County, Site, Lon, Lat) %>% 
  reframe()

# There are still some duplicated sites names due to different coordinates 
# for the same site name. Checked some of these coordinates on GIS and they are
# quite far (1 - 5 km) from each other. Therefore treating these as separate sites. 
# Might have to change names later but for now coordinates are identifiers

# Joining county_site to pop_table1 to have long format data
pop_table_sites <- left_join(county_site, pop_table1, by = "County", 
                             relationship = "many-to-many")

# Joining this with eeev data
dat_join <- left_join(pop_table_sites, positive_eeev, 
                      by = c("County", "Site", "Lon", "Lat", "week", "year"), 
                      relationship = "many-to-many")

# Create column to record testing
dat_join <- dat_join %>%
  group_by(Site, year, week, Lon, Lat) %>%
  mutate(testing = as.integer(
    any(birds_susc > 0, na.rm = TRUE) | any(n_positive_eeev > 0, na.rm = TRUE)
  )) %>%
  ungroup()

# Create column for number of positive with NA = no testing, 0 = testing but no positive, number of positives
dat_join <- dat_join %>% 
  mutate(n_positive_EEEV = case_when(
    testing == 0 ~ NA_real_,
    is.na(n_positive_eeev) ~ 0,
    TRUE ~ n_positive_eeev))

# Create a column of wnv_status where NA = no testing, 0 = no positives, 1 = at least 1 positive
dat_join <- dat_join %>% 
  mutate(EEEV_status = ifelse(n_positive_EEEV > 0, 1, n_positive_EEEV))

# Select useful columns
dat_join <- dat_join %>%
  select(
    County, Site, Lon, Lat,           
    year, month, week, day,           
    rep_date,                         
    birds_susc,                       
    testing,                          
    n_positive_EEEV,                  
    EEEV_status)

# Cleaning col names
new_col_names <- colnames(dat_join) %>% 
  str_to_lower()
setnames(dat_join, new_col_names)

# Arrange data based on week and year
dat <- arrange(dat_join, week, year)

# Sanity check
dat %>% filter(testing == 0 & !is.na(n_positive_eeev))
dat %>% filter(testing == 1 & is.na(n_positive_eeev))#all good

# Clear unnecessary objects
rm(list = setdiff(ls(), "dat"))

# Examine county level sampling effort/week
sam_eff_county <- dat %>% 
  group_by(county) %>% 
  summarise(sam_eff = sum(testing)) %>% 
  arrange(sam_eff)#all counties are sampled

# Convert lon & lat to numeric
dat$lon <- as.numeric(dat$lon)
dat$lat <- as.numeric(dat$lat)


# Checking which sites have the same name but different coordinates
dup_sites <- dat %>%
  group_by(county, site) %>%
  filter(n_distinct(lon, lat) > 1) %>%
  distinct(site, lon, lat) %>% 
  ungroup() %>%
  arrange(site)

## Reconciling names
# Create an object of unique site locations per county
county_site <- dat %>% 
  group_by(county, site, lon, lat) %>% 
  summarise() %>% 
  mutate(site = ifelse(is.na(site), county, site)) %>% 
  ungroup()

# Creating a unique ID column
county_site <- county_site %>% 
  arrange(county, lon, lat) %>% 
  group_by(county) %>% 
  mutate(site_id = paste0(county, "_", row_number())) %>%
  ungroup()

# Adding unique ID's to filtered data
dat_clean <- left_join(dat, county_site, by = c("county", "site", "lon", "lat"), 
                      relationship = "many-to-many")

# Reorder columns
dat_clean <- dat_clean[, c(1, 14, 3:13)]

# Removing sites without coordinates
dat_clean <- dat_clean %>% 
  filter(!is.na(lat)) %>% 
  filter(!is.na(lon))

# Save data
fwrite(dat_clean, file = "data/compiled_eeev_data.csv", row.names = FALSE)