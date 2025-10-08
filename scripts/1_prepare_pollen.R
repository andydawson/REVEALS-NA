library(plyr)
library(ggplot2)
library(raster)

source('scripts/config.R')

##################################################################################################################################################
## load pollen data for NA
##################################################################################################################################################

pollen_north_america = read.csv(paste0('data/pollen_north_america_LC6K_v', pollen_raw_version, '.csv'), 
                                header=TRUE)
sites_north_america = read.csv(paste0('data/sites_north_america_LC6K_v', pollen_raw_version, '.csv'), 
                               header=TRUE)

pollen_north_america = pollen_north_america[which(pollen_north_america$keep == 1),]

pollen_trans = pollen_north_america
pollen_trans = data.frame(lat = sites_north_america$latitude[match(pollen_trans$dataset_id, sites_north_america$datasetid)], 
                          long = sites_north_america$longitude[match(pollen_trans$dataset_id, sites_north_america$datasetid)], pollen_trans)
colnames(pollen_trans)[which(colnames(pollen_trans) == 'age')] = 'age_calBP'
pollen_trans$age_calBP = as.numeric(pollen_trans$age_calBP)

taxa = colnames(pollen_trans)[15:ncol(pollen_trans)]

##################################################################################################################################################
## bin pollen in time intervals
##################################################################################################################################################

pollen_trans$slice_bin = cut(pollen_trans$age_calBP, breaks_H*1000, labels=FALSE)
colnames(pollen_trans) = tolower(colnames(pollen_trans))
pollen_trans = pollen_trans[which(!is.na(pollen_trans$slice_bin)), ]

# sum samples within a time bin for each site 
pollen_bin = ddply(pollen_trans, c('dataset_id', 'lat', 'long', 'sitename', 'age_calbp', 'slice_bin'),  
                   function(x) colSums(x[tolower(taxa)]))

##################################################################################################################################################
## lookup lakesizes
##################################################################################################################################################

site_meta = data.frame(stid = pollen_trans$siteid,
                       dsid = pollen_trans$dataset_id,
                       sitename = pollen_trans$sitename,
                       lat = pollen_trans$lat,
                       long = pollen_trans$long)

site_meta = site_meta[!duplicated(site_meta),]
site_meta$dep_env = sites_north_america[match(site_meta$dsid, sites_north_america$datasetid), 'depositionalenvironment']


# four data files with lake sizes
lake_sizes = read.csv('data/area_lakes_1.4.csv')
usa_areas = read.csv('data/lake_areas_usa_v1.0.csv')
can_areas = read.csv('data/ca_lakes.csv')
wisest_areas = read.csv('data/lakes-NA-sites.csv')
wisest_areas$lake_size..ha. = as.numeric(wisest_areas$lake_size..ha.)

site_meta$edited = lake_sizes$edited[match(site_meta$stid, lake_sizes$stid)]
site_meta$notes = lake_sizes$notes[match(site_meta$stid, lake_sizes$stid)]

site_meta$lake_size = lake_sizes$area[match(site_meta$stid, lake_sizes$stid)]
sum(!is.na((site_meta$lake_size)))/length(site_meta$lake_size)

for (i in 1:nrow(site_meta)){
  #print(i)
  if (is.na(site_meta$lake_size[i])){
    idx_match = match(site_meta$stid[i], lake_sizes$stid)
    if (length(idx_match)>1){
      print(i)
    }
    if (!is.na(match(site_meta$stid[i], usa_areas$siteid))){
      site_meta$lake_size[i] = usa_areas$areaha[match(site_meta$stid[i], usa_areas$siteid)]
    # } else if (!is.na(match(site_meta$stid[i], can_areas$siteid))){
    #   site_meta$lake_size[i] = can_areas$lake_area_ha[match(site_meta$stid[i], can_areas$siteid)]
    #   
    } else if (!is.na(match(site_meta$stid[i], can_areas$siteid))){
      site_meta$lake_size[i] = can_areas$lake_area_ha[match(site_meta$stid[i], can_areas$siteid)]
      
    } else if (!is.na(match(site_meta$stid[i], wisest_areas$stid))){
      site_meta$lake_size[i] = wisest_areas$lake_size..ha.[match(site_meta$stid[i], wisest_areas$stid)]
      
    }
  }
}

sum(!is.na((site_meta$lake_size)))/length(site_meta$lake_size)
sum(is.na((site_meta$lake_size)))
length(site_meta$lake_size)

# write.csv(site_meta, 'data/lakes-NA-sites-compiled.csv', row.names=FALSE)


lake_sizes$dsid = site_meta[match(lake_sizes$stid, site_meta$stid), 'dsid']

pollen_bin <- data.frame(lake_size=rep(NA), pollen_bin)
pollen_bin$lake_size = lake_sizes$area[match(pollen_bin$dataset, lake_sizes$dsid)]
pollen_bin$lake_size_bool = !is.na(pollen_bin$lake_size)

# # missing area
# length(unique(pollen_bin$dataset[which(is.na(pollen_bin$lake_size))]))
# saveRDS(unique(pollen_bin$dataset[which(is.na(pollen_bin$lake_size))]), 'data/missing_lake_sizes.RDS')
# # with area
# length(unique(pollen_bin$dataset[which(!is.na(pollen_bin$lake_size))]))
# # zero area
# length(unique(pollen_bin$dataset[which(pollen_bin$lake_size==0)]))

# in HA; convert to radius in m
# pi * r * r
pollen_bin$lake_size = sqrt(pollen_bin$lake_size*0.01 / pi)*1000

##################################################################################################################################################
## make the 1 degree by 1 degree grid for NA
##################################################################################################################################################

source('scripts/make_grid.R')
grid_NA = make_grid(pollen_bin, 
                    coord_fun = ~ long + lat, 
                    projection = '+init=epsg:4326', 
                    resolution = 1)
saveRDS(grid_NA, 'data/grid_NA.RDS')

grid_NA_area = area(grid_NA)
cell_area = raster::extract(grid_NA_area, pollen_bin[,c('long', 'lat')])
cell_id = raster::extract(grid_NA, pollen_bin[,c('long', 'lat')])
pollen_bin = data.frame(cell_id, pollen_bin)

# saveRDS(pollen_bin, paste0('data/pollen_bin_LGM_v', version, '.RDS'))
saveRDS(pollen_bin, paste0('data/pollen_bin_v', pollen_prepared_version, '.RDS'))

