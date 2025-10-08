library(raster)
library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(elevatr)
library(sp)

source('scripts/config.R')

veg_pred = readRDS('data/veg_pred_LC6k.RDS')
grid_NA = readRDS('data/grid_NA.RDS')

coords   = raster::xyFromCell(grid_NA, veg_pred$cell_id)
veg_pred = cbind(coords, veg_pred)

veg_grid = aggregate(mediansim ~ taxon + ages + cell_id + x + y, veg_pred, sum)
veg_cast = dcast(veg_grid, cell_id + x + y + ages ~ taxon, value.var='mediansim')
veg_cast[,5:ncol(veg_cast)] = t(apply(veg_cast[,5:ncol(veg_cast)], 1, function(x) x/sum(x)))

veg_grid = melt(veg_cast, id.vars=c('cell_id', 'x', 'y', 'ages'))

lat_long_coords = data.frame(veg_grid[,c('x', 'y')])
ele_get = get_elev_point(lat_long_coords , proj_WGS84, src = 'aws')
veg_grid$elev = ele_get$elevation

nslices_plot = length(bin_label_H)
veg_grid = veg_grid[which(veg_grid$ages %in% bin_label_H),]

breaks = c(0, 1e-5, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1)
n_breaks = length(breaks) - 1
veg_grid$value_bin = cut(veg_grid$value, breaks, labels=FALSE)

ice = readRDS('data/map_data/ice/glacier_shapefiles_21-1k.RDS')
ice_years = seq(1, 21)*1000

veg_grid_ice = data.frame(ice = rep(NA), veg_grid)

ages = unique(veg_grid$ages)

ice_fort = data.frame(matrix(NA, nrow=0, ncol=10))
for (i in 1:length(ages)){
  print(i)
  idx_age = which(veg_grid_ice$ages == ages[i])
  coords  = SpatialPoints(veg_grid_ice[idx_age,c('x', 'y')], proj4string=CRS(proj_WGS84))
  
  if (i %in% c(1, 2, 3)){
    veg_grid_ice[idx_age, 'ice'] = NA
  } else {
    idx_ice_match = which.min(abs(ages[i] - ice_years)) 
    proj4string(ice[[idx_ice_match]]) = proj_WGS84
    
    ice_status = over(coords, ice[[idx_ice_match]])
    ice_fort_age = fortify(ice[[idx_ice_match]])
    ice_fort = rbind(ice_fort, 
                     data.frame(ice_fort_age, 
                                ice_year = rep(ice_years[idx_ice_match], nrow(ice_fort_age)), 
                                ages = rep(ages[i], nrow(ice_fort_age))))
    
    veg_grid_ice[idx_age, 'ice'] = ice_status[,which(substr(colnames(ice_status), 1,4)=='SYMB')]
  }
}

veg_grid = veg_grid_ice
write.csv(veg_grid, paste0('data/veg_pred_grid_v', land_cover_version, '.csv'), row.names=FALSE)


taxon2pft = read.csv('data/taxon2LCT_translation_v2.csv')
veg_grid$LCT = taxon2pft$LCT[match(veg_grid$variable, tolower(taxon2pft$taxon))]
veg_lct = aggregate(value ~ cell_id + x + y + elev + ages + LCT, veg_grid, sum, na.rm=TRUE)
write.csv(veg_lct, paste0('data/veg_pred_grid_lct_v', land_cover_version, '.csv'), row.names=FALSE)

veg_lct_wide = veg_lct %>% 
  group_by(cell_id, x, y, ages, elev) %>%
  pivot_wider(names_from = LCT, values_from = value) %>%
  arrange(cell_id, x, y, ages)
write.csv(veg_lct_wide, paste0('data/veg_pred_grid_lct_wide_v', land_cover_version, '.csv'), row.names=FALSE)

veg_lct$value_bin = cut(veg_lct$value, breaks, labels=FALSE)
veg_lct[which(veg_lct$value==0),'value_bin'] = 1
veg_lct[which(veg_lct$value>0.99999999999999),'value_bin'] = 11
write.csv(veg_lct, paste0('data/veg_pred_grid_lct_binned_v', land_cover_version, '.csv'), row.names=FALSE)
