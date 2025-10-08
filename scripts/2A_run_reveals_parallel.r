library(raster)
library(sp)
# #library(DISQOVER)
# library(reshape2)
# library(ggplot2)
library(plyr)
library(dplyr)
library(tidyverse)
# #library(neotoma)
# library(maptools)
# library(fields)
# library(maps)
# library(rgdal)
# library(truncnorm)

#devtools::install_github("PalEON-Project/stepps-cal")
# library(stepps)

# # load package
# load('DISQOVER/R/sysdata.rda')

source('scripts/utils/REVEALSinR_main.R')

source('scripts/config.R')

args = commandArgs(trailingOnly=TRUE)
id_to_run = as.integer(args[1])

##################################################################################################################################################
## read in pollen data for NA
##################################################################################################################################################

if (LC6K){
  pollen_trans = readRDS(paste0('data/pollen_bin_v', pollen_prepared_version, '.RDS'))
} else {
  pollen_trans = readRDS(paste0('data/pollen_bin_LGM_v', version, '.RDS'))
}

pollen_trans = pollen_trans[, -ncol(pollen_trans)]
taxa = colnames(pollen_trans)[8:ncol(pollen_trans)]
colnames(pollen_trans)[which(colnames(pollen_trans) == 'age_calbp')] = 'age'


pollen_trans$slice_bin = cut(pollen_trans$age, breaks_H*1000, labels=FALSE)
colnames(pollen_trans) = tolower(colnames(pollen_trans))
pollen_trans = pollen_trans[which(!is.na(pollen_trans$slice_bin)), ]

# sum samples within a time bin for each site 
pollen_bin = ddply(pollen_trans, .(dataset_id, lat, long, sitename, lake_size, slice_bin),  
                   function(x) colSums(x[tolower(taxa)]))

# some summary stats
length(unique(pollen_bin$dataset_id))
table(pollen_bin$slice_bin) / length(unique(pollen_bin$dataset_id))

##################################################################################################################################################
## fill in missing lake sizes
##################################################################################################################################################

# use 50 hectares are medium lake size
# assigned to all lakes with missing lake size
medium_radius = sqrt(50*0.01 / pi)*1000

idx_na = which(is.na(pollen_bin$lake_size))

pollen_bin$lake_size[idx_na] = medium_radius

##################################################################################################################################################
## sort some taxa stuff
##################################################################################################################################################
taxa_trans = read.csv('data/taxa_latin2common_filled.csv', stringsAsFactors = FALSE)

colnames(pollen_bin)[8:ncol(pollen_bin)] = taxa_trans$common[match(colnames(pollen_bin[,8:ncol(pollen_bin)]),taxa_trans$latin)]
colnames(pollen_bin)[is.na(colnames(pollen_bin))] <- 'unknown'

pollen_bin = pollen_bin %>%
  pivot_longer(cols=colnames(pollen_bin)[8:ncol(pollen_bin)]) %>%
  group_by(cell_id, dataset_id, lat, long, sitename, lake_size, slice_bin, name) %>% 
  dplyr::summarize(value=sum(value), .groups='keep') %>%
  pivot_wider(names_from = name, values_from = value)

matched = taxa_trans$common[match(colnames(pollen_bin[,8:ncol(pollen_bin)]),taxa_trans$common)]
taxa_common = sort(matched[which(!is.na(matched))])

pollen_bin = data.frame(pollen_bin[,1:7], pollen_bin[,taxa_common])

saveRDS(pollen_bin, paste0('data/pollen_bin_taxa_v', pollen_prepared_version, '.RDS'))

##################################################################################################################################################
## read in and prep ppes and svs
##################################################################################################################################################

N_taxa = length(taxa_common)
pars = data.frame(species = taxa_common, 
                  PPEs = numeric(N_taxa), 
                  PPE.errors = numeric(N_taxa), 
                  fallspeed = numeric(N_taxa)) 
pars[,2:4] = NA

ppes_NH = read.csv('data/RPP_FS_wieczorek.csv', stringsAsFactors=FALSE)

for (k in 1:N_taxa){
  
  taxon = taxa_common[k]
  idx_row = which(pars$species == taxon)
  
  if (is.na(pars$PPEs[idx_row])){
    
    idx_match = match(taxon, tolower(ppes_NH$taxon))
    
    if (taxon == 'pea') {
      
      idx_match = match('legume', tolower(ppes_NH$taxon))
      
    }
    
    if ((!is.na(ppes_NH[idx_match, "RPP.v2..America."])) & (!is.na(ppes_NH[idx_match, "SE..America."]))){
      
      pars[idx_row, 'PPEs'] = ppes_NH[idx_match, "RPP.v2..America."]
      pars[idx_row, 'PPE.errors'] = ppes_NH[idx_match,  "SE..America." ]
      
    } else {
      
      print(taxon)
      pars[idx_row, 'PPEs'] = ppes_NH[idx_match, "ppe"]
      pars[idx_row, 'PPE.errors'] = ppes_NH[idx_match,  "error" ]
      
    }
    
  }
  
  if (is.na(pars$fallspeed[idx_row])){
    
    idx_match = match(taxon, tolower(ppes_NH$taxon))
    
    if (taxon == 'pea') {
      
      idx_match = match('legume', tolower(ppes_NH$taxon))
      
    }
    
    if (taxon == 'larch'){
      
      next
      
    }
    
    if (!is.na(ppes_NH[idx_match,  "FS..m.s...America."])){
      
      pars[idx_row, 'fallspeed'] = ppes_NH[idx_match,  "FS..m.s...America."]
      
    } else {
      
      pars[idx_row, 'fallspeed'] = ppes_NH[idx_match, "fallspeed"]
      
    }
    
  }
  
}


ppes_NA = readRDS('data/PPEs_grass.RDS')

idx_na = which(is.na(pars$PPEs))

idx_match = match(pars[idx_na, 'species'], tolower(ppes_NA$taxon))

pars$PPEs[idx_na] = ppes_NA[idx_match, 'ppe.grass']
pars$PPE.errors[idx_na] = ppes_NA[idx_match, 'error.grass']


svs_NA = read.csv('data/svs_LC6K.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)

svs_NA$taxon = tolower(svs_NA$taxon)
svs_NA = svs_NA[which(svs_NA$use == 'y'),]

svs_NA %>% dplyr::group_by(taxon) %>% dplyr::summarize(sv = median(sv))
svs_NA = svs_NA %>% dplyr::group_by(taxon) %>% dplyr::summarize(sv = mean(sv))

idx_na = which(is.na(pars$fallspeed))
idx_match = match(pars[idx_na, 'species'], tolower(svs_NA$taxon))

pars$fallspeed[idx_na] = svs_NA[idx_match, 'sv']$sv

pars$species = tolower(pars$species)
all(match(taxa_common, pars$species))

reveals_inputs = pars
rownames(reveals_inputs) = NULL

write.csv(reveals_inputs, "data/REVEALS_input/params.csv", row.names=FALSE)

lct_trans = read.csv('data/taxon2LCT_translation_v2.csv', header=TRUE)
lct_trans$LCT = factor(lct_trans$LCT, 
                       levels = c('ET', 'ST', 'OL'), 
                       labels = c('ETS', 'STS', 'OVL'))

REVEALS_input_table = data.frame(species_latin = str_to_title(taxa_trans$latin[match(reveals_inputs$species, taxa_trans$common)]),
                                 reveals_inputs,
                                 LCT = lct_trans$LCT[match(reveals_inputs$species, lct_trans$taxon)])
REVEALS_input_table$PPEs = round(REVEALS_input_table$PPEs, digits = 2)
REVEALS_input_table$PPE.errors = round(REVEALS_input_table$PPE.errors, digits = 2)
REVEALS_input_table$fallspeed = round(REVEALS_input_table$fallspeed, digits = 3)

colnames(REVEALS_input_table) = c('Taxon Latin', 'Taxon Common', 'PPE (unitless)', 'PPE SE', 'FS (cm/s)', 'LCT')

write.csv(REVEALS_input_table, 'data/REVEALS_input_REVEALS_input_table.csv', row.names = FALSE)

##################################################################################################################################################
## run reveals
##################################################################################################################################################

ids = unique(pollen_bin$dataset)
pol_dat = pollen_bin

print(length(ids))


veg_pred = data.frame(id=numeric(0), 
                      long=numeric(0), 
                      lat=numeric(0), 
                      cell_id = numeric(0),
                      basin_radius=numeric(0),
                      ages=numeric(0), 
                      taxon=character(0), 
                      meansim=numeric(0), 
                      mediansim=numeric(0), 
                      q10sim=numeric(0), 
                      q90sim=numeric(0), 
                      sdsim=numeric(0))




id = ids[id_to_run]

counts_site = data.frame(ages=bin_label_H[pol_dat[which(pol_dat$dataset_id == id),'slice_bin']], 
                         pol_dat[which(pol_dat$dataset_id == id),which(colnames(pol_dat) %in% tolower(taxa_common))])
colnames(counts_site) = c('ages', taxa_common)

basin_radius = pol_dat[which(pol_dat$dataset_id == id), c('lake_size')][1]

coords_site = data.frame(pol_dat[which(pol_dat$dataset_id == id), c('dataset_id', 'long', 'lat', 'cell_id')][1,], basin_radius)
rownames(coords_site) = NULL

# prob don't need this step - fix later
write.csv(counts_site, 'data/REVEALS_input/reveals_input.csv', row.names=FALSE)

# cycle through and estimate background veg
# with english csv files
dbasin = 2*round(basin_radius)

site_veg_out = REVEALSinR(pollenFile = "data/REVEALS_input/reveals_input.csv",
                          pf         = "data/REVEALS_input/params.csv",
                          filetype   = "csv",
                          dwm        = "gpm neutral",
                          tBasin     = "lake",
                          dBasin     = 2*round(basin_radius), # diameter!
                          regionCutoff = 100000,
                          repeats      = 1000)

veg = melt(site_veg_out, 
           id.vars=c('Pollen.file', 'Parameter.file', 'Distance.weighting', 'Basin.type', 'ages'))

veg$taxon = unlist(lapply(veg$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[1]))
veg$type  = unlist(lapply(veg$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[2]))

veg_cast = dcast(veg, ages + taxon ~ type)
veg_pred = rbind(veg_pred, data.frame(coords_site, veg_cast))

if (LC6k){
  fname = sprintf('data/cache/veg_pred_LC6K_%06d.RDS', id_to_run)
} else {
  fname = sprintf('data/cache/veg_pred_LGM_%06d.RDS', id_to_run)
}
saveRDS(veg_pred, fname)
