library(ggplot2)
library(rgdal)
library(dplyr)
library(terra)
library(sf)

source('scripts/config.R')

map_lat_long = readRDS('data/map_data/NA_map_lat_long.RDS')

##################################################################################################################################################
## raw pollen data
##################################################################################################################################################

pollen_raw = read.csv(paste0('data/pollen_north_america_LC6K_v', pollen_raw_version, '.csv'), header=TRUE)
length(unique(pollen_raw$dataset_id))

# pollen_north_america = pollen_north_america[which(pollen_north_america$keep == 1),]

##################################################################################################################################################
## prepared pollen data
##################################################################################################################################################

pollen_bin = readRDS(paste0('data/pollen_bin_v', pollen_prepared_version, '.RDS'))
pollen_bin = pollen_bin[which(pollen_bin$long < -50),]

##################################################################################################################################################
## record summary
##################################################################################################################################################

length(unique(pollen_bin$dataset_id))

coords = pollen_bin[,c('cell_id', 'long', 'lat')]
coords = coords[!duplicated(coords),]

coords_sp = SpatialPoints(coords[, c('long', 'lat')])
crs(coords_sp) = proj_WGS84

coords_equal = spTransform(coords_sp, proj_equal)
coords_equal_df = coordinates(coords_equal)
colnames(coords_equal_df) = c('x_equal', 'y_equal')

coords_equal_df = cbind(coords, coords_equal_df)

pollen_bin_equal_coords = coords_equal_df[match(pollen_bin$cell_id, coords_equal_df$cell_id), c('x_equal', 'y_equal')]
pollen_bin_equal = data.frame(pollen_bin_equal_coords, pollen_bin)

###############################################################################################################
## map showing record locations and case study regions
###############################################################################################################


# 
# coords_ENA = data.frame(x=c(-100, -52, -100, -52), y=c(48, 48, 67, 67))
# coords_ENA_sp = SpatialPoints(coords_ENA)
# crs(coords_ENA_sp) = proj_WGS84
# 
# foo = st_bbox(coords_ENA_sp)
# 
# coords_ENA_equal = spTransform(coords_ENA_sp, proj_equal)
# 
# 
# poly_ENA = coords_ENA_sp %>% 
#   st_as_sf(coords = c("x", "y"), 
#            crs = proj_WGS84) %>% 
#   st_bbox() %>% 
#   st_as_sfc()
# poly_ENA_equal = st_transform(poly_ENA, proj_equal)
# 
# ####
# 
# y_ENA = c(rep(48, 50), seq(48,67, length=50), rep(67, 50), seq(48,67, length=50))
# x_ENA = c(seq(-100,-52, length=50), rep(-100, 50), seq(-100,-52, length=50), rep(-52, 50))
# 
# coords_ENA = data.frame(x=x_ENA, y=y_ENA)
# coords_ENA_sp = SpatialPoints(coords_ENA)
# crs(coords_ENA_sp) = proj_WGS84
# 
# # SpatialPolygons(as.matrix(coords_ENA))
# 
# # foo = st_bbox(coords_ENA_sp)
# 
# coords_ENA_equal = spTransform(coords_ENA_sp, proj_equal)
# 
# 
# poly_ENA = coords_ENA_equal %>% 
#   st_as_sf(coords = c("x", "y"), 
#            crs = proj_equal) %>% 
#   summarise(geometry = st_combine(geometry)) %>%
#   st_cast("POLYGON")
#   # st_bbox() %>% 
#   # st_as_sfc()
# # poly_ENA_equal = st_transform(poly, proj_equal)
# 
# north <- 72.5; south <- 59; east <- -52; west <- -45
# 
# extent <- st_sf(a = 1:2, crs = proj_WGS84,
#                 geom = st_sfc(st_point(c(xhi, ylo)), 
#                               st_point(c(xlo, yhi)))) |>
#   st_transform(crs = proj_equal) |>
#   st_bbox()

## CASE: ENA
xlo = -100
xhi = -52
ylo = 48
yhi = 67
xlim = c(xlo, xhi)
ylim = c(ylo, yhi)

boxes_ENA = data.frame(yhi = 67, ylo = 48, xlo = -100, xhi = -52, id="ECAN")
# boxes_ENA = transform(boxes_ENA, laby=(yhi + ylo)/2, labx=(xhi + xlo)/2)
boxes_ENA = transform(boxes_ENA, laby=yhi-3, labx=xhi-6)

# CASE: HEMLOCK
xlo = -90
xhi = -60 
ylo = 37
yhi = 48
xlim = c(xlo, xhi)
ylim = c(ylo, yhi)

boxes_HEM = data.frame(yhi = 48, ylo = 37, xlo = -90, xhi = -60, id="NEUS/SEC")
boxes_HEM = transform(boxes_HEM, laby=ylo + 2, labx=xhi-7)

# CASE: WCAN
xlo = -165
xhi = -105
ylo = 45
yhi = 75
xlim = c(xlo, xhi)
ylim = c(ylo, yhi)

xmid = -141

boxes_WCAN = data.frame(yhi = 75, ylo = 45, xlo = -165, xhi = -105, id="WCAN/AK")
boxes_WCAN = transform(boxes_WCAN, laby= ylo + 3, labx= xlo + 10)

# CASE: WNA
xlo = -130
xhi = -115
ylo = 27
yhi = 60
xlim = c(xlo, xhi)
ylim = c(ylo, yhi)

ymid = 43

boxes_WNA = data.frame(yhi = yhi, ylo = ylo, xlo = xlo, xhi = xhi, id="PCCS")
boxes_WNA = transform(boxes_WNA, laby=ylo + 3, labx=xlo + 5)

###############################################################################################################
## show regions on a map
###############################################################################################################
ylim = c(10, 10, 75, 75) 
xlim = c(-166, -50, -166, -50) 


# -178.2,6.6,-49.0,83.3

coord_limits = data.frame(x=xlim, y=ylim)

coord_limits = SpatialPoints(coord_limits[,c('x', 'y')], 
                             proj4string=CRS(proj_WGS84))
coord_limits_equal =  spTransform(coord_limits, 
                                  proj_equal)
coord_limits_equal_df = data.frame(coordinates(coord_limits_equal))
colnames(coord_limits_equal_df) = c('x', 'y')

xlim_equal = c(min(coord_limits_equal_df$x), max(coord_limits_equal_df$x))
ylim_equal = c(min(coord_limits_equal_df$y), max(coord_limits_equal_df$y))

xlim_equal = c(-3500000, 3200000)
ylim_equal = c(-5500000, 800000)

map_equal = spTransform(map_lat_long, 
                        proj_equal)


ggplot()+
  geom_polygon(data=map_lat_long, aes(long,lat, group = group), color="darkgrey", fill="grey", alpha=0.4) +
  geom_point(data=coords, aes(x=long, y=lat), colour="navy", size=1.7, alpha=0.7) +
  geom_rect(data=boxes_ENA, aes(xmin=xhi, xmax=xlo, ymin=ylo, ymax=yhi), 
            color="dodgerblue", fill="transparent", lwd=1) + 
  geom_text(data=boxes_ENA, aes(x=labx, y=laby, label=id), color="dodgerblue", size=4) + 
  geom_rect(data=boxes_HEM, aes(xmin=xhi, xmax=xlo, ymin=ylo, ymax=yhi), 
            color="darkorange", fill="transparent", lwd=1) + 
  geom_text(data=boxes_HEM, aes(x=labx, y=laby, label=id), color="darkorange", size=4) + 
  geom_rect(data=boxes_WCAN, aes(xmin=xhi, xmax=xlo, ymin=ylo, ymax=yhi), 
            color="deeppink", fill="transparent", lwd=1) + 
  geom_text(data=boxes_WCAN, aes(x=labx, y=laby, label=id), color="deeppink", size=4) + 
  geom_segment(aes(x = xmid, y = ylo, xend = xmid, yend = yhi, colour = "segment"), 
               data = boxes_WCAN, lwd = 1, linetype = 4, color="deeppink") +
  geom_rect(data=boxes_WNA, aes(xmin=xhi, xmax=xlo, ymin=ylo, ymax=yhi), 
            color="purple", fill="transparent", lwd=1) + 
  geom_text(data=boxes_WNA, aes(x=labx, y=laby, label=id), color="purple", size=4) + 
  geom_segment(aes(x = xlo, y = ymid, xend = xhi, yend = ymid, colour = "segment"), 
               data = boxes_WNA, lwd = 1, linetype = 4, color="purple") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14)) +
  coord_fixed(xlim=xlim_NA, ylim = ylim_NA)
ggsave(paste0('figures/map_records_and_regions.png'))
ggsave(paste0('figures/map_records_and_regions.pdf'))

# ENA

y_ENA = c(rep(48, 50), seq(48,67, length=50), rep(67, 50), seq(67, 48, length=50))
x_ENA = c(seq(-100,-52, length=50),  rep(-52, 50), seq(-52,-100, length=50), rep(-100, 50))

coords_ENA = data.frame(x=x_ENA, y=y_ENA)
coords_ENA_sp = SpatialPoints(coords_ENA)
crs(coords_ENA_sp) = proj_WGS84

coords_ENA_equal = spTransform(coords_ENA_sp, proj_equal)

poly_ENA = coords_ENA_equal %>% 
  st_as_sf(coords = c("x", "y"), 
           crs = proj_equal) %>% 
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

# CASE: HEMLOCK
# boxes_HEM = data.frame(yhi = 48, ylo = 37, xlo = -90, xhi = -60, id="NEUS/SEC")
# # boxes_HEM = transform(boxes_HEM, laby=(yhi + ylo)/2, labx=(xhi + xlo)/2)
# boxes_HEM = transform(boxes_HEM, laby=ylo + 2, labx=xhi-7)

y_HEM = c(rep(37, 50), seq(37,48, length=50), rep(48, 50), seq(48, 37, length=50))
x_HEM = c(seq(-90,-60, length=50),  rep(-60, 50), seq(-60,-90, length=50), rep(-90, 50))

coords_HEM = data.frame(x=x_HEM, y=y_HEM)
coords_HEM_sp = SpatialPoints(coords_HEM)
crs(coords_HEM_sp) = proj_WGS84

coords_HEM_equal = spTransform(coords_HEM_sp, proj_equal)

poly_HEM = coords_HEM_equal %>% 
  st_as_sf(coords = c("x", "y"), 
           crs = proj_equal) %>% 
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")



# CASE: WCAN
xlo = -165
xhi = -105
ylo = 45
yhi = 75
xlim = c(xlo, xhi)
ylim = c(ylo, yhi)

xmid = -141
# 
# boxes_WCAN = data.frame(yhi = 75, ylo = 45, xlo = -165, xhi = -105, id="WCAN/AK")
# # boxes_WCAN = transform(boxes_WCAN, laby=(yhi + ylo)/2, labx=(xhi + xlo)/2)
# boxes_WCAN = transform(boxes_WCAN, laby= ylo + 3, labx= xlo + 10)


y_WCAN = c(rep(45, 50), seq(45,75, length=50), rep(75, 50), seq(75, 45, length=50))
x_WCAN = c(seq(-165,-105, length=50),  rep(-105, 50), seq(-105,-165, length=50), rep(-165, 50))

coords_WCAN = data.frame(x=x_WCAN, y=y_WCAN)
coords_WCAN_sp = SpatialPoints(coords_WCAN)
crs(coords_WCAN_sp) = proj_WGS84

coords_WCAN_equal = spTransform(coords_WCAN_sp, proj_equal)

poly_WCAN = coords_WCAN_equal %>% 
  st_as_sf(coords = c("x", "y"), 
           crs = proj_equal) %>% 
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")



# CASE: WNA
xlo = -130
xhi = -115
ylo = 27
yhi = 60
xlim = c(xlo, xhi)
ylim = c(ylo, yhi)

ymid = 43

boxes_WNA = data.frame(yhi = yhi, ylo = ylo, xlo = xlo, xhi = xhi, id="PCCS")
# boxes_WNA = transform(boxes_WNA, laby=(yhi + ylo)/2, labx=(xhi + xlo)/2)

boxes_WNA = transform(boxes_WNA, laby=ylo + 3, labx=xlo + 5)


y_WNA = c(rep(27, 50), seq(27,60, length=50), rep(60, 50), seq(60, 27, length=50))
x_WNA = c(seq(-130,-115, length=50),  rep(-115, 50), seq(-115,-130, length=50), rep(-130, 50))

coords_WNA = data.frame(x=x_WNA, y=y_WNA)
coords_WNA_sp = SpatialPoints(coords_WNA)
crs(coords_WNA_sp) = proj_WGS84

coords_WNA_equal = spTransform(coords_WNA_sp, proj_equal)

poly_WNA = coords_WNA_equal %>% 
  st_as_sf(coords = c("x", "y"), 
           crs = proj_equal) %>% 
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")


lab_cases = data.frame(x=c(-58, -53, -155, -125), 
                       y=c(64, 39, 48, 30), 
                       case = c("ECAN", "NEUS/SEC", "WCAN/AK", "PCCS"),
                       colour=c("dodgerblue", "darkorange", "deeppink", "purple"))

lab_cases = data.frame(x=c(-58, -68, -150, -125), 
                       y=c(62, 34, 50, 30), 
                       case = c("ECAN", "NEUS/SEC", "WCAN/AK", "PCCS"),
                       colour=c("dodgerblue", "darkorange", "deeppink", "purple"))
lab_cases_sp = SpatialPoints(lab_cases[,c('x', 'y')])
crs(lab_cases_sp) = proj_WGS84
lab_cases_equal = spTransform(lab_cases_sp, proj_equal)
lab_cases_equal = data.frame(coordinates(lab_cases_equal), lab_cases[,c('case', 'colour')])

ggplot()+
  geom_polygon(data=map_equal, aes(long,lat, group = group), linewidth=0.3, color="grey74", fill="grey", alpha=0.4) +
  geom_point(data=coords_equal_df, aes(x=x_equal, y=y_equal), colour="navy", size=1.5, alpha=0.4) +
  geom_sf(data=poly_ENA, color=alpha("dodgerblue", 0.7), fill=NA, lwd=1) +
  # geom_text(data=lab_ENA_equal, aes(x=x, y=y), label='ECAN', color="dodgerblue", size=4) + 
  geom_sf(data=poly_HEM, color=alpha("darkorange", 0.7), fill=NA, lwd=1) +
  geom_sf(data=poly_WCAN, color=alpha("deeppink", 0.7), fill=NA, lwd=1) +
  geom_sf(data=poly_WNA, color=alpha("purple", 0.7), fill=NA, lwd=1) +
  geom_text(data=lab_cases_equal, aes(x=x, y=y, label=case, colour=colour), size=4) + 
  scale_colour_manual(values=c("dodgerblue", "darkorange", "deeppink", "purple"), 
                      labels=c("ECAN", "NEUS/SEC", "WCAN/AK", "PCCS"),
                      breaks=c("dodgerblue", "darkorange", "deeppink", "purple"),
                      name = "Case study") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        legend.position = "none") +
  coord_sf(crs = proj_equal, xlim=xlim_equal, ylim=ylim_equal)
ggsave(paste0('figures/map_records_and_regions_equal.png'))
ggsave(paste0('figures/map_records_and_regions_equal.pdf'))

##################################################################################################################################################
## figure 1: record locations for each time slices
##################################################################################################################################################

pollen_bin_equal_sub = pollen_bin_equal[which(bin_label_H[pollen_bin_equal$slice_bin] %in% ages_sub_H),]
pollen_bin_equal_sub$age = bin_label_H[pollen_bin_equal_sub$slice_bin]

pollen_bin_equal_sub$age = factor(pollen_bin_equal_sub$age, 
                                  levels=ages_sub_H, 
                                  labels=ages_sub_label_H)

p = ggplot()+
  geom_polygon(data=map_equal, aes(long,lat, group = group), linewidth=0.3, color="grey74", fill="grey", alpha=0.4) +
  geom_point(data=pollen_bin_equal_sub, aes(x=x_equal, y=y_equal), colour="navy", size=1.5, alpha=0.4) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        legend.position = "none") +
  coord_sf(crs = proj_equal, xlim=xlim_equal, ylim=ylim_equal) +
  facet_wrap(~age)
print(p)
ggsave(paste0('figures/map_records_ages_subset.png'), width=8, height=10)
ggsave(paste0('figures/map_records_ages_subset.pdf'), width=8, height=10)

# ##################################################################################################################################################
# ## figure 1: record locations
# ##################################################################################################################################################
# 
# p <- ggplot(data=coords) +
#   geom_polygon(data=na_fort, aes(x=long, y=lat, group=group), color="grey54", fill="grey") +
#   geom_point(aes(x=long, y=lat), colour="navy", size=1.7, alpha=0.7) +
#   theme_bw() + 
#   theme(axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank()) + 
#   # coord_cartesian(xlim = longlimits, ylim = latlimits) +
#   coord_fixed(xlim = longlimits, ylim = latlimits)
# print(p)
# ggsave('figures/figure1_record_coords_map.pdf')
# ggsave('figures/figure1_record_coords_map.png')
# 
# 
# coords = pollen_bin[,c('lat', 'long')]
# coords = coords[!duplicated(coords),]
# 
# p <- ggplot(data=coords) +
#   geom_polygon(data=na_fort, aes(x=long, y=lat, group=group), color="grey54", fill="grey") +
#   geom_point(aes(x=long, y=lat), colour="navy", size=1.7, alpha=0.7) +
#   theme_bw() + 
#   theme(axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank()) + 
#   # coord_cartesian(xlim = longlimits, ylim = latlimits) +
#   coord_fixed(xlim = longlimits, ylim = latlimits) 
# print(p)
# ggsave('figures/figure1_record_coords_map.pdf')
# ggsave('figures/figure1_record_coords_map.png')

##################################################################################################################################################
## figure 2: record temporal extent
##################################################################################################################################################

pollen_record_extent = pollen_bin %>%
  group_by(dataset_id) %>%
  summarize(age_old = max(age_calbp)/1000, age_young = min(age_calbp)/1000) %>% 
  arrange(desc(age_old))
pollen_record_extent$stat_id = seq(1, nrow(pollen_record_extent))

ggplot(data = pollen_record_extent) +
  geom_linerange(aes(y=stat_id, xmin=age_old, xmax=age_young), colour="grey50", lwd=0.2, alpha=0.75) +
  theme_bw(14) +
  scale_x_reverse(breaks=c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(breaks=seq(0, 1500, by=250)) +
  xlab('Time (ka)') +
  ylab('Record number') 
ggsave('figures/pollen_record_temporal_extent.pdf')
ggsave('figures/pollen_record_temporal_extent.png')

##################################################################################################################################################
## Supplementary table 1: record metadata
##################################################################################################################################################
library(neotoma2)

pollen_sites = pollen_bin[!duplicated(pollen_bin$dataset_id),]
dsids = unique(pollen_sites$dataset_id)

meta_df = pollen_sites[c('sitename', 'dataset_id', 'lat', 'long', 'lake_size')]
meta_df = data.frame(meta_df,
                     elevation = rep(NA),
                     citation = rep(NA),
                     DOI = rep(NA))

N_datasets = nrow(meta_df)

for (i in 1:N_datasets){
  this_dataset = get_datasets(x = meta_df$dataset_id[i])   
  this_DOI = doi(this_dataset)$doi
  this_citation = cite_data(this_dataset)
  this_author = strsplit(this_citation$citation, '[.]')[[1]][1]
  
  this_author = strsplit(this_citation$citation, '; pollen dataset')[[1]][1]
  print(this_citation$citation)
  # bar = get_publications(datasetid = meta_df$dataset_id[i])
  
  meta_df[i, 'elevation'] = this_dataset$altitude
  meta_df[i, 'citation'] = this_author
  meta_df[i, 'DOI'] = this_DOI
  
}

colnames(meta_df) = c('Site name', 'Dataset ID', 'Latitude', 
                      'Longitude', 'Lake size', 'Elevation (m)', 'Citation', 'DOI')

meta_df[, c('Latitude', 'Longitude', 'Lake size')] = round(meta_df[, c('Latitude', 'Longitude')], 2)

write.csv(meta_df, 'data/pollen_record_metadata_table.csv', row.names = FALSE)
