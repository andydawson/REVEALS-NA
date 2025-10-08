# USA Contiguous albers equal area
# proj_out <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
# WGS84
proj_WGS84 <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
# updated projection
proj_equal <- '+proj=aea +lat_1=33.33 +lat_2=66.66 +lat_0=70 +lon_0=-100 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

LC6K = TRUE

ylim_NA = c(12, 82) 
xlim_NA = c(-166, -50) 

latlimits = c(10, 80)
longlimits = c(-165, -50)
bounding_box = c(-165, 10, -50, 80)

pollen_raw_version = "8.0"
pollen_prepared_version = "9.0"
land_cover_version = "9.0"

breaks_LGM = c(-74, 0.1, 0.35, 0.7, seq(1.5, 21.5, by=1))
bin_number_LGM = seq(1, 24)
bin_label_LGM = c(50, 200, 500, seq(1000, 21000, by=1000))

breaks_H = c(-74, 0.1, 0.35, 0.7, seq(1.2, 11.7, by=0.5))
bin_number_H = seq(1, 25)
bin_label_H = c(50, 200, seq(500, 11500, by=500))

ages_sub_H = c(50, 500, 2000, 4000, 6000, 8000, 10000)
ages_sub_label_H = c('0.05 ka', '0.5 ka', '2 ka', '4 ka', '6 ka', '8 ka', '10 ka')



