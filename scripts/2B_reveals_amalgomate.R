version = "8.0"
LC6k = TRUE

# lc6k
files <- list.files("data/cache", pattern="veg_pred_LC6k_[0-9]+.RDS", full.names=TRUE)
dfs <- lapply(files, readRDS)
df <- do.call(rbind, dfs)
saveRDS(df, file="data/cache/veg_pred_LC6k.RDS")
saveRDS(df, file="data/veg_pred_LC6k.RDS")

# # LGM
# files <- list.files("data/cache", pattern="veg_pred_LGM_[0-9]+.RDS", full.names=TRUE)
# dfs <- lapply(files, readRDS)
# df <- do.call(rbind, dfs)
# saveRDS(df, file=paste0("data/cache/veg_pred_LGM_", version, ".RDS"))

