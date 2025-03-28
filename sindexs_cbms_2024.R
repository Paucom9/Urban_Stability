library(dplyr)

df<- read.csv("C:/Users/pcolom/OneDrive/Escriptori/Urban stability/sites_24_pau.csv", sep = ";")

head(df)

cbms_sites<- read.csv("C:/Users/pcolom/OneDrive/Escriptori/Urban stability/transect_sites_CBMS.csv", sep = ";")

head(cbms_sites)

sp_names<- read.csv("C:/Users/pcolom/OneDrive/Escriptori/Urban stability/names_cbms.csv", sep = ";")

head(sp_names)
sp_names<- sp_names[, c("IDesp", "sp_latin")]

colnames(sp_names) <- c("SPECIES", "SPECIES_name")

#####

cbms_sites$scheme <- "CBMS"

df_merged <- left_join(df, cbms_sites, by = c("SITE_ID", "scheme"))

head(df_merged)

# Delete cases of sp * site * year with SINDEX = 0
head(df_merged)

df_merged <- df_merged[df_merged$SINDEX != 0, ]


# Include species name


df_merged2 <- left_join(df_merged, sp_names, by = "SPECIES")

unique(df_merged2$SPECIES_name)

write.csv(df_merged2, "C:/Users/pcolom/OneDrive/Escriptori/Urban stability/CBMS_sindex_2024.csv")

# Standarize SINDEX to 1000m of transect length

df_merged$SINDEX_standarized <- (df_merged$SINDEX / df_merged$transect_length) * 1000

hist(df_merged$SINDEX_standarized)



