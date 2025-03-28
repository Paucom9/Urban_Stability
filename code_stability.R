library(dplyr)

df<- read.csv("C:/Users/pcolom/OneDrive/Escriptori/Urban stability/s_index_18_24_ubms_cbms.csv", sep = ";")

head(df)

df$SITE_ID <- as.factor(df$SITE_ID)

cbms_sites<- read.csv("C:/Users/pcolom/OneDrive/Escriptori/Urban stability/transect_sites_CBMS.csv", sep = ";")

head(cbms_sites)

cbms_sites$SITE_ID <- as.factor(cbms_sites$SITE_ID)

cbms_sites$SCHEME <- "CBMS"

df_merged <- left_join(df, cbms_sites, by = c("SITE_ID", "SCHEME"))

head(df_merged)

#Standarize SINDEX of CBMS sites
df_merged$SINDEX[df_merged$SCHEME == "CBMS"] <- 
  (df_merged$SINDEX[df_merged$SCHEME == "CBMS"] / df_merged$transect_length[df_merged$SCHEME == "CBMS"]) * 1000


# Calculate population stability

inverse_CV <- function(x) {
  mean(x, na.rm = TRUE) / sd(x, na.rm = TRUE)
}


df_pop_stb <- df_merged %>%
  group_by(SPECIES, SITE_ID) %>%
  summarise(pop_stability = mean(SINDEX, na.rm = TRUE) / sd(SINDEX, na.rm = TRUE),
            .groups = "drop")


hist(df_pop_stb$pop_stability)
summary(df_pop_stb$pop_stability)
