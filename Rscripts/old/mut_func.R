
#########################
#
# libraries
#
########################

library(tidyverse)
library(data.table)
library(cowplot)
library(stringr)

library(fst) 
library(arrow)
library(scales)
library(lavaan)
library(duckdb) 
library(DBI)
library(readr)

setwd("/Users/gvbarroso/Data")

num_chrom <- 22

########################
#
# analyzing mutation rates within / outside phastcons
#
########################

# get 1kb predictions of B
maps <- c("roulette", "carlson", "gnomad")
bs_1kb <- vector("list", length(maps))

for(m in 1:length(maps)) {
  
  map <- maps[m]
  print(map)
  
  dir <- paste("phastcons_1kb_B_mpp", map, sep="/")
  bs_1kb_map <- vector("list", 22)
  
  pb <- txtProgressBar(min=0, max=22, style=3)
  for(c in 1:22) {
    
    setTxtProgressBar(pb, c)
    
    bs <- fread(paste0(dir, "/B_raw_YRI_chr", c, ".csv.gz")) %>% 
      dplyr::select(., c(chrom, chromStart, chromEnd, B, num_sites))
    
    bs$elem <- "phastcons"
    bs$map <- map
    
    bs_1kb_map[[c]] <- bs
  }
  close(pb)
  
  bs_1kb[[m]] <- data.table::rbindlist(bs_1kb_map)
}

bs_1kb <- data.table::rbindlist(bs_1kb)

# getting mutation rates in phastcons
df_maps <- vector("list", length=length(maps))
pos <- 1
for(m in maps) {
  
  print(m)
  bs <- filter(bs_1kb, map==m) %>% dplyr::select(., -num_sites)
  
  chr_maps <- vector("list", 22)
  mut_map <- vector("list", 22)
  pb <- txtProgressBar(min=0, max=22, style=3)
  for(c in 1:22) {
    
    setTxtProgressBar(pb, c)
    
    mut_map[[c]] <- fread(paste("bgs_lmr/human_data/data/mutation_tables/", m,
                                "/", m, "_tbl_chr", c, ".csv.gz", sep="")) 
    
    phast <- vector("list", 12) # 12 phastcons bins
    for(l in seq(0, 55, 5)) {
      
      u <- l + 5
      
      tmp <- fread(paste("bgs_lmr/human_data/data/annotation_tables/", m, "/phastcons/top", l, "-", u,
                         "/phastcons_top", l, "-", u, "_tbl_chr", c, ".csv.gz", sep="")) 
      
      tmp$phast_bin <- u / 5
      tmp$chrom <- c
      tmp$chrom <- as.integer(tmp$chrom)
      
      tmp[, del_sites := NULL]
      tmp[, scale := NULL]
      
      phast[[l %/% 5 + 1]] <- tmp
    }
    
    phast <- data.table::rbindlist(phast) %>% setDT()
    chr_maps[[c]] <- phast
  }
  close(pb)
  
  mut_map <- data.table::rbindlist(mut_map)
  chr_maps <- data.table::rbindlist(chr_maps)
  
  mut_phast_w <- pivot_wider(chr_maps, values_from=c(scale_masked, del_sites_masked), names_from=phast_bin, names_prefix="phast_") %>% setDT()
  mut_phast_w[, del_sites_sum := rowSums(.SD, na.rm=TRUE), .SDcols=paste0("del_sites_masked_phast_", 1:12)]
  mut_phast_w <- mut_phast_w[del_sites_sum > 0,]
  
  dt <- merge(mut_phast_w, mut_map, by=c("chrom", "chromStart", "chromEnd"))
  dt <- dt[!is.na(avg_mut_masked) & avg_mut_masked > 0 & num_sites_masked > 0]
  dt <- merge(dt, bs, by=c("chrom", "chromStart", "chromEnd"))  
  
  sum_mut_bin <- dt[["avg_mut_masked"]] * dt[["num_sites_masked"]]
  sum_mut_del <- vector("list", 12)
  for(i in 1:12) {
    sum_mut_del[[i]] <- dt[["avg_mut_masked"]] * dt[[paste0("scale_masked_phast_", i)]] * dt[[paste0("del_sites_masked_phast_", i)]] %>% as.data.frame() %>% setDT()
  }
  
  sum_mut_del <- data.table::cbindlist(sum_mut_del)
  names(sum_mut_del) <- paste0("phast_U_bin_", 1:12)
  sum_mut_del$phast_U <- rowSums(sum_mut_del, na.rm=T)
  
  avg_mut_out <- (sum_mut_bin - sum_mut_del$phast_U) / (dt$num_sites_masked - dt$del_sites_sum)
  avg_mut_out[!is.finite(avg_mut_out)] <- NA
  avg_mut_in <- vector("list", 12)
  for(i in 1:12) {
    avg_mut_in[[i]] <- dt[["avg_mut_masked"]] * dt[[paste0("scale_masked_phast_", i)]] %>% as.data.frame() %>% setDT()
  }
  
  avg_mut_in <- data.table::cbindlist(avg_mut_in)
  names(avg_mut_in) <- paste0("phast_mut_bin_", 1:12)
  
  dt[, Outside := avg_mut_out]
  dt <- cbind.data.frame(dt, avg_mut_in) %>% dplyr::select(., -c(num_sites, elem, starts_with("scale")))
  
  df_maps[[pos]] <- dt
  pos <- pos + 1

  tmp <- cbind.data.frame(avg_mut_in, avg_mut_out)
  tmp <- pivot_longer(tmp, cols=c(avg_mut_out, starts_with("phast_mut"))) %>% filter(., !is.na(value))
  
  summary(filter(tmp, name=="avg_mut_out")$value)
          
  p <- ggplot(data=filter(tmp, value < 2e-8), aes(x=value, fill=name)) +
    geom_density(alpha=0.1) +
    scale_fill_viridis_d(labels=c(phast_mut_bin_1="1",
                                  phast_mut_bin_2="2",
                                  phast_mut_bin_3="3",
                                  phast_mut_bin_4="4",
                                  phast_mut_bin_5="5",
                                  phast_mut_bin_6="6",
                                  phast_mut_bin_7="7",
                                  phast_mut_bin_8="8",
                                  phast_mut_bin_9="9",
                                  phast_mut_bin_10="10",
                                  phast_mut_bin_11="11",
                                  phast_mut_bin_12="12", 
                                  avg_mut_out="Outside")) +
    labs(x=expression(mu), y="Density", fill=NULL) +
    theme_bw() +
    guides(fill=guide_legend(nrow=1)) +
    theme(panel.grid.minor=element_blank(),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16),
          axis.text.y=element_text(hjust=1),
          legend.position="bottom")
  
  save_plot(paste0("~/Desktop/phast_mu_dens_", m, ".pdf"), p, base_height=6, base_width=8)
}

df_maps <- data.table::rbindlist(df_maps)
fwrite(df_maps, "phast_muts.csv.gz")

# reload
df_maps <- fread("phast_muts.csv.gz")
df_maps <- dplyr::select(df_maps, -c(chrom, chromStart, chromEnd, avg_mut))

summarize_phast_bins <- function(dt, phast_prefix="phast_mut_bin_", bins=1:12) {
  
  res <- lapply(bins, function(i) {
    col <- paste0(phast_prefix, i)
    tmp <- dt[!is.na(get(col))]
    
    list(phast_bin=i,
         median_muts_in=median(tmp[[col]], na.rm=TRUE),
         median_muts_out=median(tmp[["Outside"]], na.rm=TRUE),
         se_muts_in=sd(tmp[[col]], na.rm=TRUE) / sqrt(nrow(tmp)),
         se_muts_out=sd(tmp[["Outside"]], na.rm=TRUE) / sqrt(nrow(tmp)),
         median_Bs=median(tmp[["B"]], na.rm=TRUE))
    })
  
  data.table::rbindlist(res)
}

tbl_gnomad <- summarize_phast_bins(filter(df_maps, map=="gnomad"))
tbl_gnomad$map <- "gnomad"

tbl_carlson <- summarize_phast_bins(filter(df_maps, map=="carlson"))
tbl_carlson$map <- "carlson"

tbl_roulette <- summarize_phast_bins(filter(df_maps, map=="roulette"))
tbl_roulette$map <- "roulette"

df1 <- rbind.data.frame(tbl_gnomad, tbl_carlson, tbl_roulette) 
df1$control_gnomad <- FALSE

df2 <- df1 %>% group_by(map) %>%
  mutate(median_ratio=median_muts_in / median_muts_out,
         se_ratio=median_ratio * sqrt((se_muts_in / median_muts_in)^2 + (se_muts_out / median_muts_out)^2))

df3_all <- df2 %>% group_by(map) %>%
  pivot_longer(cols=c(median_ratio, se_ratio),
               names_to=c(".value", "metric"),
               names_pattern="(median|se)_(.*)") %>% setDT()
df3_all$elem <- "phast_all"

# selecting 1 kb windows with only one phastcons class (5-percentile bin) represented
x_list <- lapply(1:12, function(i) {
  df <- filter(df_maps, del_sites_sum == .data[[paste0("del_sites_masked_phast_", i)]]) %>%
    dplyr::select(!!paste0("del_sites_masked_phast_", i),
                  !!paste0("phast_mut_bin_", i),
                  Outside, B, map) %>%
    mutate(phast_bin=i) %>% setDT() 
  
  setnames(df, old=c(paste0("phast_mut_bin_", i), paste0("del_sites_masked_phast_", i)),
           new=c("phast_mut", "del_sites_masked_phast"))
  
  df
})

dt_isolated <- data.table::rbindlist(x_list)

summary_mut_df <- vector("list", 3)
for(m in 1:length(maps)) {
  
  map_dt <- dplyr::filter(dt_isolated, map==maps[m], !is.na(Outside))

  df1 <- map_dt %>%
    group_by(phast_bin) %>%
    summarise(median_muts_in=median(phast_mut, na.rm=TRUE),
              median_muts_out=median(Outside, na.rm=TRUE),
              se_muts_in=sd(phast_mut, na.rm=TRUE) / sqrt(nrow(map_dt)),
              se_muts_out=sd(Outside, na.rm=TRUE) / sqrt(nrow(map_dt)),
              median_B=median(B, na.rm=TRUE),
              map=maps[m],
              .groups="drop") %>% collect()
  
  df2 <- df1 %>% 
    mutate(median_ratio=median_muts_in / median_muts_out,
           se_ratio=median_ratio * sqrt((se_muts_in / median_muts_in)^2 + (se_muts_out / median_muts_out)^2))
  
  summary_mut_df[[m]] <- df2
}

df2 <- data.table::rbindlist(summary_mut_df)

df3 <- df2 %>%
  pivot_longer(cols=c(median_ratio, se_ratio),
               names_to=c(".value", "metric"),
               names_pattern="(median|se)_(.*)")

d1 <- ggplot(df3, aes(x = phast_bin, y = median, color = median_B)) +
  facet_wrap(~map) + theme_classic() + 
  geom_point(size=3) + geom_line() +
  geom_errorbar(aes(ymin = median - se, ymax = median + se),
                width = 0.2, linewidth = 0.7) +
  labs(x="phastCons bin", y=expression(paste(mu, " ratio"))) +
  scale_x_continuous(breaks=1:12) +
  scale_y_continuous(breaks=pretty_breaks(), limits=c(min(df3$median) - df3$se[which.min(df3$median)], 1)) +
  scale_color_viridis_c(option="C", direction=1, name="B",
                        breaks=c(round(min(df3$median_B)+0.01, 2), 
                                 round(max(df3$median_B)-0.01, 2))) +
  theme(strip.text=element_text(size = 16),
        axis.title=element_text(size = 20),
        axis.text=element_text(size = 16),
        legend.text=element_text(size = 16), 
        legend.title=element_text(size = 16, margin = margin(b = 20)), 
        legend.position="bottom")

# looking at 1-bp elements
df_1bp <- filter(dt_isolated, del_sites_masked_phast==1)

summary_mut_df <- vector("list", 3)
for(m in 1:length(maps)) {
  
  map_dt <- dplyr::filter(df_1bp, map==maps[m], !is.na(Outside))
  
  df1 <- map_dt %>%
    group_by(phast_bin) %>%
    summarise(median_muts_in=median(phast_mut, na.rm=TRUE),
              median_muts_out=median(Outside, na.rm=TRUE),
              se_muts_in=sd(phast_mut, na.rm=TRUE) / sqrt(nrow(map_dt)),
              se_muts_out=sd(Outside, na.rm=TRUE) / sqrt(nrow(map_dt)),
              median_Bs=median(B, na.rm=TRUE),
              map=maps[m],
              .groups="drop") %>% collect()
  
  df2 <- df1 %>% 
    mutate(median_ratio=median_muts_in / median_muts_out,
           se_ratio=median_ratio * sqrt((se_muts_in / median_muts_in)^2 + (se_muts_out / median_muts_out)^2))
  
  summary_mut_df[[m]] <- df2
}

df2 <- data.table::rbindlist(summary_mut_df)

df3_1bp <- df2 %>%
  pivot_longer(cols=c(median_ratio, se_ratio),
               names_to=c(".value", "metric"),
               names_pattern="(median|se)_(.*)") %>% setDT()
df3_1bp$control_gnomad <- FALSE
df3_1bp$elem <- "phast_1bp"

df4_A <- rbind.data.frame(df3_all, df3_1bp)

########################
#
# controlling for sequence composition by dividing roulette and carlson by gnomad
#
########################

chr_maps <- vector("list", 22)
pb <- txtProgressBar(min=0, max=22, style=3)
for(c in 1:22) {
  
  setTxtProgressBar(pb, c)
  
  carlson <- fread(paste("bgs_lmr/human_data/data/mutation_tables/carlson/carlson_tbl_chr", c, ".csv.gz", sep="")) %>% setDT()
  roulette <- fread(paste("bgs_lmr/human_data/data/mutation_tables/roulette/roulette_tbl_chr", c, ".csv.gz", sep="")) %>% setDT()
  gnomad <- fread(paste("bgs_lmr/human_data/data/mutation_tables/gnomad/gnomad_tbl_chr", c, ".csv.gz", sep="")) %>% setDT()
  
  roulette[, rel_mut := avg_mut / gnomad$avg_mut]
  carlson[, rel_mut := avg_mut / gnomad$avg_mut]
  
  roulette[, map := "roulette"]
  carlson[, map := "carlson"]
  
  carlson_roulette <- list(carlson, roulette)
  
  tbl <- vector("list", length(carlson_roulette) * 12) # 12 phastcons bin
  pos <- 1 
  for(i in 1:length(carlson_roulette)) {
    
    mut_map <- carlson_roulette[[i]]
    m <- unique(mut_map$map)
    
    phast <- vector("list", 12) # 12 phastcons bins
    for(l in seq(0, 55, 5)) {
      
      u <- l + 5
      
      tmp <- fread(paste("bgs_lmr/human_data/data/annotation_tables/", m, "/phastcons/top", l, "-", u,
                         "/phastcons_top", l, "-", u, "_tbl_chr", c, ".csv.gz", sep="")) 
      
      tmp$phast_bin <- u / 5
      tmp$chrom <- c
      tmp$chrom <- as.integer(tmp$chrom)
      
      tmp[, del_sites := NULL]
      tmp[, scale := NULL]
      
      phast[[l %/% 5 + 1]] <- tmp
    }
    
    phast <- data.table::rbindlist(phast) %>% setDT()

    mut_phast_w <- pivot_wider(phast, values_from=c(scale_masked, del_sites_masked), names_from=phast_bin, names_prefix="phast_") %>% setDT()
    mut_phast_w[, del_sites_sum := rowSums(.SD, na.rm=TRUE), .SDcols=paste0("del_sites_masked_phast_", 1:12)]
    mut_phast_w <- mut_phast_w[del_sites_sum > 0,]
    
    dt <- merge(mut_phast_w, mut_map, by=c("chrom", "chromStart", "chromEnd"))
    dt <- dt[!is.na(avg_mut_masked) & avg_mut_masked > 0 & num_sites_masked > 0]

    sum_mut_bin <- dt[["rel_mut"]] * dt[["num_sites_masked"]]
    sum_mut_del <- vector("list", 12)
    for(i in 1:12) {
      sum_mut_del[[i]] <- dt[["rel_mut"]] * dt[[paste0("scale_masked_phast_", i)]] * dt[[paste0("del_sites_masked_phast_", i)]] %>% as.data.frame() %>% setDT()
    }
    
    sum_mut_del <- data.table::cbindlist(sum_mut_del)
    names(sum_mut_del) <- paste0("phast_U_bin_", 1:12)
    sum_mut_del$phast_U <- rowSums(sum_mut_del, na.rm=T)
    
    avg_mut_out <- (sum_mut_bin - sum_mut_del$phast_U) / (dt$num_sites_masked - dt$del_sites_sum)
    avg_mut_out[!is.finite(avg_mut_out)] <- NA
    avg_mut_in <- vector("list", 12)
    for(i in 1:12) {
      avg_mut_in[[i]] <- dt[["rel_mut"]] * dt[[paste0("scale_masked_phast_", i)]] %>% as.data.frame() %>% setDT()
    }
    
    avg_mut_in <- data.table::cbindlist(avg_mut_in)
    names(avg_mut_in) <- paste0("phast_mut_bin_", 1:12)
    
    dt[, Outside := avg_mut_out]
    dt <- cbind.data.frame(dt, avg_mut_in) %>% dplyr::select(., -c(num_sites, avg_mut, starts_with("scale")))
    
    tbl[[pos]] <- dt
    pos <- pos + 1
  }
  
  df_maps <- data.table::rbindlist(tbl)
  chr_maps[[c]] <- df_maps
}
close(pb)

chr_maps <- data.table::rbindlist(chr_maps)
chr_maps <- dplyr::select(chr_maps, -c(chrom, chromStart, chromEnd))

summarize_phast_bins <- function(dt, phast_prefix="phast_mut_bin_", bins=1:12) {
  
  res <- lapply(bins, function(i) {
    
    col <- paste0(phast_prefix, i)
    tmp <- dt[!is.na(get(col))]
    
    list(phast_bin=i,
         median_muts_in=median(tmp[[col]], na.rm=TRUE),
         median_muts_out=median(tmp[["Outside"]], na.rm=TRUE),
         se_muts_in=sd(tmp[[col]], na.rm=TRUE) / sqrt(nrow(tmp)),
         se_muts_out=sd(tmp[["Outside"]], na.rm=TRUE) / sqrt(nrow(tmp)))})
  
  data.table::rbindlist(res)
}

tbl_carlson <- summarize_phast_bins(filter(chr_maps, map=="carlson"))
tbl_carlson$map <- "carlson"

tbl_roulette <- summarize_phast_bins(filter(chr_maps, map=="roulette"))
tbl_roulette$map <- "roulette"

df1 <- rbind.data.frame(tbl_carlson, tbl_roulette) 
df1$control_gnomad <- TRUE

df2 <- df1 %>% group_by(map) %>%
  mutate(median_ratio=median_muts_in / median_muts_out,
         se_ratio=median_ratio * sqrt((se_muts_in / median_muts_in)^2 + (se_muts_out / median_muts_out)^2))

df3_all <- df2 %>% group_by(map) %>%
  pivot_longer(cols=c(median_ratio, se_ratio),
               names_to=c(".value", "metric"),
               names_pattern="(median|se)_(.*)") %>% setDT()
df3_all$elem <- "phast_all"

# selecting 1 kb windows with only one phastcons class (5-percentile bin) represented
x_list <- lapply(1:12, function(i) {
  df <- filter(chr_maps, del_sites_sum == .data[[paste0("del_sites_masked_phast_", i)]]) %>%
    dplyr::select(!!paste0("del_sites_masked_phast_", i),
                  !!paste0("phast_mut_bin_", i),
                  Outside, map) %>%
    mutate(phast_bin=i) %>% setDT() 
  
  setnames(df, old=c(paste0("phast_mut_bin_", i), paste0("del_sites_masked_phast_", i)),
           new=c("phast_mut", "del_sites_masked_phast"))
  
  df
})

dt_isolated <- data.table::rbindlist(x_list)
df_1bp <- filter(dt_isolated, del_sites_masked_phast==1)

summary_mut_df <- vector("list", 3)
for(m in 1:length(maps)) {
  
  map_dt <- dplyr::filter(df_1bp, map==maps[m], !is.na(Outside))
  
  df1 <- map_dt %>%
    group_by(phast_bin) %>%
    summarise(median_muts_in=median(phast_mut, na.rm=TRUE),
              median_muts_out=median(Outside, na.rm=TRUE),
              se_muts_in=sd(phast_mut, na.rm=TRUE) / sqrt(nrow(map_dt)),
              se_muts_out=sd(Outside, na.rm=TRUE) / sqrt(nrow(map_dt)),
              map=maps[m],
              .groups="drop") %>% collect()
  
  df2 <- df1 %>% 
    mutate(median_ratio=median_muts_in / median_muts_out,
           se_ratio=median_ratio * sqrt((se_muts_in / median_muts_in)^2 + (se_muts_out / median_muts_out)^2))
  
  summary_mut_df[[m]] <- df2
}

df2 <- data.table::rbindlist(summary_mut_df)

df3_1bp <- df2 %>%
  pivot_longer(cols=c(median_ratio, se_ratio),
               names_to=c(".value", "metric"),
               names_pattern="(median|se)_(.*)")
df3_1bp$elem <- "phast_1bp"
df3_1bp$control_gnomad <- TRUE

df4_B <- rbind.data.frame(df3_all, df3_1bp)
df4_phast <- rbind.data.frame(dplyr::select(df4_A, -median_Bs), df4_B)

########################
#
# analyzing mutation rates within / outside genes split-exons
#
########################

# exons split
df_maps <- vector("list", length=length(maps) * 22)
pos <- 1
for(m in maps) {
  
  print(m)

  chr_maps <- vector("list", 22)
  pb <- txtProgressBar(min=0, max=22, style=3)
  for(c in 1:22) {
    
    setTxtProgressBar(pb, c)
    
    mut_map <- fread(paste("bgs_lmr/human_data/data/mutation_tables/", m,
                            "/", m, "_tbl_chr", c, ".csv.gz", sep="")) 
        
    deciles <- vector("list", 11)
    for(d in 1:11) {
      tmp <- fread(paste("bgs_lmr/human_data/data/annotation_tables/", m, "/functional/exons_lof_deciles", 
                         "/d", d, "/", "exons_lof_d", d, "_tbl_chr", c, ".csv.gz", sep=""))
        
      tmp$elem <- paste0("decile_", d)
      deciles[[d]] <- tmp
    }
    
    mut_elem <- data.table::rbindlist(deciles)
    mut_elem_w <- pivot_wider(mut_elem, values_from=c(scale_masked, del_sites_masked), names_from=elem, names_prefix="elem_") %>% setDT()
    mut_elem_w[
        , del_sites_sum := rowSums(.SD, na.rm = TRUE),
        .SDcols = names(mut_elem_w)[startsWith(names(mut_elem_w), "del_sites_masked_elem")]
    ]
    mut_elem_w <- mut_elem_w[del_sites_sum > 0,]

    dt <- merge(mut_elem_w, mut_map, by=c("chrom", "chromStart", "chromEnd"))
    dt <- dt[!is.na(avg_mut_masked) & avg_mut_masked > 0 & num_sites_masked > 0]

    sum_mut_bin <- dt[["avg_mut_masked"]] * dt[["num_sites_masked"]]

    scale_cols <- grep("scale_masked_elem", names(dt), value = TRUE)
    del_cols <- grep("del_sites_masked_elem", names(dt), value = TRUE)
    sum_mut_del <- as.data.table(
      Map(function(s, d) dt$avg_mut_masked * dt[[s]] * dt[[d]], scale_cols, del_cols)
    )
    setnames(sum_mut_del, paste0("sum_mut_del_", seq_along(scale_cols)))
    sum_mut_del$elem_U <- rowSums(sum_mut_del, na.rm=T)
    
    avg_mut_out <- (sum_mut_bin - sum_mut_del$elem_U) / (dt$num_sites_masked - dt$del_sites_sum)
    avg_mut_out[!is.finite(avg_mut_out)] <- NA
    
    avg_mut_in <- as.data.table(
      Map(function(s, d) dt$avg_mut_masked * dt[[s]], scale_cols)
    )
    setnames(avg_mut_in, paste0("avg_mut_del_", seq_along(scale_cols)))
    
    dt[, Outside := avg_mut_out]
    dt <- cbind.data.frame(dt, avg_mut_in) %>% dplyr::select(., -c(num_sites, del_sites, starts_with("scale")))
    dt$map <- m

    df_maps[[pos]] <- dt
    pos <- pos + 1
  }
  close(pb)
}

df_exons <- data.table::rbindlist(df_maps)

summarize_exon_deciles <- function(dt, exon_prefix="avg_mut_del_", d=1:11) {
  
  res <- lapply(d, function(i) {
    
    col <- paste0(exon_prefix, i)
    tmp <- dt[!is.na(get(col))]
    
    list(decile=i,
         median_muts_in=median(tmp[[col]], na.rm=TRUE),
         median_muts_out=median(tmp[["Outside"]], na.rm=TRUE),
         se_muts_in=sd(tmp[[col]], na.rm=TRUE) / sqrt(nrow(tmp)),
         se_muts_out=sd(tmp[["Outside"]], na.rm=TRUE) / sqrt(nrow(tmp)))})
  
  data.table::rbindlist(res)
}

exons_gnomad <- summarize_exon_deciles(filter(df_exons, map=="gnomad"))
exons_gnomad$map <- "gnomad"

exons_carlson <- summarize_exon_deciles(filter(df_exons, map=="carlson"))
exons_carlson$map <- "carlson"

exons_roulette <- summarize_exon_deciles(filter(df_exons, map=="roulette"))
exons_roulette$map <- "roulette"

df_exons <- rbind.data.frame(exons_gnomad, exons_carlson, exons_roulette) 
df_exons$control_gnomad <- FALSE

df2 <- df_exons %>% group_by(map) %>%
  mutate(median_ratio=median_muts_in / median_muts_out,
         se_ratio=median_ratio * sqrt((se_muts_in / median_muts_in)^2 + (se_muts_out / median_muts_out)^2))

df_exons <- df2 %>% group_by(map) %>%
  pivot_longer(cols=c(median_ratio, se_ratio),
               names_to=c(".value", "metric"),
               names_pattern="(median|se)_(.*)") %>% setDT()
df_exons$elem <- paste0("decile_", df_exons$decile)
df_exons[, decile := NULL]

# promoters
df_maps <- vector("list", length=length(maps) * 22)
pos <- 1
for(m in maps) {
  
  print(m)
  
  chr_maps <- vector("list", 22)
  pb <- txtProgressBar(min=0, max=22, style=3)
  for(c in 1:22) {
    
    setTxtProgressBar(pb, c)
    
    mut_map <- fread(paste("bgs_lmr/human_data/data/mutation_tables/", m,
                           "/", m, "_tbl_chr", c, ".csv.gz", sep="")) 
    
    mut_elem <- fread(paste("bgs_lmr/human_data/data/annotation_tables/", m, "/functional/promoters/", 
                      "/promoters_tbl_chr", c, ".csv.gz", sep=""))
      
    mut_elem$elem <- "promoter"
    mut_elem$chrom <- c
    mut_elem$chrom <- as.integer(mut_elem$chrom)
    mut_elem <- mut_elem[del_sites_masked > 0,]
    
    dt <- merge(mut_elem, mut_map, by=c("chrom", "chromStart", "chromEnd"))
    dt <- dt[!is.na(avg_mut_masked) & avg_mut_masked > 0 & num_sites_masked > 0]
    
    sum_mut_bin <- dt[["avg_mut_masked"]] * dt[["num_sites_masked"]]
    avg_mut_in <- dt[["avg_mut_masked"]] * dt[["scale_masked"]]
    avg_mut_out <- (sum_mut_bin - avg_mut_in * dt[["del_sites_masked"]]) / (dt$num_sites_masked - dt$del_sites_masked)
    avg_mut_out[!is.finite(avg_mut_out)] <- NA
    
    dt[, Outside := avg_mut_out]
    dt[, Within := avg_mut_in]
    dt <- dplyr::select(dt, -c(num_sites, del_sites, starts_with("scale")))
    dt$map <- m
    
    df_maps[[pos]] <- dt
    pos <- pos + 1
  }
  close(pb)
}

df_promoters <- data.table::rbindlist(df_maps)

df_promoters <- df_promoters %>%
  group_by(map) %>%
  summarise(median_muts_in=median(Within, na.rm=TRUE),
            median_muts_out=median(Outside, na.rm=TRUE),
            se_muts_in=sd(Within, na.rm=TRUE) / sqrt(n()),
            se_muts_out=sd(Outside, na.rm=TRUE) / sqrt(n()),
            .groups="drop") %>%
  collect()

df2 <- df_promoters %>% group_by(map) %>%
  mutate(median_ratio=median_muts_in / median_muts_out,
         se_ratio=median_ratio * sqrt((se_muts_in / median_muts_in)^2 + (se_muts_out / median_muts_out)^2))

df_promoters <- df2 %>% group_by(map) %>%
  pivot_longer(cols=c(median_ratio, se_ratio),
               names_to=c(".value", "metric"),
               names_pattern="(median|se)_(.*)") %>% setDT()
df_promoters$elem <- "promoter"
df_promoters$control_gnomad <- FALSE

# enhancers
df_maps <- vector("list", length=length(maps) * 22)
pos <- 1
for(m in maps) {
  
  print(m)
  
  chr_maps <- vector("list", 22)
  pb <- txtProgressBar(min=0, max=22, style=3)
  for(c in 1:22) {
    
    setTxtProgressBar(pb, c)
    
    mut_map <- fread(paste("bgs_lmr/human_data/data/mutation_tables/", m,
                           "/", m, "_tbl_chr", c, ".csv.gz", sep="")) 
    
    mut_elem <- fread(paste("bgs_lmr/human_data/data/annotation_tables/", m, "/functional/enhancers/", 
                            "/enhancers_tbl_chr", c, ".csv.gz", sep=""))
    
    mut_elem$elem <- "enhancer"
    mut_elem$chrom <- c
    mut_elem$chrom <- as.integer(mut_elem$chrom)
    mut_elem <- mut_elem[del_sites_masked > 0,]
    
    dt <- merge(mut_elem, mut_map, by=c("chrom", "chromStart", "chromEnd"))
    dt <- dt[!is.na(avg_mut_masked) & avg_mut_masked > 0 & num_sites_masked > 0]
    
    sum_mut_bin <- dt[["avg_mut_masked"]] * dt[["num_sites_masked"]]
    avg_mut_in <- dt[["avg_mut_masked"]] * dt[["scale_masked"]]
    avg_mut_out <- (sum_mut_bin - avg_mut_in * dt[["del_sites_masked"]]) / (dt$num_sites_masked - dt$del_sites_masked)
    avg_mut_out[!is.finite(avg_mut_out)] <- NA
    
    dt[, Outside := avg_mut_out]
    dt[, Within := avg_mut_in]
    dt <- dplyr::select(dt, -c(num_sites, del_sites, starts_with("scale")))
    dt$map <- m
    
    df_maps[[pos]] <- dt
    pos <- pos + 1
  }
  close(pb)
}

df_enhancers <- data.table::rbindlist(df_maps)

df_enhancers <- df_enhancers %>%
  group_by(map) %>%
  summarise(median_muts_in=median(Within, na.rm=TRUE),
            median_muts_out=median(Outside, na.rm=TRUE),
            se_muts_in=sd(Within, na.rm=TRUE) / sqrt(n()),
            se_muts_out=sd(Outside, na.rm=TRUE) / sqrt(n()),
            .groups="drop") %>%
  collect()

df2 <- df_enhancers %>% group_by(map) %>%
  mutate(median_ratio=median_muts_in / median_muts_out,
         se_ratio=median_ratio * sqrt((se_muts_in / median_muts_in)^2 + (se_muts_out / median_muts_out)^2))

df_enhancers <- df2 %>% group_by(map) %>%
  pivot_longer(cols=c(median_ratio, se_ratio),
               names_to=c(".value", "metric"),
               names_pattern="(median|se)_(.*)") %>% setDT()
df_enhancers$elem <- "enhancer"
df_enhancers$control_gnomad <- FALSE

df_genes_A <- rbind.data.frame(df_exons, df_promoters, df_enhancers)

########################
#
# controlling for sequence composition by dividing roulette and carlson by gnomad
#
########################

# exons split
chr_maps <- vector("list", 22*2)
pb <- txtProgressBar(min=0, max=22, style=3)
for(c in 1:22) {
  
  setTxtProgressBar(pb, c)
  
  carlson <- fread(paste("bgs_lmr/human_data/data/mutation_tables/carlson/carlson_tbl_chr", c, ".csv.gz", sep="")) %>% setDT()
  roulette <- fread(paste("bgs_lmr/human_data/data/mutation_tables/roulette/roulette_tbl_chr", c, ".csv.gz", sep="")) %>% setDT()
  gnomad <- fread(paste("bgs_lmr/human_data/data/mutation_tables/gnomad/gnomad_tbl_chr", c, ".csv.gz", sep="")) %>% setDT()
  
  roulette[, rel_mut := avg_mut / gnomad$avg_mut]
  carlson[, rel_mut := avg_mut / gnomad$avg_mut]
  
  roulette[, map := "roulette"]
  carlson[, map := "carlson"]
  
  carlson_roulette <- list(carlson, roulette)
  
  tbl <- vector("list", length(carlson_roulette) * 11) # 11 deciles
  pos <- 1 
  for(i in 1:length(carlson_roulette)) {
    
    mut_map <- carlson_roulette[[i]]
    m <- unique(mut_map$map)
    
    deciles <- vector("list", 11)
    for(d in 1:11) {
      tmp <- fread(paste("bgs_lmr/human_data/data/annotation_tables/", m, "/functional/exons_lof_deciles", 
                         "/d", d, "/", "exons_lof_d", d, "_tbl_chr", c, ".csv.gz", sep=""))
      
      tmp$elem <- paste0("decile_", d)
      deciles[[d]] <- tmp
    }
    
    mut_elem <- data.table::rbindlist(deciles)
    mut_elem_w <- pivot_wider(mut_elem, values_from=c(scale_masked, del_sites_masked), names_from=elem, names_prefix="elem_") %>% setDT()
    mut_elem_w[
      , del_sites_sum := rowSums(.SD, na.rm = TRUE),
      .SDcols = names(mut_elem_w)[startsWith(names(mut_elem_w), "del_sites_masked_elem")]
    ]
    mut_elem_w <- mut_elem_w[del_sites_sum > 0,]
    
    dt <- merge(mut_elem_w, mut_map, by=c("chrom", "chromStart", "chromEnd"))
    dt <- dt[!is.na(avg_mut_masked) & avg_mut_masked > 0 & num_sites_masked > 0]
    
    sum_mut_bin <- dt[["rel_mut"]] * dt[["num_sites_masked"]]
    
    scale_cols <- grep("scale_masked_elem", names(dt), value = TRUE)
    del_cols <- grep("del_sites_masked_elem", names(dt), value = TRUE)
    sum_mut_del <- as.data.table(
      Map(function(s, d) dt$rel_mut * dt[[s]] * dt[[d]],
          scale_cols, del_cols)
    )
    setnames(sum_mut_del, paste0("sum_mut_del_", seq_along(scale_cols)))
    sum_mut_del$elem_U <- rowSums(sum_mut_del, na.rm=T)
    
    avg_mut_out <- (sum_mut_bin - sum_mut_del$elem_U) / (dt$num_sites_masked - dt$del_sites_sum)
    avg_mut_out[!is.finite(avg_mut_out)] <- NA
    
    avg_mut_in <- as.data.table(
      Map(function(s, d) dt$rel_mut * dt[[s]], scale_cols)
    )
    setnames(avg_mut_in, paste0("avg_mut_del_", seq_along(scale_cols)))
    
    dt[, Outside := avg_mut_out]
    dt <- cbind.data.frame(dt, avg_mut_in) %>% dplyr::select(., -c(num_sites, del_sites, starts_with("scale")))
    dt$map <- m
    
    chr_maps[[pos]] <- dt
    pos <- pos + 1
  }
}
close(pb)

df_exons <- data.table::rbindlist(chr_maps)

summarize_exon_deciles <- function(dt, exon_prefix="avg_mut_del_", d=1:11) {
  
  res <- lapply(d, function(i) {
    
    col <- paste0(exon_prefix, i)
    tmp <- dt[!is.na(get(col))]
    
    list(decile=i,
         median_muts_in=median(tmp[[col]], na.rm=TRUE),
         median_muts_out=median(tmp[["Outside"]], na.rm=TRUE),
         se_muts_in=sd(tmp[[col]], na.rm=TRUE) / sqrt(nrow(tmp)),
         se_muts_out=sd(tmp[["Outside"]], na.rm=TRUE) / sqrt(nrow(tmp)))})
  
  data.table::rbindlist(res)
}

exons_carlson <- summarize_exon_deciles(filter(df_exons, map=="carlson"))
exons_carlson$map <- "carlson"

exons_roulette <- summarize_exon_deciles(filter(df_exons, map=="roulette"))
exons_roulette$map <- "roulette"

df_exons <- rbind.data.frame(exons_carlson, exons_roulette) 
df_exons$control_gnomad <- TRUE

df2 <- df_exons %>% group_by(map) %>%
  mutate(median_ratio=median_muts_in / median_muts_out,
         se_ratio=median_ratio * sqrt((se_muts_in / median_muts_in)^2 + (se_muts_out / median_muts_out)^2))

df_exons <- df2 %>% group_by(map) %>%
  pivot_longer(cols=c(median_ratio, se_ratio),
               names_to=c(".value", "metric"),
               names_pattern="(median|se)_(.*)") %>% setDT()
df_exons$elem <- paste0("decile_", df_exons$decile)
df_exons[, decile := NULL]

# promoters
chr_maps <- vector("list", 22*2)
pos <- 1
pb <- txtProgressBar(min=0, max=22, style=3)
for(c in 1:22) {
  
  setTxtProgressBar(pb, c)
  
  carlson <- fread(paste("bgs_lmr/human_data/data/mutation_tables/carlson/carlson_tbl_chr", c, ".csv.gz", sep="")) %>% setDT()
  roulette <- fread(paste("bgs_lmr/human_data/data/mutation_tables/roulette/roulette_tbl_chr", c, ".csv.gz", sep="")) %>% setDT()
  gnomad <- fread(paste("bgs_lmr/human_data/data/mutation_tables/gnomad/gnomad_tbl_chr", c, ".csv.gz", sep="")) %>% setDT()
  
  roulette[, rel_mut := avg_mut / gnomad$avg_mut]
  carlson[, rel_mut := avg_mut / gnomad$avg_mut]
  
  roulette[, map := "roulette"]
  carlson[, map := "carlson"]
  
  carlson_roulette <- list(carlson, roulette)
  
  for(i in 1:length(carlson_roulette)) {
    
    mut_map <- carlson_roulette[[i]]
    m <- unique(mut_map$map)
    
    mut_elem <- fread(paste("bgs_lmr/human_data/data/annotation_tables/", m, "/functional/promoters/", 
                            "/promoters_tbl_chr", c, ".csv.gz", sep=""))
    
    mut_elem$elem <- "promoter"
    mut_elem$chrom <- c
    mut_elem$chrom <- as.integer(mut_elem$chrom)
    mut_elem <- mut_elem[del_sites_masked > 0,]
    
    dt <- merge(mut_elem, mut_map, by=c("chrom", "chromStart", "chromEnd"))
    dt <- dt[!is.na(avg_mut_masked) & avg_mut_masked > 0 & num_sites_masked > 0]
    
    sum_mut_bin <- dt[["rel_mut"]] * dt[["num_sites_masked"]]
    avg_mut_in <- dt[["rel_mut"]] * dt[["scale_masked"]]
    avg_mut_out <- (sum_mut_bin - avg_mut_in * dt[["del_sites_masked"]]) / (dt$num_sites_masked - dt$del_sites_masked)
    avg_mut_out[!is.finite(avg_mut_out)] <- NA
    
    dt[, Outside := avg_mut_out]
    dt[, Within := avg_mut_in]
    dt <- dplyr::select(dt, -c(num_sites, del_sites, starts_with("scale")))
    dt$map <- m
    
    chr_maps[[pos]] <- dt
    pos <- pos + 1
  }
}
close(pb)

df_promoters <- data.table::rbindlist(chr_maps)

df_promoters <- df_promoters %>%
  group_by(map) %>%
  summarise(median_muts_in=median(Within, na.rm=TRUE),
            median_muts_out=median(Outside, na.rm=TRUE),
            se_muts_in=sd(Within, na.rm=TRUE) / sqrt(n()),
            se_muts_out=sd(Outside, na.rm=TRUE) / sqrt(n()),
            .groups="drop") %>%
  collect()

df2 <- df_promoters %>% group_by(map) %>%
  mutate(median_ratio=median_muts_in / median_muts_out,
         se_ratio=median_ratio * sqrt((se_muts_in / median_muts_in)^2 + (se_muts_out / median_muts_out)^2))

df_promoters <- df2 %>% group_by(map) %>%
  pivot_longer(cols=c(median_ratio, se_ratio),
               names_to=c(".value", "metric"),
               names_pattern="(median|se)_(.*)") %>% setDT()
df_promoters$elem <- "promoter"
df_promoters$control_gnomad <- TRUE

# enhancers
chr_maps <- vector("list", 22*2)
pos <- 1
pb <- txtProgressBar(min=0, max=22, style=3)
for(c in 1:22) {
  
  setTxtProgressBar(pb, c)
  
  carlson <- fread(paste("bgs_lmr/human_data/data/mutation_tables/carlson/carlson_tbl_chr", c, ".csv.gz", sep="")) %>% setDT()
  roulette <- fread(paste("bgs_lmr/human_data/data/mutation_tables/roulette/roulette_tbl_chr", c, ".csv.gz", sep="")) %>% setDT()
  gnomad <- fread(paste("bgs_lmr/human_data/data/mutation_tables/gnomad/gnomad_tbl_chr", c, ".csv.gz", sep="")) %>% setDT()
  
  roulette[, rel_mut := avg_mut / gnomad$avg_mut]
  carlson[, rel_mut := avg_mut / gnomad$avg_mut]
  
  roulette[, map := "roulette"]
  carlson[, map := "carlson"]
  
  carlson_roulette <- list(carlson, roulette)
  
  for(i in 1:length(carlson_roulette)) {
    
    mut_map <- carlson_roulette[[i]]
    m <- unique(mut_map$map)
    
    mut_elem <- fread(paste("bgs_lmr/human_data/data/annotation_tables/", m, "/functional/enhancers/", 
                            "/enhancers_tbl_chr", c, ".csv.gz", sep=""))
    
    mut_elem$elem <- "enhancer"
    mut_elem$chrom <- c
    mut_elem$chrom <- as.integer(mut_elem$chrom)
    mut_elem <- mut_elem[del_sites_masked > 0,]
    
    dt <- merge(mut_elem, mut_map, by=c("chrom", "chromStart", "chromEnd"))
    dt <- dt[!is.na(avg_mut_masked) & avg_mut_masked > 0 & num_sites_masked > 0]
    
    sum_mut_bin <- dt[["rel_mut"]] * dt[["num_sites_masked"]]
    avg_mut_in <- dt[["rel_mut"]] * dt[["scale_masked"]]
    avg_mut_out <- (sum_mut_bin - avg_mut_in * dt[["del_sites_masked"]]) / (dt$num_sites_masked - dt$del_sites_masked)
    avg_mut_out[!is.finite(avg_mut_out)] <- NA
    
    dt[, Outside := avg_mut_out]
    dt[, Within := avg_mut_in]
    dt <- dplyr::select(dt, -c(num_sites, del_sites, starts_with("scale")))
    dt$map <- m
    
    chr_maps[[pos]] <- dt
    pos <- pos + 1
  }
}
close(pb)

df_enhancers <- data.table::rbindlist(chr_maps)

df_enhancers <- df_enhancers %>%
  group_by(map) %>%
  summarise(median_muts_in=median(Within, na.rm=TRUE),
            median_muts_out=median(Outside, na.rm=TRUE),
            se_muts_in=sd(Within, na.rm=TRUE) / sqrt(n()),
            se_muts_out=sd(Outside, na.rm=TRUE) / sqrt(n()),
            .groups="drop") %>%
  collect()

df2 <- df_enhancers %>% group_by(map) %>%
  mutate(median_ratio=median_muts_in / median_muts_out,
         se_ratio=median_ratio * sqrt((se_muts_in / median_muts_in)^2 + (se_muts_out / median_muts_out)^2))

df_enhancers <- df2 %>% group_by(map) %>%
  pivot_longer(cols=c(median_ratio, se_ratio),
               names_to=c(".value", "metric"),
               names_pattern="(median|se)_(.*)") %>% setDT()
df_enhancers$elem <- "enhancer"
df_enhancers$control_gnomad <- TRUE

df_genes_B <- rbind.data.frame(df_exons, df_promoters, df_enhancers) %>% setDT()

df_genes <- rbind.data.frame(df_genes_A, df_genes_B)
df_genes$group <- paste0(df_genes$map, as.numeric(str_detect(df_genes$elem, "decile")))

df_genes[, elem := factor(elem, levels = c(paste0("decile_", 1:11), "enhancer", "promoter"))]

########################
#
# plotting
#
########################

p1 <- ggplot(data=filter(df4_phast, elem=="phast_all"), aes(x = phast_bin, y = median, color = map)) +
  geom_point(size=2.5) + geom_line(aes(linetype=control_gnomad)) +
  geom_errorbar(aes(ymin = median - se, ymax = median + se),
                width = 0.2, linewidth = 0.7) +
  labs(x="phastCons Bin", y=expression(paste(mu, " ratio"))) +
  scale_x_continuous(breaks=1:12) +
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_color_manual(values=c("cyan3", "green4", "brown1"), name=NULL) +
  theme_classic() + 
  theme(axis.title=element_text(size = 20),
        axis.text=element_text(size = 16),
        legend.text=element_text(size = 16), 
        legend.title=element_text(size = 16), 
        legend.position="none")

p2 <- ggplot(data=filter(df4_phast, elem=="phast_1bp"), aes(x = phast_bin, y = median, color = map)) +
  geom_point(size=2.5) + geom_line(aes(linetype=control_gnomad)) +
  geom_errorbar(aes(ymin = median - se, ymax = median + se),
                width = 0.2, linewidth = 0.7) +
  labs(x="phastCons Bin", y=NULL) +
  scale_x_continuous(breaks=1:12) +
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_color_manual(values=c("cyan3", "green4", "brown1"), name=NULL) +
  theme_classic() + 
  theme(axis.title=element_text(size = 20),
        axis.text=element_text(size = 16),
        legend.text=element_text(size = 16), 
        legend.title=element_text(size = 16), 
        legend.position="bottom")

legend <- cowplot::get_legend(p2)
p2 <- p2 + theme(legend.position="none")

p3 <- ggplot(df_genes, aes(x=elem, y=median)) +
  geom_point(aes(color=map, group=interaction(map, control_gnomad, group)), size=2.5) +
  geom_line(aes(color=map, linetype=control_gnomad, group=interaction(map, control_gnomad, group)),
            linewidth=0.8) +
  geom_errorbar(aes(ymin=median - se,
                    ymax=median + se, 
                    color=map),
                width=0.2, linewidth=0.7) +
  labs(x="Genic element", y=NULL) +
  scale_x_discrete(labels = c(1:11, "e", "p")) +
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_color_manual(values=c("cyan3", "green4", "brown1"), name=NULL) +
  theme_classic() + 
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=16),
        legend.text=element_text(size=16), 
        legend.title=element_text(size=16), 
        legend.position="none")

psingle <- plot_grid(plot_grid(p1, p2, p3, nrow=1,
                               rel_widths=c(1.2, 1, 1)), 
                     legend, nrow=2, rel_heights=c(1, 0.1))

empty <- ggplot() + theme_void()
pfull <- plot_grid(psingle, empty, d1, nrow = 3,
                   labels = c("A", "", "B"),
                   rel_heights = c(1, 0.05, 1))

save_plot("~/Desktop/single_mut_ratios_elems.pdf", pfull, base_height=10, base_width=17)        
save_plot("~/Desktop/single_mut_ratios_elems.png", pfull, base_height=10, base_width=17)  

#########################
#
# checking mutation rates in bins without any phastcons elements 
# TODO adapt
#
#########################

ds <- open_dataset("phast_wide")
big_dt <- ds |> collect() 
big_dt <- setDT(big_dt)

big_dt$Far <- big_dt$avg_mut
big_dt$Far[big_dt$del_sites > 0] <- NA
big_dt[, del_sites := NULL]
big_dt[, avg_mut := NULL]
big_dt[, Phast := NULL]

tbl <- big_dt[, .SD[sample(.N, .N * 0.25)], by = .(Outside, Far, map)]

tbl <- pivot_longer(tbl, cols=c("Within", "Outside", "Far"),
                    names_to="location", values_to="mut_rate")

wilcox_results_phast <- tbl %>%
  group_by(map, phast_bin, location) %>%
  collect() %>%                 
  group_by(map, phast_bin) %>%
  summarise(
    within  = list(mut_rate[location == "Within"]),
    outside = list(mut_rate[location == "Outside"]),
    far     = list(mut_rate[location == "Far"]),
    .groups = "drop") %>%
  mutate(
    p_within_outside = wilcox.test(within[[1]], outside[[1]])$p.value,
    p_within_far     = wilcox.test(within[[1]], far[[1]])$p.value,
    p_outside_far    = wilcox.test(outside[[1]], far[[1]])$p.value) %>%
  mutate(
    stars_within_outside = case_when(
      #p_within_outside < 0.001 / 36 ~ "***",
      #p_within_outside < 0.01 / 36 ~ "**",
      p_within_outside < 0.05 / 36 ~ "*",
      TRUE ~ "ns"),
    stars_outside_far = case_when(
      #p_outside_far < 0.001 / 36 ~ "***",
      #p_outside_far < 0.01  / 36 ~ "**",
      p_outside_far < 0.05  / 36 ~ "*",
      TRUE ~ "ns"))

wilcox_for_labels <- wilcox_results_phast %>%
  select(map, phast_bin, stars_within_outside, stars_outside_far) %>%
  tidyr::pivot_longer( cols = starts_with("stars_"), names_to = "comparison", values_to = "p_stars" )

box_stats_phast <- tbl %>%
  group_by(map, phast_bin, location) %>%
  summarise(q1 = quantile(mut_rate, 0.25, na.rm = TRUE),
            med = median(mut_rate, na.rm = TRUE),
            q3 = quantile(mut_rate, 0.75, na.rm = TRUE),
            iqr = q3 - q1,
            lower_whisker = max(min(mut_rate, na.rm = TRUE), q1 - 1.5 * iqr),
            upper_whisker = min(max(mut_rate, na.rm = TRUE), q3 + 1.5 * iqr)) %>%
  collect()

box_stats_phast <- box_stats_phast %>%
  mutate(map = as.character(map), phast_bin = as.character(phast_bin))

wilcox_for_labels <- wilcox_for_labels %>%
  mutate(map = as.character(map), phast_bin = as.character(phast_bin))

cap_df <- box_stats_phast %>%
  group_by(phast_bin) %>%
  summarise(cap95 = quantile(upper_whisker, 0.95, na.rm = TRUE), .groups = "drop")

label_df_phast <- box_stats_phast %>%
  group_by(map, phast_bin) %>%
  summarise(max_whisk = max(upper_whisker, na.rm = TRUE), .groups = "drop" ) %>%
  left_join(wilcox_for_labels, by = c("map", "phast_bin")) %>%
  left_join(cap_df, by = "phast_bin") %>%
  mutate( p_stars = replace_na(p_stars, "ns"),
          y_raw = max_whisk * 1.05, y = pmin(y_raw, cap95 * 1.05) )

label_df_phast <- label_df_phast %>%
  mutate(x_offset = case_when(comparison == "stars_within_outside" ~ -0.15,
                              comparison == "stars_outside_far" ~ +0.15, TRUE ~ 0))

box_stats_phast$map_num <- as.numeric(factor(box_stats_phast$map))
label_df_phast$map_num <- as.numeric(factor(label_df_phast$map))

y_min <- min(box_stats_phast$lower_whisker, na.rm = TRUE)
y_max_data <- max(box_stats_phast$upper_whisker, na.rm = TRUE)
y_max_label <- max(label_df_phast$y, na.rm = TRUE)
y_top <- max(y_max_data, y_max_label) * 1.025

box_stats_phast <- box_stats_phast %>%
  mutate(phast_bin = as.numeric(as.character(phast_bin)))
label_df_phast <- label_df_phast %>%
  mutate(phast_bin = as.numeric(as.character(phast_bin)))
wilcox_results_phast <- wilcox_results_phast %>%
  mutate(phast_bin = as.numeric(as.character(phast_bin)))

levels_order <- sort(unique(box_stats_phast$phast_bin))
box_stats_phast <- box_stats_phast %>%
  mutate(phast_bin = factor(phast_bin, levels = levels_order, ordered = TRUE))
label_df_phast  <- label_df_phast  %>%
  mutate(phast_bin = factor(phast_bin,  levels = levels_order, ordered = TRUE))
wilcox_results_phast <- wilcox_results_phast %>%
  mutate(phast_bin = factor(phast_bin, levels = levels_order, ordered = TRUE))

b3 <- ggplot(box_stats_phast, aes(x = map_num, fill = location)) +
  coord_cartesian(ylim = c(y_min, y_top)) +
  geom_boxplot(
    aes(ymin = lower_whisker, lower = q1, middle = med,
        upper = q3, ymax = upper_whisker,
        group = interaction(map, location)),
    stat = "identity",
    position = position_dodge(width = 0.75),
    alpha = 0.6,
    width = 0.6) +
  facet_wrap(~phast_bin, scales = "free_y", ncol = 3) +
  geom_text(data = label_df_phast,
            aes(x = map_num + x_offset, y = y, label = p_stars),
            inherit.aes = FALSE,
            size = 6,
            fontface = "bold",
            vjust = 0) +
  labs(x = NULL, y = expression(mu)) +
  scale_fill_manual(values = c("cyan3", "brown1", "green4"), name = NULL) +
  scale_x_continuous(
    breaks = unique(box_stats_phast$map_num),
    labels = unique(box_stats_phast$map)) +
  theme_classic() +
  theme(strip.text = element_text(face = "bold"),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = "bottom")

save_plot("~/Desktop/phast_mut_diff_extended.pdf", b3, base_height=12, base_width=10)




###########################
#
# phastcons... again
#
###########################

phast <- vector("list", 12) # 12 phastcons bins
for(l in seq(0, 55, 5)) {
  
  u <- l + 5
  
  tmp <- fread(paste("bgs_lmr/human_data/data/annotations/phastcons/top", l, "-", u,
                     "/phastcons_top", l, "-", u, "_chr1.bed.gz", sep="")) 
  
  tmp <- tmp[chromEnd <= last_pos,] # trimming as a function of the other data
  
  tmp$phast_bin <- u / 5
  tmp$chrom <- 1
  tmp$chrom <- as.integer(tmp$chrom)

  phast[[l %/% 5 + 1]] <- tmp
}

phast <- data.table::rbindlist(phast) %>% setDT()

# preparing tables for genomic ranges manipulation
setnames(phast, c("chromStart", "chromEnd"), c("start", "end"))
phast[, end := end -1]

dat[, start := pos]
dat[, end := pos]

setkey(dat, chrom, start, end)
setkey(phast, chrom, start, end)

dat <- foverlaps(dat, phast)
dat[, c("start", "end", "i.start", "i.end") := NULL]

dat[is.na(phast_bin), phast_bin := 15] # non-phastcons sites

# stacked plot x=phast_bin, colors are triplets
plot_df <- dat[, .N, by=.(triplet, phast_bin)]
plot_df[, prop := N / sum(N), by=phast_bin]

plot_df[, phast_bin := factor(phast_bin, levels=sort(unique(phast_bin)))]
trip_levels <- sort(unique(plot_df$triplet))   
plot_df[, triplet := factor(triplet, levels = trip_levels)]

all_grid <- CJ(phast_bin = levels(plot_df$phast_bin),
               triplet  = levels(plot_df$triplet),
               unique = TRUE)
plot_full <- merge(all_grid, plot_df, by = c("phast_bin", "triplet"), all.x = TRUE)
plot_full[is.na(N), N := 0]
plot_full[is.na(prop), prop := 0]

setorder(plot_full, phast_bin, triplet)

plot_full[, ymax := cumsum(prop), by = phast_bin]
plot_full[, ymin := shift(ymax, fill = 0), by = phast_bin]

width <- 0.9
halfw <- width / 2

cpg_tris <- c("ACG", "CGC", "TCG", "CGA", "CGC", "CGG", "CGT")

s2 <- ggplot(plot_df, aes(x=phast_bin, y=prop, fill=triplet)) +
  theme_classic() + geom_col() +
  # geom_rect(
  #   data = plot_full[triplet %in% cpg_tris & prop > 0],
  #   aes(
  #     xmin = as.integer(phast_bin) - halfw,
  #     xmax = as.integer(phast_bin) + halfw,
  #     ymin = ymin,
  #     ymax = ymax
  #   ),
  #   inherit.aes = FALSE,
  #   fill = NA,
  #   color = "black",
  #   linewidth = 0.6
  # ) +
  scale_fill_viridis_d(option="C", direction=1, guide=guide_legend(nrow=4), name=NULL) +
  scale_x_discrete(breaks=c(1:nbins, nbins + 3), labels=c(as.character(1:nbins), "Outside")) +
  labs(x="PhastCons Bin", y="Proportion", fill="Trinucleotide") +
  theme(axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12),
        legend.position="bottom",
        legend.box="horizontal")
save_plot("~/Desktop/mut_rates/phastCons_triplets.pdf", s2, base_height=7, base_width=10)

# repeat plots w/ phastcons
dat[, mean_roulette := mean(roulette, na.rm=T), by=.(phast_bin)]
dat[, mean_carlson := mean(carlson, na.rm=T), by=.(phast_bin)]
dat[, mean_gnomad := mean(gnomad, na.rm=T), by=.(phast_bin)]
dat[, se_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(phast_bin)]
dat[, se_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(phast_bin)]
dat[, se_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(phast_bin)]

unique_vals <- dat[, .(phast_bin, mean_roulette, mean_carlson, mean_gnomad, se_gnomad, se_carlson, se_roulette)]
unique_vals <- unique(unique_vals)

dt_m <- unique_vals %>%
  pivot_longer(cols=c(mean_roulette, mean_carlson, mean_gnomad, se_roulette, se_carlson, se_gnomad),
               names_to=c(".value", "source"), names_pattern="(mean|se)_(.*)")

p2 <- ggplot(dt_m, aes(x=phast_bin, y=mean, color=source, group=paste0(source, phast_bin < 15))) +
  geom_line(linewidth=1) + geom_point(size=3) + 
  geom_errorbar(aes(ymin=mean - se, ymax=mean + se), width=0.2, linewidth=0.6) +
  scale_x_continuous(breaks=c(1:12, 15), labels=c(as.character(1:12), "Outside")) + theme_classic() + 
  scale_color_manual(values=c("cyan3", "brown1", "green4"), name=NULL) +
  labs(x="PhastCons bin", y="Mean Rate", color="Map",
       title="Mutation rates per phastCons bin") +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        strip.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position="bottom")
save_plot("~/Desktop/mut_rates/phastcons_mut_rates.pdf", p2, base_height=6, base_width=10)

# computing ratios of mutation rates "within" and "outside" phastcons and benegas elements
dat[, c("mean_roulette", "mean_carlson", "mean_gnomad", "se_roulette", "se_carlson", "se_gnomad") := NULL]

## phastCons
dat[, mean_roulette_phast := mean(roulette, na.rm=T), by=.(phast_bin, triplet)]
dat[, mean_carlson_phast := mean(carlson, na.rm=T), by=.(phast_bin, triplet)]
dat[, mean_gnomad_phast := mean(gnomad, na.rm=T), by=.(phast_bin, triplet)]
dat[, se_roulette_phast := sd(roulette, na.rm=T) / sqrt(.N), by=.(phast_bin, triplet)]
dat[, se_carlson_phast := sd(carlson, na.rm=T) / sqrt(.N), by=.(phast_bin, triplet)]
dat[, se_gnomad_phast := sd(gnomad, na.rm=T) / sqrt(.N), by=.(phast_bin, triplet)]

denoms <- dat[phast_bin == 15, .(den_roulette=mean(mean_roulette_phast, na.rm=T),
                                    den_carlson=mean(mean_carlson_phast, na.rm=T),
                                    den_gnomad=mean(mean_gnomad_phast, na.rm=T)), by=triplet]

nums <- dat[phast_bin %in% 1:12, .(num_roulette=mean(mean_roulette_phast, na.rm=T),
                                      num_carlson=mean(mean_carlson_phast, na.rm=T),
                                      num_gnomad=mean(mean_gnomad_phast, na.rm=T)), by=.(triplet, phast_bin)]

ratios_phast <- merge(nums, denoms, by="triplet", all.x = TRUE)
setorder(ratios_phast, triplet, phast_bin)

ratios_phast[, `:=`(ratio_roulette=num_roulette / den_roulette,
                    ratio_carlson=num_carlson / den_carlson,
                    ratio_gnomad=num_gnomad / den_gnomad)]

ratios_phast[, c("num_roulette", "num_carlson", "num_gnomad", "den_roulette", "den_carlson", "den_gnomad") := NULL]

ratios_phast_m <- pivot_longer(ratios_phast, cols=starts_with("ratio_"),
                               names_to="map", values_to="ratio") %>% setDT()

ratios_phast_m[, map := factor(sub("^ratio_", "", map), levels = c("roulette", "carlson", "gnomad"))]

p3 <- ggplot(filter(ratios_phast_m, map!="gnomad"),
             aes(x=triplet, y=ratio, color=map)) + 
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_manual(values=c("brown1", "cyan3"), name=NULL) +
  labs(x=NULL, y="Ratio of mutation rates", color="Map", 
       title="Ratios of mutation rates across phastCons bins and triplets") +
  theme(axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.box="horizontal")
save_plot("~/Desktop/mut_rates/ratios_phast_triplet.pdf", p3, base_height=7, base_width=12)

p3b <- ggplot(filter(ratios_phast_m, map!="gnomad"),
              aes(x=triplet, y=ratio, color=phast_bin)) + 
  facet_wrap(~map, ncol=1) + 
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Bin", breaks=c(1, 12)) +
  labs(x=NULL, y="Ratio of mutation rates", color="Map", 
       title="Ratios of mutation rates across Benegas bins and triplets") +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title=element_text(size=18),
        legend.box="horizontal")

## benegas
dat[, mean_roulette_benegas := mean(roulette, na.rm=T), by=.(benegas_bin, triplet)]
dat[, mean_carlson_benegas := mean(carlson, na.rm=T), by=.(benegas_bin, triplet)]
dat[, mean_gnomad_benegas := mean(gnomad, na.rm=T), by=.(benegas_bin, triplet)]
dat[, se_roulette_benegas := sd(roulette, na.rm=T) / sqrt(.N), by=.(benegas_bin, triplet)]
dat[, se_carlson_benegas := sd(carlson, na.rm=T) / sqrt(.N), by=.(benegas_bin, triplet)]
dat[, se_gnomad_benegas := sd(gnomad, na.rm=T) / sqrt(.N), by=.(benegas_bin, triplet)]

denoms <- dat[benegas_bin == 15, .(den_roulette=mean(mean_roulette_benegas, na.rm=T),
                                    den_carlson=mean(mean_carlson_benegas, na.rm=T),
                                    den_gnomad=mean(mean_gnomad_benegas, na.rm=T)), by=triplet]

nums <- dat[benegas_bin %in% 1:12, .(num_roulette=mean(mean_roulette_benegas, na.rm=T),
                                      num_carlson=mean(mean_carlson_benegas, na.rm=T),
                                      num_gnomad=mean(mean_gnomad_benegas, na.rm=T)), by=.(triplet, benegas_bin)]

ratios_benegas <- merge(nums, denoms, by="triplet", all.x = TRUE)
setorder(ratios_benegas, triplet, benegas_bin)

ratios_benegas[, `:=`(ratio_roulette=num_roulette / den_roulette,
                     ratio_carlson=num_carlson / den_carlson,
                     ratio_gnomad=num_gnomad / den_gnomad)]

ratios_benegas[, c("num_roulette", "num_carlson", "num_gnomad", "den_roulette", "den_carlson", "den_gnomad") := NULL]
ratios_benegas_m <- pivot_longer(ratios_benegas, cols=starts_with("ratio_"), names_to="map", values_to="ratio") %>% setDT()

ratios_benegas_m[, map := factor(sub("^ratio_", "", map),
                               levels = c("roulette", "carlson", "gnomad"))]

p4 <- ggplot(filter(ratios_benegas_m, map!="gnomad"),
             aes(x=triplet, y=ratio, color=map)) + 
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_manual(values=c("brown1", "cyan3"), name=NULL) +
  labs(x=NULL, y="Ratio of mutation rates", color="Map", 
       title="Ratios of mutation rates across Benegas bins and triplets") +
  theme(axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.box="horizontal")
save_plot("~/Desktop/mut_rates/ratios_benegas_triplet.pdf", p4, base_height=7, base_width=12)

p4b <- ggplot(filter(ratios_benegas_m, map!="gnomad"),
             aes(x=triplet, y=ratio, color=benegas_bin)) + 
  facet_wrap(~map, ncol=1) + 
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Bin", breaks=c(1, 12)) +
  labs(x=NULL, y="Ratio of mutation rates", color="Map", 
       title="Ratios of mutation rates across Benegas bins and triplets") +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title=element_text(size=18),
        legend.box="horizontal")

# mutation rates per trinucleotide
dat[, `:=`(mean_roulette_triplets=mean(roulette, na.rm=T),
           mean_carlson_triplets=mean(carlson, na.rm=T),
           mean_gnomad_triplets=mean(gnomad, na.rm=T)), by=triplet]

tmp <- dat[!duplicated(triplet),] %>% 
  pivot_longer(., cols=ends_with("_triplets"), names_to="map", values_to="rate")

p5 <- ggplot(tmp, aes(x=triplet, y=rate, color=map)) + 
  geom_point(size=2.5) + theme_classic() + 
  scale_color_manual(values=c("cyan3", "green4", "brown1"), 
                     labels=c(mean_roulette_triplets="roulette", 
                              mean_carlson_triplets="carlson",
                              mean_gnomad_triplets="gnomad"), name=NULL) +
  labs(x=NULL, y=expression(mu), color="Map", 
       title="Mutation rates across triplets") +
  theme(axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, vjust=1, hjust=0.5),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.box="horizontal")
save_plot("~/Desktop/mut_rates/rates_triplets.pdf", p5, base_height=6, base_width=10)
