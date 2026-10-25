# =============================================================================
# Project 3: Cropping Cycles from Space - analysis in R
# Sentinel-2 NDVI/NDWI and CHIRPS rainfall, cropland of Bojonegoro, East Java
# -----------------------------------------------------------------------------
# Self-initiated portfolio project using public data.
# Author: Utik Tri Wulan Cahya, PhD
#
# Input (exported from Google Earth Engine with P3_cropping_cycles_GEE.js):
#   data/P3_timeseries_10day.csv   10-day NDVI, NDWI, rainfall, cloud-free share
#   data/P3_ndvi_monthly_2023.tif  monthly NDVI composites for 2023 (60 m)
#   data/P3_boundary.geojson       Bojonegoro Regency boundary
#
# Run from the project folder:  Rscript P3_cropping_cycles_analysis.R
# Input files are read from data/ if that folder exists, otherwise from the
# project folder itself; results are written to output/.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(terra)
  library(sf)
})

IN <- if (dir.exists("data")) "data" else "."
inp <- function(f) file.path(IN, f)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
sink("output/analysis_output.txt", split = TRUE)
section <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

# --- Settings ----------------------------------------------------------------
MIN_VALID   <- 0.5    # periods with < 50% cloud-free cropland are treated as unreliable
SPLINE_DF   <- 50     # flexibility of the smoothing spline (6 years of 10-day data)
MIN_PROM    <- 0.03   # minimum NDVI rise above the surrounding troughs to count as a crop cycle
AMP_FRAC    <- 0.20   # start/end of season = 20% of the trough-to-peak amplitude
MIN_SEP_D   <- 60     # minimum days between two crop-cycle peaks
ORIGIN      <- as.Date("2019-01-01")

# Neutral ink and a colour-vision-deficiency validated categorical palette
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e5e1"; muted <- "#b9b8b2"
col_ndvi <- "#2f7a3a"; col_rain <- "#2a78d6"
season_pal <- c("Wet season (MT-I)" = "#2a78d6", "Dry season I (MT-II)" = "#eb6834",
                "Dry season II (MT-III)" = "#1baf7a")

theme_pub <- function(base = 11) {
  theme_minimal(base_size = base, base_family = "DejaVu Sans") +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
          axis.title = element_text(colour = ink2),
          plot.title = element_text(face = "bold", size = base + 3),
          plot.subtitle = element_text(colour = ink2),
          plot.caption = element_text(colour = ink2, size = base - 2, hjust = 0),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = grid_col, linewidth = 0.35),
          legend.position = "top", legend.justification = "left",
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title.position = "plot", plot.caption.position = "plot")
}

# =============================================================================
section("1. DATA AND QUALITY CONTROL")
# =============================================================================
d <- read.csv(inp("P3_timeseries_10day.csv"))
d$date  <- as.Date(d$date_mid)
d$t     <- as.numeric(d$date - ORIGIN)
d$good  <- d$valid_share >= MIN_VALID
d$year  <- as.integer(format(d$date, "%Y"))
d$month <- as.integer(format(d$date, "%m"))

cat("10-day periods:", nrow(d), "(", min(d$date_start), "to", max(d$date_start), ")\n")
cat("Periods with < 50% cloud-free cropland (excluded from smoothing):",
    sum(!d$good), sprintf("(%.0f%%)", 100 * mean(!d$good)), "\n")
cat("Unreliable periods by month (all years):\n")
print(table(factor(month.abb[d$month[!d$good]], levels = month.abb)))

bnd <- st_read(inp("P3_boundary.geojson"), quiet = TRUE)
area_km2 <- as.numeric(st_area(st_transform(bnd, 32749))) / 1e6
r <- rast(inp("P3_ndvi_monthly_2023.tif"))
crop_mask <- !is.na(sum(!is.na(r)) > 0) & (sum(!is.na(r)) > 0)
crop_km2 <- global(crop_mask, "sum", na.rm = TRUE)[1, 1] * prod(res(r)) / 1e6
cat(sprintf("Regency area: %.0f km2; cropland with NDVI data (ESA WorldCover 2021): %.0f km2 (%.0f%%)\n",
            area_km2, crop_km2, 100 * crop_km2 / area_km2))

# =============================================================================
section("2. SMOOTHING AND CROP-CYCLE DETECTION")
# =============================================================================
g  <- d[d$good, ]
ss <- smooth.spline(g$t, g$NDVI, w = g$valid_share, df = SPLINE_DF)
sw <- smooth.spline(g$t, g$NDWI, w = g$valid_share, df = SPLINE_DF)
fit_rmse <- sqrt(mean((g$NDVI - predict(ss, g$t)$y)^2))
cat(sprintf("Weighted smoothing spline (df = %d): RMSE against reliable observations = %.3f NDVI\n",
            SPLINE_DF, fit_rmse))

tt <- seq(0, max(d$t), by = 1)
sm <- data.frame(t = tt, date = ORIGIN + tt,
                 NDVI = predict(ss, tt)$y, NDWI = predict(sw, tt)$y)

# Local maxima and minima of the smoothed curve
dy <- diff(sm$NDVI)
imax <- which(diff(sign(dy)) < 0) + 1
imin <- which(diff(sign(dy)) > 0) + 1

# Keep peaks with enough prominence above the troughs on either side
peaks <- data.frame()
for (p in imax) {
  left_min  <- if (any(imin < p)) max(imin[imin < p]) else 1
  right_min <- if (any(imin > p)) min(imin[imin > p]) else nrow(sm)
  prom <- sm$NDVI[p] - max(sm$NDVI[left_min], sm$NDVI[right_min])
  peaks <- rbind(peaks, data.frame(i_peak = p, i_start = left_min, i_end = right_min, prom = prom))
}
peaks <- peaks[peaks$prom >= MIN_PROM, ]
# Enforce minimum separation (keep the higher peak)
peaks <- peaks[order(-sm$NDVI[peaks$i_peak]), ]
keep <- c()
for (k in seq_len(nrow(peaks))) {
  if (all(abs(peaks$i_peak[k] - peaks$i_peak[keep]) >= MIN_SEP_D)) keep <- c(keep, k)
}
peaks <- peaks[sort(keep), ]
peaks <- peaks[order(peaks$i_peak), ]

season_of <- function(m) ifelse(m %in% c(12, 1, 2, 3, 4), "Wet season (MT-I)",
                         ifelse(m %in% 5:8, "Dry season I (MT-II)", "Dry season II (MT-III)"))
# Start of season (SOS) / end of season (EOS): first/last day the curve is
# above 20% of the trough-to-peak amplitude (threshold method, as in TIMESAT)
sos <- eos <- integer(nrow(peaks))
for (k in seq_len(nrow(peaks))) {
  a <- peaks$i_start[k]; p <- peaks$i_peak[k]; b <- peaks$i_end[k]
  thr_up <- sm$NDVI[a] + AMP_FRAC * (sm$NDVI[p] - sm$NDVI[a])
  thr_dn <- sm$NDVI[b] + AMP_FRAC * (sm$NDVI[p] - sm$NDVI[b])
  sos[k] <- a - 1 + which(sm$NDVI[a:p] >= thr_up)[1]
  eos[k] <- p - 1 + max(which(sm$NDVI[p:b] >= thr_dn))
}
seasons <- data.frame(
  trough_before = sm$date[peaks$i_start],
  start      = sm$date[sos],
  peak       = sm$date[peaks$i_peak],
  end        = sm$date[eos],
  trough_after = sm$date[peaks$i_end],
  peak_NDVI  = round(sm$NDVI[peaks$i_peak], 3),
  prominence = round(peaks$prom, 3)
)
seasons$season <- season_of(as.integer(format(seasons$peak, "%m")))
# Agricultural year runs from September to August and is assigned by the
# start date of each cycle (a cycle that starts in July still belongs to the
# agricultural year in progress, even if it peaks in September)
sm_ <- as.integer(format(seasons$start, "%m")); sy <- as.integer(format(seasons$start, "%Y"))
seasons$crop_year <- ifelse(sm_ >= 9, paste0(sy, "/", sy + 1), paste0(sy - 1, "/", sy))
seasons$green_up_to_peak_days <- as.integer(seasons$peak - seasons$start)
seasons$cycle_length_days     <- as.integer(seasons$end - seasons$start)
# Edge seasons (truncated by the start or end of the record) are flagged
seasons$complete <- peaks$i_start > 1 & peaks$i_end < nrow(sm)
cat("Detected crop cycles (start/end = 20% of NDVI amplitude before/after each peak):\n")
print(seasons, row.names = FALSE)
write.csv(seasons, "output/tables/detected_crop_cycles.csv", row.names = FALSE)

cy <- seasons[seasons$crop_year %in% c("2019/2020", "2020/2021", "2021/2022", "2022/2023", "2023/2024"), ]
cat("\nCrop cycles per agricultural year (Sep-Aug), complete years only:\n")
print(table(cy$crop_year))

cat("\nMedian timing by season type (all detected cycles):\n")
for (s in names(season_pal)) {
  x <- seasons[seasons$season == s, ]
  if (nrow(x) == 0) next
  cat(sprintf("  %-24s n = %d | peak around %s | median green-up to peak %d days | median peak NDVI %.2f\n",
              s, nrow(x), format(as.Date(paste0("2001-", format(x$peak, "%m-%d")))[order(x$peak)][ceiling(nrow(x)/2)], "%d %b"),
              as.integer(median(x$green_up_to_peak_days)), median(x$peak_NDVI)))
}

# =============================================================================
section("3. RAINFALL AND THE START OF THE WET-SEASON CROP")
# =============================================================================
d$rain_mm <- pmax(d$rain_mm, 0)
annual <- aggregate(rain_mm ~ year, d, sum)
cat("Annual rainfall (CHIRPS, regency mean, mm; approximate from 10-day sums):\n")
print(round(annual))

mon <- aggregate(rain_mm ~ month, d, function(x) sum(x) / length(unique(d$year)))
mon$rain_mm <- round(mon$rain_mm)
mon$month_name <- month.abb[mon$month]
write.csv(mon, "output/tables/mean_monthly_rainfall.csv", row.names = FALSE)
cat("\nMean monthly rainfall 2019-2024 (mm):\n"); print(mon[, c("month_name", "rain_mm")], row.names = FALSE)

# Rainy-season onset: first 10-day period from 1 September with >= 50 mm,
# where the following 20 days also receive >= 50 mm (avoids false starts)
onset <- data.frame()
for (y in 2019:2024) {
  x <- d[d$date >= as.Date(paste0(y, "-09-01")) & d$date < as.Date(paste0(y + 1, "-03-01")), ]
  if (nrow(x) < 3) next
  hit <- NA
  for (i in seq_len(nrow(x) - 2)) {
    if (x$rain_mm[i] >= 50 && sum(x$rain_mm[(i + 1):(i + 2)]) >= 50) { hit <- i; break }
  }
  if (!is.na(hit)) onset <- rbind(onset, data.frame(season_start_year = y, rain_onset = as.Date(x$date_start[hit])))
}
# Wet-season crop green-up that follows each onset
onset$wet_green_up <- as.Date(NA); onset$lag_days <- NA_integer_
wet <- seasons[seasons$season == "Wet season (MT-I)", ]
for (i in seq_len(nrow(onset))) {
  cand <- wet$start[wet$start >= onset$rain_onset[i] - 30 & wet$start <= onset$rain_onset[i] + 120]
  if (length(cand)) { onset$wet_green_up[i] <- min(cand); onset$lag_days[i] <- as.integer(min(cand) - onset$rain_onset[i]) }
}
cat("\nRainy-season onset and green-up of the wet-season crop:\n"); print(onset, row.names = FALSE)
write.csv(onset, "output/tables/rain_onset_vs_green_up.csv", row.names = FALSE)

cat("\nLowest smoothed NDVI in each dry season (Aug-Dec):\n")
dry_min <- do.call(rbind, lapply(2019:2024, function(y) {
  x <- sm[sm$date >= as.Date(paste0(y, "-08-01")) & sm$date <= as.Date(paste0(y, "-12-31")), ]
  data.frame(year = y, min_NDVI = round(min(x$NDVI), 3), date = x$date[which.min(x$NDVI)])
}))
print(dry_min, row.names = FALSE)
write.csv(dry_min, "output/tables/dry_season_minimum_ndvi.csv", row.names = FALSE)

# =============================================================================
section("4. FIGURES")
# =============================================================================
rain_col <- "#8a8984"
seasons_c <- seasons[!is.na(seasons$start), ]

# --- Figure 1: NDVI time series with crop cycles, and rainfall ---------------
d$quality <- ifelse(d$good, "Reliable (>= 50% cloud-free)", "Unreliable (cloudy)")
xlim_all <- c(as.Date("2019-01-01"), as.Date("2025-01-01"))
p_top <- ggplot() +
  geom_rect(data = seasons_c, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = season), alpha = 0.16) +
  geom_point(data = d[!d$good, ], aes(date, NDVI), shape = 21, colour = muted, fill = "white", size = 1.6) +
  geom_point(data = d[d$good, ], aes(date, NDVI), colour = ink2, size = 1.3, alpha = 0.75) +
  geom_line(data = sm, aes(date, NDVI), colour = col_ndvi, linewidth = 0.9) +
  geom_point(data = seasons_c, aes(peak, peak_NDVI, colour = season), size = 2.8, show.legend = FALSE) +
  scale_fill_manual(values = season_pal, breaks = names(season_pal), name = NULL) +
  scale_colour_manual(values = season_pal, guide = "none") +
  scale_x_date(limits = xlim_all, date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
  labs(y = "NDVI (cropland mean)", x = NULL,
       subtitle = "Grey dots: 10-day Sentinel-2 composites (hollow = cloudy, excluded). Green line: weighted smoothing spline. Shading: detected crop cycles.") +
  theme_pub() + theme(legend.position = "top")
p_bot <- ggplot(d, aes(as.Date(date_start) + 5, rain_mm)) +
  geom_col(fill = rain_col, width = 8) +
  scale_x_date(limits = xlim_all, date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
  labs(y = "Rainfall\n(mm / 10 days)", x = NULL) + theme_pub()
fig1 <- (p_top / p_bot) + plot_layout(heights = c(3, 1.2)) +
  plot_annotation(
    title = "Bojonegoro's cropland greens up twice a year, following the rains",
    caption = "Data: Copernicus Sentinel-2 (Cloud Score+ masked), CHIRPS rainfall, ESA WorldCover 2021 cropland, geoBoundaries. Self-initiated portfolio project using public data.",
    theme = theme(plot.title = element_text(face = "bold", size = 15, family = "DejaVu Sans", colour = ink),
                  plot.caption = element_text(colour = ink2, size = 8.5, hjust = 0, family = "DejaVu Sans")))
ggsave("output/figures/fig1_ndvi_rainfall_timeseries.png", fig1, width = 11, height = 6.4, dpi = 300, bg = "white")

# --- Figure 2: crop calendar by agricultural year ----------------------------
ag_start <- function(cyear) as.Date(paste0(substr(cyear, 1, 4), "-09-01"))
cal <- seasons_c
cal$x0 <- as.numeric(cal$start - ag_start(cal$crop_year))
cal$x1 <- as.numeric(cal$end   - ag_start(cal$crop_year))
cal$xp <- as.numeric(cal$peak  - ag_start(cal$crop_year))
on <- onset
on$crop_year <- paste0(on$season_start_year, "/", on$season_start_year + 1)
on$x <- as.numeric(on$rain_onset - ag_start(on$crop_year))
years_order <- sort(unique(cal$crop_year))
on <- on[on$crop_year %in% years_order, ]
cal$crop_year <- factor(cal$crop_year, levels = rev(years_order))
on$crop_year  <- factor(on$crop_year,  levels = rev(years_order))
month_breaks <- as.numeric(seq(as.Date("2001-09-01"), as.Date("2002-11-01"), by = "month") - as.Date("2001-09-01"))
fig2 <- ggplot(cal) +
  geom_segment(aes(x = x0, xend = x1, y = crop_year, yend = crop_year, colour = season),
               linewidth = 6, lineend = "round") +
  geom_point(aes(x = xp, y = crop_year), shape = 21, fill = "white", colour = ink, size = 2.4, stroke = 0.8) +
  geom_point(data = on, aes(x = x, y = crop_year), shape = 25, fill = ink, colour = ink, size = 2.6,
             position = position_nudge(y = 0.33)) +
  scale_colour_manual(values = season_pal, breaks = names(season_pal), name = NULL) +
  scale_x_continuous(breaks = month_breaks, labels = c(month.abb[c(9:12, 1:8)], "Sep", "Oct", "Nov"),
                     limits = c(0, 440), expand = c(0.01, 0)) +
  geom_vline(xintercept = 365, colour = ink2, linewidth = 0.4, linetype = "dashed") +
  labs(title = "The wet-season crop starts between late November and early December every year",
       subtitle = "Bars: start to end of each detected crop cycle (20% of NDVI amplitude). White dot: NDVI peak. Black triangle: rainy-season onset.",
       x = "Agricultural year (September to August)", y = NULL,
       caption = paste0("2018/2019 begins in January 2019 (start of the record). Dashed line: end of the agricultural year.\n",
                        "Rainy-season onset: first 10-day period from 1 September with >= 50 mm, followed by >= 50 mm in the next 20 days.")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 3))) +
  theme_pub() + theme(panel.grid.major.y = element_blank(), legend.text = element_text(margin = margin(r = 12)))
ggsave("output/figures/fig2_crop_calendar.png", fig2, width = 10.5, height = 4.8, dpi = 300, bg = "white")

# --- Figure 3: monthly NDVI maps, 2023 ---------------------------------------
r2 <- aggregate(r, fact = 2, fun = "mean", na.rm = TRUE)          # 120 m for display
crop_any <- sum(!is.na(r2)) > 0
bnd_utm <- st_transform(bnd, 32749)
mdf <- do.call(rbind, lapply(1:12, function(i) {
  x <- as.data.frame(c(r2[[i]], crop_any), xy = TRUE, na.rm = FALSE)
  names(x)[3:4] <- c("NDVI", "crop")
  x <- x[!is.na(x$crop) & x$crop == 1, ]
  x$month <- factor(month.name[i], levels = month.name)
  x
}))
mdf$NDVI <- pmin(pmax(mdf$NDVI, 0), 0.9)
cloud_share <- tapply(is.na(mdf$NDVI), mdf$month, mean)
cat("\nShare of cropland without cloud-free data in the monthly 2023 maps:\n"); print(round(cloud_share, 3))
fig3 <- ggplot() +
  geom_raster(data = mdf[is.na(mdf$NDVI), ], aes(x, y), fill = "#d9d8d3") +
  geom_raster(data = mdf[!is.na(mdf$NDVI), ], aes(x, y, fill = NDVI)) +
  geom_sf(data = bnd_utm, fill = NA, colour = ink, linewidth = 0.25) +
  facet_wrap(~ month, ncol = 4) +
  scale_fill_gradientn(colours = c("#f3efe0", "#d9e3b0", "#9cc37a", "#4f9a4f", "#1f6b33", "#0b3d1d"),
                       limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3), name = "NDVI") +
  coord_sf(datum = NA, expand = FALSE) +
  labs(title = "Monthly NDVI of Bojonegoro's cropland, 2023",
       subtitle = "Sentinel-2 monthly median composites at 60 m (shown at 120 m). Light grey: cropland with no cloud-free observation that month.",
       caption = "Crops are green in January and April-May and brown from September to November 2023; fields in the northern river valley stay green longer.") +
  theme_void(base_family = "DejaVu Sans") +
  theme(plot.title = element_text(face = "bold", size = 15, colour = ink),
        plot.subtitle = element_text(colour = ink2, size = 10), plot.caption = element_text(colour = ink2, size = 9, hjust = 0),
        strip.text = element_text(colour = ink, size = 10, margin = margin(b = 3)),
        legend.position = "right", plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 10, 10, 10), plot.title.position = "plot", plot.caption.position = "plot")
ggsave("output/figures/fig3_monthly_ndvi_maps_2023.png", fig3, width = 10, height = 6.0, dpi = 250, bg = "white")

# --- Figure 4: dry years and dry-season greenness ----------------------------
yr <- merge(annual, dry_min, by = "year")
yr$label <- paste0(yr$year)
fig4 <- ggplot(yr, aes(rain_mm, min_NDVI)) +
  geom_point(size = 3.2, colour = col_ndvi) +
  geom_text(aes(label = label), nudge_y = 0.012, size = 3.6, colour = ink) +
  scale_x_continuous(limits = c(1500, 2900)) + scale_y_continuous(limits = c(0.22, 0.45)) +
  labs(title = "In drier years, cropland turned browner in the dry season",
       subtitle = "Lowest NDVI from August to December vs annual rainfall (6 years; Spearman rho = 0.94)",
       x = "Annual rainfall (mm, CHIRPS regency mean)", y = "Dry-season minimum NDVI",
       caption = "2019 and 2023 were dry years (2019: strong positive Indian Ocean Dipole; 2023: El Nino).") +
  theme_pub()
ggsave("output/figures/fig4_dry_years_ndvi.png", fig4, width = 8.2, height = 5, dpi = 300, bg = "white")
sp <- cor.test(yr$rain_mm, yr$min_NDVI, method = "spearman", exact = TRUE)
cat(sprintf("\nSpearman correlation, annual rainfall vs dry-season minimum NDVI (n = 6): rho = %.2f, p = %.3f\n",
            sp$estimate, sp$p.value))

# --- Cover image for Upwork (16:9) -------------------------------------------
cover <- (p_top + labs(subtitle = NULL) + theme(legend.position = "top")) / p_bot +
  plot_layout(heights = c(3, 1.1)) +
  plot_annotation(
    title = "Cropping Cycles from Space: Sentinel-2 NDVI and Rainfall",
    subtitle = "Google Earth Engine + R | Bojonegoro cropland, East Java, 2019-2024 | crop-cycle detection and rainfall onset",
    caption = "Self-initiated portfolio project using public data (Copernicus Sentinel-2, CHIRPS, ESA WorldCover)",
    theme = theme(plot.title = element_text(face = "bold", size = 17, family = "DejaVu Sans", colour = ink),
                  plot.subtitle = element_text(size = 12, family = "DejaVu Sans", colour = ink2),
                  plot.caption = element_text(colour = ink2, size = 9, hjust = 0, family = "DejaVu Sans"),
                  plot.margin = margin(16, 22, 10, 22)))
ggsave("output/figures/cover_upwork.png", cover, width = 10, height = 5.625, dpi = 200, bg = "white")

write.csv(data.frame(date = sm$date, NDVI_smoothed = round(sm$NDVI, 4), NDWI_smoothed = round(sm$NDWI, 4))[seq(1, nrow(sm), 5), ],
          "output/tables/smoothed_series_5day.csv", row.names = FALSE)

section("SESSION INFO")
print(sessionInfo())
sink()
message("Done.")
