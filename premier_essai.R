#20/10/2026

# Ce script est le premier essai de prise en main de données satellitaires provenant de différentes
# sources pour reconstruire l'héritage environnemental des lagons de Nouvell-Calédonie, de la 
# Réunion et de Mayotte.

# Je vais tout d'abord télécharger des fichiers provenant de Marine Copernicus

# chargement des libraries

library(ncdf4)
library(tidync)
library(dplyr)
library(terra)
library(heatwaveR)
library(lubridate)
library(ggplot2)

# Marine Copernicus -------------------------------------------------------

## produit OSTIA -------------------------------------------------------

# J'ai télécharger via Copernicus Marine la variable "Sea surface foundation temperature"
# provenant du produit OSTIA, entre le 17/01/2024 et le 31/12/2025

df <- tidync("METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_1792504939433.nc") |>
  hyper_tibble()

# Convertir la date en "vraie" date
df$time <- as.Date(df$time)

# Convertir les coordonnées en nombres (l'arrondi évite les écarts de décimales)
df$longitude <- round(as.numeric(df$longitude), 3)
df$latitude  <- round(as.numeric(df$latitude), 3)

# la colonne analysed_sst est en Kelvin et pas en °C donc on doit la convertir
# on retir 273,15 à cette colonne

df$sst_celsius <- df$analysed_sst - 273.15

# Calculer la moyenne en °C pour chaque jour

moy_jour <- df |>
  group_by(time) |>
  summarise(sst_moy = mean(sst_celsius, na.rm = TRUE))

# Visualisation

ggplot(moy_jour, aes(x = time, y = sst_moy)) +
  geom_line(color = "steelblue") +
  labs(title = "SST moyenne journalière",
       x = NULL, y = "Température (°C)") +
  theme_minimal()

# Analyse spatiale

# Le but ici est de créer des cartes représentant la température moyenne mensuelle
# pour chaue pixel

# créer une colonne année
df$annee <- as.numeric(format(as.Date(df$time), "%Y"))

# créer une colonne mois
df$mois <- as.numeric(format(df$time, "%m"))

# créer une colonne jour
df$day <- as.numeric(format(df$time, "%d"))

# créer une colonne jour de l'année
df$doy <- yday(df$time)

# créer des cartes de températures de surface de l'eau par mois pour chaue
# pixel

sst_mois <- df |>
  group_by(longitude, latitude, mois) |>
  summarise(sst = mean(sst_celsius, na.rm = TRUE), .groups = "drop") |>
  mutate(mois = factor(mois, levels = 1:12,
                       labels = c("janv", "févr", "mars", "avr", "mai", "juin",
                                  "juil", "août", "sept", "oct", "nov", "déc")))
# Visualisation

ggplot(sst_mois, aes(x = longitude, y = latitude, fill = sst)) +
  geom_raster() +
  facet_wrap(~ mois, ncol = 4) +
  scale_fill_viridis_c(option = "inferno", name = "SST (°C)") +
  coord_quickmap() +
  labs(title = "Température de surface moyenne par mois",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

## produit GLORYS12 -------------------------------------------------------

df <- tidync("data/cmems_mod_glo_phy_my_0.083deg_P1D-m_1792520915883.nc") |>
  hyper_tibble()

# Convertir la date en "vraie" date
df$time <- as.Date(df$time)

# Convertir les coordonnées en nombres (l'arrondi évite les écarts de décimales)
df$longitude <- round(as.numeric(df$longitude), 3)
df$latitude  <- round(as.numeric(df$latitude), 3)

# Calculer la moyenne en °C pour chaque jour

moy_jour <- df |>
  group_by(time) |>
  summarise(sst_moy = mean(thetao, na.rm = TRUE))

# Visualisation

ggplot(moy_jour, aes(x = time, y = sst_moy)) +
  geom_line(color = "steelblue") +
  labs(title = "SST moyenne journalière",
       x = NULL, y = "Température (°C)") +
  theme_minimal()

# Analyse spatiale

# Le but ici est de créer des cartes représentant la température moyenne mensuelle
# pour chaue pixel

# créer une colonne année
df$annee <- as.numeric(format(as.Date(df$time), "%Y"))

# créer une colonne mois
df$mois <- as.numeric(format(df$time, "%m"))

# créer une colonne jour
df$day <- as.numeric(format(df$time, "%d"))

# créer une colonne jour de l'année
df$doy <- yday(df$time)

# créer des cartes de températures de surface de l'eau par mois pour chaue
# pixel

sst_mois <- df |>
  group_by(longitude, latitude, mois) |>
  summarise(sst = mean(thetao, na.rm = TRUE), .groups = "drop") |>
  mutate(mois = factor(mois, levels = 1:12,
                       labels = c("janv", "févr", "mars", "avr", "mai", "juin",
                                  "juil", "août", "sept", "oct", "nov", "déc")))
# Visualisation

ggplot(sst_mois, aes(x = longitude, y = latitude, fill = sst)) +
  geom_raster() +
  facet_wrap(~ mois, ncol = 4) +
  scale_fill_viridis_c(option = "inferno", name = "SST (°C)") +
  coord_quickmap() +
  labs(title = "Température de surface moyenne par mois",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

### Climatologie mensuelle -------------------------------------------------------

# Pour calculer une climatologie on doit définir une période de référence
# Normalement cette climatologie de référence se calcule sur 30 ans, ici 
# on a ue 25 ans de données, donc tant pis. 

clim_mens <- df |> 
  group_by(mois) |> 
  summarise(clim = mean(thetao, na.rm = TRUE),
            sd = sd(thetao, na.rm = TRUE))

# Visualisation

ggplot(clim_mens, aes(mois, clim)) +
  geom_ribbon(aes(ymin = clim - sd, ymax = clim + sd), alpha = 0.2) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = 1:12) +
  labs(title = "Climatologie mensuelle de la SST",
       x = "Mois", y = "SST (°C)") +
  theme_minimal()

### Anomalie mensuelle -------------------------------------------------------

anomalie <- df |> 
  group_by(annee, mois) |> 
  summarise(temp = mean(thetao), .groups = "drop") |> 
  left_join(clim_mens, by = "mois") |> 
  mutate(anomalie = temp - clim,
         date = as.Date(paste(annee, mois, 1, sep = "-")))

# Visualisation

ggplot(anomalie, aes(date, anomalie, fill = anomalie > 0)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue")) +
  labs(x = NULL, y = "Anomalie (°C)") +
  theme_minimal()

