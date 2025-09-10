setwd("C:/Users/Nathan/Documents/Practicum/herbicide_soil_DOBBINS/")
# install.packages("wesanderson")

# Load required libraries
library(tidyverse)
library(readxl)
library(wesanderson)
library(ggplot2)

# File path
file <- "data/CO2_Resp_RAW_Data.xlsx"

# Sheets that contain respiration data
sheets <- c("Lab Exp (Day 0) respiration", "Lab Exp. (Day 1)", "Lab Exp (Day 3)", 
            "Lab Exp (Day 7)", "Lab Exp (Day 14)", "Lab Exp (Day 21)", "Lab Exp (Day 28)")

# Read and combine all sheets, track the day
all_data <- sheets %>%
  map_dfr(~ {
    df <- read_excel(file, sheet = .x)
    
    df <- df %>%
      mutate(
        `Start Time` = if ("Start Time" %in% names(df)) as.character(`Start Time`) else `Start Time`,
        `End Time` = if ("End Time" %in% names(df)) as.character(`End Time`) else `End Time`,
        `Time Difference` = if ("Time Difference" %in% names(df)) as.character(`Time Difference`) else `Time Difference`,
        Day = .x
      )
    
    return(df)
  })

### === 2-4 D Panel === ###
filtered <- all_data %>%
  filter(grepl("2-4 D|Control", Trial, ignore.case = TRUE)) %>%
  filter(!grepl("outlier", Trial, ignore.case = TRUE))

filtered <- filtered %>%
  mutate(Concentration = case_when(
    grepl("1/8", Trial) ~ "12.5%",
    grepl("1/4", Trial) ~ "25%",
    grepl("1/2", Trial) ~ "50%",
    grepl("Rec", Trial) ~ "100%",
    grepl("Control", Trial) ~ "Control",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Concentration)) %>%
  mutate(Sample_ID = str_trim(Trial))

cumulative <- filtered %>%
  group_by(Sample_ID, Concentration) %>%
  summarise(Cumulative_CO2 = sum(`Adjusted ppm`, na.rm = TRUE), .groups = "drop")

cumulative$Concentration <- factor(cumulative$Concentration, 
                                   levels = c("Control", "12.5%", "25%", "50%", "100%"))

ggplot(cumulative, aes(x = Concentration, y = Cumulative_CO2, fill = Concentration)) +
  geom_boxplot(color = "black", width = 0.6, outlier.shape = NA) +
  labs(
    title = expression("2-4 D"),
    x = "Concentration",
    y = expression("Cumulative CO"[2]*" ("*mu*"g/gdw/hr)")
  ) +
  scale_fill_manual(values = c("white", rep("lightblue", 4))) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
### === 2-4 D Bar Chart === ###
p_24d_bar <- ggplot(cumulative, aes(x = Concentration, y = Cumulative_CO2, fill = Concentration)) +
  geom_bar(stat = "summary", fun = mean, width = 0.6, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, color = "black") +
  labs(
    title = "2,4-D",
    x = "% of Recommended Dose",
    y = expression("Cumulative CO"[2]*" ("*mu*"g/gdw/hr)")
  ) +
  scale_fill_manual(values = c("Control" = "white", 
                               "12.5%" = "#E3F2FD", 
                               "25%" = "#BBDEFB", 
                               "50%" = "#64B5F6", 
                               "100%" = "#1976D2")) +
  theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size=19),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),   
    legend.position = "none"
  )

print(p_24d_bar)
ggsave("Figures/CO2_24D_barplot.png", plot = p_24d_bar, width = 8, height = 6, dpi = 300, bg = "white")

### === Glyphosate Panel === ###
glypho_filtered <- all_data %>%
  filter(grepl("Glyphosate|Control", Trial, ignore.case = TRUE)) %>%
  filter(!grepl("outlier", Trial, ignore.case = TRUE))

glypho_filtered <- glypho_filtered %>%
  mutate(Concentration = case_when(
    grepl("1/8", Trial) ~ "12.5%",
    grepl("1/4", Trial) ~ "25%",
    grepl("1/2", Trial) ~ "50%",
    grepl("Rec", Trial) ~ "100%",
    grepl("Control", Trial) ~ "Control",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Concentration)) %>%
  mutate(Sample_ID = str_trim(Trial))

cumulative_glypho <- glypho_filtered %>%
  group_by(Sample_ID, Concentration) %>%
  summarise(Cumulative_CO2 = sum(`Adjusted ppm`, na.rm = TRUE), .groups = "drop")

cumulative_glypho$Concentration <- factor(cumulative_glypho$Concentration, 
                                          levels = c("Control", "12.5%", "25%", "50%", "100%"))

ggplot(cumulative_glypho, aes(x = Concentration, y = Cumulative_CO2, fill = Concentration)) +
  geom_boxplot(color = "black", width = 0.6, outlier.shape = NA) +
  labs(
    title = expression("Glyphosate"),
    x = "Concentration",
    y = expression("Cumulative CO"[2]*" ("*mu*"g/gdw/hr)")
  ) +
  scale_fill_manual(values = c("white", rep("lightblue", 4))) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# 2-4 D
ggplot(subset(cumulative, Herbicide == "2-4 D"), aes(x = Concentration, y = Cumulative_CO2, fill = Concentration)) +
  geom_bar(stat = "summary", fun = mean) +
  labs(
    title = expression("2-4 D"),
    x = "Concentration",
    y = expression("Cumulative CO"[2]*" ("*mu*"g/gdw/hr)")
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Glyphosate
ggplot(subset(cumulative, Herbicide == "Glyphosate"), aes(x = Concentration, y = Cumulative_CO2, fill = Concentration)) +
  geom_bar(stat = "summary", fun = mean) +
  labs(
    title = "Glyphosate",
    x = "Concentration",
    y = expression("Cumulative CO"[2]*" ("*mu*"g/gdw/hr)")
  ) +
  theme_minimal() +
  theme(legend.position = "none")

### === Glyphosate Bar Chart === ###

p_glyph_bar <- ggplot(cumulative_glypho, aes(x = Concentration, y = Cumulative_CO2, fill = Concentration)) +
  geom_bar(stat = "summary", fun = mean, width = 0.6, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, color = "black") +
  labs(
    title = "Glyphosate",
    x = "% of Recommended Dose",
    y = expression("Cumulative CO"[2]*" ("*mu*"g/gdw/hr)")
  ) +
  scale_fill_manual(values = c("Control" = "white", 
                               "12.5%" = "#E3F2FD", 
                               "25%" = "#BBDEFB", 
                               "50%" = "#64B5F6", 
                               "100%" = "#1976D2")) +
  theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size=19),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),   
    legend.position = "none"
  )

print(p_glyph_bar)
ggsave("Figures/CO2_Glyphosate_barplot.png", plot = p_glyph_bar, width = 8, height = 6, dpi = 300, bg = "white")