## Plot case data and diversity through time for state-level SARS

# Check and load packages, install missing. 
packages = c("ggplot2", "tidyverse", "Cairo", "PopGenome", "zoo")

lapply(packages, FUN = function(i) {
  if (!require(i, character.only = T)) {
    install.packages(i, dependencies = T, repos = "http://cran.us.r-project.org")
    library(i, character.only = T)
  }
})

#Read in modified metadata, raw case data, period case counts, and set name for output file. 
args <- commandArgs(trailingOnly = T)
all_metadata <- read_csv(args[1])
cases <- read_csv(args[2])
period_case_counts <- read_csv(args[3])
prefix = args[4]


## Plot variant turnover from metadata
variants_dates <- all_metadata %>% select(date, GISAID_clade, who_clades)
variants_dates$date <- as.Date(variants_dates$date)
begin <- unique(sort(variants_dates$date))[2]
end <- tail(unique(sort(variants_dates$date)), 1)
all_periods <- all_metadata %>% select(time_period) %>% drop_na()
all_periods <- sort(unique(all_periods$time_period))
all_periods <- data_frame(all_periods)
colnames(all_periods) <- "period"

variant_colors <- setNames(c("navyblue", "green4", "yellow", 
  "red", "purple", "orange"), unique(variants_dates$who_clades)) 
subset_colors <- setNames(c("red", "purple", "orange"), c("Delta", "Other", "Omicron"))

whole_counts_variants <- ggplot(data = variants_dates, aes(x = date, fill = who_clades)) + 
  geom_area(stat = "bin", binwidth = 7, alpha = .5) + 
  ggtitle("Covid variant replacement through time in Montana") +
  theme(plot.title = element_text(hjust = .5)) + 
  labs(fill = "CoV clades") +
  xlab("Date of sample") + 
  ylab("Total counts") + 
  scale_x_date(date_breaks = "1 week", date_labels = "%Y-%m",
               limits = c(as.Date(begin, "%Y-%m-%d"), as.Date(end, "%Y-%m-%d")), 
               expand = c(0,0)) +
  theme_void() + 
  theme(axis.text.x = element_text(angle = 90, size = 7), plot.title = element_text(hjust = .5))

frequencies_all <- ggplot(data = variants_dates, aes(x = date)) + 
  geom_histogram(aes(fill = who_clades), position = "fill", stat = "bin", binwidth = 7, alpha = .5) + 
  ggtitle("Covid variant replacement through time in Montana") +
  theme(plot.title = element_text(hjust = .5)) + 
  labs(fill = "CoV clades") +
  xlab("Date of sample") + 
  ylab("Frequency of variant") + 
  scale_x_date(date_breaks = "1 week", date_labels = "%Y-%m",
               limits = c(as.Date(begin, "%Y-%m-%d"), as.Date(end, "%Y-%m-%d")), 
               expand = c(0,0)) +
  scale_y_continuous(breaks = c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1)) +
  theme(axis.text.x = element_text(angle = 90, size = 7), axis.text.y = element_text(),
        plot.title = element_text(hjust = .5), panel.grid = element_line(color = "white", size =.5)) +
  scale_fill_manual(values = variant_colors)

frequencies_omicron_weekly <- ggplot(data = variants_dates, aes(x = date)) + 
  geom_histogram(aes(fill = who_clades), position = "fill", stat = "bin", binwidth = 7, alpha = .5) + 
  ggtitle("Covid variant replacement through time in Montana") +
  theme(plot.title = element_text(hjust = .5)) + 
  labs(fill = "CoV clades") +
  xlab("Date of sample") + 
  ylab("Frequency of variant") + 
  scale_x_date(date_breaks = "1 week", date_labels = "%Y-%m",
               limits = c(as.Date("2021-12-01", "%Y-%m-%d"), as.Date(end, "%Y-%m-%d")), 
               expand = c(0,0)) +
  scale_y_continuous(breaks = c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1)) +
  theme(axis.text.x = element_text(angle = 90, size = 7), axis.text.y = element_text(),
        plot.title = element_text(hjust = .5), panel.grid = element_line(color = "white", size =.5)) + 
  scale_fill_manual(values = subset_colors)

frequencies_omicron_daily <- ggplot(data = variants_dates, aes(x = date)) + 
  geom_histogram(aes(fill = who_clades), position = "fill", stat = "bin", binwidth = 1, alpha = .5) + 
  ggtitle("Covid variant replacement through time in Montana") +
  theme(plot.title = element_text(hjust = .5)) + 
  labs(fill = "CoV clades") +
  xlab("Date of sample") + 
  ylab("Frequency of variant") + 
  scale_x_date(date_breaks = "1 day", date_labels = "%Y-%m",
               limits = c(as.Date("2021-12-01", "%Y-%m-%d"), as.Date(end, "%Y-%m-%d")), 
               expand = c(0,0)) +
  scale_y_continuous(breaks = c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1)) +
  #theme_void() + 
  theme(axis.text.x = element_text(angle = 90, size = 7), axis.text.y = element_text(),
        plot.title = element_text(hjust = .5), panel.grid = element_line(color = "white", size =.5)) + 
  scale_fill_manual(values = subset_colors)


frequencies_all
ggsave(filename = "variant_frequency_all_weekly.pdf", width = 10, height = 5)
frequencies_omicron_weekly
ggsave(filename = "variant_frequency_omicron_weekly.pdf", width = 10, height = 5)
frequencies_omicron_daily
ggsave(filename = "variant_frequency_omicron_daily.pdf", width = 10, height = 5)


## Nucleotide diversity by week and total case counts
## Calc diversity through time from aligned fastas for each period. 

#Plot things. 

first_date = str_split(period_case_counts$period[1], "_")[[1]][1]
end_date = str_split(tail(period_case_counts$period,1), "_")[[1]][1]

all_fastas <- readData("results/aligned_fastas")


all_fastas <- F_ST.stats(all_fastas)
all_fastas <- neutrality.stats(all_fastas)
all_fastas <- linkage.stats(all_fastas)
all_fastas <- sweeps.stats(all_fastas)
all_fastas <- detail.stats(all_fastas)

spike_regions <- c(21563:25384)
regions_list <- list()

for (i in 1:length(all_fastas@region.names)) {
  regions_list[i] <- list(spike_regions)
}

spike_only <- splitting.data(all_fastas, positions = regions_list, whole.data = F, type = 2)

spike_only <- F_ST.stats(spike_only)
spike_only <- neutrality.stats(spike_only)
spike_only <- linkage.stats(spike_only)
spike_only <- sweeps.stats(spike_only)
spike_only <- detail.stats(spike_only)


nuc_diversity <- all_fastas@Pi
rows <- rownames(nuc_diversity)
rows <- str_remove_all(rows, ".fasta")
nuc_diversity <- rownames_to_column(data_frame(nuc_diversity))
nuc_diversity$rowname <- rows
colnames(nuc_diversity) <- c("period", "Pi")

spike_diversity <- spike_only@Pi
spike_diversity <- rownames_to_column(data_frame(spike_diversity))
spike_diversity$rowname <- rows
colnames(spike_diversity) <- c("period", "SpikePi")

nuc_diversity <- nuc_diversity %>% arrange(desc(period)) %>% 
  mutate(rolling_pi=rollmean(Pi, k = 2, fill = NA, align = "left")) %>% arrange(period)
nuc_diversity_all <- all_periods %>% full_join(nuc_diversity)
nuc_diversity_all <- nuc_diversity_all %>% full_join(period_case_counts)

tajimas <- all_fastas@Tajima.D
tajimas <- rownames_to_column(data_frame(tajimas))
tajimas$rowname <- rows
colnames(tajimas) <- c("period", "TajimaD")

nuc_diversity_all <- nuc_diversity_all %>% full_join(tajimas)

nuc_diversity_all <- nuc_diversity_all %>% mutate(count_divide = count/1000)

nuc_diversity_all <- nuc_diversity_all %>% full_join(spike_diversity)

nuc_diversity_all <- nuc_diversity_all %>% mutate(Pi_norm = Pi/max(Pi, na.rm = T))
nuc_diversity_all <- nuc_diversity_all %>% mutate(Spike_norm = SpikePi/max(SpikePi, na.rm = T))
nuc_diversity_all <- nuc_diversity_all %>% mutate(count_norm = count/max(count, na.rm = T))

diversity_plot_moving <- ggplot(data = nuc_diversity_all, aes(x = period, y = rolling_pi, group = 1)) +
  geom_smooth(stat = "smooth", position = "identity", show.legend = T, color = "blue", span = 0.2) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
              axis.text.x = element_blank(), plot.title = element_text(hjust=0.5), 
              panel.background = element_rect(fill = "white", colour = "grey50")) +
  ylab("Rolling average diversity") + 
  ggtitle("Weekly nucleotide diversity and positive COVID cases in Montana")

diversity_plot <- ggplot(data = nuc_diversity_all, aes(x = period, y = Pi, group = 1)) +
  geom_smooth(stat = "smooth", position = "identity", show.legend = T, color = "blue", span = 0.2) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), plot.title = element_text(hjust=0.5), 
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  ylab("Nucleotide diversity") + 
  ggtitle("Weekly nucleotide diversity and positive COVID cases in Montana")

tajimas_plot <- ggplot(data = nuc_diversity_all, aes(x = period, y = TajimaD, group = 1))+
  geom_smooth(stat = "smooth", position = "identity", show.legend = T, color = "blue", span = 0.2) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), plot.title = element_text(hjust=0.5), 
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  ylab("Tajima's D") + 
  ggtitle("Weekly Tajima's D and positive COVID cases in Montana")


weekly_case_plot <- ggplot(data = nuc_diversity_all, aes(x = period, y = count, group = 1)) + 
  geom_area(stat = "identity", fill = "blue", alpha = 0.4) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) + 
  ylab("Case count") + xlab("Week") + 
  scale_x_discrete(breaks = c("2020-03-27_2020-04-03", "2020-11-06_2020-11-13", 
                  "2021-09-24_2021-10-01", "2022-01-14_2022-01-21","2022-05-20_2022-05-27"), 
                  labels = c("March, 2020", "November, 2020", "September, 2021", "January, 2022", "May, 2022"), 
                  expand = c(0,0))


full_plot <- ggplot(data = nuc_diversity_all, aes(x = period, group = 1)) +
  geom_smooth(aes(y = Pi, colour = "Genome-wide"), stat = "smooth", position = "identity", show.legend = T, span = 0.2, alpha = 0.2) +
  geom_smooth(aes(y = SpikePi, colour = "Spike protein"), stat = "smooth", position = "identity", show.legend = T, span = 0.2, alpha = 0.2,
              linetype = "dashed") +
  scale_color_manual(name = "", values = c("blue", "red", "green")) +
  geom_area(aes(y = count_divide, fill = count_divide), stat = "identity", fill = "deepskyblue1", alpha = 0.3) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  scale_x_discrete(breaks = c("2020-01-01_2020-01-08", "2020-11-18_2020-11-25", "2021-11-17_2021-11-24", "2022-01-19_2022-01-26",
                              "2022-05-25_2022-06-01"), labels = c("January, 2020", "November, 2020", "November, 2021", "January, 2022", "May, 2022"), 
                   expand = c(0,0)) + 
  labs(y = expression("Nucleotide diversity ("~pi~")")) +
  scale_y_continuous(sec.axis = sec_axis(trans=~.*1000, name = "Weekly cases"))

nuc_plot_name <- paste0("results/nuc_diversity", prefix, "_nucleotide_diversity.pdf")
ggsave(nuc_plot_name, height = 8, width = 11, units = "in")
