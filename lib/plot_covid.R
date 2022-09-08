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

all_fastas <- readData("results/aligned_fastas",include.unknown=T)


all_fastas <- F_ST.stats(all_fastas,mode="nucleotide")
#all_fastas <- neutrality.stats(all_fastas)
#all_fastas <- linkage.stats(all_fastas)
#all_fastas <- sweeps.stats(all_fastas)
#all_fastas <- detail.stats(all_fastas)

gene_names = c("leader", "nsp2", "nsp3", "nsp4", "proteinase", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "RNA_RNA_polymerase", "helicase", "exonuclease", "endoRNAase", "ribose_methyltransferase", "spike", "ORF3a", "envelopes", "membrane", "ORF6", "ORF7a", "ORF7b", "ORF8", "nucleocapsid", "ORF10")

genes = list(266:805, 806:2719, 2720:8554, 8555:10054, 10055:10972, 10973:11842, 11843:12091, 12092:12685, 12686:13024, 13025:13441, c(13442:13468,13468:16236), 16237:18039, 18040:19620, 19621:20658, 20659:21552, 21563:25384, 25393:26220, 26245:26472, 26523:27191, 27202:27387, 27394:27759, 27756:27887, 27894:28259, 28274:29553, 29558:29674)

gene_output = vector("list",length(genes))
diversity = vector("list",length(genes))

for (i in 1:length(genes)) {
  regions_temp = list()
  foo = as.integer(unlist(genes[i]))
  for (j in 1:length(all_fastas@region.names)) {
    regions_temp[j] = list(foo)
  }
  out = splitting.data(all_fastas, positions = regions_temp, whole.data = F, type = 2)
  out = F_ST.stats(out,mode="nucleotide")
  gene_output[i] = out
}


nuc_diversity <- all_fastas@nuc.diversity.within
rows <- rownames(nuc_diversity)
rows <- str_remove_all(rows, ".fasta")
nuc_diversity <- rownames_to_column(data_frame(nuc_diversity))
nuc_diversity$rowname <- rows
colnames(nuc_diversity) <- c("period", "Pi")

for (i in 1:length(genes)) {
  foo = gene_output[[i]]@nuc.diversity.within
  foo = rownames_to_column(data_frame(foo))
  foo$rowname = rows
  colnames(foo) = c("period", paste0(gene_names[[i]],"Pi"))
  diversity[i] = list(foo)
}

nuc_diversity <- nuc_diversity %>% arrange(desc(period)) %>% 
  mutate(rolling_pi=rollmean(Pi, k = 2, fill = NA, align = "left")) %>% arrange(period)
nuc_diversity_all <- all_periods %>% full_join(nuc_diversity)
nuc_diversity_all <- nuc_diversity_all %>% full_join(period_case_counts)

#tajimas <- all_fastas@Tajima.D
#tajimas <- rownames_to_column(data_frame(tajimas))
#tajimas$rowname <- rows
#colnames(tajimas) <- c("period", "TajimaD")

#nuc_diversity_all <- nuc_diversity_all %>% full_join(tajimas)

nuc_diversity_all <- nuc_diversity_all %>% mutate(count_divide = count/2000)

for (i in 1:length(genes)) {
  nuc_diversity_all = nuc_diversity_all %>% full_join(diversity[[i]])
}

nuc_diversity_all <- nuc_diversity_all %>% mutate(Pi_norm = Pi/max(Pi, na.rm = T))
for (i in 1:length(gene_names)) {
  foo = paste0(gene_names[i],"Pi")
  colname = paste0(gene_names[i],"_norm")
  bar = which(colnames(nuc_diversity_all) == foo)
  nuc_diversity_all <- nuc_diversity_all %>% mutate(colname = nuc_diversity_all[[bar]]/max(nuc_diversity_all[[bar]], na.rm = T))
}
nuc_diversity_all <- nuc_diversity_all %>% mutate(count_norm = count/max(count, na.rm = T))

#diversity_plot_moving <- ggplot(data = nuc_diversity_all, aes(x = period, y = rolling_pi, group = 1)) +
#  geom_smooth(stat = "smooth", position = "identity", show.legend = T, color = "blue", span = 0.2) +
#  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
#              axis.text.x = element_blank(), plot.title = element_text(hjust=0.5), 
#              panel.background = element_rect(fill = "white", colour = "grey50")) +
#  ylab("Rolling average diversity") + 
#  ggtitle("Weekly nucleotide diversity and positive COVID cases in Montana")

#diversity_plot <- ggplot(data = nuc_diversity_all, aes(x = period, y = Pi, group = 1)) +
#  geom_smooth(stat = "smooth", position = "identity", show.legend = T, color = "blue", span = 0.2) +
#  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
#        axis.text.x = element_blank(), plot.title = element_text(hjust=0.5), 
#        panel.background = element_rect(fill = "white", colour = "grey50")) +
#  ylab("Nucleotide diversity") + 
#  ggtitle("Weekly nucleotide diversity and positive COVID cases in Montana")

#tajimas_plot <- ggplot(data = nuc_diversity_all, aes(x = period, y = TajimaD, group = 1))+
#  geom_smooth(stat = "smooth", position = "identity", show.legend = T, color = "blue", span = 0.2) +
#  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
#        axis.text.x = element_blank(), plot.title = element_text(hjust=0.5), 
#        panel.background = element_rect(fill = "white", colour = "grey50")) +
#  ylab("Tajima's D") + 
#  ggtitle("Weekly Tajima's D and positive COVID cases in Montana")


#weekly_case_plot <- ggplot(data = nuc_diversity_all, aes(x = period, y = count, group = 1)) + 
#  geom_area(stat = "identity", fill = "blue", alpha = 0.4) +
#  theme(panel.background = element_rect(fill = "white", colour = "grey50")) + 
#  ylab("Case count") + xlab("Week") + 
#  scale_x_discrete(breaks = c("2020-03-27_2020-04-03", "2020-11-06_2020-11-13", 
#                  "2021-09-24_2021-10-01", "2022-01-14_2022-01-21","2022-05-20_2022-05-27"), 
#                  labels = c("March, 2020", "November, 2020", "September, 2021", "January, 2022", "May, 2022"), 
#                  expand = c(0,0))

#Hack because R is the worst.
`+.uneval` <- function(a,b) {
    `class<-`(modifyList(a,b), "uneval")
}


full_plot <- ggplot(data = nuc_diversity_all, aes(x = period, group = 1)) + geom_smooth(aes(y = Pi, colour = "Genome-wide"), stat = "smooth", position = "identity", show.legend = T, span = 0.2, alpha = 0.2)
for (i in 1:length(genes)) {
  full_plot = full_plot + geom_smooth(aes_string(y=paste0(gene_names[i],"Pi")) + aes(color = !! gene_names[i]), stat = "smooth", position = "identity", show.legend = T, span = 0.2, alpha = 0.2)
}
full_plot = full_plot +  geom_area(aes(y = count_divide, fill = count_divide), stat = "identity", fill = "deepskyblue1", alpha = 0.3) + theme(panel.background = element_rect(fill = "white", colour = "grey50")) + scale_x_discrete(breaks = c("2020-01-01_2020-01-08", "2020-02-26_2020-03-04", "2020-04-29_2020-05-06", "2020-07-01_2020-07-08", "2020-08-26_2020-09-02", "2020-11-04_2020-11-11", "2020-12-30_2021-01-06", "2021-02-24_2021-03-03", "2021-04-28_2021-05-05", "2021-06-30_2021-07-07", "2021-08-25_2021-09-01", "2021-10-27_2021-11-03", "2021-12-29_2022-01-05", "2022-02-23_2022-03-02", "2022-04-27_2022-05-04","2022-06-29_2022-07-06"), labels = c("1/20","3/20","5/20","7/20","9/20","11/20","1/21","3/21","5/21","7/21","9/21","11/21","1/22","3/22","5/22","7/22"), expand = c(0,0)) + labs(title = prefix, x = "Time", y = expression("Nucleotide diversity ("~pi~")")) + scale_y_continuous(limits=c(0,40), sec.axis = sec_axis(trans=~.*2000, name = "Weekly cases"))

nuc_plot_name <- paste0("results/", prefix, "_all_genes.pdf")
ggsave(nuc_plot_name, height = 8, width = 11, units = "in")

full_plot <- ggplot(data = nuc_diversity_all, aes(x = period, group = 1)) + geom_smooth(aes(y = Pi, colour = "Genome-wide"), stat = "smooth", position = "identity", show.legend = T, span = 0.2, alpha = 0.2)
for (i in 16) {
  full_plot = full_plot + geom_smooth(aes_string(y=paste0(gene_names[i],"Pi")) + aes(color = !! gene_names[i]), stat = "smooth", position = "identity", show.legend = T, span = 0.2, alpha = 0.2)
}
full_plot = full_plot +  geom_area(aes(y = count_divide, fill = count_divide), stat = "identity", fill = "deepskyblue1", alpha = 0.3) + theme(panel.background = element_rect(fill = "white", colour = "grey50")) + scale_x_discrete(breaks = c("2020-01-01_2020-01-08", "2020-02-26_2020-03-04", "2020-04-29_2020-05-06", "2020-07-01_2020-07-08", "2020-08-26_2020-09-02", "2020-11-04_2020-11-11", "2020-12-30_2021-01-06", "2021-02-24_2021-03-03", "2021-04-28_2021-05-05", "2021-06-30_2021-07-07", "2021-08-25_2021-09-01", "2021-10-27_2021-11-03", "2021-12-29_2022-01-05", "2022-02-23_2022-03-02", "2022-04-27_2022-05-04", "2022-06-29_2022-07-06"), labels = c("1/20","3/20","5/20","7/20","9/20","11/20","1/21","3/21","5/21","7/21","9/21","11/21","1/22","3/22","5/22","7/22"), expand = c(0,0)) + labs(title = prefix, x = "Time", y = expression("Nucleotide diversity ("~pi~")")) + scale_y_continuous(limits=c(0,40), sec.axis = sec_axis(trans=~.*2000, name = "Weekly cases"))

nuc_plot_name <- paste0("results/", prefix, "_spike_only.pdf")
ggsave(nuc_plot_name, height = 8, width = 11, units = "in")

write.table(nuc_diversity_all,paste0(prefix,"_nuc_diversity.csv"),sep=",")
