
# Subset to parvimonas in faeces
ps <- subset_samples(rar_physeq, origin=="faeces")
ps <- subset_taxa(ps, Genus %in% c("g__Parvimonas", "g__Parvimonas_NA"))
tax_table(ps)

data_phylo_long <- ps %>%
  tax_glom(taxrank = "Genus", bad_empty = c()) %>%   # Aglomerate to taxnomic range
  psmelt()    

# Proporciones
genus.crc.abund <- subset(data_phylo_long, sample.type=="crc")$Abundance
100*sum(genus.crc.abund>=2) / length(genus.crc.abund)

genus.noncrc.abund <- subset(data_phylo_long, sample.type=="non-crc")$Abundance
100*sum(genus.noncrc.abund>=2) / length(genus.noncrc.abund)

# data_phylo_long <- subset(data_phylo_long, Abundance>=1)

# Edad frente a cantidad
library(broom)
temp <- data_phylo_long
temp$sample.type <- gsub("non-crc", "Non CRC", temp$sample.type)
temp$sample.type <- gsub("crc", "CRC", temp$sample.type)
levels(temp$sex)[levels(temp$sex)=='female'] <- 'Female'
levels(temp$sex)[levels(temp$sex)=='male'] <- 'Male'
temp$Sex <- temp$sex



model <- lm(Abundance~sample.type+age+Sex, data=temp)
summary(model)

model.diag.metrics <- broom::augment(model)
p1 <- ggplot(model.diag.metrics, aes(age, Abundance, color=Sex, fill=Sex)) +
  facet_wrap(~sample.type) + 
  geom_point() +
  stat_smooth(method = NULL, se = TRUE, alpha = 0.2, linetype=2) +
  # geom_segment(aes(xend = age, yend = .fitted), color = "red", size = 0.3) +
  ggtitle("Abundance of Parvimonas depending of age, sex and patient group", 
          glue("No significant correlation found between the abundance of Parvimonas and sex or age.")
  ) + xlab("Age") +
  theme_bw() + theme(plot.title = element_text(size = 15, color = '#333333', face="bold"))

p2 <- ggplot(temp, aes(x = age, color = Sex, fill = Sex)) +
  geom_density(position = "identity", alpha = 0.3) +
  facet_wrap(~sample.type) + 
  ggtitle("Density distribution of samples", 
          glue("Density distribution of samples depending of age, sex and patient group")
  ) + xlab("Age") + ylab("Density") +
  theme_bw() + theme(plot.title = element_text(size = 15, color = '#333333', face="bold"))

p <- p1 + theme(legend.position="none") + p2  + plot_layout(guides = "collect", nrow=2)
ggsave(glue("test.png"), p, dpi="retina", unit="px", height=3500, width=3500)



# ANOVA
temp2 <- subset(temp, sample.type=="CRC")
temp2$group <- paste(temp2$sample.type, temp2$sex)
ggplot(temp2, aes(x=group, y=Abundance, color=group)) + geom_boxplot()
res.aov <- aov(Abundance ~ sex + age_group, data = temp2)
summary(res.aov)
shapiro.test(x = residuals(object = res.aov))

kruskal.test(Abundance ~ group, data = temp2)
