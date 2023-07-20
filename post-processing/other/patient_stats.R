f = "/home/usuario/Proyectos/Results/patients.tsv"
df <- read.csv(f, sep="\t")

categorical.vars <- c("sample.type","sex","bristol_scale",
                      "has_oral_disease", "has_other_disease", 
                      "other_disease", "exercises", "had_lengthy_antibiotic_treatment", 
                      "consumes_alcohol", "consumes_tobacco", "consumes_caffeine", 
                      "sleep_disorder", "environment")
numerical.vars <- c("age","height","weight")

# Transform columns to factors
df[,categorical.vars] <- lapply(df[,categorical.vars], factor)

# Numerical variables t-test
temp <- subset(df, sex=="hombre")[c("sample.type", numerical.vars)]
res <- lapply(temp[,-1], function(x) t.test(x ~ temp$sample.type)); res
do.call(rbind, res)[,c(1,3)]

temp <- subset(df, sex=="mujer")[c("sample.type", numerical.vars)]
res <- lapply(temp[,-1], function(x) t.test(x ~ temp$sample.type)); res
do.call(rbind, res)[,c(1,3)]

# Numerical variables medians
df %>% dplyr::group_by(sample.type, sex) %>% dplyr::summarise(median_age=median(na.omit(age)), median_height=median(na.omit(height)), median_weight=median(na.omit(weight)))
df %>% dplyr::group_by(sex) %>% dplyr::summarise(median_age=median(na.omit(age)), median_height=median(na.omit(height)), median_weight=median(na.omit(weight)))


# Categorical variables chisq test
temp <- df[categorical.vars]
res <- lapply(temp[,-1], function(x) chisq.test(temp[,1], x, simulate.p.value = TRUE)); res
do.call(rbind, res)[,c(1,3)]






