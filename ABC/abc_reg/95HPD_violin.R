##计算95HPD

library(ggplot2)
library(dplyr)
library(gridExtra)
library(coda)


current_directory <- setwd("/fast3/group_crf/home/b20gengjlin3/xipho_simul/abc_s/APS/abc_reg")
matching_files <- list.files(current_directory, pattern = ".tangent.post.gz")
oF1 <- gzfile(matching_files[1])
oF1 <- gzfile(matching_files[2])
oF1 <- gzfile(matching_files[3])
oF1 <- gzfile(matching_files[4])
oF1 <- gzfile(matching_files[5])


dat1 <- read.table(oF1, header = F)
dat1$Source <- "APS"

colnames(dat1) <- c('gen', 'popsize1', 'popsize2', "Ne","s","pop")
dat1$l <- dat1$popsize1 + dat1$popsize2
dat1$lambda <- dat1$popsize2 / dat1$l

dat1_new <- dat1[,c("gen","lambda","Ne","s","pop")]

library(ggplot2)
# 直方图
ggplot(dat1_new, aes(x = s)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of s", x = "s", y = "Frequency") +
  theme_minimal() +
  xlim(-1, 1)  

# 密度图
ggplot(dat1_new, aes(x = s)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density Plot of s", x = "s", y = "Density") +
  theme_minimal() + 
  xlim(-1, 1)

###扫描峰值


# 计算密度
dens <- density(dat1_new$s)

# 找到峰值
peak <- dens$x[which.max(dens$y)]

# 绘制图形并标注峰值
ggplot(dat1_new, aes(x = s)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = peak, color = "red", linetype = "dashed") +
  annotate("text", x = peak, y = max(dens$y), 
           label = paste("Peak at", round(peak, 3)), 
           vjust = -1, color = "red") +
  labs(title = "Density Plot of s", x = "s", y = "Density") +
  theme_minimal() + 
  xlim(-1, 1)


###
#lambda
###
# 直方图
ggplot(dat1_new, aes(x = lambda)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of lambda", x = "lambda", y = "Frequency") +
  theme_minimal()

# 密度图
ggplot(dat1_new, aes(x = lambda)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density Plot of lambda", x = "lambda", y = "Density") +
  theme_minimal()

# 计算密度
dens <- density(dat1_new$lambda)

# 找到峰值
peak <- dens$x[which.max(dens$y)]

# 绘制图形并标注峰值
ggplot(dat1_new, aes(x = lambda)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = peak, color = "red", linetype = "dashed") +
  annotate("text", x = peak, y = max(dens$y), 
           label = paste("Peak at", round(peak, 3)), 
           vjust = -1, color = "red") +
  labs(title = "Density Plot of s", x = "s", y = "Density") +
  theme_minimal() + 
  xlim(-1, 1)
####
#gen
####
# 直方图
ggplot(dat1_new, aes(x = gen)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of T", x = "T", y = "Frequency") +
  theme_minimal()

# 密度图
ggplot(dat1_new, aes(x = gen)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density Plot of T", x = "T", y = "Density") +
  theme_minimal()

####
#Ne
####
# 直方图
ggplot(dat1_new, aes(x = Ne)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of Ne", x = "Ne", y = "Frequency") +
  theme_minimal()

# 密度图
ggplot(dat1_new, aes(x = Ne)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density Plot of Ne", x = "Ne", y = "Density") +
  theme_minimal()


















current_directory <- setwd("/fast3/group_crf/home/b20gengjlin3/xipho_simul/abc_s/APS/abc_reg/123")
matching_files <- list.files(current_directory, pattern = ".tangent.post.gz")
oF1 <- gzfile(matching_files[1])


dat1 <- read.table(oF1, header = F)
dat1$Source <- "APS"

colnames(dat1) <- c('gen', 'popsize1', 'popsize2', "Ne","s","pop")
dat1$l <- dat1$popsize1 + dat1$popsize2
dat1$lambda <- dat1$popsize2 / dat1$l

dat1_new <- dat1[,c("gen","lambda","Ne","s","pop")]

library(ggplot2)
# 直方图
ggplot(dat1_new, aes(x = s)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of s", x = "s", y = "Frequency") +
  theme_minimal()

# 密度图
ggplot(dat1_new, aes(x = s)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density Plot of s", x = "s", y = "Density") +
  theme_minimal()
# 直方图
ggplot(dat1_new, aes(x = lambda)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of lambda", x = "lambda", y = "Frequency") +
  theme_minimal()

# 密度图
ggplot(dat1_new, aes(x = lambda)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density Plot of lambda", x = "lambda", y = "Density") +
  theme_minimal()
####
####
# 直方图
ggplot(dat1_new, aes(x = gen)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of T", x = "T", y = "Frequency") +
  theme_minimal()

# 密度图
ggplot(dat1_new, aes(x = gen)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density Plot of T", x = "T", y = "Density") +
  theme_minimal()
####
####
# 直方图
ggplot(dat1_new, aes(x = Ne)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of Ne", x = "Ne", y = "Frequency") +
  theme_minimal()

# 密度图
ggplot(dat1_new, aes(x = Ne)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density Plot of Ne", x = "Ne", y = "Density") +
  theme_minimal()


