
wilcox.test(duplets$X4.521625163826997751e.02, mu = 0)
shapiro.test(duplets$X4.521625163826997751e.02)
?Sys.glob
?lapply
duplets <- lapply(Sys.glob("160938??_duplets_score.csv"), read.csv, header=FALSE)
View(duplets)
length(duplets)

#Import the observed duplets scores and counts
duplets_s1 <- read.csv("16093801_duplets_score.csv", header=FALSE)
duplets_s3 <- read.csv("16093803_duplets_score.csv", header=FALSE)
duplets_s4 <- read.csv("16093804_duplets_score.csv", header=FALSE)
duplets_s5 <- read.csv("16093805_duplets_score.csv", header=FALSE)
duplets_s6 <- read.csv("16093806_duplets_score.csv", header=FALSE)
duplets_s7 <- read.csv("16093807_duplets_score.csv", header=FALSE)
duplets_s8 <- read.csv("16093808_duplets_score.csv", header=FALSE)
duplets_s9 <- read.csv("16093809_duplets_score.csv", header=FALSE)
duplets_s10 <- read.csv("16093810_duplets_score.csv", header=FALSE)
duplets_s11 <- read.csv("16093811_duplets_score.csv", header=FALSE)
duplets_s12 <- read.csv("16093812_duplets_score.csv", header=FALSE)


?hist()

dev.cur()
dev.off(1)

h1 <- hist(duplets_s1$V1, main = "Duplets Score S1", xlab="Doublet score",  ylab = "Cell count")
h3 <- hist(duplets_s3$V1, main = "Duplets Score S3", xlab="Doublet score",  ylab = "Cell count")
h4 <- hist(duplets_s4$V1, main = "Duplets Score S4", xlab="Doublet score",  ylab = "Cell count")
h5 <- hist(duplets_s5$V1, main = "Duplets Score S5", xlab="Doublet score",  ylab = "Cell count")
h6 <- hist(duplets_s6$V1, main = "Duplets Score S6", xlab="Doublet score",  ylab = "Cell count")
h7 <- hist(duplets_s7$V1, main = "Duplets Score S7", xlab="Doublet score",  ylab = "Cell count")
h8 <- hist(duplets_s8$V1, main = "Duplets Score S8", xlab="Doublet score",  ylab = "Cell count")
h9 <- hist(duplets_s9$V1, main = "Duplets Score S9", xlab="Doublet score",  ylab = "Cell count")
h10 <- hist(duplets_s10$V1, main = "Duplets Score S10", xlab="Doublet score",  ylab = "Cell count")
h11 <- hist(duplets_s11$V1, main = "Duplets Score S11", xlab="Doublet score",  ylab = "Cell count")
h12 <- hist(duplets_s12$V1, main = "Duplets Score S12", xlab="Doublet score",  ylab = "Cell count")
plot(rnorm(50), rnorm(50))

#Import the simulated duplets scores

sim_duplets <- lapply(Sys.glob("160938??_sim_duplets_score.csv"), read.csv, header=FALSE,)
length(sim_duplets)
sim_duplets
