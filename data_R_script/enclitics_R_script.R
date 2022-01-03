#============================================================================================
# Author: Robert Daugs

# Publication: English modal enclitic constructions: A diachronic, usage-based study of 'd
#              and 'll. Cognitive Linguistics

# R version 4.1.2 (2021-11-01)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)

# last updated: Dec 31, 2021
#============================================================================================

## DISCLAIMER

# The script becomes computationally quite expensive. R 64-bit strongly recommended and 16GB
# Ram memory.

#============================================================================================

## SETUP

# clear environment
rm(list = ls())

# set working directory
setwd("~/data_R_script")

# load/write required functions
source("vnc_sd.R") # VNC analysis, standard deviation as distance measure, modified version
                   # (cf. Gries & Hilpert 2012)

logw0 <- function(x) {
  ifelse(x == 0, 0, log(x))
} # function to deal with log(0)s

# load packages
library(collostructions) # package not (yet) available on CRAN, see https://sfla.ch
library(data.table)
library(plyr)
library(dplyr)
library(tidyverse)
library(zipfR)

#============================================================================================

## LOAD AND PREPARE DATA

# read files
input.d.FIC <- read.delim("d_COHA_FIC.txt",
                          header = TRUE)
input.would.FIC <- read.delim("would_COHA_FIC.txt",
                              header = TRUE)
input.ll.FIC <- read.delim("ll_COHA_FIC.txt",
                           header = TRUE)
input.will.FIC <- read.delim("will_COHA_FIC.txt",
                             header = TRUE)

# order data chronologically by year
input.d.FIC <- input.d.FIC[order(input.d.FIC$YEAR),]
input.would.FIC <- input.would.FIC[order(input.would.FIC$YEAR),]
input.ll.FIC <- input.ll.FIC[order(input.ll.FIC$YEAR),]
input.will.FIC <- input.will.FIC[order(input.will.FIC$YEAR),]

# corpus size in words, per decade in FIC
COHA.1830s.FIC <- c(7357832); COHA.1840s.FIC <- c(8565934)	
COHA.1850s.FIC <- c(8920782); COHA.1860s.FIC <- c(9061293)	
COHA.1870s.FIC <- c(10038777); COHA.1880s.FIC <- c(11005103)
COHA.1890s.FIC <- c(10980588); COHA.1900s.FIC <- c(11859939);
COHA.1910s.FIC <- c(11758911); COHA.1920s.FIC <- c(12368670);
COHA.1930s.FIC <- c(11696524); COHA.1940s.FIC <- c(11777039);
COHA.1950s.FIC <- c(11886489); COHA.1960s.FIC <- c(11448002);
COHA.1970s.FIC <- c(11482403); COHA.1980s.FIC <- c(11975709);
COHA.1990s.FIC <- c(13083540); COHA.2000s.FIC <- c(14365534)

CORP.SIZE.FIC <- c(COHA.1830s.FIC, COHA.1840s.FIC, COHA.1850s.FIC, 
                   COHA.1860s.FIC, COHA.1870s.FIC, COHA.1880s.FIC, 
                   COHA.1890s.FIC, COHA.1900s.FIC, COHA.1910s.FIC, 
                   COHA.1920s.FIC, COHA.1930s.FIC, COHA.1940s.FIC, 
                   COHA.1950s.FIC, COHA.1960s.FIC, COHA.1970s.FIC,
                   COHA.1980s.FIC, COHA.1990s.FIC, COHA.2000s.FIC)

# create frequency tables
# raw
freq.d.by.dec.FIC <- aggregate(cbind(count = CXN) ~ DECADE, 
                               data = input.d.FIC,
                               function(x){NROW(x)})
freq.would.by.dec.FIC <- aggregate(cbind(count = CXN) ~ DECADE, 
                                   data = input.would.FIC,
                                   function(x){NROW(x)})
freq.ll.by.dec.FIC <- aggregate(cbind(count = CXN) ~ DECADE, 
                                data = input.ll.FIC,
                                function(x){NROW(x)})
freq.will.by.dec.FIC <- aggregate(cbind(count = CXN) ~ DECADE, 
                                  data = input.will.FIC,
                                  function(x){NROW(x)})
freq.all.by.dec.FIC <- cbind(freq.d.by.dec.FIC$count, freq.would.by.dec.FIC$count,
                             freq.ll.by.dec.FIC$count, freq.will.by.dec.FIC$count)
colnames(freq.all.by.dec.FIC) <- c("SUBJ.d.V", "SUBJ.would.V",
                                   "SUBJ.ll.V", "SUBJ.will.V")

# pmw
list.freq.all.by.dec.pmw.FIC <- list()
for (b in 1:18){
  
  list.freq.all.by.dec.pmw.FIC[[b]] <- lapply(freq.all.by.dec.FIC[b,],
                                              function(x) x/CORP.SIZE.FIC[[b]]*1000000)
  
}
freq.all.by.dec.pmw.FIC <- data.frame(matrix(unlist(list.freq.all.by.dec.pmw.FIC),
                                             nrow = length(list.freq.all.by.dec.pmw.FIC),
                                             byrow = T))
colnames(freq.all.by.dec.pmw.FIC) <- c("SUBJ.d.V", "SUBJ.would.V", 
                                       "SUBJ.ll.V", "SUBJ.will.V")
TIME <- c(seq(1830, 2000, by = 10)) 
freq.all.by.dec.pmw.FIC <- cbind(TIME, freq.all.by.dec.pmw.FIC)

#============================================================================================

## FREQUENCY PROFILES OF ENCLITIC PATTERNS AND FULL FORMS

# identify trends in frequency data
cor.test(freq.all.by.dec.pmw.FIC$SUBJ.d.V, TIME, method = "kendall")
cor.test(freq.all.by.dec.pmw.FIC$SUBJ.ll.V, TIME, method = "kendall")
cor.test(freq.all.by.dec.pmw.FIC$SUBJ.would.V, TIME, method = "kendall")
cor.test(freq.all.by.dec.pmw.FIC$SUBJ.will.V, TIME, method = "kendall")

# create tables for identifying developmental stages in usage frequency (pmw)
vnc.d.FIC <- data.frame(as.factor(freq.all.by.dec.pmw.FIC$SUBJ.d.V), 
                        as.factor(TIME))
colnames(vnc.d.FIC) <- c("freq", "time")
write.table(vnc.d.FIC, file = "vnc_d_FIC.txt", row.names = F, quote = F, 
            sep = "\t")

vnc.ll.FIC <- data.frame(freq.all.by.dec.pmw.FIC$SUBJ.ll.V, TIME)
colnames(vnc.ll.FIC) <- c("freq", "time")
write.table(vnc.ll.FIC, file = "vnc_ll_FIC.txt", row.names = F, quote = F, 
            sep = "\t")

# implement VNC algorithm
vnc_sd("vnc_d_FIC.txt")
vnc_sd("vnc_ll_FIC.txt")

# define stages based on first point of inflection in scree plots
d.mean.p1.FIC <- mean(freq.all.by.dec.pmw.FIC$SUBJ.d.V[1:8]) # 1830-1900
d.mean.p2.FIC <- mean(freq.all.by.dec.pmw.FIC$SUBJ.d.V[9:18]) # 1910-2000
ll.mean.p1.FIC <- mean(freq.all.by.dec.pmw.FIC$SUBJ.ll.V[1:7]) # 1830-1890
ll.mean.p2.FIC <- mean(freq.all.by.dec.pmw.FIC$SUBJ.ll.V[8:18]) # 1900-2000

# plot developments in usage frequency (pmw)
png("Fig1.png", height = 5, width = 10, un = "in", res = 300) # saves plots in working dir.
par(mfrow = c(1, 2))
par(mar = c(4.1, 4.1, 1.1, 1.1))
par("cex" = 1.2)

plot(TIME, freq.all.by.dec.pmw.FIC$SUBJ.d.V, col = "grey50", 
     type = "b", ylim = c(0, 1700), frame.plot = F, xaxt = "n", 
     xlab = NA, ylab = "Token frequency (pmw)", pch = 15)
axis(side = 1, at = seq(1830, 2000, by=10))
segments(1830, d.mean.p1.FIC, 1900, d.mean.p1.FIC, lwd = 5, 
         col = gray(.5, alpha = .5))
segments(1910, d.mean.p2.FIC, 2000, d.mean.p2.FIC, lwd = 5, 
         col = gray(.5, alpha = .5))
lines(TIME, freq.all.by.dec.pmw.FIC$SUBJ.would.V, col = "darkgrey")
text(1855, 1400, labels = "SUBJ would V", cex = .8, col = "darkgrey")
text(1855, 400, labels = "[[SUBJ'd] V]", cex = .8)
 
plot(TIME, freq.all.by.dec.pmw.FIC$SUBJ.ll.V, col= "grey30",
     type = "b", ylim = c(0, 1700), frame.plot = F, xaxt = "n",
     xlab = NA, ylab = NA, pch = 16)
axis(side = 1, at = seq(1830, 2000, by=10))
segments(1830, ll.mean.p1.FIC, 1890, ll.mean.p1.FIC, lwd = 5, 
         col = gray(.5, alpha = .5))
segments(1900, ll.mean.p2.FIC, 2000, ll.mean.p2.FIC, lwd = 5, 
         col = gray(.5, alpha = .5))
lines(TIME, freq.all.by.dec.pmw.FIC$SUBJ.will.V, col="darkgrey")
text(1855, 1400, labels = "SUBJ will V", cex = .8, col="darkgrey")
text(1855, 250, labels = "[[SUBJ'll] V]", cex = .8)

mtext("Time", side = 1, outer = TRUE,  
      cex = 1.2, line = -1.2, adj = .53)
dev.off() # to display plots in R, select lines between png() and dev.off()

#============================================================================================

## CHANGES IN SUBJECT-HOST VARIABILITY

# prepare data, merge rare SUBJs into larger categories
input.d.FIC$SUBJ <- as.factor(input.d.FIC$SUBJ)
subj.d.FIC <- table(input.d.FIC$SUBJ, input.d.FIC$DECADE)
PRON.d <- as.vector(levels(input.d.FIC$SUBJ))
subj.d.FIC.matrix <- as.data.frame.matrix(subj.d.FIC, row.names = F)
subj.d.FIC.df <- cbind(PRON.d, subj.d.FIC.matrix)
colnames(subj.d.FIC.df)[1] <- "PRON"
subj.d.FIC.df.clean <- subj.d.FIC.df %>%
  mutate(PRON = ifelse(grepl("body", PRON),
                       paste("Xbody"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("thing", PRON), 
                       paste("Xthing"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("one", PRON),
                       paste("Xone"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("^him$|^her$|^me$|them|us", PRON), 
                       paste("OBJ.CASE"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("hims|plenty|those|wh(i|oe)",
                             PRON), 
                       paste("RESIDUE"),
                       as.character(PRON))) %>%
  group_by(PRON) %>%
  summarize_all(sum)

input.ll.FIC$SUBJ <- as.factor(input.ll.FIC$SUBJ)
subj.ll.FIC <- table(input.ll.FIC$SUBJ, input.ll.FIC$DECADE)
PRON.ll <- as.vector(levels(input.ll.FIC$SUBJ))
subj.ll.FIC.matrix <- as.data.frame.matrix(subj.ll.FIC, row.names = F)
subj.ll.FIC.df <- cbind(PRON.ll, subj.ll.FIC.matrix)
colnames(subj.ll.FIC.df)[1] <- "PRON"
subj.ll.FIC.df.clean <- subj.ll.FIC.df %>%
  mutate(PRON = ifelse(grepl("body", PRON), 
                       paste("Xbody"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("thing", PRON), 
                       paste("Xthing"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("one", PRON), 
                       paste("Xone"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("^him$|^her$|^me$|them|us", PRON), 
                       paste("OBJ.CASE"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("self|our|^his$|th.se|mine|which", PRON), 
                       paste("RESIDUE"),
                       as.character(PRON))) %>%
  group_by(PRON) %>%
  summarize_all(sum)

input.would.FIC$SUBJ <- as.factor(input.would.FIC$SUBJ)
subj.would.FIC <- table(input.would.FIC$SUBJ, input.would.FIC$DECADE)
PRON.would <- as.vector(levels(input.would.FIC$SUBJ))
subj.would.FIC.matrix <- as.data.frame.matrix(subj.would.FIC, row.names = F)
subj.would.FIC.df <- cbind(PRON.would, subj.would.FIC.matrix)
colnames(subj.would.FIC.df)[1] <- "PRON"
subj.would.FIC.df.clean <- subj.would.FIC.df %>%
  mutate(PRON = ifelse(grepl("body", PRON),
                       paste("Xbody"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("thing", PRON), 
                       paste("Xthing"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("one", PRON),
                       paste("Xone"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("^him$|^her$|^me$|them$|us", 
                             PRON), 
                       paste("OBJ.CASE"),
                       as.character(PRON))) %>%
  mutate(PRON = 
           ifelse(grepl("th.se|mine|^his$|rs$|sel(f|v)|plenty|wh(i|oe|om|os)|lots|our|other",
                        PRON),
                  paste("RESIDUE"),
                  as.character(PRON))) %>%
  group_by(PRON) %>%
  summarize_all(sum)

input.will.FIC$SUBJ <- as.factor(input.will.FIC$SUBJ)
subj.will.FIC <- table(input.will.FIC$SUBJ, input.will.FIC$DECADE)
PRON.will <- as.vector(levels(input.will.FIC$SUBJ))
subj.will.FIC.matrix <- as.data.frame.matrix(subj.will.FIC, row.names = F)
subj.will.FIC.df <- cbind(PRON.will, subj.will.FIC.matrix)
colnames(subj.will.FIC.df)[1] <- "PRON"
subj.will.FIC.df.clean <- subj.will.FIC.df %>%
  mutate(PRON = ifelse(grepl("body", PRON),
                       paste("Xbody"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("thing", PRON), 
                       paste("Xthing"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("one", PRON),
                       paste("Xone"),
                       as.character(PRON))) %>%
  mutate(PRON = ifelse(grepl("^him$|^her$|^me$|them|us", PRON), 
                       paste("OBJ.CASE"),
                       as.character(PRON))) %>%
  mutate(PRON = 
           ifelse(grepl("mine|rs$|^his$|sel(f|v)|plenty|th.se|wh(i|oe|om|os)|lots|our|other",
                        PRON),
                  paste("RESIDUE"),
                  as.character(PRON))) %>%
  group_by(PRON) %>%
  summarize_all(sum)

# calculate normalized entropy (Hnorm) (Gries 2021)
perc.list.d.FIC <- list() # SUBJ'd
for (m in 2:19){
  
  perc.list.d.FIC[m] <- subj.d.FIC.df.clean[m]/sum(subj.d.FIC.df.clean[m])
  
}

hnorm.list.d.FIC <- list()
for (p in 2:19){
  
  hnorm.list.d.FIC[p] <- -sum(perc.list.d.FIC[[p]]*logw0(perc.list.d.FIC[[p]]))/
    logw0(length(perc.list.d.FIC[[p]]))
  
}
hnorm.list.d.vec.FIC <- unlist(hnorm.list.d.FIC)

perc.list.ll.FIC <- list() # SUBJ'll
for (m in 2:19){
  
  perc.list.ll.FIC[m] <- subj.ll.FIC.df.clean[m]/sum(subj.ll.FIC.df.clean[m])
  
}
hnorm.list.ll.FIC <- list()
for (p in 2:19){
  
  hnorm.list.ll.FIC[p] <- -sum(perc.list.ll.FIC[[p]]*logw0(perc.list.ll.FIC[[p]]))/
    logw0(length(perc.list.ll.FIC[[p]]))
  
}
hnorm.list.ll.vec.FIC <- unlist(hnorm.list.ll.FIC)

perc.list.would.FIC <- list() # SUBJ would
for (m in 2:19){
  
  perc.list.would.FIC[m] <- subj.would.FIC.df.clean[m]/sum(subj.would.FIC.df.clean[m])
  
}

hnorm.list.would.FIC <- list()
for (p in 2:19){
  
  hnorm.list.would.FIC[p] <- -sum(perc.list.would.FIC[[p]]*logw0(perc.list.would.FIC[[p]]))/
    logw0(length(perc.list.would.FIC[[p]]))
  
}
hnorm.list.would.vec.FIC <- unlist(hnorm.list.would.FIC)

perc.list.will.FIC <- list() # SUBJ will
for (m in 2:19){
  
  perc.list.will.FIC[m] <- subj.will.FIC.df.clean[m]/sum(subj.will.FIC.df.clean[m])
  
}

hnorm.list.will.FIC <- list()
for (p in 2:19){
  
  hnorm.list.will.FIC[p] <- -sum(perc.list.will.FIC[[p]]*logw0(perc.list.will.FIC[[p]]))/
    logw0(length(perc.list.will.FIC[[p]]))
  
}
hnorm.list.will.vec.FIC <- unlist(hnorm.list.will.FIC)

# identify trends in Hnorm data
cor.test(TIME, hnorm.list.d.vec.FIC, method = "kendall")
cor.test(TIME, hnorm.list.ll.vec.FIC, method = "kendall")
cor.test(TIME, hnorm.list.would.vec.FIC, method = "kendall")
cor.test(TIME, hnorm.list.will.vec.FIC, method = "kendall")

# create tables for identifying developmental stages in Hnorm
vnc.d.hnorm.FIC <- data.frame(as.factor(hnorm.list.d.vec.FIC), 
                              as.factor(TIME))
colnames(vnc.d.hnorm.FIC) <- c("h.rel", "time")
write.table(vnc.d.hnorm.FIC, file = "vnc_d_hnorm_FIC.txt", row.names = F,
            quote = F, sep = "\t")

vnc.ll.hnorm.FIC <- data.frame(as.factor(hnorm.list.ll.vec.FIC), 
                               as.factor(TIME))
colnames(vnc.ll.hnorm.FIC) <- c("h.rel", "time")
write.table(vnc.ll.hnorm.FIC, file = "vnc_ll_hnorm_FIC.txt", row.names = F,
            quote = F, sep = "\t")

# implement VNC algorithm
vnc_sd("vnc_d_hnorm_FIC.txt")
vnc_sd("vnc_ll_hnorm_FIC.txt")

# define stages based on first point of inflection in scree plots
d.mean.p1.h.norm.FIC <- mean(hnorm.list.d.vec.FIC[1:10]) # 1830-1920
d.mean.p2.h.norm.FIC <- mean(hnorm.list.d.vec.FIC[11:18]) # 1930-2000
ll.mean.p1.h.norm.FIC <- mean(hnorm.list.ll.vec.FIC[1:2]) # 1830-1840
ll.mean.p2.h.norm.FIC <- mean(hnorm.list.ll.vec.FIC[3:18]) # 1850-2000

# plot developments in Hnorm
png("Fig2.png", height = 5, width = 10, un = "in", res = 300)
par(mfrow = c(1, 2))
par(mar = c(4.1, 4.1, 1.1, 1.1))
par("cex" = 1.2)

plot(TIME, hnorm.list.d.vec.FIC, type = "b", pch = 15, 
     ylim = c(.45, .9), xlab = NA, col = "grey50",
     ylab = expression("Normalized entropy H"[" norm"]),
     frame.plot = F, xaxt="n")
axis(side = 1, at = seq(1830, 2000, by = 10))
segments(1830, d.mean.p1.h.norm.FIC, 1920, d.mean.p1.h.norm.FIC, lwd = 5,
         col = gray(.5, alpha = .5))
segments(1930, d.mean.p2.h.norm.FIC, 2000, d.mean.p2.h.norm.FIC, lwd = 5,
         col = gray(.5, alpha = .5))
lines(TIME, hnorm.list.would.vec.FIC, col="darkgrey")
text(1980, .78, labels = "SUBJ would V", cex = .8, col = "darkgrey")
text(1980, .58, labels = "[[SUBJ'd] V]", cex = .8)

plot(TIME, hnorm.list.ll.vec.FIC, type = "b", pch = 16,
     ylim = c(.45, .9), xlab = NA, ylab = NA, col = "grey30",
     frame.plot = F, xaxt = "n")
axis(side = 1, at = seq(1830, 2000, by = 10))
segments(1830, ll.mean.p1.h.norm.FIC, 1840, ll.mean.p1.h.norm.FIC, lwd = 5,
         col = gray(.5, alpha = .5))
segments(1850, ll.mean.p2.h.norm.FIC, 2000, ll.mean.p2.h.norm.FIC, lwd = 5,
         col = gray(.5, alpha = .5))
lines(TIME, hnorm.list.will.vec.FIC, col = "darkgrey")
text(1980, .78, labels = "SUBJ will V", cex = .8, col = "darkgrey")
text(1980, .52, labels = "[[SUBJ'll] V]", cex = .8)

mtext("Time", side = 1, outer = TRUE,  
      cex = 1.2, line = -1.2, adj = .53)
dev.off()

# prepare data to assess relative frequency changes for selected pronouns
subj.d.long.select.FIC <- subj.d.FIC.df %>%
  gather(DECADE, FREQ, -PRON) %>%
  filter(DECADE != "TOTAL") %>%
  mutate(PRON = ifelse(PRON %in% c("she", "he"), paste("(s)he"),
                       as.character(PRON))) %>%
  filter(PRON == "I" | PRON == "you" | PRON == "(s)he" | 
           PRON == "we" | PRON == "they" | PRON == "it" |
           PRON == "that"| PRON == "there" | PRON == "who") %>%
  group_by(DECADE, PRON) %>%
  summarize_all(sum) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(DECADE = gsub("s", "", DECADE)) %>%
  mutate(DECADE = as.integer(DECADE)) %>%
  mutate(FREQ.PMW = (FREQ/rep(CORP.SIZE.FIC, each = 9)*1000000)) %>%
  mutate(PRON = paste(PRON, "_d", sep = "")) %>%
  select(-3)

subj.would.long.select.FIC <- subj.would.FIC.df %>%
  gather(DECADE, FREQ, -PRON) %>%
  mutate(PRON = ifelse(PRON %in% c("she", "he"), paste("(s)he"),
                       as.character(PRON))) %>%
  filter(PRON == "I" | PRON == "you" | PRON == "(s)he" | 
           PRON == "we" | PRON == "they" | PRON == "it" |
           PRON == "that" | PRON == "there" | PRON == "who") %>%
  group_by(DECADE, PRON) %>%
  summarize_all(sum) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(DECADE = gsub("s", "", DECADE)) %>%
  mutate(DECADE = as.integer(DECADE)) %>%
  mutate(FREQ.PMW = (FREQ/rep(CORP.SIZE.FIC, each = 9)*1000000)) %>%
  mutate(PRON = paste(PRON, "_would", sep = "")) %>%
  select(-3)

subj.d.would.select.FIC <- rbind(subj.d.long.select.FIC, 
                                 subj.would.long.select.FIC)

subj.ll.long.select.FIC <- subj.ll.FIC.df %>%
  gather(DECADE, FREQ, -PRON) %>%
  filter(DECADE != "TOTAL") %>%
  mutate(PRON = ifelse(PRON %in% c("she", "he"), paste("(s)he"),
                       as.character(PRON))) %>%
  filter(PRON == "I" | PRON == "you" | PRON == "(s)he" | 
           PRON == "we" | PRON == "they" | PRON == "it" |
           PRON == "that"| PRON == "there" | PRON == "who") %>%
  group_by(DECADE, PRON) %>%
  summarize_all(sum) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(DECADE = gsub("s", "", DECADE)) %>%
  mutate(DECADE = as.integer(DECADE)) %>%
  mutate(FREQ.PMW = (FREQ/rep(CORP.SIZE.FIC, each = 9)*1000000)) %>%
  mutate(PRON = paste(PRON, "_ll", sep = "")) %>%
  select(-3)

subj.will.long.select.FIC <- subj.will.FIC.df %>%
  gather(DECADE, FREQ, -PRON) %>%
  mutate(PRON = ifelse(PRON %in% c("she", "he"), paste("(s)he"),
                       as.character(PRON))) %>%
  filter(PRON == "I" | PRON == "you" | PRON == "(s)he" | 
           PRON == "we" | PRON == "they" | PRON == "it" |
           PRON == "that"| PRON == "there" | PRON == "who") %>%
  group_by(DECADE, PRON) %>%
  summarize_all(sum) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(DECADE = gsub("s", "", DECADE)) %>%
  mutate(DECADE = as.integer(DECADE)) %>%
  mutate(FREQ.PMW = (FREQ/rep(CORP.SIZE.FIC, each = 9)*1000000)) %>%
  mutate(PRON = paste(PRON, "_will", sep = "")) %>%
  select(-3)

subj.ll.will.select.FIC <- rbind(subj.ll.long.select.FIC, 
                                 subj.will.long.select.FIC)

# prepare data for spineplot
d.would.subj.I.FIC <- subj.d.would.select.FIC %>%
  filter(PRON %in% c("I_d", "I_would")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
d.would.subj.I.FIC$PRON <- relevel(d.would.subj.I.FIC$PRON, "I_would")

ll.will.subj.I.FIC <- subj.ll.will.select.FIC %>%
  filter(PRON %in% c("I_ll", "I_will")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
ll.will.subj.I.FIC$PRON <- relevel(ll.will.subj.I.FIC$PRON, "I_will")

d.would.subj.you.FIC <- subj.d.would.select.FIC %>%
  filter(PRON %in% c("you_d", "you_would")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
d.would.subj.you.FIC$PRON <- relevel(d.would.subj.you.FIC$PRON, "you_would")

ll.will.subj.you.FIC <- subj.ll.will.select.FIC %>%
  filter(PRON %in% c("you_ll", "you_will")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
ll.will.subj.you.FIC$PRON <- relevel(ll.will.subj.you.FIC$PRON, "you_will")

d.would.subj.we.FIC <- subj.d.would.select.FIC %>%
  filter(PRON %in% c("we_d", "we_would")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
d.would.subj.we.FIC$PRON <- relevel(d.would.subj.we.FIC$PRON, "we_would")

ll.will.subj.we.FIC <- subj.ll.will.select.FIC %>%
  filter(PRON %in% c("we_ll", "we_will")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
ll.will.subj.we.FIC$PRON <- relevel(ll.will.subj.we.FIC$PRON, "we_will")

d.would.subj.s.he.FIC <- subj.d.would.select.FIC %>%
  filter(PRON %in% c("(s)he_d", "(s)he_would")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
d.would.subj.s.he.FIC$PRON <- relevel(d.would.subj.s.he.FIC$PRON, "(s)he_would")

ll.will.subj.s.he.FIC <- subj.ll.will.select.FIC %>%
  filter(PRON %in% c("(s)he_ll", "(s)he_will")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
ll.will.subj.s.he.FIC$PRON <- relevel(ll.will.subj.s.he.FIC$PRON, "(s)he_will")

d.would.subj.they.FIC <- subj.d.would.select.FIC %>%
  filter(PRON %in% c("they_d", "they_would")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
d.would.subj.they.FIC$PRON <- relevel(d.would.subj.they.FIC$PRON, "they_would")

ll.will.subj.they.FIC <- subj.ll.will.select.FIC %>%
  filter(PRON %in% c("they_ll", "they_will")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
ll.will.subj.they.FIC$PRON <- relevel(ll.will.subj.they.FIC$PRON, "they_will")

d.would.subj.it.FIC <- subj.d.would.select.FIC %>%
  filter(PRON %in% c("it_d", "it_would")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
d.would.subj.it.FIC$PRON <- relevel(d.would.subj.it.FIC$PRON, "it_would")

ll.will.subj.it.FIC <- subj.ll.will.select.FIC %>%
  filter(PRON %in% c("it_ll", "it_will")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
ll.will.subj.it.FIC$PRON <- relevel(ll.will.subj.it.FIC$PRON, "it_will")

d.would.subj.that.FIC <- subj.d.would.select.FIC %>%
  filter(PRON %in% c("that_d", "that_would")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
d.would.subj.that.FIC$PRON <- relevel(d.would.subj.that.FIC$PRON, "that_would")

ll.will.subj.that.FIC <- subj.ll.will.select.FIC %>%
  filter(PRON %in% c("that_ll", "that_will")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
ll.will.subj.that.FIC$PRON <- relevel(ll.will.subj.that.FIC$PRON, "that_will")

d.would.subj.there.FIC <- subj.d.would.select.FIC %>%
  filter(PRON %in% c("there_d", "there_would")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
d.would.subj.there.FIC$PRON <- relevel(d.would.subj.there.FIC$PRON, "there_would")

ll.will.subj.there.FIC <- subj.ll.will.select.FIC %>%
  filter(PRON %in% c("there_ll", "there_will")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
ll.will.subj.there.FIC$PRON <- relevel(ll.will.subj.there.FIC$PRON, "there_will")

d.would.subj.who.FIC <- subj.d.would.select.FIC %>%
  filter(PRON %in% c("who_d", "who_would")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
d.would.subj.who.FIC$PRON <- relevel(d.would.subj.who.FIC$PRON, "who_would")

ll.will.subj.who.FIC <- subj.ll.will.select.FIC %>%
  filter(PRON %in% c("who_ll", "who_will")) %>%
  droplevels() %>%
  mutate(PRON = as.factor(PRON))
ll.will.subj.who.FIC$PRON <- relevel(ll.will.subj.who.FIC$PRON, "who_will")

# plot relative frequency changes for above selection
png("Fig3.png", height = 10.5, width = 8.5, un = "in", res = 300)
par(mfrow = c(6, 3))
par(mar = c(2.1, 2.1, 1.1, 2.1))
par("cex" = 1.2)

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = d.would.subj.I.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey50", "grey90"))
legend("bottomright", legend = c("I would V", 
                                 "I'd V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey50"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = d.would.subj.you.FIC),
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey50", "grey90"))
legend("bottomright", legend = c("you would V", 
                                 "you'd V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey50"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = d.would.subj.we.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey50", "grey90"))
legend("bottomright", legend = c("we would V", 
                                 "we'd V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey50"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = d.would.subj.s.he.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey50", "grey90"))
legend("topleft", legend = c("(s)he would V", 
                             "(s)he'd V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey50"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = d.would.subj.they.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey50", "grey90"))
legend("topleft", legend = c("they would V", 
                             "they'd V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey50"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = d.would.subj.it.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey50", "grey90"))
legend("topleft", legend = c("it would V", 
                             "it'd V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey50"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = d.would.subj.that.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey50", "grey90"))
legend("topleft", legend = c("that would V", 
                             "that'd V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey50"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = d.would.subj.there.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey50", "grey90"))
legend("topleft", legend = c("there would V", 
                             "there'd V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey50"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = d.would.subj.who.FIC),
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey50", "grey90"))
legend("topleft", legend = c("who would V", 
                             "who'd V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey50"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = ll.will.subj.I.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey30", "grey90"))
legend("bottomright", legend = c("I will V", 
                                 "I'll V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey30"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = ll.will.subj.you.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey30", "grey90"))
legend("bottomright", legend = c("you will V", 
                                 "you'll V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey30"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = ll.will.subj.we.FIC),
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey30", "grey90"))
legend("bottomright", legend = c("we will V", 
                                 "we'll V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey30"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = ll.will.subj.s.he.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey30", "grey90"))
legend("bottomright", legend = c("(s)he will V", 
                                 "(s)he'll V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey30"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = ll.will.subj.they.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey30", "grey90"))
legend("bottomright", legend = c("they will V", 
                                 "they'll V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey30"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = ll.will.subj.it.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey30", "grey90"))
legend("topleft", legend = c("it will V", 
                             "it'll V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey30"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = ll.will.subj.that.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey30", "grey90"))
legend("topleft", legend = c("that will V", 
                             "that'll V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey30"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = ll.will.subj.there.FIC), 
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey30", "grey90"))
legend("topleft", legend = c("there will V", 
                             "there'll V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey30"))

spineplot(xtabs(FREQ.PMW ~ DECADE + PRON, data = ll.will.subj.who.FIC),
          off = 0, xlab = NA, ylab = NA, 
          yaxlabels = NA, col = c("grey30", "grey90"))
legend("topleft", legend = c("who will V", 
                             "who'll V"),
       pch = c(22,22), cex = .72,
       pt.bg = c("grey90", "grey30"))
dev.off()

#============================================================================================

## CHANGES IN THE PRODUCTIVITY IN THE VERB SLOT

# Note: Results will not be 100% identical to publication, given that the LNRE models are
# fitted to random samples. Trends should be fairly stable.

# create subsets
subsets.d.FIC <- split(input.d.FIC, input.d.FIC$DECADE)
subsets.ll.FIC <- split(input.ll.FIC, input.ll.FIC$DECADE)
subsets.would.FIC <- split(input.would.FIC, input.would.FIC$DECADE)
subsets.will.FIC <- split(input.will.FIC, input.will.FIC$DECADE)

# create vector lists of verbs for each decade 
vec.list.subsets.d.FIC <- list()
for (g in 1:18){
  
  vec.list.subsets.d.FIC[[g]] <- sapply(subsets.d.FIC[[g]]$VERB,
                                        function(x) as.vector(x))
}

vec.list.subsets.ll.FIC <- list()
for (g in 1:18){
  
  vec.list.subsets.ll.FIC[[g]] <- sapply(subsets.ll.FIC[[g]]$VERB,
                                         function(x) as.vector(x))
}

vec.list.subsets.would.FIC <- list()
for (g in 1:18){
  
  vec.list.subsets.would.FIC[[g]] <- sapply(subsets.would.FIC[[g]]$VERB,
                                            function(x) as.vector(x))
}

vec.list.subsets.will.FIC <- list()
for (g in 1:18){
  
  vec.list.subsets.will.FIC[[g]] <- sapply(subsets.will.FIC[[g]]$VERB,
                                           function(x) as.vector(x))
}

# merge periods into larger clusters
vec.d.clust.1.FIC <- c(unlist(vec.list.subsets.d.FIC[1:4])) # N<10,000! less reliable
vec.d.clust.2.FIC <- c(unlist(vec.list.subsets.d.FIC[5:6])) # N<10,000! less reliable
vec.d.clust.3.FIC <- c(unlist(vec.list.subsets.d.FIC[7:8])) # N<10,000! less reliable
vec.d.clust.4.FIC <- c(unlist(vec.list.subsets.d.FIC[9:10]))
vec.d.clust.5.FIC <- c(unlist(vec.list.subsets.d.FIC[11:12]))
vec.d.clust.6.FIC <- c(unlist(vec.list.subsets.d.FIC[13:14]))
vec.d.clust.7.FIC <- c(unlist(vec.list.subsets.d.FIC[15:16]))
vec.d.clust.8.FIC <- c(unlist(vec.list.subsets.d.FIC[17:18]))

vec.ll.clust.1.FIC <- c(unlist(vec.list.subsets.ll.FIC[1:4]))
vec.ll.clust.2.FIC <- c(unlist(vec.list.subsets.ll.FIC[5:6]))
vec.ll.clust.3.FIC <- c(unlist(vec.list.subsets.ll.FIC[7:8]))
vec.ll.clust.4.FIC <- c(unlist(vec.list.subsets.ll.FIC[9:10]))
vec.ll.clust.5.FIC <- c(unlist(vec.list.subsets.ll.FIC[11:12]))
vec.ll.clust.6.FIC <- c(unlist(vec.list.subsets.ll.FIC[13:14]))
vec.ll.clust.7.FIC <- c(unlist(vec.list.subsets.ll.FIC[15:16]))
vec.ll.clust.8.FIC <- c(unlist(vec.list.subsets.ll.FIC[17:18]))

vec.would.clust.1.FIC <- c(unlist(vec.list.subsets.would.FIC[1:4]))
vec.would.clust.2.FIC <- c(unlist(vec.list.subsets.would.FIC[5:6]))
vec.would.clust.3.FIC <- c(unlist(vec.list.subsets.would.FIC[7:8]))
vec.would.clust.4.FIC <- c(unlist(vec.list.subsets.would.FIC[9:10]))
vec.would.clust.5.FIC <- c(unlist(vec.list.subsets.would.FIC[11:12]))
vec.would.clust.6.FIC <- c(unlist(vec.list.subsets.would.FIC[13:14]))
vec.would.clust.7.FIC <- c(unlist(vec.list.subsets.would.FIC[15:16]))
vec.would.clust.8.FIC <- c(unlist(vec.list.subsets.would.FIC[17:18]))

vec.will.clust.1.FIC <- c(unlist(vec.list.subsets.will.FIC[1:4]))
vec.will.clust.2.FIC <- c(unlist(vec.list.subsets.will.FIC[5:6]))
vec.will.clust.3.FIC <- c(unlist(vec.list.subsets.will.FIC[7:8]))
vec.will.clust.4.FIC <- c(unlist(vec.list.subsets.will.FIC[9:10]))
vec.will.clust.5.FIC <- c(unlist(vec.list.subsets.will.FIC[11:12]))
vec.will.clust.6.FIC <- c(unlist(vec.list.subsets.will.FIC[13:14]))
vec.will.clust.7.FIC <- c(unlist(vec.list.subsets.will.FIC[15:16]))
vec.will.clust.8.FIC <- c(unlist(vec.list.subsets.will.FIC[17:18]))

# create frequency spectrum objects (spc)
spc.d.clust.1.FIC <- vec2spc(vec.d.clust.1.FIC)
spc.d.clust.2.FIC <- vec2spc(vec.d.clust.2.FIC)
spc.d.clust.3.FIC <- vec2spc(vec.d.clust.3.FIC)
spc.d.clust.4.FIC <- vec2spc(vec.d.clust.4.FIC)
spc.d.clust.5.FIC <- vec2spc(vec.d.clust.5.FIC)
spc.d.clust.6.FIC <- vec2spc(vec.d.clust.6.FIC)
spc.d.clust.7.FIC <- vec2spc(vec.d.clust.7.FIC)
spc.d.clust.8.FIC <- vec2spc(vec.d.clust.8.FIC)

spc.ll.clust.1.FIC <- vec2spc(vec.ll.clust.1.FIC)
spc.ll.clust.2.FIC <- vec2spc(vec.ll.clust.2.FIC)
spc.ll.clust.3.FIC <- vec2spc(vec.ll.clust.3.FIC)
spc.ll.clust.4.FIC <- vec2spc(vec.ll.clust.4.FIC)
spc.ll.clust.5.FIC <- vec2spc(vec.ll.clust.5.FIC)
spc.ll.clust.6.FIC <- vec2spc(vec.ll.clust.6.FIC)
spc.ll.clust.7.FIC <- vec2spc(vec.ll.clust.7.FIC)
spc.ll.clust.8.FIC <- vec2spc(vec.ll.clust.8.FIC)

spc.would.clust.1.FIC <- vec2spc(vec.would.clust.1.FIC)
spc.would.clust.2.FIC <- vec2spc(vec.would.clust.2.FIC)
spc.would.clust.3.FIC <- vec2spc(vec.would.clust.3.FIC)
spc.would.clust.4.FIC <- vec2spc(vec.would.clust.4.FIC)
spc.would.clust.5.FIC <- vec2spc(vec.would.clust.5.FIC)
spc.would.clust.6.FIC <- vec2spc(vec.would.clust.6.FIC)
spc.would.clust.7.FIC <- vec2spc(vec.would.clust.7.FIC)
spc.would.clust.8.FIC <- vec2spc(vec.would.clust.8.FIC)

spc.will.clust.1.FIC <- vec2spc(vec.will.clust.1.FIC)
spc.will.clust.2.FIC <- vec2spc(vec.will.clust.2.FIC)
spc.will.clust.3.FIC <- vec2spc(vec.will.clust.3.FIC)
spc.will.clust.4.FIC <- vec2spc(vec.will.clust.4.FIC)
spc.will.clust.5.FIC <- vec2spc(vec.will.clust.5.FIC)
spc.will.clust.6.FIC <- vec2spc(vec.will.clust.6.FIC)
spc.will.clust.7.FIC <- vec2spc(vec.will.clust.7.FIC)
spc.will.clust.8.FIC <- vec2spc(vec.will.clust.8.FIC)

# create lists for bootstrapped potential productivity values
boot.list.d.clust.1.FIC <-
  boot.list.d.clust.2.FIC <- 
  boot.list.d.clust.3.FIC <- 
  boot.list.d.clust.4.FIC <- 
  boot.list.d.clust.5.FIC <- 
  boot.list.d.clust.6.FIC <- 
  boot.list.d.clust.7.FIC <- 
  boot.list.d.clust.8.FIC <- list()

boot.list.ll.clust.1.FIC <- 
  boot.list.ll.clust.2.FIC <- 
  boot.list.ll.clust.3.FIC <- 
  boot.list.ll.clust.4.FIC <- 
  boot.list.ll.clust.5.FIC <- 
  boot.list.ll.clust.6.FIC <- 
  boot.list.ll.clust.7.FIC <- 
  boot.list.ll.clust.8.FIC <- list()

boot.list.would.clust.1.FIC <- 
  boot.list.would.clust.2.FIC <- 
  boot.list.would.clust.3.FIC <- 
  boot.list.would.clust.4.FIC <- 
  boot.list.would.clust.5.FIC <- 
  boot.list.would.clust.6.FIC <- 
  boot.list.would.clust.7.FIC <- 
  boot.list.would.clust.8.FIC <- list()

boot.list.will.clust.1.FIC <- 
  boot.list.will.clust.2.FIC <- 
  boot.list.will.clust.3.FIC <- 
  boot.list.will.clust.4.FIC <- 
  boot.list.will.clust.5.FIC <- 
  boot.list.will.clust.6.FIC <- 
  boot.list.will.clust.7.FIC <- 
  boot.list.will.clust.8.FIC <- list()

# *'d*; fit LNRE models to 1,000 parametric bootstrap samples for each cluster
boot.lnre.d.clust.1.FIC <- lnre("fzm", spc.d.clust.1.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.d.clust.2.FIC <- lnre("fzm", spc.d.clust.2.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.d.clust.3.FIC <- lnre("fzm", spc.d.clust.3.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.d.clust.4.FIC <- lnre("fzm", spc.d.clust.4.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.d.clust.5.FIC <- lnre("fzm", spc.d.clust.5.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.d.clust.6.FIC <- lnre("fzm", spc.d.clust.6.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.d.clust.7.FIC <- lnre("fzm", spc.d.clust.7.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.d.clust.8.FIC <- lnre("fzm", spc.d.clust.8.FIC, m.max = 1, bootstrap = 1000)

# obtain hapax-growth curves from models and calculate PP-values, then remove large files
extr.hgc.d.1 <- lapply(boot.lnre.d.clust.1.FIC$bootstrap,
                       function(x) lnre.vgc(x, 1:200000,
                                            m.max = 1))
for (j in seq_along(extr.hgc.d.1)){
  
  boot.list.d.clust.1.FIC[j] <- lapply(extr.hgc.d.1[[j]]$V1[200000],
                                       function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.d.2 <- lapply(boot.lnre.d.clust.2.FIC$bootstrap,
                       function(x) lnre.vgc(x, 1:200000,
                                            m.max = 1))
for (j in seq_along(extr.hgc.d.2)){
  
  boot.list.d.clust.2.FIC[j] <- lapply(extr.hgc.d.2[[j]]$V1[200000],
                                       function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.d.3 <- lapply(boot.lnre.d.clust.3.FIC$bootstrap,
                       function(x) lnre.vgc(x, 1:200000,
                                            m.max = 1))
for (j in seq_along(extr.hgc.d.3)){
  
  boot.list.d.clust.3.FIC[j] <- lapply(extr.hgc.d.3[[j]]$V1[200000],
                                       function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.d.4 <- lapply(boot.lnre.d.clust.4.FIC$bootstrap,
                       function(x) lnre.vgc(x, 1:200000,
                                            m.max = 1))
for (j in seq_along(extr.hgc.d.4)){
  
  boot.list.d.clust.4.FIC[j] <- lapply(extr.hgc.d.4[[j]]$V1[200000],
                                       function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.d.5 <- lapply(boot.lnre.d.clust.5.FIC$bootstrap,
                       function(x) lnre.vgc(x, 1:200000,
                                            m.max = 1))
for (j in seq_along(extr.hgc.d.5)){
  
  boot.list.d.clust.5.FIC[j] <- lapply(extr.hgc.d.5[[j]]$V1[200000],
                                       function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.d.6 <- lapply(boot.lnre.d.clust.6.FIC$bootstrap,
                       function(x) lnre.vgc(x, 1:200000,
                                            m.max = 1))
for (j in seq_along(extr.hgc.d.6)){
  
  boot.list.d.clust.6.FIC[j] <- lapply(extr.hgc.d.6[[j]]$V1[200000],
                                       function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.d.7 <- lapply(boot.lnre.d.clust.7.FIC$bootstrap,
                       function(x) lnre.vgc(x, 1:200000,
                                            m.max = 1))
for (j in seq_along(extr.hgc.d.7)){
  
  boot.list.d.clust.7.FIC[j] <- lapply(extr.hgc.d.7[[j]]$V1[200000],
                                       function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.d.8 <- lapply(boot.lnre.d.clust.8.FIC$bootstrap,
                       function(x) lnre.vgc(x, 1:200000,
                                            m.max = 1))
for (j in seq_along(extr.hgc.d.8)){
  
  boot.list.d.clust.8.FIC[j] <- lapply(extr.hgc.d.8[[j]]$V1[200000],
                                       function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

# repeat procedure for *'ll*
boot.lnre.ll.clust.1.FIC <- lnre("fzm", spc.ll.clust.1.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.ll.clust.2.FIC <- lnre("fzm", spc.ll.clust.2.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.ll.clust.3.FIC <- lnre("fzm", spc.ll.clust.3.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.ll.clust.4.FIC <- lnre("fzm", spc.ll.clust.4.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.ll.clust.5.FIC <- lnre("fzm", spc.ll.clust.5.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.ll.clust.6.FIC <- lnre("fzm", spc.ll.clust.6.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.ll.clust.7.FIC <- lnre("fzm", spc.ll.clust.7.FIC, m.max = 1, bootstrap = 1000)
boot.lnre.ll.clust.8.FIC <- lnre("fzm", spc.ll.clust.8.FIC, m.max = 1, bootstrap = 1000)

# hapax-growth curves
extr.hgc.ll.1 <- lapply(boot.lnre.ll.clust.1.FIC$bootstrap,
                        function(x) lnre.vgc(x, 1:200000,
                                             m.max = 1))
for (j in seq_along(extr.hgc.ll.1)){
  
  boot.list.ll.clust.1.FIC[j] <- lapply(extr.hgc.ll.1[[j]]$V1[200000],
                                        function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.ll.2 <- lapply(boot.lnre.ll.clust.2.FIC$bootstrap,
                        function(x) lnre.vgc(x, 1:200000,
                                             m.max = 1))
for (j in seq_along(extr.hgc.ll.2)){
  
  boot.list.ll.clust.2.FIC[j] <- lapply(extr.hgc.ll.2[[j]]$V1[200000],
                                        function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.ll.3 <- lapply(boot.lnre.ll.clust.3.FIC$bootstrap,
                        function(x) lnre.vgc(x, 1:200000,
                                             m.max = 1))
for (j in seq_along(extr.hgc.ll.3)){
  
  boot.list.ll.clust.3.FIC[j] <- lapply(extr.hgc.ll.3[[j]]$V1[200000],
                                        function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.ll.4 <- lapply(boot.lnre.ll.clust.4.FIC$bootstrap,
                        function(x) lnre.vgc(x, 1:200000,
                                             m.max = 1))
for (j in seq_along(extr.hgc.ll.4)){
  
  boot.list.ll.clust.4.FIC[j] <- lapply(extr.hgc.ll.4[[j]]$V1[200000],
                                        function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.ll.5 <- lapply(boot.lnre.ll.clust.5.FIC$bootstrap,
                        function(x) lnre.vgc(x, 1:200000,
                                             m.max = 1))
for (j in seq_along(extr.hgc.ll.5)){
  
  boot.list.ll.clust.5.FIC[j] <- lapply(extr.hgc.ll.5[[j]]$V1[200000],
                                        function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.ll.6 <- lapply(boot.lnre.ll.clust.6.FIC$bootstrap,
                        function(x) lnre.vgc(x, 1:200000,
                                             m.max = 1))
for (j in seq_along(extr.hgc.ll.6)){
  
  boot.list.ll.clust.6.FIC[j] <- lapply(extr.hgc.ll.6[[j]]$V1[200000],
                                        function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.ll.7 <- lapply(boot.lnre.ll.clust.7.FIC$bootstrap,
                        function(x) lnre.vgc(x, 1:200000,
                                             m.max = 1))
for (j in seq_along(extr.hgc.ll.7)){
  
  boot.list.ll.clust.7.FIC[j] <- lapply(extr.hgc.ll.7[[j]]$V1[200000],
                                        function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.ll.8 <- lapply(boot.lnre.ll.clust.8.FIC$bootstrap,
                        function(x) lnre.vgc(x, 1:200000,
                                             m.max = 1))
for (j in seq_along(extr.hgc.ll.8)){
  
  boot.list.ll.clust.8.FIC[j] <- lapply(extr.hgc.ll.8[[j]]$V1[200000],
                                        function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

# repeat procedure for *would*
boot.lnre.would.clust.1.FIC <- lnre("fzm", spc.would.clust.1.FIC, m.max = 1,
                                     bootstrap = 1000)
boot.lnre.would.clust.2.FIC <- lnre("fzm", spc.would.clust.2.FIC, m.max = 1, 
                                     bootstrap = 1000)
boot.lnre.would.clust.3.FIC <- lnre("fzm", spc.would.clust.3.FIC, m.max = 1,
                                     bootstrap = 1000)
boot.lnre.would.clust.4.FIC <- lnre("fzm", spc.would.clust.4.FIC, m.max = 1,
                                     bootstrap = 1000)
boot.lnre.would.clust.5.FIC <- lnre("fzm", spc.would.clust.5.FIC, m.max = 1, 
                                     bootstrap = 1000)
boot.lnre.would.clust.6.FIC <- lnre("fzm", spc.would.clust.6.FIC, m.max = 1,
                                     bootstrap = 1000)
boot.lnre.would.clust.7.FIC <- lnre("fzm", spc.would.clust.7.FIC, m.max = 1,
                                     bootstrap = 1000)
boot.lnre.would.clust.8.FIC <- lnre("fzm", spc.would.clust.8.FIC, m.max = 1, 
                                     bootstrap = 1000)

# hapax-growth curves
extr.hgc.would.1 <- lapply(boot.lnre.would.clust.1.FIC$bootstrap,
                           function(x) lnre.vgc(x, 1:200000,
                                                m.max = 1))
for (j in seq_along(extr.hgc.would.1)){
  
  boot.list.would.clust.1.FIC[j] <- lapply(extr.hgc.would.1[[j]]$V1[200000],
                                           function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.would.2 <- lapply(boot.lnre.would.clust.2.FIC$bootstrap,
                           function(x) lnre.vgc(x, 1:200000,
                                                m.max = 1))
for (j in seq_along(extr.hgc.would.2)){
  
  boot.list.would.clust.2.FIC[j] <- lapply(extr.hgc.would.2[[j]]$V1[200000],
                                           function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.would.3 <- lapply(boot.lnre.would.clust.3.FIC$bootstrap,
                           function(x) lnre.vgc(x, 1:200000,
                                                m.max = 1))
for (j in seq_along(extr.hgc.would.3)){
  
  boot.list.would.clust.3.FIC[j] <- lapply(extr.hgc.would.3[[j]]$V1[200000],
                                           function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.would.4 <- lapply(boot.lnre.would.clust.4.FIC$bootstrap,
                           function(x) lnre.vgc(x, 1:200000,
                                                m.max = 1))
for (j in seq_along(extr.hgc.would.4)){
  
  boot.list.would.clust.4.FIC[j] <- lapply(extr.hgc.would.4[[j]]$V1[200000],
                                           function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.would.5 <- lapply(boot.lnre.would.clust.5.FIC$bootstrap,
                           function(x) lnre.vgc(x, 1:200000,
                                                m.max = 1))
for (j in seq_along(extr.hgc.would.5)){
  
  boot.list.would.clust.5.FIC[j] <- lapply(extr.hgc.would.5[[j]]$V1[200000],
                                           function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.would.6 <- lapply(boot.lnre.would.clust.6.FIC$bootstrap,
                           function(x) lnre.vgc(x, 1:200000,
                                                m.max = 1))
for (j in seq_along(extr.hgc.would.6)){
  
  boot.list.would.clust.6.FIC[j] <- lapply(extr.hgc.would.6[[j]]$V1[200000],
                                           function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.would.7 <- lapply(boot.lnre.would.clust.7.FIC$bootstrap,
                           function(x) lnre.vgc(x, 1:200000,
                                                m.max = 1))
for (j in seq_along(extr.hgc.would.7)){
  
  boot.list.would.clust.7.FIC[j] <- lapply(extr.hgc.would.7[[j]]$V1[200000],
                                           function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.would.8 <- lapply(boot.lnre.would.clust.8.FIC$bootstrap,
                           function(x) lnre.vgc(x, 1:200000,
                                                m.max = 1))
for (j in seq_along(extr.hgc.would.8)){
  
  boot.list.would.clust.8.FIC[j] <- lapply(extr.hgc.would.8[[j]]$V1[200000],
                                           function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

# repeat procedure for *will*
boot.lnre.will.clust.1.FIC <- lnre("fzm", spc.will.clust.1.FIC, m.max = 1,
                                    bootstrap = 1000)
boot.lnre.will.clust.2.FIC <- lnre("fzm", spc.will.clust.2.FIC, m.max = 1,
                                    bootstrap = 1000)
boot.lnre.will.clust.3.FIC <- lnre("fzm", spc.will.clust.3.FIC, m.max = 1,
                                    bootstrap = 1000)
boot.lnre.will.clust.4.FIC <- lnre("fzm", spc.will.clust.4.FIC, m.max = 1,
                                    bootstrap = 1000)
boot.lnre.will.clust.5.FIC <- lnre("fzm", spc.will.clust.5.FIC, m.max = 1,
                                    bootstrap = 1000)
boot.lnre.will.clust.6.FIC <- lnre("fzm", spc.will.clust.6.FIC, m.max = 1,
                                    bootstrap = 1000)
boot.lnre.will.clust.7.FIC <- lnre("fzm", spc.will.clust.7.FIC, m.max = 1,
                                    bootstrap = 1000)
boot.lnre.will.clust.8.FIC <- lnre("fzm", spc.will.clust.8.FIC, m.max = 1,
                                    bootstrap = 1000)

# hapax-growth curves
extr.hgc.will.1 <- lapply(boot.lnre.will.clust.1.FIC$bootstrap,
                          function(x) lnre.vgc(x, 1:200000,
                                               m.max = 1))
for (j in seq_along(extr.hgc.will.1)){
  
  boot.list.will.clust.1.FIC[j] <- lapply(extr.hgc.will.1[[j]]$V1[200000],
                                          function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.will.2 <- lapply(boot.lnre.will.clust.2.FIC$bootstrap,
                          function(x) lnre.vgc(x, 1:200000,
                                               m.max = 1))
for (j in seq_along(extr.hgc.will.2)){
  
  boot.list.will.clust.2.FIC[j] <- lapply(extr.hgc.will.2[[j]]$V1[200000],
                                          function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.will.3 <- lapply(boot.lnre.will.clust.3.FIC$bootstrap,
                          function(x) lnre.vgc(x, 1:200000,
                                               m.max = 1))
for (j in seq_along(extr.hgc.will.3)){
  
  boot.list.will.clust.3.FIC[j] <- lapply(extr.hgc.will.3[[j]]$V1[200000],
                                          function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.will.4 <- lapply(boot.lnre.will.clust.4.FIC$bootstrap,
                          function(x) lnre.vgc(x, 1:200000,
                                               m.max = 1))
for (j in seq_along(extr.hgc.will.4)){
  
  boot.list.will.clust.4.FIC[j] <- lapply(extr.hgc.will.4[[j]]$V1[200000],
                                          function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.will.5 <- lapply(boot.lnre.will.clust.5.FIC$bootstrap,
                          function(x) lnre.vgc(x, 1:200000,
                                               m.max = 1))
for (j in seq_along(extr.hgc.will.5)){
  
  boot.list.will.clust.5.FIC[j] <- lapply(extr.hgc.will.5[[j]]$V1[200000],
                                          function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.will.6 <- lapply(boot.lnre.will.clust.6.FIC$bootstrap,
                          function(x) lnre.vgc(x, 1:200000,
                                               m.max = 1))
for (j in seq_along(extr.hgc.will.6)){
  
  boot.list.will.clust.6.FIC[j] <- lapply(extr.hgc.will.6[[j]]$V1[200000],
                                          function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.will.7 <- lapply(boot.lnre.will.clust.7.FIC$bootstrap,
                          function(x) lnre.vgc(x, 1:200000,
                                               m.max = 1))
for (j in seq_along(extr.hgc.will.7)){
  
  boot.list.will.clust.7.FIC[j] <- lapply(extr.hgc.will.7[[j]]$V1[200000],
                                          function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

extr.hgc.will.8 <- lapply(boot.lnre.will.clust.8.FIC$bootstrap,
                          function(x) lnre.vgc(x, 1:200000,
                                               m.max = 1))
for (j in seq_along(extr.hgc.will.8)){
  
  boot.list.will.clust.8.FIC[j] <- lapply(extr.hgc.will.8[[j]]$V1[200000],
                                          function(x) x/200000)
}
rm(list = ls(pattern = "^extr"))

# plot developments in extrapolated potential productivity
png("Fig4.png", height = 10, width = 10, un = "in", res = 300)
par(mfrow = c(2, 1))
par(mar = c(4.1, 5.1, 2.1, 0.6)) 
par("cex" = 1.2)

boxplot(unlist(boot.list.d.clust.1.FIC), 
        unlist(boot.list.would.clust.1.FIC), NA, 
        unlist(boot.list.d.clust.2.FIC), 
        unlist(boot.list.would.clust.2.FIC), NA, 
        unlist(boot.list.d.clust.3.FIC), 
        unlist(boot.list.would.clust.3.FIC), NA, 
        unlist(boot.list.d.clust.4.FIC), 
        unlist(boot.list.would.clust.4.FIC), NA,
        unlist(boot.list.d.clust.5.FIC), 
        unlist(boot.list.would.clust.5.FIC), NA, 
        unlist(boot.list.d.clust.6.FIC), 
        unlist(boot.list.would.clust.6.FIC), NA, 
        unlist(boot.list.d.clust.7.FIC), 
        unlist(boot.list.would.clust.7.FIC), NA, 
        unlist(boot.list.d.clust.8.FIC), 
        unlist(boot.list.would.clust.8.FIC),
        outline = F, notch = T, 
        col = c("grey90", "grey90", NA, "grey90", "grey90", NA,
                "grey90", "grey90", NA, "grey50", "grey90", NA,
                "grey50", "grey90", NA, "grey50", "grey90", NA,
                "grey50", "grey90", NA, "grey50", "grey90"),
        border = c("grey70", 1, NA, "grey70", 1, NA, "grey70", 1, NA,
                   1, 1, NA, 1, 1, NA, 1, 1, NA, 1, 1, NA, 1, 1),
        ylim = c(0, .0125), frame.plot = F, xaxt = "n",
        ylab = NA)

axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 19.5, 22.5),
     labels = c("1830s-\n1860s", "1870s-\n1880s", 
                "1890s-\n1900s", "1910s-\n1920s",
                "1930s-\n1940s", "1950s-\n1960s", 
                "1970s-\n1980s", "1990s-\n2000s"), lwd = 0)
legend("top", legend = c("[[SUBJ'd] V]", "SUBJ would V"),
       bty = "n", pch = 22, ncol = 2, pt.bg = c("grey50", "grey90"), cex = .9)

points(1:23, c(mean(unlist(boot.list.d.clust.1.FIC)),
               mean(unlist(boot.list.would.clust.1.FIC)), NA,
               mean(unlist(boot.list.d.clust.2.FIC)),
               mean(unlist(boot.list.would.clust.2.FIC)), NA,
               mean(unlist(boot.list.d.clust.3.FIC)), 
               mean(unlist(boot.list.would.clust.3.FIC)), NA,
               mean(unlist(boot.list.d.clust.4.FIC)),
               mean(unlist(boot.list.would.clust.4.FIC)), NA,
               mean(unlist(boot.list.d.clust.5.FIC)),
               mean(unlist(boot.list.would.clust.5.FIC)), NA,
               mean(unlist(boot.list.d.clust.6.FIC)),
               mean(unlist(boot.list.would.clust.6.FIC)), NA,
               mean(unlist(boot.list.d.clust.7.FIC)),
               mean(unlist(boot.list.would.clust.7.FIC)), NA,
               mean(unlist(boot.list.d.clust.8.FIC)),
               mean(unlist(boot.list.would.clust.8.FIC))),
       col = c("grey70", 1, NA, "grey70", 1, NA, "grey70", 1, NA,
               0, 1, NA, 0, 1, NA, 0, 1, NA, 0, 1, NA, 0, 1),  pch = 3)

boxplot(unlist(boot.list.ll.clust.1.FIC), 
        unlist(boot.list.will.clust.1.FIC), NA, 
        unlist(boot.list.ll.clust.2.FIC), 
        unlist(boot.list.will.clust.2.FIC), NA, 
        unlist(boot.list.ll.clust.3.FIC), 
        unlist(boot.list.will.clust.3.FIC), NA, 
        unlist(boot.list.ll.clust.4.FIC), 
        unlist(boot.list.will.clust.4.FIC), NA,
        unlist(boot.list.ll.clust.5.FIC), 
        unlist(boot.list.will.clust.5.FIC), NA, 
        unlist(boot.list.ll.clust.6.FIC), 
        unlist(boot.list.will.clust.6.FIC), NA, 
        unlist(boot.list.ll.clust.7.FIC), 
        unlist(boot.list.will.clust.7.FIC), NA, 
        unlist(boot.list.ll.clust.8.FIC), 
        unlist(boot.list.will.clust.8.FIC),
        outline = F, notch = T, 
        col = c("grey30", "grey90", NA, "grey30", "grey90", NA,
                "grey30", "grey90", NA, "grey30", "grey90", NA,
                "grey30", "grey90", NA, "grey30", "grey90", NA,
                "grey30", "grey90", NA, "grey30", "grey90"),
        ylim = c(0,.0125), frame.plot = F, xaxt = "n",
        ylab = NA)

axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 19.5, 22.5),
     labels = c("1830s-\n1860s", "1870s-\n1880s", 
                "1890s-\n1900s", "1910s-\n1920s",
                "1930s-\n1940s", "1950s-\n1960s", 
                "1970s-\n1980s", "1990s-\n2000s"), lwd = 0)
legend("top", legend = c("[[SUBJ'll] V]", "SUBJ will V"),
       bty = "n", pch = 22, ncol = 2, pt.bg = c("grey30", "grey90"), cex = .9)

points(1:23, c(mean(unlist(boot.list.ll.clust.1.FIC)),
               mean(unlist(boot.list.will.clust.1.FIC)), NA,
               mean(unlist(boot.list.ll.clust.2.FIC)),
               mean(unlist(boot.list.will.clust.2.FIC)), NA,
               mean(unlist(boot.list.ll.clust.3.FIC)), 
               mean(unlist(boot.list.will.clust.3.FIC)), NA,
               mean(unlist(boot.list.ll.clust.4.FIC)),
               mean(unlist(boot.list.will.clust.4.FIC)), NA,
               mean(unlist(boot.list.ll.clust.5.FIC)),
               mean(unlist(boot.list.will.clust.5.FIC)), NA,
               mean(unlist(boot.list.ll.clust.6.FIC)),
               mean(unlist(boot.list.will.clust.6.FIC)), NA,
               mean(unlist(boot.list.ll.clust.7.FIC)),
               mean(unlist(boot.list.will.clust.7.FIC)), NA,
               mean(unlist(boot.list.ll.clust.8.FIC)),
               mean(unlist(boot.list.will.clust.8.FIC))),
       col = c(0, 1, NA, 0, 1, NA, 0, 1, NA, 0, 1, NA, 0, 1, NA, 
               0, 1, NA, 0, 1, NA, 0, 1), pch = 3)

mtext(expression("P"["extr"]*" (N"["extr"]*"=200,000)"), side = 2, outer = TRUE,  
      cex = 1.3, line = -2.2, adj = 0.53)
dev.off()

# assess trends
median.PP.d <- c(median(unlist(boot.list.d.clust.1.FIC)),
                 median(unlist(boot.list.d.clust.2.FIC)),
                 median(unlist(boot.list.d.clust.3.FIC)),
                 median(unlist(boot.list.d.clust.4.FIC)),
                 median(unlist(boot.list.d.clust.5.FIC)),
                 median(unlist(boot.list.d.clust.6.FIC)),
                 median(unlist(boot.list.d.clust.7.FIC)),
                 median(unlist(boot.list.d.clust.8.FIC)))

median.PP.ll <- c(median(unlist(boot.list.ll.clust.1.FIC)),
                  median(unlist(boot.list.ll.clust.2.FIC)),
                  median(unlist(boot.list.ll.clust.3.FIC)),
                  median(unlist(boot.list.ll.clust.4.FIC)),
                  median(unlist(boot.list.ll.clust.5.FIC)),
                  median(unlist(boot.list.ll.clust.6.FIC)),
                  median(unlist(boot.list.ll.clust.7.FIC)),
                  median(unlist(boot.list.ll.clust.8.FIC)))

median.PP.would <- c(median(unlist(boot.list.would.clust.1.FIC)),
                     median(unlist(boot.list.would.clust.2.FIC)),
                     median(unlist(boot.list.would.clust.3.FIC)),
                     median(unlist(boot.list.would.clust.4.FIC)),
                     median(unlist(boot.list.would.clust.5.FIC)),
                     median(unlist(boot.list.would.clust.6.FIC)),
                     median(unlist(boot.list.would.clust.7.FIC)),
                     median(unlist(boot.list.would.clust.8.FIC)))

median.PP.will <- c(median(unlist(boot.list.will.clust.1.FIC)),
                    median(unlist(boot.list.will.clust.2.FIC)),
                    median(unlist(boot.list.will.clust.3.FIC)),
                    median(unlist(boot.list.will.clust.4.FIC)),
                    median(unlist(boot.list.will.clust.5.FIC)),
                    median(unlist(boot.list.will.clust.6.FIC)),
                    median(unlist(boot.list.will.clust.7.FIC)),
                    median(unlist(boot.list.will.clust.8.FIC)))

cor.test(median.PP.d, 1:8, method = "kendall")
cor.test(median.PP.ll, 1:8, method = "kendall")
cor.test(median.PP.would, 1:8, method = "kendall")
cor.test(median.PP.will, 1:8, method = "kendall")

#============================================================================================

## COLLOSTRUCTIONAL ANALYSIS (CA)

# merge data for enclitics and full forms
input.d.would <- merge.data.table(input.d.FIC, input.would.FIC, all = T)
input.ll.will <- merge.data.table(input.ll.FIC, input.will.FIC, all = T)

# create three equidistant time periods (EDP); select most common pronouns to reduce noise; 
# run CA
output.ca.d.would.3EDP <- input.d.would %>%
  mutate(VERB = ifelse(VERB %in% "'ve",
                       paste("have"),
                       as.character(VERB))) %>%
  mutate(TIME_period = ifelse(DECADE %in% c("1830s",
                                            "1840s",
                                            "1850s",
                                            "1860s",
                                            "1870s",
                                            "1880s"),
                              "1830s-1880s",
                              ifelse(DECADE %in% c("1890s",
                                                   "1900s",
                                                   "1910s",
                                                   "1920s",
                                                   "1930s",
                                                   "1940s"),
                                     "1890s-1940s","1950s-2000s"))) %>%
  filter(SUBJ %in% c("I", "you", "he", "she", "it",
                     "we", "they", "there", "that", 
                     "what", "who", "this", "which")) %>%
  droplevels() %>%
  select(4:7) %>%
  collex.covar.mult() # configuration: CXN x SUBJ x VERB X TIME, coll.strength: t-score

output.ca.ll.will.3EDP <- input.ll.will %>%
  mutate(TIME_period = ifelse(DECADE %in% c("1830s",
                                            "1840s",
                                            "1850s",
                                            "1860s",
                                            "1870s",
                                            "1880s"),
                              "1830s-1880s",
                              ifelse(DECADE %in% c("1890s",
                                                   "1900s",
                                                   "1910s",
                                                   "1920s",
                                                   "1930s",
                                                   "1940s"),
                                     "1890s-1940s","1950s-2000s"))) %>%
  filter(SUBJ %in% c("I", "you", "he", "she", "it",
                     "we", "they", "there", "that", 
                     "what", "who", "this", "which")) %>%
  droplevels() %>%
  select(4:7) %>%
  collex.covar.mult()

# inspect data
View(output.ca.d.would.3EDP)
View(output.ca.ll.will.3EDP)

# prepare data for Cleveland dot plots; filtered for SUBJ and periods; ordered by descending 
# t-score; for *'d* and *would*
p1.I.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "I" & CXN == "'d" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.I.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "I" & CXN == "would" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.I.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "I" & CXN == "'d" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.I.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "I" & CXN == "would" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.I.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "I" & CXN == "'d" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.I.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "I" & CXN == "would" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.we.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "we" & CXN == "'d" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.we.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "we" & CXN == "would" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.we.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "we" & CXN == "'d" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.we.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "we" & CXN == "would" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.we.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "we" & CXN == "'d" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.we.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "we" & CXN == "would" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.you.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "you" & CXN == "'d" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.you.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "you" & CXN == "would" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.you.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "you" & CXN == "'d" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.you.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "you" & CXN == "would" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.you.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "you" & CXN == "'d" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.you.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "you" & CXN == "would" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.they.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "they" & CXN == "'d" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.they.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "they" & CXN == "would" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.they.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "they" & CXN == "'d" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.they.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "they" & CXN == "would" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.they.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "they" & CXN == "'d" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.they.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "they" & CXN == "would" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.it.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "it" & CXN == "'d" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.it.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "it" & CXN == "would" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.it.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "it" & CXN == "'d" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.it.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "it" & CXN == "would" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.it.d <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "it" & CXN == "'d" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.it.would <- output.ca.d.would.3EDP %>%
  filter(SUBJ == "it" & CXN == "would" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))

# for *'ll* and *will*
p1.I.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "I" & CXN == "'ll" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.I.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "I" & CXN == "will" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.I.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "I" & CXN == "'ll" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.I.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "I" & CXN == "will" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.I.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "I" & CXN == "'ll" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.I.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "I" & CXN == "will" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.we.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "we" & CXN == "'ll" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.we.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "we" & CXN == "will" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.we.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "we" & CXN == "'ll" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.we.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "we" & CXN == "will" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.we.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "we" & CXN == "'ll" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.we.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "we" & CXN == "will" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.you.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "you" & CXN == "'ll" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.you.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "you" & CXN == "will" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.you.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "you" & CXN == "'ll" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.you.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "you" & CXN == "will" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.you.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "you" & CXN == "'ll" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.you.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "you" & CXN == "will" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.they.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "they" & CXN == "'ll" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.they.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "they" & CXN == "will" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.they.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "they" & CXN == "'ll" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.they.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "they" & CXN == "will" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.they.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "they" & CXN == "'ll" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.they.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "they" & CXN == "will" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.it.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "it" & CXN == "'ll" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p1.it.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "it" & CXN == "will" & TIME_period == "1830s-1880s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.it.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "it" & CXN == "'ll" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p2.it.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "it" & CXN == "will" & TIME_period == "1890s-1940s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.it.ll <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "it" & CXN == "'ll" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))
p3.it.will <- output.ca.ll.will.3EDP %>%
  filter(SUBJ == "it" & CXN == "will" & TIME_period == "1950s-2000s") %>%
  droplevels() %>%
  arrange(desc(T))

# plot results

# Note: For earlier periods, it seemed sensible to add a frequency threshold in some cases,
# otherwise there would have been some obscure combinations. 

# Also, for some lines, you'll get "Error in xy.coords(x, y) : 'x' and 'y' lengths differ". 
# This can be ignored, as the combinations in question are genuinely not in the data set and
# can therefore not produce any results.   
png("Fig5.png", height = 14, width = 19, un = "in", res = 300)
par(mar = c(3.1, 4.1, 2.1, .6))
par(mfrow = c(3, 5))
par("cex" = 1.3)

# Since R 4.x, the first dotchart always comes with graphical errors (e.g. missing symbols,
# incomplete abline) (bug???); the following two lines take care of that
plot.new()
par(new=T)

# *'d/would*
# period 1, 1830s-1880s
dotchart(rev(p1.I.d[p1.I.d$OBS >= 10,][1:15,]$T), 
         labels = rev(p1.I.d[p1.I.d$OBS >= 10,][1:15,]$VERB), 
         main = "I'd V vs I would V", pch = 16, xlab = NA, offset = .45,
         ylab = "1830s-1880s", cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p1.I.would[p1.I.would$VERB == "like",]$T), 15 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "give",]$T), 14 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "go",]$T), 13 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "try",]$T), 12 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "let",]$T), 11 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "swear",]$T), 10 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "hang",]$T), 9 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "thank",]$T), 8 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "die",]$T), 7 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "sell",]$T), 6 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "advise",]$T), 5 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "trust",]$T), 4 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "marry",]$T), 3 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "put",]$T), 2 , pch = 2)
points(rev(p1.I.would[p1.I.would$VERB == "stay",]$T), 1 , pch = 2)

dotchart(rev(p1.we.d[p1.we.d$OBS >= 4,][1:15,]$T), 
         labels = rev(p1.we.d[p1.we.d$OBS >= 4,][1:15,]$VERB), 
         main = "we'd V vs we would V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p1.we.would[p1.we.would$VERB == "meet",]$T),15 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "put",]$T), 14 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "come",]$T), 13 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "try",]$T), 12 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "let",]$T), 11 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "see",]$T), 10 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "get",]$T), 9 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "go",]$T), 8 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "find",]$T), 7 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "give",]$T), 6 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "make",]$T), 5 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "take",]$T), 4 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "like",]$T), 3 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "do",]$T), 2 , pch = 2)
points(rev(p1.we.would[p1.we.would$VERB == "have",]$T), 1 , pch = 2)

dotchart(rev(p1.you.d[p1.you.d$OBS >= 9,][1:15,]$T), 
         labels = rev(p1.you.d[p1.you.d$OBS >= 9,][1:15,]$VERB), 
         main = "you'd V vs you would V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p1.you.would[p1.you.would$VERB == "like",]$T), 15 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "think",]$T), 14 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "let",]$T), 13 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "want",]$T), 12 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "see",]$T), 11 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "find",]$T), 10 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "see",]$T), 9 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "wish",]$T), 8 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "call",]$T), 7 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "make",]$T), 6 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "come",]$T), 5 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "get",]$T), 4 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "run",]$T), 3 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "know",]$T), 2 , pch = 2)
points(rev(p1.you.would[p1.you.would$VERB == "tell",]$T), 1 , pch = 2)

dotchart(rev(p1.they.d[p1.they.d$OBS >= 7,][1:15,]$T), 
         labels = rev(p1.they.d[p1.they.d$OBS >= 7,][1:15,]$VERB), 
         main = "they'd V vs they would V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col= "darkgrey")
points(rev(p1.they.would[p1.they.would$VERB == "send",]$T), 15 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "let",]$T), 14 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "want",]$T), 13 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "find",]$T), 12 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "think",]$T), 11 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "know",]$T), 10 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "take",]$T), 9 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "get",]$T), 8 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "make",]$T), 7 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "give",]$T), 6 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "do",]$T), 5 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "come",]$T), 4 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "go",]$T), 3 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "say",]$T), 2 , pch = 2)
points(rev(p1.they.would[p1.they.would$VERB == "like",]$T), 1 , pch = 2)

dotchart(rev(p1.it.would[1:15,]$T), 
         labels = rev(p1.it.would[1:15,]$VERB), 
         main = "it'd V vs it would V", pch = 2, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(-75,15 , pch = 16) # actual score: -186.28; manually set to -75 > avoids long tail
points(rev(p1.it.d[p1.it.d$VERB == "seem",]$T), 14 , pch = 16) 
points(rev(p1.it.d[p1.it.d$VERB == "have",]$T), 13 , pch = 16) 
points(rev(p1.it.d[p1.it.d$VERB == "make",]$T), 12 , pch = 16)
points(rev(p1.it.d[p1.it.d$VERB == "appear",]$T), 11 , pch = 16)
points(rev(p1.it.d[p1.it.d$VERB == "take",]$T), 10 , pch = 16)
points(rev(p1.it.d[p1.it.d$VERB == "require",]$T), 9 , pch = 16)
points(rev(p1.it.d[p1.it.d$VERB == "break",]$T), 8 , pch = 16)
points(rev(p1.it.d[p1.it.d$VERB == "kill",]$T), 7 , pch = 16)
points(rev(p1.it.d[p1.it.d$VERB == "afford",]$T), 6 , pch = 16)
points(rev(p1.it.d[p1.it.d$VERB == "please",]$T), 5 , pch = 16)
points(rev(p1.it.d[p1.it.d$VERB == "cost",]$T), 4 , pch = 16)
points(rev(p1.it.d[p1.it.d$VERB == "do",]$T), 3 , pch = 16)
points(rev(p1.it.d[p1.it.d$VERB == "burst",]$T), 2 , pch = 16)
points(rev(p1.it.d[p1.it.d$VERB == "save",]$T), 1 , pch = 16)

# period 2, 1890s-1940s
dotchart(rev(p2.I.d[1:15,]$T), 
         labels = rev(p2.I.d[1:15,]$VERB), offset = .5, 
         main = "I'd V vs I would V", pch = 16, xlab = NA, 
         ylab = "1890s-1940s", cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p2.I.would[p2.I.would$VERB == "like",]$T),15 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "give",]$T),14 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "hate",]$T),13 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "love",]$T),12 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "say",]$T),11 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "go",]$T),10 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "do",]$T),9 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "have",]$T),8 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "get",]$T),7 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "feel",]$T),6 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "take",]$T),5 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "let",]$T),4 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "tell",]$T),3 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "know",]$T),2 , pch = 2)
points(rev(p2.I.would[p2.I.would$VERB == "put",]$T),1 , pch = 2)

dotchart(rev(p2.we.d[p2.we.d$OBS >= 9,][1:15,]$T), 
         labels = rev(p2.we.d[p2.we.d$OBS >= 9,][1:15,]$VERB), 
         main = "we'd V vs we would V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p2.we.would[p2.we.would$VERB == "have",]$T),15 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "get",]$T),14 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "go",]$T),13 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "like",]$T),12 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "be",]$T),11 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "find",]$T),10 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "come",]$T),9 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "love",]$T),8 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "do",]$T),7 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "see",]$T),6 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "drop",]$T),5 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "run",]$T),4 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "lose",]$T),3 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "need",]$T),2 , pch = 2)
points(rev(p2.we.would[p2.we.would$VERB == "hear",]$T),1 , pch = 2)

dotchart(rev(p2.you.d[1:15,]$T), 
         labels = rev(p2.you.d[1:15,]$VERB), 
         main = "you'd V vs you would V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p2.you.would[p2.you.would$VERB == "think",]$T),15 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "like",]$T),14 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "be",]$T),13 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "come",]$T),12 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "let",]$T),11 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "get",]$T),10 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "call",]$T),9 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "know",]$T),8 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "want",]$T),7 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "find",]$T),6 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "do",]$T),5 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "tell",]$T),4 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "have",]$T),3 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "see",]$T),2 , pch = 2)
points(rev(p2.you.would[p2.you.would$VERB == "go",]$T),1 , pch = 2)

dotchart(rev(p2.they.d[p2.they.d$OBS >= 9,][1:15,]$T), 
         labels = rev(p2.they.d[p2.they.d$OBS >= 9,][1:15,]$VERB), 
         main = "they'd V vs they would V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p2.they.would[p2.they.would$VERB == "get",]$T),15 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "put",]$T),14 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "let",]$T),13 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "think",]$T),12 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "know",]$T),11 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "come",]$T),10 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "kill",]$T),9 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "laugh",]$T),8 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "shoot",]$T),7 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "throw",]$T),6 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "want",]$T),5 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "find",]$T),4 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "start",]$T),3 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "stay",]$T),2 , pch = 2)
points(rev(p2.they.would[p2.they.would$VERB == "pay",]$T),1 , pch = 2)

dotchart(rev(p2.it.would[1:15,]$T), 
         labels = rev(p2.it.would[1:15,]$VERB), 
         main = "it'd V vs  it would V", pch = 2, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p2.it.d[p2.it.d$VERB == "be",]$T),15 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "seem",]$T),14 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "mean",]$T),13 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "take",]$T),12 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "make",]$T),11 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "appear",]$T),10 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "break",]$T),9 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "cost",]$T),8 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "have",]$T),7 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "help",]$T),6 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "hurt",]$T),5 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "serve",]$T),4 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "require",]$T),3 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "please",]$T),2 , pch = 16)
points(rev(p2.it.d[p2.it.d$VERB == "spoil",]$T),1 , pch = 16)

# period 3, 1950s-2000s
dotchart(rev(p3.I.d[1:15,]$T), 
         labels = rev(p3.I.d[1:15,]$VERB), offset = .35,
         main = "I'd V vs I would V", pch = 16, xlab = NA, 
         ylab = "1950s-2000s", cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p3.I.would[p3.I.would$VERB == "like",]$T),15 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "say",]$T),14 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "love",]$T),13 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "get",]$T),12 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "hate",]$T),11 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "go",]$T),10 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "appreciate",]$T),9 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "do",]$T),8 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "tell",]$T),7 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "feel",]$T),6 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "come",]$T),5 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "guess",]$T),4 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "give",]$T),3 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "have",]$T),2 , pch = 2)
points(rev(p3.I.would[p3.I.would$VERB == "call",]$T),1 , pch = 2)

dotchart(rev(p3.we.d[1:15,]$T), 
         labels = rev(p3.we.d[1:15,]$VERB), 
         main = "we'd V vs we would V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p3.we.would[p3.we.would$VERB == "have",]$T),15 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "be",]$T),14 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "get",]$T),13 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "like",]$T),12 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "go",]$T),11 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "need",]$T),10 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "sit",]$T),9 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "love",]$T),8 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "make",]$T),7 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "come",]$T),6 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "talk",]$T),5 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "eat",]$T),4 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "do",]$T),3 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "know",]$T),2 , pch = 2)
points(rev(p3.we.would[p3.we.would$VERB == "see",]$T),1 , pch = 2)

dotchart(rev(p3.you.d[1:15,]$T), 
         labels = rev(p3.you.d[1:15,]$VERB), 
         main = "you'd V vs you would V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p3.you.would[p3.you.would$VERB == "think",]$T),15 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "like",]$T),14 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "be",]$T),13 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "want",]$T),12 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "call",]$T),11 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "get",]$T),10 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "know",]$T),9 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "expect",]$T),8 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "do",]$T),7 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "see",]$T),6 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "come",]$T),5 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "care",]$T),4 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "say",]$T),3 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "let",]$T),2 , pch = 2)
points(rev(p3.you.would[p3.you.would$VERB == "have",]$T),1 , pch = 2)

dotchart(rev(p3.they.d[1:15,]$T), 
         labels = rev(p3.they.d[1:15,]$VERB), 
         main = "they'd V vs they would V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p3.they.would[p3.they.would$VERB == "come",]$T),15 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "put",]$T),14 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "get",]$T),13 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "let",]$T),12 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "find",]$T),11 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "know",]$T),10 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "run",]$T),9 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "want",]$T),8 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "kill",]$T),7 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "need",]$T),6 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "go",]$T),5 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "eat",]$T),4 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "talk",]$T),3 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "become",]$T),2 , pch = 2)
points(rev(p3.they.would[p3.they.would$VERB == "pick",]$T),1 , pch = 2)

dotchart(rev(p3.it.would[1:15,]$T), 
         labels = rev(p3.it.would[1:15,]$VERB), 
         main = "it'd V vs it would V", pch = 2, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p3.it.d[p3.it.d$VERB == "be",]$T),15 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "take",]$T),14 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "mean",]$T),13 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "seem",]$T),12 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "help",]$T),11 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "cost",]$T),10 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "hurt",]$T),9 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "work",]$T),8 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "serve",]$T),7 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "end",]$T),6 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "appear",]$T),5 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "make",]$T),4 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "happen",]$T),3 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "sound",]$T),2 , pch = 16)
points(rev(p3.it.d[p3.it.d$VERB == "please",]$T),1 , pch = 16)

mtext("Collostructional strength (t-score)", side = 1, outer = TRUE,
      cex = 1.4, line = -1, adj = 0.53)
dev.off()

# *'ll/will*
png("Fig6.png", height = 14, width = 19, un = "in", res = 300)
par(mar = c(3.1, 4.1, 2.1, .6))
par(mfrow = c(3, 5))
par("cex" = 1.3)
plot.new()
par(new=T)

# period 1, 1830s-1880s 
dotchart(rev(p1.I.ll[1:15,]$T), offset = .35,
         labels = rev(p1.I.ll[1:15,]$VERB),
         main = "I'll V vs I will V", pch = 16, xlab = NA, 
         ylab = "1830s-1880s", cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p1.I.will[p1.I.will$VERB == "tell",]$T),15 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "warrant",]$T),14 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "swear",]$T),13 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "give",]$T),12 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "go",]$T),11 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "try",]$T),10 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "teach",]$T),9 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "engage",]$T),8 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "bet",]$T),7 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "risk",]$T),6 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "lick",]$T),5 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "show",]$T),4 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "wager",]$T),3 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "fix",]$T),2 , pch = 2)
points(rev(p1.I.will[p1.I.will$VERB == "bear",]$T),1 , pch = 2)

dotchart(rev(p1.we.ll[p1.we.ll$OBS >= 9,][1:15,]$T), 
         labels = rev(p1.we.ll[p1.we.ll$OBS >= 9,][1:15,]$VERB), 
         main = "we'll V vs we will V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p1.we.will[p1.we.will$VERB == "see",]$T),15 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "talk",]$T),14 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "settle",]$T),13 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "seek",]$T),12 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "meet",]$T),11 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "hunt",]$T),10 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "drink",]$T),9 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "have",]$T),8 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "sing",]$T),7 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "fight",]$T),6 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "fix",]$T),5 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "begin",]$T),4 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "set",]$T),3 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "pull",]$T),2 , pch = 2)
points(rev(p1.we.will[p1.we.will$VERB == "join",]$T),1 , pch = 2)

dotchart(rev(p1.you.ll[p1.you.ll$OBS >= 9,][1:15,]$T), 
         labels = rev(p1.you.ll[p1.you.ll$OBS >= 9,][1:15,]$VERB), 
         main = "you'll V vs you will V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p1.you.will[p1.you.will$VERB == "find",]$T),15 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "excuse",]$T),14 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "catch",]$T),13 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "see",]$T),12 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "like",]$T),11 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "allow",]$T),10 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "oblige",]$T),9 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "promise",]$T),8 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "feel",]$T),7 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "believe",]$T),6 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "spoil",]$T),5 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "forgive",]$T),4 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "lose",]$T),3 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "learn",]$T),2 , pch = 2)
points(rev(p1.you.will[p1.you.will$VERB == "hear",]$T),1 , pch = 2)

dotchart(rev(p1.they.ll[p1.they.ll$OBS >= 9,][1:15,]$T), 
         labels = rev(p1.they.ll[p1.they.ll$OBS >= 9,][1:15,]$VERB), 
         main = "they'll V vs they ll V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p1.they.will[p1.they.will$VERB == "hang",]$T),15 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "want",]$T),14 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "think",]$T),13 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "kill",]$T),12 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "hear",]$T),11 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "let",]$T),10 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "keep",]$T),9 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "say",]$T),8 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "send",]$T),7 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "know",]$T),6 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "come",]$T),5 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "find",]$T),4 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "do",]$T),3 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "make",]$T),2 , pch = 2)
points(rev(p1.they.will[p1.they.will$VERB == "give",]$T),1 , pch = 2)

dotchart(rev(p1.it.will[1:15,]$T), offset = .15,
         labels = rev(p1.it.will[1:15,]$VERB), 
         main = "it'll V vs it will V", pch = 2, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p1.it.ll[p1.it.ll$VERB == "be",]$T),15 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "do",]$T),14 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "make",]$T),13 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "take",]$T),12 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "cost",]$T),11 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "come",]$T),10 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "seem",]$T),9 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "serve",]$T),8 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "save",]$T),7 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "bring",]$T),6 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "prove",]$T),5 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "afford",]$T),4 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "require",]$T),3 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "end",]$T),2 , pch = 16)
points(rev(p1.it.ll[p1.it.ll$VERB == "please",]$T),1 , pch = 16)

# period 2, 1890s-1940s
dotchart(rev(p2.I.ll[1:15,]$T), 
         labels = rev(p2.I.ll[1:15,]$VERB), offset = .5,
         main = "I'll V vs I will V", pch = 16, xlab = NA, 
         ylab = "1890s-1940s", cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p2.I.will[p2.I.will$VERB == "tell",]$T),15 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "go",]$T),14 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "bet",]$T),13 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "take",]$T),12 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "show",]$T),11 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "get",]$T),10 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "give",]$T),9 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "try",]$T),8 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "send",]$T),7 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "do",]$T),6 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "put",]$T),5 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "fix",]$T),4 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "say",]$T),3 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "ask",]$T),2 , pch = 2)
points(rev(p2.I.will[p2.I.will$VERB == "wait",]$T),1 , pch = 2)

dotchart(rev(p2.we.ll[1:15,]$T), 
         labels = rev(p2.we.ll[1:15,]$VERB), 
         main = "we'll V vs we will V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p2.we.will[p2.we.will$VERB == "have",]$T),15 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "go",]$T),14 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "get",]$T),13 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "see",]$T),12 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "talk",]$T),11 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "start",]$T),10 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "make",]$T),9 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "try",]$T),8 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "take",]$T),7 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "put",]$T),6 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "fix",]$T),5 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "wait",]$T),4 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "let",]$T),3 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "eat",]$T),2 , pch = 2)
points(rev(p2.we.will[p2.we.will$VERB == "begin",]$T),1 , pch = 2)

dotchart(rev(p2.you.ll[1:15,]$T), 
         labels = rev(p2.you.ll[1:15,]$VERB), 
         main = "you'll V vs you will V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p2.you.will[p2.you.will$VERB == "have",]$T),15 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "find",]$T),14 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "get",]$T),13 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "see",]$T),12 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "excuse",]$T),11 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "like",]$T),10 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "know",]$T),9 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "feel",]$T),8 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "want",]$T),7 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "hear",]$T),6 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "need",]$T),5 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "let",]$T),4 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "pardon",]$T),3 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "forgive",]$T),2 , pch = 2)
points(rev(p2.you.will[p2.you.will$VERB == "understand",]$T),1 , pch = 2)

dotchart(rev(p2.they.ll[1:15,]$T), 
         labels = rev(p2.they.ll[1:15,]$VERB), 
         main = "they'll V vs they will V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p2.they.will[p2.they.will$VERB == "think",]$T),15 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "kill",]$T),14 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "get",]$T),13 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "hang",]$T),12 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "fit",]$T),11 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "know",]$T),10 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "hate",]$T),9 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "steal",]$T),8 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "want",]$T),7 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "let",]$T),6 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "believe",]$T),5 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "come",]$T),4 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "burn",]$T),3 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "rush",]$T),2 , pch = 2)
points(rev(p2.they.will[p2.they.will$VERB == "be",]$T),1 , pch = 2)

dotchart(rev(p2.it.will[1:15,]$T), 
         labels = rev(p2.it.will[1:15,]$VERB), 
         main = "it'll V vs it will V", pch = 2, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p2.it.ll[p2.it.ll$VERB == "be",]$T),15 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "make",]$T),14 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "take",]$T),13 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "mean",]$T),12 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "do",]$T),11 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "help",]$T),10 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "cost",]$T),9 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "save",]$T),8 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "seem",]$T),7 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "please",]$T),6 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "give",]$T),5 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "come",]$T),4 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "hurt",]$T),3 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "end",]$T),2 , pch = 16)
points(rev(p2.it.ll[p2.it.ll$VERB == "serve",]$T),1 , pch = 16)

# period 3, 1950s-2000s
dotchart(rev(p3.I.ll[1:15,]$T), 
         labels = rev(p3.I.ll[1:15,]$VERB), offset = .5,
         main = "I'll V vs I will V", pch = 16, xlab = NA, 
         ylab = "1950s-2000s", cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p3.I.will[p3.I.will$VERB == "tell",]$T),15 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "bet",]$T),14 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "get",]$T),13 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "call",]$T),12 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "take",]$T),11 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "try",]$T),10 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "show",]$T),9 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "give",]$T),8 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "go",]$T),7 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "see",]$T),6 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "buy",]$T),5 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "kill",]$T),4 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "put",]$T),3 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "let",]$T),2 , pch = 2)
points(rev(p3.I.will[p3.I.will$VERB == "wait",]$T),1 , pch = 2)

dotchart(rev(p3.we.ll[1:15,]$T), 
         labels = rev(p3.we.ll[1:15,]$VERB), 
         main = "we'll V vs we will V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p3.we.will[p3.we.will$VERB == "have",]$T),15 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "get",]$T),14 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "see",]$T),13 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "talk",]$T),12 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "go",]$T),11 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "need",]$T),10 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "start",]$T),9 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "find",]$T),8 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "take",]$T),7 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "meet",]$T),6 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "put",]$T),5 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "work",]$T),4 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "figure",]$T),3 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "call",]$T),2 , pch = 2)
points(rev(p3.we.will[p3.we.will$VERB == "do",]$T),1 , pch = 2)

dotchart(rev(p3.you.ll[1:15,]$T), 
         labels = rev(p3.you.ll[1:15,]$VERB), 
         main = "you'll V vs you will V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p3.you.will[p3.you.will$VERB == "have",]$T),15 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "see",]$T),14 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "get",]$T),13 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "find",]$T),12 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "need",]$T),11 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "like",]$T),10 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "excuse",]$T),9 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "want",]$T),8 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "feel",]$T),7 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "know",]$T),6 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "learn",]$T),5 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "be",]$T),4 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "notice",]$T),3 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "understand",]$T),2 , pch = 2)
points(rev(p3.you.will[p3.you.will$VERB == "pardon",]$T),1 , pch = 2)

dotchart(rev(p3.they.ll[1:15,]$T), 
         labels = rev(p3.they.ll[1:15,]$VERB), 
         main = "they'll V vs they will V", pch = 16, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p3.they.will[p3.they.will$VERB == "be",]$T),15 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "kill",]$T),14 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "think",]$T),13 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "know",]$T),12 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "come",]$T),11 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "want",]$T),10 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "have",]$T),9 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "get",]$T),8 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "start",]$T),7 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "believe",]$T),6 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "use",]$T),5 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "eat",]$T),4 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "stop",]$T),3 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "put",]$T),2 , pch = 2)
points(rev(p3.they.will[p3.they.will$VERB == "say",]$T),1 , pch = 2)

dotchart(rev(p3.it.will[1:15,]$T), 
         labels = rev(p3.it.will[1:15,]$VERB), 
         main = "it'll V vs it will V", pch = 2, xlab = NA, 
         cex = 1.4, cex.main = 1, xlim = c(-75,75))
abline(v = 0, lty = 2, col = "darkgrey")
points(rev(p3.it.ll[p3.it.ll$VERB == "be",]$T),15 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "take",]$T),14 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "work",]$T),13 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "cost",]$T),12 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "help",]$T),11 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "happen",]$T),10 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "mean",]$T),9 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "make",]$T),8 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "pass",]$T),7 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "end",]$T),6 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "hurt",]$T),5 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "serve",]$T),4 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "seem",]$T),3 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "last",]$T),2 , pch = 16)
points(rev(p3.it.ll[p3.it.ll$VERB == "heal",]$T),1 , pch = 16)

mtext("Collostructional strength (t-score)", side = 1, outer = TRUE,  
      cex = 1.4, line = -1, adj = 0.53)
dev.off()

#============================================================================================ 

## END

#============================================================================================ 