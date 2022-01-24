#============================================================================================
# Update on script "enclitics_R_script.R"
# Jan 23, 2022

# fZM-modeling with random number seeds

# The original script did not use random number seeds for the fZM-models - certainly, not the
# best practice. The current script now includes them and, although the numbers are not 
# identical down to the fourth or fifth decimal (which is doubtful to be linguistically
# relevant anyway), it produces the exact same trends as reported in the published paper. The
# overall findings thus appear to be quite robust.

#============================================================================================
## CHANGES IN THE PRODUCTIVITY IN THE VERB SLOT (FZM-MODELING WITH RANDOM NUMBER SEEDS)

# run "enclitics_R_script.R" first up to line 912

# *'d*; fit LNRE models to 1,000 parametric bootstrap samples for each cluster
set.seed(1934) # ensures that results are reproducible
boot.lnre.d.clust.1.FIC <- lnre("fzm", spc.d.clust.1.FIC, m.max = 1, bootstrap = 1000)
set.seed(5487)
boot.lnre.d.clust.2.FIC <- lnre("fzm", spc.d.clust.2.FIC, m.max = 1, bootstrap = 1000)
set.seed(2341)
boot.lnre.d.clust.3.FIC <- lnre("fzm", spc.d.clust.3.FIC, m.max = 1, bootstrap = 1000)
set.seed(7809)
boot.lnre.d.clust.4.FIC <- lnre("fzm", spc.d.clust.4.FIC, m.max = 1, bootstrap = 1000)
set.seed(9998)
boot.lnre.d.clust.5.FIC <- lnre("fzm", spc.d.clust.5.FIC, m.max = 1, bootstrap = 1000)
set.seed(4973)
boot.lnre.d.clust.6.FIC <- lnre("fzm", spc.d.clust.6.FIC, m.max = 1, bootstrap = 1000)
set.seed(1331)
boot.lnre.d.clust.7.FIC <- lnre("fzm", spc.d.clust.7.FIC, m.max = 1, bootstrap = 1000)
set.seed(1029)
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
set.seed(3203)
boot.lnre.ll.clust.1.FIC <- lnre("fzm", spc.ll.clust.1.FIC, m.max = 1, bootstrap = 1000)
set.seed(1998)
boot.lnre.ll.clust.2.FIC <- lnre("fzm", spc.ll.clust.2.FIC, m.max = 1, bootstrap = 1000)
set.seed(9821)
boot.lnre.ll.clust.3.FIC <- lnre("fzm", spc.ll.clust.3.FIC, m.max = 1, bootstrap = 1000)
set.seed(3476)
boot.lnre.ll.clust.4.FIC <- lnre("fzm", spc.ll.clust.4.FIC, m.max = 1, bootstrap = 1000)
set.seed(2999)
boot.lnre.ll.clust.5.FIC <- lnre("fzm", spc.ll.clust.5.FIC, m.max = 1, bootstrap = 1000)
set.seed(2989)
boot.lnre.ll.clust.6.FIC <- lnre("fzm", spc.ll.clust.6.FIC, m.max = 1, bootstrap = 1000)
set.seed(9911)
boot.lnre.ll.clust.7.FIC <- lnre("fzm", spc.ll.clust.7.FIC, m.max = 1, bootstrap = 1000)
set.seed(8080)
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
set.seed(4445)
boot.lnre.would.clust.1.FIC <- lnre("fzm", spc.would.clust.1.FIC, m.max = 1,
                                    bootstrap = 1000)
set.seed(1922)
boot.lnre.would.clust.2.FIC <- lnre("fzm", spc.would.clust.2.FIC, m.max = 1, 
                                    bootstrap = 1000)
set.seed(3546)
boot.lnre.would.clust.3.FIC <- lnre("fzm", spc.would.clust.3.FIC, m.max = 1,
                                    bootstrap = 1000)
set.seed(1234)
boot.lnre.would.clust.4.FIC <- lnre("fzm", spc.would.clust.4.FIC, m.max = 1,
                                    bootstrap = 1000)
set.seed(9804)
boot.lnre.would.clust.5.FIC <- lnre("fzm", spc.would.clust.5.FIC, m.max = 1, 
                                    bootstrap = 1000)
set.seed(1000)
boot.lnre.would.clust.6.FIC <- lnre("fzm", spc.would.clust.6.FIC, m.max = 1,
                                    bootstrap = 1000)
set.seed(2999)
boot.lnre.would.clust.7.FIC <- lnre("fzm", spc.would.clust.7.FIC, m.max = 1,
                                    bootstrap = 1000)
set.seed(2855)
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
set.seed(2222)
boot.lnre.will.clust.1.FIC <- lnre("fzm", spc.will.clust.1.FIC, m.max = 1,
                                   bootstrap = 1000)
set.seed(7766)
boot.lnre.will.clust.2.FIC <- lnre("fzm", spc.will.clust.2.FIC, m.max = 1,
                                   bootstrap = 1000)
set.seed(8080)
boot.lnre.will.clust.3.FIC <- lnre("fzm", spc.will.clust.3.FIC, m.max = 1,
                                   bootstrap = 1000)
set.seed(8111)
boot.lnre.will.clust.4.FIC <- lnre("fzm", spc.will.clust.4.FIC, m.max = 1,
                                   bootstrap = 1000)
set.seed(1117)
boot.lnre.will.clust.5.FIC <- lnre("fzm", spc.will.clust.5.FIC, m.max = 1,
                                   bootstrap = 1000)
set.seed(6100)
boot.lnre.will.clust.6.FIC <- lnre("fzm", spc.will.clust.6.FIC, m.max = 1,
                                   bootstrap = 1000)
set.seed(5004)
boot.lnre.will.clust.7.FIC <- lnre("fzm", spc.will.clust.7.FIC, m.max = 1,
                                   bootstrap = 1000)
set.seed(4809)
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
