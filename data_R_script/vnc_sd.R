vnc_sd <- function (filename = NULL, distance.measure = "sd") 
{
  {
  user.input <- read.delim(filename)
  stopifnot(ncol(user.input) == 2)
  stopifnot(all(apply(user.input, 2, is.vector)))
  input <- user.input[, 1]
  years <- user.input[, 2]
  names(input) <- years
  }
  data.collector <- list()
  data.collector[["0"]] <- input
  position.collector <- list()
  position.collector[[1]] <- 0
  overall.distance <- 0
  number.of.steps <- length(input) - 1
  for (i in 1:number.of.steps) {
    cat(i/number.of.steps, "\n", sep = "")
    difference.checker <- numeric()
    for (j in 1:(length(unique(names(input))) - 1)) {
      first.name <- unique(names(input))[j]
      second.name <- unique(names(input))[(j + 1)]
      pooled.sample <- input[names(input) %in% c(first.name, 
                                                 second.name)]
      difference.checker[j] <- ifelse(sum(pooled.sample) == 
                                        0, 0, sd(pooled.sample))
    }
    pos.to.be.merged <- which.min(difference.checker)
    distance <- min(difference.checker)
    overall.distance <- overall.distance + distance
    lower.name <- unique(names(input))[pos.to.be.merged]
    higher.name <- unique(names(input))[(pos.to.be.merged + 
                                           1)]
    matches <- names(input) %in% c(lower.name, higher.name)
    new.mean.age <- round(mean(as.numeric(names(input)[names(input) %in% 
                                                         c(lower.name, higher.name)])), 4)
    position.collector[[(i + 1)]] <- which(names(input) == 
                                             lower.name | names(input) == higher.name)
    names(input)[names(input) %in% c(lower.name, higher.name)] <- as.character(new.mean.age)
    data.collector[[(i + 1)]] <- input
    names(data.collector)[(i + 1)] <- distance
  }
  par(mfrow=c(1,2))
  plot(0, mean(years), xlim = range(years), ylim = c(0, 1.1 * 
                                                       sum(as.numeric(names(data.collector)))), 
       xlab = "Time", ylab = "Distance in summed standard deviations", 
       type = "n", axes = FALSE, main="VNC dendrogram")
  axis(1, at = years)
  axis(2)
  cur.y <- rep(0, length(data.collector[[1]]))
  for (k in 1:(length(data.collector) - 1)) {
    cur.x <- as.numeric(unique(names(data.collector[[k]])))
    arrows(cur.x, cur.y, cur.x, cur.y + as.numeric(unique(names(data.collector[(k + 
                                                                                  1)]))), col = "grey", length = 0)
    cur.y <- rep(unique(cur.y) + as.numeric(unique(names(data.collector[(k + 
                                                                           1)]))), length(cur.x) - 1)
    left <- min(position.collector[[(k + 1)]])
    right <- max(position.collector[[(k + 1)]])
    lower.x <- as.numeric(names(data.collector[[k]])[left])
    higher.x <- as.numeric(names(data.collector[[k]])[right])
    arrows(lower.x, unique(cur.y), higher.x, unique(cur.y), 
           length = 0, col = "grey")
  }
  #superimpose frequency development on dendrogram
  par(new=T)
  plot(years, input, type="b", axes=F, 
       xlab = NA, ylab = NA)
  #scree plot
  plot(rev(names(data.collector)) ~ c(1:length(years)), 
       xlab = "Clusters", ylab = "Distance in standard deviations", 
       type = "c", main="Scree plot")
  grid()
  text(c(1:length(years))[-length(years)], as.numeric(rev(names(data.collector)))[-length(years)], 
       labels = round(as.numeric(rev(names(data.collector))), 
                      1)[-length(years)], cex = 1)

}
