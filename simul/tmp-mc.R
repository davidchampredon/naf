library(snowfall)

n.cpu <- 64

f <- function(i){
  set.seed(i)
  x <- data.frame(a=i, b=runif(1e6))
  return(x)
}


sfInit(parallel = (n.cpu>1), cpu = n.cpu)
x.list <- list()
x.list <- sfSapply(1:n.cpu, f, simplify= FALSE)
sfStop()
x.all <- do.call('rbind.data.frame',x.list)

head(x.all)
object.size(x.all)/1e6
