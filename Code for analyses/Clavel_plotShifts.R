

#This code is to make the plotshift output 
#Run this before running 'Clavel_models_shifts_jumps.R'

library(fields)

plotShifts <-
  function(phylo, chain, burnin=1000, ...){
    #require(fields)
    args <- list(...)
    # options
    if(is.null(args[["fun"]])) args$fun <- mean
    if(is.null(args[["show.tip.label"]])) args$show.tip.label <- TRUE
    if(is.null(args[["horizontal"]])) args$horizontal <- TRUE
    if(is.null(args[["color"]])) args$color <- c("blue", "red")
    if(is.null(args[["scale"]])) args$scale <- FALSE
    if(is.null(args[["log"]])) args$log <- FALSE
    if(is.null(args[["palette"]])) args$palette <- FALSE
    if(is.null(args[["main"]])) args$main <- NULL
    if(is.null(args[["cex"]])) args$cex <- 0.8
    if(is.null(args[["width"]])) args$width <- 1
    
    if(inherits(chain, "mcmc")){
      tot <- nrow(chain)
      if(burnin>tot) stop("Error! the burnin value is higher than the chain length")
      chain <- chain[c(burnin:tot),-1]
      meanRate <- apply(chain, 2, args$fun)
      if(args$log==TRUE) meanRate <- log(meanRate)
    }else{
      meanRate <- chain
      if(args$log==TRUE) meanRate <- log(chain)
    }
    
    # check the order of the tree; prunning algorithm use "postorder"
    # if(attr(phylo,"order")!="postorder") phylo <- reorder.phylo(phylo, "postorder")
    
    # colors mapping
    if(any(args$palette==FALSE)){
      Colors = colorRampPalette(args$color)( 100 )
    }else{
      Colors = args$palette
    }
    # 0 index induce error I scale it between 1 and 100
    linScale <- function(x, from, to) round( (x - min(x)) / max(x - min(x)) * (to - from) + from)
    col <- linScale(meanRate, from=1, to=100)
    
    if(args$scale==TRUE){
      phylo$edge.length <- phylo$edge.length*meanRate
    }
    plot(phylo, edge.color = Colors[col], show.tip.label = args$show.tip.label, main = args$main, cex = args$cex, edge.width = args$width)
    
    image.plot(z = as.matrix(meanRate),col = Colors,
               legend.only = T, horizontal = args$horizontal)
    
  }
