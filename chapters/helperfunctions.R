# helper functions
prior.prop <- function(M, I, pred.sd, sd.scale){
  ### Prior probability of all positive
  sd <- rep(pred.sd, M) # sd of beta (common)
  #s2_country <- MCMCpack::rinvgamma(M, 3, sd.scale) #variance of country deviation from beta 
  s2_country <- rnorm(M, 0, sd.scale)^2
  mu <- rnorm(M, 0, sd)
  res <- exp(pnorm(0, mu, sqrt(s2_country), lower.tail = F, log.p = T) * I)
  return(mean(res))
}

getBFs <- function(bridge0,bridge1,bridge2,m1,m2,param,I){
  bfs=1:4
  bfs[1] <- bridge0
  # post prob common effect
  post <- posterior::as_draws_df(m1, variable=paste0("b_",param))
  PostProbOne <- mean(post[,1] > 0)
  bfs[2] <- bridge1 + log(PostProbOne) - log(0.5)
  # post prob varying effect + fixed effect 
  post <- data.frame(ranef(m2, summary=FALSE, pars = param)$site[,,1] + 
                       fixef(m2, summary = FALSE, pars = param)[,1])
  PriorProb <- prior.prop(M=1e6, I = I, pred.sd = 0.5, sd.scale = 1)
  PostProb <-  mean(apply(post > 0, 1, mean) == 1)
  bfs[3] <- bridge2 + log(PostProb) - log(PriorProb)
  bfs[4] <- bridge2
  return(bfs = bfs)
}
specify_decimalBF = function(x, k=2, short=TRUE){
  if(x==Inf){
    y=printnum(Inf)
  } else if(x>99999){
    y=sfsmisc::pretty10exp(x,digits = k, lab.type = "latex", lab.sep = "times")
    if(short) y=paste0("$",strsplit(y," ")[[1]][3])
  } else if((x<=99999&x>99)|x==1) {
    y=sprintf("%.0f",x)
  } else if(x<.01&x>=.0005){
    y=sprintf("%.3f",x)
  } else if(x<=99&x>=.01&x!=1){
    y=sprintf(paste0("%.",k,"f"),x)
  } else if (x<.0005&x>0){
    y=sfsmisc::pretty10exp(x,digits = k, lab.type = "latex", lab.sep = "times")
    if(short) y=strsplit(y," ")[[1]][3]
  } else if (x<=0){
    y=sprintf(paste0("%.",k,"f"),x)
  }
  return(y)
}

writeBF <- function(bfs, short=TRUE, log=FALSE){
  if(log) {bfs <- max(bfs)-bfs} else {bfs <- 1/exp(bfs-max(bfs))}
  n <- sapply(bfs,function(x) specify_decimalBF(x,2))
  if(short){
    out = sapply(n, function(x) strsplit(x," ")[[1]][3])
    out = ifelse(is.na(out),n,paste0("$",out))
    n = out
  }
  bf <- paste("1-to-",n,sep="")
  return(bf)
}
writeBF01 <- function(bfs){
  bfs <- exp(bfs-max(bfs))
  return(bfs)
}
writeCI <- function(x,digits=2, percentage=FALSE){
  y = paste0(printnum(x[1],digits=digits), 
             " [", printnum(x[2],digits=digits), ", ", 
             printnum(x[3],digits=digits),"]")
  if(percentage) y = paste0(printnum(x[1]*100,digits=digits-2), 
                            "% [", printnum(x[2]*100,digits=digits-2), "%, ", 
                            printnum(x[3]*100,digits=digits-2),"%]")
  return(y)
}
printBFab <- function(bfa, bfb){
  bfab <- 1/exp(bfb-bfa)
  bfab <- specify_decimalBF(bfab,2)
  return(bfab)
}
post_pred <- function(model, effect, interaction=F){
  post <- as.data.frame(as_draws_df(model))
  a_sim <- rnorm(length(post$Intercept), post$Intercept, post$sd_site__Intercept)
  b_sim <- rnorm(length(post$Intercept), post[,paste0("b_",effect)], post[,paste0("sd_site__",effect)])
  if(interaction){
    c_sim <- rnorm(length(post$Intercept), post[,paste0("b_",strsplit(effect,split=":")[[1]][1])], 
                   post[,paste0("sd_site__",strsplit(effect,split=":")[[1]][1])])
    d_sim <- rnorm(length(post$Intercept), post[,paste0("b_",strsplit(effect,split=":")[[1]][2])], 
                   post[,paste0("sd_site__",strsplit(effect,split=":")[[1]][2])])
  }
  p_link_asim <- function(x){
    logodds <- with(post, a_sim + b_sim*x)
    logodds <- with(post, a_sim + b_sim*x)
    return(inv_logit_scaled(logodds))
  }
  p_raw_sim <- sapply(c(-1/2,1/2), function(i) p_link_asim(i))
  diff <- p_raw_sim[,1]-p_raw_sim[,2]
  return(list(diff=diff,p_sim=p_raw_sim))
}

# plots
# Define custom y axis function
base_breaks_y <- function(yLimits) {
  d <- data.frame(x=-Inf, xend=-Inf, y=yLimits[1], yend=yLimits[2])
  list(ggplot2::geom_segment(data=d, ggplot2::aes(x=x, y=y, xend=xend, yend=yend), size = 0.75, inherit.aes=FALSE))
}
base_breaks_x <- function(xLimits) {
  d <- data.frame(x=xLimits[1], xend=xLimits[2], y=-Inf, yend=-Inf)
  list(ggplot2::geom_segment(data=d, ggplot2::aes(x=x, y=y, xend=xend, yend=yend), size = 0.75, inherit.aes=FALSE))
}

my_theme = function(){
    ggplot2::theme(
      text = element_text(family="",face="plain"),
      panel.grid.minor=   ggplot2::element_blank(),
      plot.title=         ggplot2::element_text(size = 14),
      panel.grid.major=   ggplot2::element_blank(),
      axis.title.x=       ggplot2::element_text(size = 12),
      axis.title.y=       ggplot2::element_text(size = 12),
      axis.text.x=        ggplot2::element_text(size = 9),
      axis.text.y=        ggplot2::element_text(size = 9),
      panel.background=   ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.background=    ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.border=       ggplot2::element_blank(),
      axis.line=          ggplot2::element_blank(),
      axis.ticks.y=       ggplot2::element_line(size = 0.5),
      axis.ticks.x=       ggplot2::element_line(size = 0.5),
      axis.ticks.length=  grid::unit(1.5, "mm"),
      plot.margin=        grid::unit(c(0.1, 0.1, 0.6, 0.6), "cm"), 
      legend.position=    "bottom",
      legend.background=  ggplot2::element_blank(),
      strip.background=   ggplot2::element_blank(),
      strip.text.x=       ggplot2::element_text(size = 11)
    )
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

