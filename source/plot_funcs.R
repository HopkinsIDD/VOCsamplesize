## ---------- Functions to plot figures ---------- ##

##' Function to plot simulated variant frequency estimates
##' with binomial and hypergeometric distributions on top
##'
##' @param simdata simulation data
##' @param ssize sample sizes to plot
##' @param varfreq input variant frequencies to plot

plot.dist.validation.var <- function(simdata,ssize,varfreq){
  
  dtemp <- simdata %>% filter(v1.in==varfreq) %>% filter(samplesize==ssize)
  
  sim_dens <- ggplot(dtemp) + geom_histogram(aes(x=v1.out,y=(after_stat(count))/sum(after_stat(count))),bins=20)
  sim_dens <- ggplot_build(sim_dens)$data[[1]] %>% select(x,y)
  bw <- sim_dens$x[2]-sim_dens$x[1]
  
  ggplot() +
    geom_bar(data=sim_dens,aes(x=x,y=y),stat="identity",fill="black",alpha=0.4) +
    geom_density(data=dtemp,aes(x=v1.out.hyper,y=(after_stat(density))*bw),adjust=5,colour="light blue",linewidth=1) +
    geom_density(data=dtemp,aes(x=v1.out.binom,y=(after_stat(density))*bw),adjust=5,colour="blue",linewidth=1) +
    geom_vline(xintercept = varfreq, linetype="dashed", color="black", linewidth=0.5) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    xlab("v1.out") +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank())
  
}

##' Function to plot heatmap of variant detection probabilities
##' given sample size and variant prevalence or time to detection
##' for a variant undergoing logistic growth
##'
##' @param df data to plot (includes one probability value per coefficient of detection ratio)
##' @param c_ratio coefficient of detection ratio to plot
##' @param c_label unique column label from df
##' @param t0 initial variant prevalence in population
##' @param r logistic growth rate

plot.ss.logistic <- function(df,c_ratio,c_label,t0,r){
  
  f <- sym(c_label)
  
  prev_breaks <- sapply(seq(0,40,10),phylosamp::varfreq_freq_logistic,t0,r)
  prev_labels <- sprintf("%.4f", prev_breaks)
  
  ggplot(df,aes(x=t,y=n,fill=!!f)) +
    geom_tile() +
    geom_contour(aes(z = !!f), color='black', linewidth=0.4, alpha=0.5, breaks=c(0.5,0.95)) +
    geom_text_contour(aes(z = !!f), breaks=c(0.5,0.95)) +
    scale_fill_gradient(low="light yellow", high="blue",name="Pr(d ≤ t)",limits=c(0,1)) +
    scale_x_continuous(limits=c(-1,41),expand = c(0,0),
                       sec.axis = sec_axis(trans = ~ phylosamp::varfreq_freq_logistic(.,t0,r),
                                           name="Prevalence",breaks=prev_breaks,labels=prev_labels)) +
    scale_y_continuous(expand = c(0.01,0)) +
    labs(x="Days",y="Sample size per day") +
    theme_minimal() +
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8),
          panel.grid.major = element_line(color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = "none")
}
