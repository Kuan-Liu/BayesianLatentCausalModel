library(xtable)

results<-read.table("simulation_results/n125_nindic10_H.txt")

tablerow <- function(estimator, pevalues, vevalues, truevalue, cilvalues=NULL, ciuvalues=NULL) {
  return(data.frame(Estimator=estimator,
                    Mean=mean(pevalues),
                    RB=100.0 * mean((pevalues - truevalue)/abs(truevalue)),
                    SD=sd(pevalues),
                    SE=mean(vevalues),
                    CP=ifelse(is.null(cilvalues),
                              100.00 * mean((truevalue > (pevalues + qnorm(0.025) * vevalues)) & (truevalue < (pevalues + qnorm(0.975) * vevalues))),
                              100.00 * mean((truevalue > cilvalues) & (truevalue < ciuvalues)))
  ))
}

tab <- NULL
tab <- rbind(tab,tablerow('u11-u00 True', results[,1], 0, mean(results[,1])))
tab <- rbind(tab,tablerow('u11-u00 REG2', results[,10], results[,13], mean(results[,1])))
tab <- rbind(tab,tablerow('u11-u00 MSM', results[,16], results[,19], mean(results[,1])))
tab <- rbind(tab,tablerow('u11-u00 TMLE', results[,28], results[,31], mean(results[,1])))
tab <- rbind(tab,tablerow('u11-u00 Bayes', results[,34], results[,37], mean(results[,1]), results[,40],results[,41]))

for (i in 1:ncol(tab))
  if (is.numeric(tab[,i]))
    tab[,i] <- round(tab[,i], 3)
print(tab)
xtable(tab)