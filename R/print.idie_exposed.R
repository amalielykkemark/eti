#'@title Print function for \code{idie_exposed}
#'@description Print function for \code{idie_exposed} a targeted minimum-loss based estimator for the Interventional
#'Disparity Indirect Effect (IDIE) among the exposed.
#'@name print.idie_exposed
#'@author Amalie Lykkemark Moller \email{amalielykkemark@@live.dk}
#'@rdname print.idie_exposed
#'@method print idie_exposed
#'@export
print.idie_exposed <- function(x,...) {
  if(identical(class(x), "idie_exposed")){

    table<-data.table(data.frame(Estimate=unlist(x$estimate),
                                 SE=unlist(x$se)))
    table[,CI95lower:=Estimate-1.96*SE]
    table[,CI95upper:=Estimate+1.96*SE]

    table[,Parameter:=c('Risk under stochastic intervention','Risk under no intervetion','Risk difference (IDIE among the exposed)')]
    table<-table[,.(Parameter,Estimate,SE,CI95lower,CI95upper)]

    print(table)
  }
}
