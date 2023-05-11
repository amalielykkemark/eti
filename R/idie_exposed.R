#'@title Targeted minimum-loss based estimator for the Interventional Disparity Indirect Effect (IDIE) among the exposed
#'@description The \code{idie_exposed} is a Targeted Minumum-loss based
#'estimator (TMLE) for the Interventional Disparity Indirect Effect (IDIE) among the exposed.
#' We consider a structure, where the mediator (Z) is an effect of the exposure (A) and a cause of the outcome (Y),
#' A -> Z -> Y. The function estimates the expected change in outcome risk among the exposed (A=1) if
#'hypothetically the exposed had the same probability of the mediator (Z)
#'as observed for similar unexposed (A=0) individuals.
#'The expected outcome risk under this hypothetical intervention is compared to the outcome risk
#'among the exposed when the distribution of the mediator was set
#'to the observed level among the exposed. Importantly, for this estimator the
#'exposure, mediator, and outcome must all be binary. The
#'underlying model for the exposure, mediator, and outcome, which are needed
#'to estimate any of the parameters, can be modeled using Super learning. Super
#'learning can be used to produce a weighted combination of candidate algorithms
#'that optimize the cross-validated loss-function or to select the single best
#'performing algorithm among the candidate algorithms, also known as the discrete
#'Super learner. One must define a library of candidate algorithms which should
#'be considered by the Super learner. If the Super learner library contains only
#'one algorithm, results will be estimated based on this algorithm alone, and
#'thus, not using Super Learning.
#'
#'@name idie_exposed
#'
#'@author Amalie Lykkemark Moller \email{amlm@@sund.ku.dk} Helene Charlotte Wiese Rytgaard, Thomas Alexander Gerds, and Christian Torp-Pedersen
#'@usage idie_exposed(data, discrete.SL=TRUE, exposure.A=NA, mediator.Z=NA, outcome.Y=NA,
#'cov.A, cov.Z, cov.Y, SL.lib.A=FALSE, SL.lib.Z=FALSE, SL.lib.Y=FALSE, iterations=10)
#'@param data  A data frame/data table with a binary exposure, a binary
#'  mediator, a binary outcome, and covariates.
#'@param discrete.SL  If \code{FALSE} TMLE will use the weighted combination of
#'  the algorithms in the super learner library that minimizes the loss
#'  function. Defaults to \code{TRUE}, in which case the discrete super learner
#'  i used, i.e. the best performing algorithm.
#'@param exposure.A  Name of the binary exposure.
#'@param mediator.Z  Name of the binary mediator, which is the
#'  target of the hypothetical intervention.
#'@param outcome.Y  Name of the binary outcome.
#'@param cov.A  A vector containing names of possible confounders which should
#'  be included in models of the exposure.
#'@param cov.Z  A vector of confounders which should be included in models of
#'  the mediator. Do not include the exposure as the function does this.
#'@param cov.Y  A vector of confounders which should be included in models of
#'  the outcome. Do not include the exposure and the mediator as the
#'  function does this.
#'@param SL.lib.A  A vector of algorithms that should be considered by the super
#'  learner when modelling the exposure. All algorithms must be specified as
#'  Super Learner objects.
#'@param SL.lib.Z  A vector of algorithms for modelling the mediator. All
#'  algorithms must be specified as Super Learner objects.
#'@param SL.lib.Y  A vector of algorithms for modelling the outcome. All
#'  algorithms must be specified as Super Learner objects.
#'@param iterations  Number of iterations for the updating step in TMLE.
#'  Defaults to 10.
#'
#'@details The structure of the data should be as follows: \item For the binary
#'exposure (\code{exposure.A}) 1 = exposed and 0 = unexposed. \item For the
#'binary mediator (\code{mediator.Z}) 1 = treatment and 0 = no
#'treatment. \item For the binary outcome (\code{outcome.Y}) 1 = event and 0 =
#'no event.
#'
#'@return The function outputs the absolute outcome risk among the exposed had
#'  their chance of the mediator been the same as for similar unexposed individuals,
#'  the absolute outcome risk among the exposed under no intervention, where the
#'  probability of the mediator is as observed (psi1), the
#'  absolute risk difference between the two, the interventional disparity indirect
#'  effect among the exposed, and standard errors for each estimate.
#'
#' @import data.table
#' @import SuperLearner
#'
#' @examples
#'library(data.table)
#'require(tmleExposed)
#'n=5000
#'set.seed(1)
#'sex <- rbinom(n,1,0.4)
#'age <- rnorm(n,65,sd=5)
#'disease <- rbinom(n,1,0.6)
#'
#'A <- rbinom(n, 1, plogis(-3+0.05*age+1*sex))
#'Z <- rbinom(n, 1, plogis(5-0.08*age+1*sex-1.2*disease-0.8*A+0.01*A*disease))
#'Y <- rbinom(n, 1, plogis(-9+0.09*age+0.5*sex+0.8*disease-1.2*Z+0.7*A))
#'
#'d <- data.table(id=1:n, exposure=as.integer(A), mediator=as.integer(Z),
#'                outcome=as.integer(Y), age, sex, disease)
#'
#'##### Define algorithms for the Super Learner library #####
#'lib = c('SL.glm','SL.step.interaction')
#'
#' #intervention: changing probability of the mediator (Z=1) among the exposed (A=1)
#' #to what it would have been had they been unexposed (A=0).
#' #target parameter: the change in outcome among the exposed (A=1) had their chance of
#' #the mediator (Z=1) been as among similar unexposed individuals (A=0).
#'
#'res<-idie_exposed(data=d,
#'                  exposure.A='exposure',
#'                  mediator.Z='mediator',
#'                  outcome.Y='outcome',
#'                  cov.A=c('sex','age'),
#'                  cov.Z =c('sex','age','disease'),
#'                  cov.Y=c('sex','age','disease'),
#'                  SL.lib.A = lib,
#'                  SL.lib.Z = lib,
#'                  SL.lib.Y = lib,
#'                  discrete.SL = TRUE)
#'
#'summary(res)
#'
#'
#'@export
idie_exposed<-function(data,
                       discrete.SL=TRUE,
                       exposure.A=NA,
                       mediator.Z=NA,
                       outcome.Y=NA,
                       cov.A,
                       cov.Z,
                       cov.Y,
                       SL.lib.A=FALSE,
                       SL.lib.Z=FALSE,
                       SL.lib.Y=FALSE,
                       iterations=10){

  requireNamespace('data.table')
  requireNamespace('SuperLearner')
  requireNamespace('riskRegression')

  #variable and input check
  if(is.na(exposure.A)|is.na(mediator.Z)|is.na(outcome.Y)) {
    stop(paste('Please speficy names for the exposure, mediator, and outcome variables'))
  }

  if(is.data.table(data)){
    dt <- data.table::copy(data)
  }else{
    dt <- data.table::as.data.table(data)
  }

  if(exists('id',dt)==FALSE){
    dt[,id:=1:.N]
  }

  #convert variable names
  if (exposure.A!='A'){
    setnames(dt,exposure.A,'A')
  }

  if (mediator.Z!='Z'){
    setnames(dt,mediator.Z,'Z')
  }

  if (outcome.Y!='Y'){
    setnames(dt,outcome.Y,'Y')
  }

  #stop if bigger than 2 instead of this
  if (length(unique(dt[,A]))!=2 | length(unique(dt[,Z]))!=2 | length(unique(dt[,Y]))!=2) {
    stop('Exposure, mediator, and outcome must be binary.')
  }

  pibar <- dt[,mean(A==1)]
  pifit <- SuperLearner::SuperLearner(Y=dt[,A], #dt[,as.numeric(A==1)],
                                      X=dt[,.SD,.SDcols=cov.A],
                                      family = binomial(),
                                      SL.library = SL.lib.A)

  gammafit <-SuperLearner::SuperLearner(Y=dt[,Z],#dt[,as.numeric(Z==1)],
                                        X=dt[,.SD,.SDcols=c(cov.Z,'A')],
                                        family = binomial(),
                                        SL.library = SL.lib.Z)

  Qfit <-SuperLearner::SuperLearner(Y=dt[,Y],#dt[,as.numeric(Y==1)],
                                    X=dt[,.SD,.SDcols=c(cov.Y,'A','Z')],
                                    family = binomial(),
                                    SL.library = SL.lib.Y)

  if (discrete.SL==FALSE){
    dt[, pihat:=predict(pifit, newdata=dt[,.SD,.SDcols=c(cov.A)], onlySL = T)$pred]

    dt[, gammahat:=predict(gammafit, newdata=copy(dt[,.SD,.SDcols=c(cov.Z,'A')]), onlySL = T)$pred]
    dt[, gammahat.a0:=predict(gammafit, newdata=copy(dt[,.SD,.SDcols=c(cov.Z)])[, A:=0], onlySL = T)$pred]
    dt[, gammahat.a1:=predict(gammafit, newdata=copy(dt[,.SD,.SDcols=c(cov.Z)])[,A:=1], onlySL = T)$pred]

    dt.full <- data.table(rbind(copy(dt)[, Z:=1],
                                copy(dt)[, Z:=0]),
                          Z.obs=c(dt[, Z], dt[, Z]))

    dt.full[, Qhat:=predict(Qfit, newdata=dt.full[,.SD,.SDcols=c(cov.Y,'A','Z')], onlySL = T)$pred]
    dt.full[, Qhat.a1:=predict(Qfit, newdata=copy(dt.full[,.SD,.SDcols=c(cov.Y,'Z')])[,A:=1], onlySL = T)$pred]
    dt.full[, Qhat.a1.z0:=predict(Qfit, newdata=copy(dt.full[,.SD,.SDcols=c(cov.Y)])[, `:=`(A=1, Z=0)], onlySL = T)$pred]
    dt.full[, Qhat.a1.z1:=predict(Qfit, newdata=copy(dt.full[,.SD,.SDcols=c(cov.Y)])[, `:=`(A=1, Z=1)], onlySL = T)$pred]
  }

  if (discrete.SL==TRUE){
    p<-pifit$libraryNames[which.max(pifit$coef)]
    pifit.discrete<-pifit$fitLibrary[[p]]
    g<-gammafit$libraryNames[which.max(gammafit$coef)]
    gammafit.discrete<-gammafit$fitLibrary[[g]]
    Q<-Qfit$libraryNames[which.max(Qfit$coef)]
    Qfit.discrete<-Qfit$fitLibrary[[Q]]
    dt[, pihat:=predict(pifit.discrete, newdata=copy(dt[,.SD,.SDcols=c(cov.A)]), type="response")]

    dt[, gammahat:=predict(gammafit.discrete, newdata=copy(dt[,.SD,.SDcols=c(cov.Z,'A')]), type="response")]
    dt[, gammahat.a0:=predict(gammafit.discrete, newdata=copy(dt[,.SD,.SDcols=c(cov.Z)])[, A:=0], type="response")]
    dt[, gammahat.a1:=predict(gammafit.discrete, newdata=copy(dt[,.SD,.SDcols=c(cov.Z)])[, A:=1], type="response")]

    dt.full <- data.table(rbind(copy(dt)[, Z:=1],
                                copy(dt)[, Z:=0]),
                          Z.obs=c(dt[, Z], dt[, Z]))

    dt.full[, Qhat:=predict(Qfit.discrete, newdata=dt.full, type="response")]
    dt.full[, Qhat.a1:=predict(Qfit.discrete, newdata=copy(dt.full)[, A:=1], type="response")]
    dt.full[, Qhat.a1.z0:=predict(Qfit.discrete, newdata=copy(dt.full)[, `:=`(A=1, Z=0)], type="response")]
    dt.full[, Qhat.a1.z1:=predict(Qfit.discrete, newdata=copy(dt.full)[, `:=`(A=1, Z=1)], type="response")]
  }

  # no intervention
  dt.full[, psi.1:=sum(Qhat.a1*(gammahat.a1*Z+(1-gammahat.a1)*(1-Z))), by="id"]
  # as among the unexposed
  dt.full[, psi.0:=sum(Qhat.a1*(gammahat.a0*Z+(1-gammahat.a0)*(1-Z))), by="id"]

  init.est0 <- dt.full[Z==Z.obs, 1/pibar*mean((A == 1)*(psi.0))]
  init.est1 <- tmle.est1 <- dt.full[Z==Z.obs & A == 1, mean(Y)]

  dt.full[Z==Z.obs,eic0:= ((A==1)/pibar*((Z*gammahat.a0+(1-Z)*(1-gammahat.a0))/
                                             (Z*gammahat.a1+(1-Z)*(1-gammahat.a1)))*
                               (Y - Qhat.a1) +
                               (A==0)/(1-pihat)*(pihat)/pibar*(Qhat.a1 - psi.0) +
                               (A==1)/pibar*(psi.0 - init.est0))]

  dt.full[Z==Z.obs,eic1:= ((A==1)/pibar*(Y - init.est1))]

  se0 <- sqrt(mean(dt.full[Z==Z.obs,eic0]^2)/nrow(dt))
  se1 <- sqrt(mean(dt.full[Z==Z.obs,eic1]^2)/nrow(dt))
  dt.full[Z==Z.obs,eic:=eic0-eic1]
  se.diff <-  sqrt(mean((dt.full[Z==Z.obs,eic])^2)/nrow(dt))

  dt.full.copy <- copy(dt.full)

  #-------- tmle
  #-- psi0;
  for (iter in 1:iterations) {  #-- iterative updating;
    #-- update Q;
    dt.full[, H.Y:=(A==1)/pibar*((Z*gammahat.a0+(1-Z)*(1-gammahat.a0))/
                                   (Z*gammahat.a1+(1-Z)*(1-gammahat.a1)))]
    eps.Y <- coef(glm(Y ~ offset(qlogis(Qhat))+1,
                            data=dt.full[Z==Z.obs],
                            weights=H.Y,
                            family=quasibinomial()))

    dt.full[, Qhat:=plogis(qlogis(Qhat) + eps.Y)]
    dt.full[, Qhat.a1:=plogis(qlogis(Qhat.a1) + eps.Y)]
    dt.full[, Qhat.a1.z0:=plogis(qlogis(Qhat.a1.z0) + eps.Y)]
    dt.full[, Qhat.a1.z1:=plogis(qlogis(Qhat.a1.z1) + eps.Y)]

    #-- update Z;
    dt.full[, H.Z.weight:=(A==0)/(1-pihat)*pihat/pibar]
    dt.full[, H.Z.covar:=(Qhat.a1.z1 - Qhat.a1.z0)]
    eps.Z <- coef(glm(Z ~ offset(qlogis(gammahat))+H.Z.covar-1,
                            weights=H.Z.weight,
                            data=dt.full[Z==Z.obs], family=quasibinomial()))

    dt.full[, gammahat:=plogis(qlogis(gammahat) + eps.Z*H.Z.covar)]
    dt.full[, gammahat.a1:=plogis(qlogis(gammahat.a1) + eps.Z*H.Z.covar)]
    dt.full[, gammahat.a0:=plogis(qlogis(gammahat.a0) + eps.Z*H.Z.covar)]

    dt.full[, psi.0:=sum(Qhat.a1*(gammahat.a0*Z+(1-gammahat.a0)*(1-Z))), by="id"]
    #updated estimate
    tmle.est0 <- psihat.0 <- dt.full[Z==Z.obs, (1/pibar)*mean((A == 1)*psi.0)]

    #-- update A;
    dt.full[, H.A:=(psi.0-tmle.est0)/pibar]

    eps.A <- coef(glm(A==1 ~ offset(qlogis(pihat))+H.A-1, data=dt.full[Z==Z.obs], family=quasibinomial()))
    dt.full[, pihat:=plogis(qlogis(pihat) + eps.A*H.A)]

    #updated estimate
    tmle.est0 <- psihat.0 <- dt.full[Z==Z.obs, (1/pibar)*mean((A == 1)*psi.0)]

    solve.eic0 <- abs(dt.full[Z==Z.obs, mean((A==1)/pibar*((Z*gammahat.a0+(1-Z)*(1-gammahat.a0))/
                                                             (Z*gammahat.a1+(1-Z)*(1-gammahat.a1)))*
                                               (Y - Qhat.a1) +
                                               (A==0)/(1-pihat)*(pihat)/pibar*(Qhat.a1 - psi.0) +
                                               (A==1)/pibar*(psi.0 - psihat.0))])

    if (solve.eic0<=se0/(log(nrow(dt))*sqrt(nrow(dt)))) break

    if(iter==iterations) {
      warning(paste('Efficient influence function for psi0 was not solved in',iterations,'iterations'))
    }
  }

  #-- psi1;

  dt.full <- copy(dt.full.copy)

  #tmle.est1 <- dt.full[Z==Z.obs & A == 1, mean(Y)]

  solve.eic1 <- abs(dt.full[Z==Z.obs, mean((A==1)/pibar*(Y - tmle.est1))]) # obs overflodig?

  # dt.full <- copy(dt.full.copy)
  #
  # for (iter in 1:iterations) {  #-- iterative updating;
  #
  #   #-- update Q;
  #   dt.full[, H.Y:=1/pibar]
  #   eps.Y <- coef(glm(Y ~ offset(qlogis(Qhat.a1))+H.Y-1, data=dt.full[Z==Z.obs & A==1], family=binomial()))
  #
  #   dt.full[, Qhat.a1:=plogis(qlogis(Qhat.a1) + eps.Y/pibar)]
  #   dt.full[, Qhat.a1.z0:=plogis(qlogis(Qhat.a1.z0) + eps.Y/pibar)]
  #   dt.full[, Qhat.a1.z1:=plogis(qlogis(Qhat.a1.z1) + eps.Y/pibar)]
  #
  #   dt.full[, psi.1:=sum(Qhat.a1*(gammahat.a1*Z+(1-gammahat.a1)*(1-Z))), by="id"]
  #
  #   tmle.est1 <- dt.full[Z==Z.obs, (1/pibar)*mean((pihat)*psi.1)]
  #
  #   #-- update A;
  #   dt.full[, H.A:=(psi.1-tmle.est1)/pibar]
  #   eps.A <- coef(glm(A==1 ~ offset(qlogis(pihat))+H.A-1, data=dt.full[Z==Z.obs], family=binomial()))
  #   dt.full[, pihat:=plogis(qlogis(pihat) + eps.A*H.A)]
  #
  #   #-- updated estimate
  #   tmle.est1 <- dt.full[Z==Z.obs, (1/pibar)*mean((A)*psi.1)]
  #
  #   solve.eic1 <- abs(dt.full[Z==Z.obs, mean((A==1)/pibar*(Y - Qhat.a1) +
  #                                              (psi.1 - tmle.est1)/pibar*((A==1) - (pihat)) +
  #                                              (pihat)/pibar*(psi.1-tmle.est1))])
  #
  #   if (solve.eic1<=se1/(log(nrow(dt))*sqrt(nrow(dt)))) break
  #   if(iter==iterations) {
  #     warning(paste('Efficient influence function for psi1 was not solved in',iterations,'iterations'))
  #   }
  # }

  # prepare output
  out<-list()
  psi.diff.tmle<-tmle.est0-tmle.est1

  out$estimate=list(psi0=tmle.est0,psi1=tmle.est1,psi=psi.diff.tmle)
  out$se=list(se0=se0,se1=se1,se.diff=se.diff)

  #CV risk for exposure, mediator, and outcome
  out$superlearner.CVrisk$A.exposure<-pifit$cvRisk
  out$superlearner.CVrisk$Z.mediator<-gammafit$cvRisk
  out$superlearner.CVrisk$Y.outcome<-Qfit$cvRisk
  #weights assigned to each algorithm in the super learner library
  if (discrete.SL==TRUE){
    out$superlearner.discrete$A.exposure<-p
    out$superlearner.discrete$Z.mediator<-g
    out$superlearner.discrete$Y.outcome<-Q
  }
  if (discrete.SL==FALSE){
    out$superlearner.weight$A.exposure<-pifit$coef
    out$superlearner.weight$Z.mediator<-gammafit$coef
    out$superlearner.weight$Y.outcome<-Qfit$coef
  }
  out$distributions=rbind(distribution.A1=dt.full[Z==Z.obs & A==1, summary(pihat)],
                          distribution.Z.a1=dt.full[Z==Z.obs & A==1, summary(gammahat.a1)],
                          distribution.Z.a0=dt.full[Z==Z.obs & A==1, summary(gammahat.a0)],
                          distribution.Y=dt.full[Z==Z.obs & A==1, summary(Qhat)],
                          distribution.Y.a1=dt.full[Z==Z.obs & A==1, summary(Qhat.a1)])

  out$output.dataset<-dt.full

  class(out)<-'idie_exposed'

  return(invisible(out))
}



