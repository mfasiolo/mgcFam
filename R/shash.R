
# Added shash
shash <- function (link = list("identity", "logeb", "identity", "slogit"), b = 0.01, a1 = -10, a2 = 1) 
{ 
  sech <- function(.x){ 1 / cosh(.x) }
  
  if(a2 < a1){ stop("a2 must be >= a1") }
  
  npar <- 4
  if (length(link) != npar) stop("shash requires 4 links specified as character strings")
  okLinks <- list("identity", "logeb", "identity", "slogit")
  stats <- list()
  param.names <- c("mu", "tau", "eps", "phi")
  for (i in c(1, 3)) { # Links for mu, eps and phi
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- make.link(link[[i]]) else 
      stop(link[[i]]," link not available for ", param.names[i]," parameter of shashlss")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
                          mu.eta=stats[[i]]$mu.eta),
                     class="family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  } 
  
  # Tau=log(sigma) uses the link: eta = log(exp(tau) - b)
  if (link[[2]] %in% okLinks[[2]]) { ## creating the logeb link
    stats[[2]] <- list()
    stats[[2]]$valideta <- function(eta) TRUE 
    stats[[2]]$link = link[[2]]
    stats[[2]]$linkfun <- eval(parse(text=paste("function(mu) log(exp(mu) - ",b,")", sep='')))
    stats[[2]]$linkinv <- eval(parse(text=paste("function(eta) log(exp(eta) +",b,")", sep='')))
    stats[[2]]$mu.eta <- eval(parse(text=
                                      paste("function(eta) { ee <- exp(eta); ee/(ee +",b,") }")))
    stats[[2]]$d2link <-  eval(parse(text=
                                       paste("function(mu) { em<-exp(mu); fr<-em/(em-",b,"); fr*(1-fr) }",sep='')))
    stats[[2]]$d3link <-  eval(parse(text=
                                       paste("function(mu) { em<-exp(mu); fr<-em/(em-",b,"); oo<-fr*(1-fr); oo-2*oo*fr }",sep='')))
    stats[[2]]$d4link <-  eval(parse(text=
                                       paste("function(mu) { em<-exp(mu); b<-",b,"; -b*em*(b^2+4*b*em+em^2)/(em-b)^4 }",sep='')))
  } else stop(link[[2]]," link not available for scale parameter of shash")
  
  # phi=log(delta) uses the link: eta = logit( (phi-a1)/a2 )
  if (link[[4]] %in% okLinks[[4]]) { ## creating the logeb link
    stats[[4]] <- list()
    stats[[4]]$valideta <- function(eta) TRUE 
    stats[[4]]$link = link[[4]]
    stats[[4]]$linkfun <- eval(parse(text=paste("function(mu) qlogis((mu - ",a1,")/", (a2-a1),")", sep='')))
    stats[[4]]$linkinv <- eval(parse(text=paste("function(eta) ", (a2-a1)," * plogis(eta) + ",a1, sep='')))
    stats[[4]]$mu.eta <- eval(parse(text=paste("function(eta) { oo <- plogis(eta); ", (a2-a1), "*oo*(1-oo) }", sep = '')))
    stats[[4]]$d2link <-  eval(parse(text=
                                       paste("function(mu) { z<-(mu-", a1,")/", (a2-a1),";  (1/(1-z)^2-1/z^2)/",(a2-a1),"^2}",sep='')))
    stats[[4]]$d3link <-  eval(parse(text=
                                       paste("function(mu) { z<-(mu-", a1,")/", (a2-a1),";  (2/(1-z)^3+2/z^3)/",(a2-a1),"^3}",sep='')))
    stats[[4]]$d4link <-  eval(parse(text=
                                       paste("function(mu) { z<-(mu-", a1,")/", (a2-a1),";  (6/(1-z)^4-6/z^4)/",(a2-a1),"^4}",sep='')))
  } else stop(link[[4]]," link not available for kurtosis parameter of shash")
  
  
  # variance <- function(mu) exp(get(".Theta"))  ##### XXX ##### Necessary?
  
  # validmu <- function(mu) all( is.finite(mu) )
  
  residuals <- function(object, type = c("deviance", "response")) {
    
    mu <-  object$fitted[ , 1]
    tau <- object$fitted[ , 2]
    eps <- object$fitted[ , 3]
    phi <- object$fitted[ , 4]
    
    sig <- exp( tau )
    del <- exp( phi )
    
    type <- match.arg(type)
    
    # raw residuals  
    rsd <- object$y - mu - sig*del*exp(0.25)*(besselK(0.25, nu = (1/del+1)/2)+besselK(0.25, nu = (1/del-1)/2))/sqrt(8*pi)
    
    if (type=="response"){ 
      return(rsd)
    }
    else { ## compute deviance residuals
      sgn <- sign(rsd)
      
      z <- (object$y - mu) / (sig*del)
      
      dTasMe <- del*asinh(z) - eps
      CC <- cosh( dTasMe )
      SS <- sinh( dTasMe )
      
      l <- -tau - 0.5*log(2*pi) + log(CC) - 0.5*log1p(z^2) - 0.5*SS^2
      
      # By putting ls to zero we are using only log-likelihood
      ls <- 0
      
      rsd <- pmax(0, 2*(ls - l))
      
      rsd <- sqrt(rsd)*sgn
    }
    rsd
  } ## residuals
  
  ll <- function(y, X, coef, wt, family, offset = NULL, deriv=0, d1b=0, d2b=0, Hp=NULL, rank=0, fh=NULL, D=NULL) {
    ## function defining the shash model log lik. 
    ## deriv: 0 - eval
    ##        1 - grad and Hess
    ##        2 - diagonal of first deriv of Hess
    ##        3 - first deriv of Hess
    ##        4 - everything.
    npar <- 4
    
    # If offset is not null or a vector of zeros, give an error
    if( !is.null(offset[[1]]) && sum(abs(offset)) )  stop("offset not still available for this family")
    
    jj <- attr(X, "lpi") ## extract linear predictor index
    
    eta <-  X[ , jj[[1]], drop=FALSE] %*% coef[jj[[1]]]
    eta1 <- X[ , jj[[2]], drop=FALSE] %*% coef[jj[[2]]]
    eta2 <- X[ , jj[[3]], drop=FALSE] %*% coef[jj[[3]]]
    eta3 <- X[ , jj[[4]], drop=FALSE] %*% coef[jj[[4]]]
    
    mu <-  family$linfo[[1]]$linkinv( eta )
    tau <- family$linfo[[2]]$linkinv( eta1 )
    eps <- family$linfo[[3]]$linkinv( eta2 )
    phi <- family$linfo[[4]]$linkinv( eta3 )
    
    sig <- exp( tau )
    del <- exp( phi )
    
    n <- length(y)
    
    z <- (y - mu) / (sig*del)
    
    dTasMe <- del*asinh(z) - eps
    g <- -dTasMe
    CC <- cosh( dTasMe )
    SS <- sinh( dTasMe )
    
    l <- sum( -tau - 0.5*log(2*pi) + log(CC) - 0.5*log1p(z^2) - 0.5*SS^2 )
    
    if (deriv>0) {
      
      zsd <- z*sig*del
      sSp1 <- sqrt(z^2+1)
      asinhZ <- asinh(z)
      
      ## First derivatives
      De <- tanh(g) - 0.5*sinh(2*g)
      Dm <- 1/(del*sig*sSp1)*(del*(De)+z/sSp1)
      Dt <- zsd*Dm - 1
      Dp <- Dt + 1 - del*asinhZ*De
      
      L1 <- cbind(Dm,Dt,De,Dp)
      
      ## the second derivatives  
      Dme <- (sech(g)^2 - cosh(2*g)) / (sig*sSp1)
      Dte <- zsd*Dme
      Dmm <- Dme/(sig*sSp1) + z*De/(sig^2*del*(z^2+1)^(3/2)) + (z^2-1)/(del*sig*del*sig*(z^2+1)^2)
      Dmt <- zsd*Dmm - Dm
      Dee <- -2*cosh(g)^2 + sech(g)^2 + 1 
      Dtt <-  zsd*Dmt
      Dep <- Dte - del*asinhZ*Dee
      Dmp <- Dmt + De/(sig*sSp1) - del*asinhZ*Dme
      Dtp <- zsd*Dmp
      Dpp <- Dtp - del*asinhZ*Dep + del*(z/sSp1-asinhZ)*De
      
      # Put them in matrix form
      L2 <- cbind(Dmm, Dmt, Dme, Dmp, Dtt, Dte ,Dtp ,Dee ,Dep ,Dpp)  
      
      ## need some link derivatives for derivative transform
      IG1 <- cbind(family$linfo[[1]]$mu.eta(eta), family$linfo[[2]]$mu.eta(eta1), 
                   family$linfo[[3]]$mu.eta(eta2), family$linfo[[4]]$mu.eta(eta3))
      G2 <- cbind(family$linfo[[1]]$d2link(mu), family$linfo[[2]]$d2link(tau), 
                  family$linfo[[3]]$d2link(eps), family$linfo[[4]]$d2link(phi))
    }
    
    L3 <- L4 <- G3 <- G4 <- 0 ## defaults
    
    if (deriv>1) {
      
      ## the third derivatives
      Deee <-  -2*(sinh(2*g)+sech(g)^2*tanh(g))
      Dmee <- Deee/(sig*sSp1)
      Dmme <- Dmee/(sig*sSp1) + z*Dee/(sig*sig*del*(z^2+1)^(3/2))
      Dmmm <- 2*z*Dme/(sig*sig*del*sSp1^3) + Dmme/(sig*sSp1) + 
        (2*z^2-1)*De/(sig^3*del^2*sSp1^5) + 2*z*(z^2-3)/((sig*del)^3*(z^2+1)^3)
      Dmmt <- zsd*Dmmm - 2*Dmm
      Dtee <- zsd*Dmee
      Dmte <- zsd*Dmme - Dme
      Dtte <- zsd*Dmte
      Dmtt <- zsd*Dmmt - Dmt
      Dttt <- zsd*Dmtt
      Dmep <- Dmte + Dee/(sig*sSp1) - del*asinhZ*Dmee
      Dtep <- zsd*Dmep
      Deep <- Dtee - del*asinhZ*Deee
      Depp <- Dtep - del*asinhZ*Deep + del*( z/sSp1-asinhZ )*Dee
      Dmmp <- Dmmt + 2*Dme/(sig*sSp1) + z*De/(del*sig*sig*sSp1^3) - del*asinhZ*Dmme
      Dmtp <- zsd*Dmmp - Dmp
      Dttp <- zsd*Dmtp
      Dmpp <- Dmtp + Dep/(sig*sSp1) + z^2*De/(sig*sSp1^3) - 
        del*asinhZ*Dmep + del*Dme*(z/sSp1 - asinhZ)
      Dtpp <- zsd*Dmpp
      Dppp <- Dtpp - del*asinhZ*Depp + del*(z/sSp1-asinhZ)*(2*Dep + De) + del*(z/sSp1)^3 * De
      
      ## Put them in matrix form
      L3 <- cbind(Dmmm,Dmmt,Dmme,Dmmp,Dmtt,Dmte,Dmtp,Dmee,Dmep,Dmpp,
                  Dttt,Dtte,Dttp,Dtee,Dtep,Dtpp,Deee,Deep,Depp,Dppp)
      
      G3 <- cbind(family$linfo[[1]]$d3link(mu), family$linfo[[2]]$d3link(tau), 
                  family$linfo[[3]]$d3link(eps), family$linfo[[4]]$d3link(phi))
    }
    
    if (deriv>3) {
      ## the fourth derivatives
      ## 35 4th derivatives:  mmmm,mmmt,mmme,mmmp,mmtt,mmte,mmtp,mmee,mmep,mmpp,
      ##                      mttt,mtte,mttp,mtee,mtep,mtpp,meee,meep,mepp,mppp,
      ##                      tttt,ttte,tttp,ttee,ttep,ttpp,teee,teep,tepp,tppp,
      ##	              eeee,eeep,eepp,eppp,pppp
      ## j2...r3
      ## The code for these is auto-generated, by auto-translation of simplified
      ## maxima expressions, furhter auto-simplification in R, and then some 
      ## further non-automatic simplification (which could be taken further) 
      m <- mu; t <- tau; p <- phi; e <- eps
      ## auto generated code...
      exp1 <- exp(1);
      aaa1 <- -t;
      aaa2 <- y-m;
      aaa3 <- exp1^p*asinh(exp1^(aaa1-p)*aaa2)-e; ## as abb7 and add1
      abb8 <- cosh(aaa3);
      abb9 <- sinh(aaa3);
      abb1 <- exp1^((-2*t)-2*p);
      abb3 <- aaa2^2;
      abb4 <- 1/exp1^t;
      abb5 <- -t-p;
      abb7 <- exp1^(2*abb5)*abb3+1
      abb6 <- 1/sqrt(abb7);
      aee5 <- aaa3 + e
      aff04 <- abb1*abb3+1;
      aff05 <- abb4^2
      aff08 <- 2*abb5;
      aff10 <- 1/abb7; ## abb6^2
      aff13 <- abb8^2; ## cosh(aaa3)^2
      aff14 <- exp1^(aaa1+aff08);
      aff15 <- abb6^3
      aff17 <- abb9^2; ##sinh(aaa3)^2
      agg15 <- 1/abb6
      agg17 <- 1/abb8;
      aii11 <- aaa3 + e
      aii12 <- aii11-abb4*aaa2*abb6;
      aii17 <- abb6^3
      ajj15 <- aaa2^3;
      ann05 <- exp1^p;
      ann06 <- asinh(exp1^abb5*aaa2);
      aoo09 <- -aaa2/(exp1^t*agg15);
      app02 <- -2*t;
      app04 <- exp1^(app02-2*p)*abb3+1;
      app08 <- exp1^(app02+aff08);
      app10 <- 1/abb7^2;
      app14 <- exp1^(aaa1+4*abb5);
      app16 <- 1/agg15^5;
      app21 <- 1/exp1^(3*t);
      aqq03 <- exp1^(app02-2*p);
      aqq05 <- aqq03*abb3+1;
      aqq27 <- 1/aff13;
      arr06 <- exp1^aff08*aaa2^2+1;
      arr07 <- 1/sqrt(arr06)^3;
      arr12 <- 1/arr06;
      ass16 <- aii11-aaa2/(exp1^t*agg15);
      ass23 <- 1/abb8;
      ass28 <- 1/aff13;
      att19 <- aaa2^4;
      avv19 <- aii11-abb4*aaa2*abb6;
      ayy14 <- -abb4*aaa2*abb6;
      ayy16 <- aii11+ayy14;
      ayy17 <- aii11+ayy14-aff14*ajj15*aii17;
      ayy24 <- ayy16^2;
      azz19 <- aaa2^5;
      bdd07 <- sqrt(exp1^aff08*aaa2^2+1);
      bdd08 <- 1/bdd07^3;
      bdd14 <- 1/bdd07;
      bdd15 <- aii11-abb4*aaa2*bdd14;
      bgg4 <- aee5-aaa2/(exp1^t*sqrt(exp1^(2*abb5)*aaa2^2+1));
      bhh13 <- -abb4*aaa2*bdd14;
      bhh14 <- ann05*ann06; ## aaa3 + e
      bii11 <- aii11+aoo09;
      bii15 <- aii11+aoo09-aff14*ajj15*aii17;
      
      
      bjj07 <- 4*abb5;
      bjj08 <- exp1^(app02+bjj07);
      bjj11 <- 1/abb7^3;
      bjj14 <- 1/exp1^(4*t);
      bjj18 <- exp1^(aaa1+6*abb5);
      bjj21 <- 1/agg15^7;
      bjj24 <- exp1^(aff08-3*t);
      bjj26 <- exp1^(aaa1+bjj07);
      j2  <-  (-(6*bjj14*app10*abb9^4)/abb8^4)-(12*bjj24*aaa2*app16*abb9^3)/abb8^3+8*bjj14*app10*aqq27*aff17+
        4*app08*app10*aqq27*aff17-15*bjj08*abb3*bjj11*aqq27*aff17-4*bjj14*app10*aff17+4*app08*app10*aff17-
        15*bjj08*abb3*bjj11*aff17-9*bjj26*aaa2*app16*abb8*abb9+24*bjj24*aaa2*app16*abb8*abb9+
        15*bjj18*ajj15*bjj21*abb8*abb9+9*bjj26*aaa2*app16*agg17*abb9+12*bjj24*aaa2*app16*agg17*abb9-
        15*bjj18*ajj15*bjj21*agg17*abb9-4*bjj14*app10*aff13+4*app08*app10*aff13-15*bjj08*abb3*bjj11*aff13-
        2*bjj14*app10-4*app08*app10+15*bjj08*abb3*bjj11+(6*exp1^((-4*t)-4*p))/app04^2-(48*exp1^((-6*t)-6*p)*abb3)/app04^3+
        (48*exp1^((-8*t)-8*p)*aaa2^4)/app04^4;
      bkk33 <- 1/abb8^3;
      bkk34 <- abb9^3;
      k2  <-  (-(6*bjj14*aaa2*app10*abb9^4)/abb8^4)+6*app21*aff15*bkk33*bkk34-12*bjj24*abb3*app16*bkk33*bkk34+
        8*bjj14*aaa2*app10*aqq27*aff17+13*app08*aaa2*app10*aqq27*aff17-15*bjj08*ajj15*bjj11*aqq27*aff17-
        4*bjj14*aaa2*app10*aff17+13*app08*aaa2*app10*aff17-15*bjj08*ajj15*bjj11*aff17-
        12*app21*aff15*abb8*abb9+3*aff14*aff15*abb8*abb9-18*bjj26*abb3*app16*abb8*abb9+24*bjj24*abb3*app16*abb8*abb9+
        15*bjj18*att19*bjj21*abb8*abb9-6*app21*aff15*agg17*abb9-3*aff14*aff15*agg17*abb9+18*bjj26*abb3*app16*agg17*abb9+
        12*bjj24*abb3*app16*agg17*abb9-15*bjj18*att19*bjj21*agg17*abb9-4*bjj14*aaa2*app10*aff13+13*app08*aaa2*app10*aff13-
        15*bjj08*ajj15*bjj11*aff13-2*bjj14*aaa2*app10-13*app08*aaa2*app10+15*bjj08*ajj15*bjj11+
        (24*exp1^((-4*t)-4*p)*aaa2)/app04^2-(72*exp1^((-6*t)-6*p)*ajj15)/app04^3+(48*exp1^((-8*t)-8*p)*aaa2^5)/app04^4;
      bll16 <- exp1^(aff08-2*t);
      l2  <-  (-(6*app21*aff15*abb9^4)/abb8^4)-(6*bll16*aaa2*app10*abb9^3)/abb8^3+8*app21*aff15*aqq27*aff17+
        aff14*aff15*aqq27*aff17-3*app14*abb3*app16*aqq27*aff17-4*app21*aff15*aff17+aff14*aff15*aff17-3*app14*abb3*app16*aff17+
        12*bll16*aaa2*app10*abb8*abb9+(6*bll16*aaa2*app10*abb9)/abb8-4*app21*aff15*aff13+aff14*aff15*aff13-
        3*app14*abb3*app16*aff13-2*app21*aff15-aff14*aff15+3*app14*abb3*app16;
      bmm34 <- 1/abb8^3;
      bmm35 <- abb9^3;
      m2  <-  (6*app21*aff15*ass16*abb9^4)/abb8^4+6*app08*aaa2*app10*ass16*bmm34*bmm35-6*bjj24*abb3*app16*bmm34*bmm35-
        8*app21*aff15*ass16*ass28*aff17-aff14*aff15*ass16*ass28*aff17+3*bjj26*abb3*app16*ass16*ass28*aff17+
        6*app08*aaa2*app10*ass28*aff17-12*bjj08*ajj15*bjj11*ass28*aff17+4*app21*aff15*ass16*aff17-aff14*aff15*ass16*aff17+
        3*bjj26*abb3*app16*ass16*aff17+6*app08*aaa2*app10*aff17-12*bjj08*ajj15*bjj11*aff17-12*app08*aaa2*app10*ass16*abb8*abb9+
        2*aff14*aff15*abb8*abb9-15*bjj26*abb3*app16*abb8*abb9+12*bjj24*abb3*app16*abb8*abb9+15*bjj18*att19*bjj21*abb8*abb9-
        6*app08*aaa2*app10*ass16*ass23*abb9-2*aff14*aff15*ass23*abb9+15*bjj26*abb3*app16*ass23*abb9+6*bjj24*abb3*app16*ass23*abb9-
        15*bjj18*att19*bjj21*ass23*abb9+4*app21*aff15*ass16*aff13-aff14*aff15*ass16*aff13+3*bjj26*abb3*app16*ass16*aff13+
        6*app08*aaa2*app10*aff13-12*bjj08*ajj15*bjj11*aff13+2*app21*aff15*ass16+aff14*aff15*ass16-3*bjj26*abb3*app16*ass16-
        6*app08*aaa2*app10+12*bjj08*ajj15*bjj11+(24*exp1^((-4*t)-4*p)*aaa2)/app04^2-(72*exp1^((-6*t)-6*p)*ajj15)/app04^3+
        (48*exp1^((-8*t)-8*p)*aaa2^5)/app04^4;
      n2  <-  (-(6*bjj14*abb3*app10*abb9^4)/abb8^4)+10*app21*aaa2*aff15*bkk33*bkk34-12*bjj24*ajj15*app16*bkk33*bkk34-
        4*aff05*aff10*aqq27*aff17+8*bjj14*abb3*app10*aqq27*aff17+19*app08*abb3*app10*aqq27*aff17-15*bjj08*att19*bjj11*aqq27*aff17-
        4*aff05*aff10*aff17-4*bjj14*abb3*app10*aff17+19*app08*abb3*app10*aff17-15*bjj08*att19*bjj11*aff17-
        20*app21*aaa2*aff15*abb8*abb9+9*aff14*aaa2*aff15*abb8*abb9-24*bjj26*ajj15*app16*abb8*abb9+24*bjj24*ajj15*app16*abb8*abb9+
        15*bjj18*azz19*bjj21*abb8*abb9-10*app21*aaa2*aff15*agg17*abb9-9*aff14*aaa2*aff15*agg17*abb9+24*bjj26*ajj15*app16*agg17*abb9+
        12*bjj24*ajj15*app16*agg17*abb9-15*bjj18*azz19*bjj21*agg17*abb9-4*aff05*aff10*aff13-4*bjj14*abb3*app10*aff13+
        19*app08*abb3*app10*aff13-15*bjj08*att19*bjj11*aff13+4*aff05*aff10-2*bjj14*abb3*app10-19*app08*abb3*app10+
        15*bjj08*att19*bjj11-(4*aqq03)/aqq05+(44*exp1^((-4*t)-4*p)*abb3)/aqq05^2-(88*exp1^((-6*t)-6*p)*att19)/aqq05^3+
        (48*exp1^((-8*t)-8*p)*aaa2^6)/aqq05^4;
      o2  <-  (-(6*app21*aaa2*aff15*abb9^4)/abb8^4)+4*aff05*aff10*bkk33*bkk34-6*bll16*abb3*app10*bkk33*bkk34+
        8*app21*aaa2*aff15*aqq27*aff17+3*aff14*aaa2*aff15*aqq27*aff17-3*app14*ajj15*app16*aqq27*aff17-
        4*app21*aaa2*aff15*aff17+3*aff14*aaa2*aff15*aff17-3*app14*ajj15*app16*aff17-8*aff05*aff10*abb8*abb9+
        12*bll16*abb3*app10*abb8*abb9-4*aff05*aff10*agg17*abb9+6*bll16*abb3*app10*agg17*abb9-4*app21*aaa2*aff15*aff13+
        3*aff14*aaa2*aff15*aff13-3*app14*ajj15*app16*aff13-2*app21*aaa2*aff15-3*aff14*aaa2*aff15+3*app14*ajj15*app16;
      p2  <-  (6*app21*aaa2*aff15*ass16*abb9^4)/abb8^4-4*aff05*aff10*ass16*bmm34*bmm35+6*app08*abb3*app10*ass16*bmm34*bmm35-
        6*bjj24*ajj15*app16*bmm34*bmm35-8*app21*aaa2*aff15*ass16*ass28*aff17-3*aff14*aaa2*aff15*ass16*ass28*aff17+
        3*bjj26*ajj15*app16*ass16*ass28*aff17+10*app08*abb3*app10*ass28*aff17-12*bjj08*att19*bjj11*ass28*aff17+
        4*app21*aaa2*aff15*ass16*aff17-3*aff14*aaa2*aff15*ass16*aff17+3*bjj26*ajj15*app16*ass16*aff17+10*app08*abb3*app10*aff17-
        12*bjj08*att19*bjj11*aff17+8*aff05*aff10*ass16*abb8*abb9-12*app08*abb3*app10*ass16*abb8*abb9+
        6*aff14*aaa2*aff15*abb8*abb9-21*bjj26*ajj15*app16*abb8*abb9+12*bjj24*ajj15*app16*abb8*abb9+15*bjj18*azz19*bjj21*abb8*abb9+
        4*aff05*aff10*ass16*ass23*abb9-6*app08*abb3*app10*ass16*ass23*abb9-
        6*aff14*aaa2*aff15*ass23*abb9+21*bjj26*ajj15*app16*ass23*abb9+6*bjj24*ajj15*app16*ass23*abb9-
        15*bjj18*azz19*bjj21*ass23*abb9+4*app21*aaa2*aff15*ass16*aff13-3*aff14*aaa2*aff15*ass16*aff13+
        3*bjj26*ajj15*app16*ass16*aff13+10*app08*abb3*app10*aff13-12*bjj08*att19*bjj11*aff13+2*app21*aaa2*aff15*ass16+
        3*aff14*aaa2*aff15*ass16-3*bjj26*ajj15*app16*ass16-10*app08*abb3*app10+12*bjj08*att19*bjj11-(4*aqq03)/aqq05+
        (44*exp1^((-4*t)-4*p)*abb3)/aqq05^2-(88*exp1^((-6*t)-6*p)*att19)/aqq05^3+(48*exp1^((-8*t)-8*p)*aaa2^6)/aqq05^4;
      q2  <-  (-(6*aff05*arr12*abb9^4)/abb8^4)-(2*aff14*aaa2*arr07*abb9^3)/abb8^3+(8*aff05*arr12*aff17)/aff13-
        4*aff05*arr12*aff17+4*aff14*aaa2*arr07*abb8*abb9+(2*aff14*aaa2*arr07*abb9)/abb8-4*aff05*arr12*aff13-2*aff05*arr12;
      r2  <-  (6*aff05*aff10*ass16*abb9^4)/abb8^4+2*aff14*aaa2*aff15*ass16*bmm34*bmm35-4*bll16*abb3*app10*bmm34*bmm35-
        8*aff05*aff10*ass16*ass28*aff17+2*aff14*aaa2*aff15*ass28*aff17-3*app14*ajj15*app16*ass28*aff17+
        4*aff05*aff10*ass16*aff17+2*aff14*aaa2*aff15*aff17-3*app14*ajj15*app16*aff17-4*aff14*aaa2*aff15*ass16*abb8*abb9+
        8*bll16*abb3*app10*abb8*abb9-2*aff14*aaa2*aff15*ass16*ass23*abb9+4*bll16*abb3*app10*ass23*abb9+
        4*aff05*aff10*ass16*aff13+2*aff14*aaa2*aff15*aff13-3*app14*ajj15*app16*aff13+2*aff05*aff10*ass16-
        2*aff14*aaa2*aff15+3*app14*ajj15*app16;
      bss21 <- 2*aff14*abb3*aff15-3*bjj26*att19*app16;
      bss23 <- -abb4*aaa2*abb6;
      bss25 <- aii11+bss23;
      bss26 <- aii11+bss23-aff14*ajj15*aff15;
      bss29 <- bss25^2;
      bss33 <- (-4*aff14*aaa2*aff15)+18*bjj26*ajj15*app16-(15*exp1^(aaa1+6*abb5)*aaa2^5)/agg15^7;
      s2  <-  (-(6*aff05*aff10*bss29*abb9^4)/abb8^4)-2*aff14*aaa2*aff15*bss29*bmm34*bmm35+2*aff05*aff10*bss26*bmm34*bmm35+
        8*app08*abb3*app10*bss25*bmm34*bmm35+8*aff05*aff10*bss29*ass28*aff17+aff14*aaa2*aff15*bss26*ass28*aff17-
        4*aff14*aaa2*aff15*bss25*ass28*aff17+6*bjj26*ajj15*app16*bss25*ass28*aff17+2*abb4*abb6*bss21*ass28*aff17-
        2*bjj08*att19*bjj11*ass28*aff17-4*aff05*aff10*bss29*aff17+aff14*aaa2*aff15*bss26*aff17-4*aff14*aaa2*aff15*bss25*aff17+
        6*bjj26*ajj15*app16*bss25*aff17+2*abb4*abb6*bss21*aff17-2*bjj08*att19*bjj11*aff17+4*aff14*aaa2*aff15*bss29*abb8*abb9-
        4*aff05*aff10*bss26*abb8*abb9-16*app08*abb3*app10*bss25*abb8*abb9-bss33*abb8*abb9+2*aff14*aaa2*aff15*bss29*ass23*abb9-
        2*aff05*aff10*bss26*ass23*abb9-8*app08*abb3*app10*bss25*ass23*abb9+bss33*ass23*abb9-4*aff05*aff10*bss29*aff13+
        aff14*aaa2*aff15*bss26*aff13-4*aff14*aaa2*aff15*bss25*aff13+6*bjj26*ajj15*app16*bss25*aff13+2*abb4*abb6*bss21*aff13-
        2*bjj08*att19*bjj11*aff13-2*aff05*aff10*bss29-aff14*aaa2*aff15*bss26+4*aff14*aaa2*aff15*bss25-6*bjj26*ajj15*app16*bss25-
        2*abb4*abb6*bss21+2*bjj08*att19*bjj11-(4*aqq03)/aqq05+(44*exp1^((-4*t)-4*p)*abb3)/aqq05^2-
        (88*exp1^((-6*t)-6*p)*att19)/aqq05^3+(48*exp1^((-8*t)-8*p)*aaa2^6)/aqq05^4;
      btt24 <- aaa2^6;
      t2  <-  (-(6*bjj14*ajj15*app10*abb9^4)/abb8^4)+12*app21*abb3*aff15*bkk33*bkk34-12*bjj24*att19*app16*bkk33*bkk34-
        7*aff05*aaa2*aff10*aqq27*aff17+8*bjj14*ajj15*app10*aqq27*aff17+22*app08*ajj15*app10*aqq27*aff17-
        15*bjj08*azz19*bjj11*aqq27*aff17-7*aff05*aaa2*aff10*aff17-4*bjj14*ajj15*app10*aff17+22*app08*ajj15*app10*aff17-
        15*bjj08*azz19*bjj11*aff17-abb4*abb6*abb8*abb9-24*app21*abb3*aff15*abb8*abb9+13*aff14*abb3*aff15*abb8*abb9-
        27*bjj26*att19*app16*abb8*abb9+24*bjj24*att19*app16*abb8*abb9+15*bjj18*btt24*bjj21*abb8*abb9+abb4*abb6*agg17*abb9-
        12*app21*abb3*aff15*agg17*abb9-13*aff14*abb3*aff15*agg17*abb9+27*bjj26*att19*app16*agg17*abb9+
        12*bjj24*att19*app16*agg17*abb9-15*bjj18*btt24*bjj21*agg17*abb9-7*aff05*aaa2*aff10*aff13-4*bjj14*ajj15*app10*aff13+
        22*app08*ajj15*app10*aff13-15*bjj08*azz19*bjj11*aff13+7*aff05*aaa2*aff10-2*bjj14*ajj15*app10-22*app08*ajj15*app10+
        15*bjj08*azz19*bjj11-(8*aqq03*aaa2)/aqq05+(56*exp1^((-4*t)-4*p)*ajj15)/aqq05^2-(96*exp1^((-6*t)-6*p)*azz19)/aqq05^3+
        (48*exp1^((-8*t)-8*p)*aaa2^7)/aqq05^4;
      u2  <-  (-(6*app21*abb3*aff15*abb9^4)/abb8^4)+6*aff05*aaa2*aff10*bkk33*bkk34-6*bll16*ajj15*app10*bkk33*bkk34-
        abb4*abb6*aqq27*aff17+8*app21*abb3*aff15*aqq27*aff17+4*aff14*abb3*aff15*aqq27*aff17-3*app14*att19*app16*aqq27*aff17-
        abb4*abb6*aff17-4*app21*abb3*aff15*aff17+4*aff14*abb3*aff15*aff17-3*app14*att19*app16*aff17-
        12*aff05*aaa2*aff10*abb8*abb9+12*bll16*ajj15*app10*abb8*abb9-6*aff05*aaa2*aff10*agg17*abb9+
        6*bll16*ajj15*app10*agg17*abb9-abb4*abb6*aff13-4*app21*abb3*aff15*aff13+4*aff14*abb3*aff15*aff13-
        3*app14*att19*app16*aff13+abb4*abb6-2*app21*abb3*aff15-4*aff14*abb3*aff15+3*app14*att19*app16;
      v2  <-  (6*app21*abb3*aff15*avv19*abb9^4)/abb8^4-6*aff05*aaa2*aff10*avv19*bmm34*bmm35+
        6*app08*ajj15*app10*avv19*bmm34*bmm35-6*bjj24*att19*app16*bmm34*bmm35+abb4*abb6*avv19*ass28*aff17-
        8*app21*abb3*aff15*avv19*ass28*aff17-4*aff14*abb3*aff15*avv19*ass28*aff17+3*bjj26*att19*app16*avv19*ass28*aff17+
        12*app08*ajj15*app10*ass28*aff17-12*bjj08*azz19*bjj11*ass28*aff17+abb4*abb6*avv19*aff17+4*app21*abb3*aff15*avv19*aff17-
        4*aff14*abb3*aff15*avv19*aff17+3*bjj26*att19*app16*avv19*aff17+12*app08*ajj15*app10*aff17-12*bjj08*azz19*bjj11*aff17+
        12*aff05*aaa2*aff10*avv19*abb8*abb9-12*app08*ajj15*app10*avv19*abb8*abb9+9*aff14*abb3*aff15*abb8*abb9-
        24*bjj26*att19*app16*abb8*abb9+12*bjj24*att19*app16*abb8*abb9+15*bjj18*btt24*bjj21*abb8*abb9+
        6*aff05*aaa2*aff10*avv19*ass23*abb9-6*app08*ajj15*app10*avv19*ass23*abb9-9*aff14*abb3*aff15*ass23*abb9+
        24*bjj26*att19*app16*ass23*abb9+6*bjj24*att19*app16*ass23*abb9-15*bjj18*btt24*bjj21*ass23*abb9+abb4*abb6*avv19*aff13+
        4*app21*abb3*aff15*avv19*aff13-4*aff14*abb3*aff15*avv19*aff13+3*bjj26*att19*app16*avv19*aff13+12*app08*ajj15*app10*aff13-
        12*bjj08*azz19*bjj11*aff13-abb4*abb6*avv19+2*app21*abb3*aff15*avv19+4*aff14*abb3*aff15*avv19-3*bjj26*att19*app16*avv19-
        12*app08*ajj15*app10+12*bjj08*azz19*bjj11-(8*aqq03*aaa2)/aqq05+(56*exp1^((-4*t)-4*p)*ajj15)/aqq05^2-
        (96*exp1^((-6*t)-6*p)*azz19)/aqq05^3+(48*exp1^((-8*t)-8*p)*aaa2^7)/aqq05^4;
      w2  <-  (-(6*aff05*aaa2*aff10*abb9^4)/abb8^4)+2*abb4*abb6*bkk33*bkk34-2*aff14*abb3*aff15*bkk33*bkk34+
        (8*aff05*aaa2*aff10*aff17)/aff13-4*aff05*aaa2*aff10*aff17-4*abb4*abb6*abb8*abb9+4*aff14*abb3*aff15*abb8*abb9-
        2*abb4*abb6*agg17*abb9+2*aff14*abb3*aff15*agg17*abb9-4*aff05*aaa2*aff10*aff13-2*aff05*aaa2*aff10;
      x2  <-  (6*aff05*aaa2*aff10*avv19*abb9^4)/abb8^4-2*abb4*abb6*avv19*bmm34*bmm35+2*aff14*abb3*aff15*avv19*bmm34*bmm35-
        4*bll16*ajj15*app10*bmm34*bmm35-8*aff05*aaa2*aff10*avv19*ass28*aff17+3*aff14*abb3*aff15*ass28*aff17-
        3*app14*att19*app16*ass28*aff17+4*aff05*aaa2*aff10*avv19*aff17+3*aff14*abb3*aff15*aff17-3*app14*att19*app16*aff17+
        4*abb4*abb6*avv19*abb8*abb9-4*aff14*abb3*aff15*avv19*abb8*abb9+8*bll16*ajj15*app10*abb8*abb9+2*abb4*abb6*avv19*ass23*abb9-
        2*aff14*abb3*aff15*avv19*ass23*abb9+4*bll16*ajj15*app10*ass23*abb9+4*aff05*aaa2*aff10*avv19*aff13+
        3*aff14*abb3*aff15*aff13-3*app14*att19*app16*aff13+2*aff05*aaa2*aff10*avv19-3*aff14*abb3*aff15+3*app14*att19*app16;
      byy24 <- 2*aff14*ajj15*aff15-3*bjj26*azz19*app16;
      byy35 <- (-6*aff14*abb3*aff15)+21*bjj26*att19*app16-(15*exp1^(aaa1+6*abb5)*aaa2^6)/agg15^7;
      y2  <-  (-(6*aff05*aaa2*aff10*bss29*abb9^4)/abb8^4)+2*abb4*abb6*bss29*bmm34*bmm35-2*aff14*abb3*aff15*bss29*bmm34*bmm35+
        2*aff05*aaa2*aff10*bss26*bmm34*bmm35+8*app08*ajj15*app10*bss25*bmm34*bmm35+8*aff05*aaa2*aff10*bss29*ass28*aff17-
        abb4*abb6*bss26*ass28*aff17+aff14*abb3*aff15*bss26*ass28*aff17-6*aff14*abb3*aff15*bss25*ass28*aff17+
        6*bjj26*att19*app16*bss25*ass28*aff17+abb4*abb6*byy24*ass28*aff17+abb4*aaa2*abb6*bss21*ass28*aff17-
        2*bjj08*azz19*bjj11*ass28*aff17-4*aff05*aaa2*aff10*bss29*aff17-abb4*abb6*bss26*aff17+aff14*abb3*aff15*bss26*aff17-
        6*aff14*abb3*aff15*bss25*aff17+6*bjj26*att19*app16*bss25*aff17+abb4*abb6*byy24*aff17+abb4*aaa2*abb6*bss21*aff17-
        2*bjj08*azz19*bjj11*aff17-4*abb4*abb6*bss29*abb8*abb9+4*aff14*abb3*aff15*bss29*abb8*abb9-
        4*aff05*aaa2*aff10*bss26*abb8*abb9-16*app08*ajj15*app10*bss25*abb8*abb9-byy35*abb8*abb9-2*abb4*abb6*bss29*ass23*abb9+
        2*aff14*abb3*aff15*bss29*ass23*abb9-2*aff05*aaa2*aff10*bss26*ass23*abb9-8*app08*ajj15*app10*bss25*ass23*abb9+
        byy35*ass23*abb9-4*aff05*aaa2*aff10*bss29*aff13-abb4*abb6*bss26*aff13+aff14*abb3*aff15*bss26*aff13-
        6*aff14*abb3*aff15*bss25*aff13+6*bjj26*att19*app16*bss25*aff13+abb4*abb6*byy24*aff13+abb4*aaa2*abb6*bss21*aff13-
        2*bjj08*azz19*bjj11*aff13-2*aff05*aaa2*aff10*bss29+abb4*abb6*bss26-aff14*abb3*aff15*bss26+6*aff14*abb3*aff15*bss25-
        6*bjj26*att19*app16*bss25-abb4*abb6*byy24-abb4*aaa2*abb6*bss21+2*bjj08*azz19*bjj11-(8*aqq03*aaa2)/aqq05+
        (56*exp1^((-4*t)-4*p)*ajj15)/aqq05^2-(96*exp1^((-6*t)-6*p)*azz19)/aqq05^3+(48*exp1^((-8*t)-8*p)*aaa2^7)/aqq05^4;
      bzz7 <- abb8^2;
      bzz9 <- abb9^2;
      z2  <-  (-(6*abb4*abb6*abb9^4)/abb8^4)+(8*abb4*abb6*bzz9)/bzz7-4*abb4*abb6*bzz9-4*abb4*abb6*bzz7-2*abb4*abb6;
      a3  <-  (6*abb4*abb6*aii12*abb9^4)/abb8^4-(2*aff14*abb3*aii17*abb9^3)/abb8^3-(8*abb4*abb6*aii12*aff17)/aff13+
        4*abb4*abb6*aii12*aff17+4*aff14*abb3*aii17*abb8*abb9+(2*aff14*abb3*aii17*abb9)/abb8+4*abb4*abb6*aii12*aff13+
        2*abb4*abb6*aii12;
      cbb09 <- 1/agg15^5;
      cbb18 <- 2*aff14*abb3*aii17-3*app14*att19*cbb09;
      cbb24 <- aii11+ayy14-aff14*aaa2^3*aii17;
      b3  <-  (-(6*abb4*abb6*ayy24*abb9^4)/abb8^4)+2*abb4*abb6*cbb24*bmm34*bmm35+4*aff14*abb3*aii17*ayy16*bmm34*bmm35+
        8*abb4*abb6*ayy24*ass28*aff17+cbb18*ass28*aff17-4*abb4*abb6*ayy24*aff17+cbb18*aff17-4*abb4*abb6*cbb24*abb8*abb9-
        8*aff14*abb3*aii17*ayy16*abb8*abb9-2*abb4*abb6*cbb24*ass23*abb9-4*aff14*abb3*aii17*ayy16*ass23*abb9-
        4*abb4*abb6*ayy24*aff13+cbb18*aff13-2*abb4*abb6*ayy24-2*aff14*abb3*aii17+3*app14*att19*cbb09;
      ccc23 <- aii11+ayy14+aff14*ajj15*aii17-3*app14*azz19*cbb09;
      ccc24 <- ayy16^3;
      ccc28 <- (-4*aff14*abb3*aii17)+18*app14*att19*cbb09-(15*exp1^(aaa1+6*abb5)*aaa2^6)/agg15^7;
      c3  <-  (6*abb4*abb6*ccc24*abb9^4)/abb8^4-6*aff14*abb3*aii17*ayy24*bmm34*bmm35-6*abb4*abb6*ayy16*ayy17*bmm34*bmm35-
        8*abb4*abb6*ccc24*ass28*aff17+abb4*abb6*ccc23*ass28*aff17+3*aff14*abb3*aii17*ayy17*ass28*aff17-
        3*cbb18*ayy16*ass28*aff17+4*abb4*abb6*ccc24*aff17+abb4*abb6*ccc23*aff17+3*aff14*abb3*aii17*ayy17*aff17-
        3*cbb18*ayy16*aff17+12*aff14*abb3*aii17*ayy24*abb8*abb9+12*abb4*abb6*ayy16*ayy17*abb8*abb9-ccc28*abb8*abb9+
        6*aff14*abb3*aii17*ayy24*ass23*abb9+6*abb4*abb6*ayy16*ayy17*ass23*abb9+ccc28*ass23*abb9+4*abb4*abb6*ccc24*aff13+
        abb4*abb6*ccc23*aff13+3*aff14*abb3*aii17*ayy17*aff13-3*cbb18*ayy16*aff13+2*abb4*abb6*ccc24-abb4*abb6*ccc23-
        3*aff14*abb3*aii17*ayy17+3*cbb18*ayy16-(8*abb1*aaa2)/aff04+(56*exp1^((-4*t)-4*p)*ajj15)/aff04^2-
        (96*exp1^((-6*t)-6*p)*azz19)/aff04^3+(48*exp1^((-8*t)-8*p)*aaa2^7)/aff04^4;
      cdd24 <- aaa2^7;
      d3  <-  (-(6*bjj14*att19*app10*abb9^4)/abb8^4)+12*app21*ajj15*aff15*bkk33*bkk34-12*bjj24*azz19*app16*bkk33*bkk34-
        7*aff05*abb3*aff10*aqq27*aff17+8*bjj14*att19*app10*aqq27*aff17+22*app08*att19*app10*aqq27*aff17-
        15*bjj08*btt24*bjj11*aqq27*aff17-7*aff05*abb3*aff10*aff17-4*bjj14*att19*app10*aff17+22*app08*att19*app10*aff17-
        15*bjj08*btt24*bjj11*aff17-abb4*aaa2*abb6*abb8*abb9-24*app21*ajj15*aff15*abb8*abb9+13*aff14*ajj15*aff15*abb8*abb9-
        27*bjj26*azz19*app16*abb8*abb9+24*bjj24*azz19*app16*abb8*abb9+15*bjj18*cdd24*bjj21*abb8*abb9+abb4*aaa2*abb6*agg17*abb9-
        12*app21*ajj15*aff15*agg17*abb9-13*aff14*ajj15*aff15*agg17*abb9+27*bjj26*azz19*app16*agg17*abb9+
        12*bjj24*azz19*app16*agg17*abb9-15*bjj18*cdd24*bjj21*agg17*abb9-7*aff05*abb3*aff10*aff13-4*bjj14*att19*app10*aff13+
        22*app08*att19*app10*aff13-15*bjj08*btt24*bjj11*aff13+7*aff05*abb3*aff10-2*bjj14*att19*app10-22*app08*att19*app10+
        15*bjj08*btt24*bjj11-(8*aqq03*abb3)/aqq05+(56*exp1^((-4*t)-4*p)*att19)/aqq05^2-(96*exp1^((-6*t)-6*p)*btt24)/aqq05^3+
        (48*exp1^((-8*t)-8*p)*aaa2^8)/aqq05^4;
      e3  <-  (-(6*app21*ajj15*aff15*abb9^4)/abb8^4)+6*aff05*abb3*aff10*bkk33*bkk34-6*bll16*att19*app10*bkk33*bkk34-
        abb4*aaa2*abb6*aqq27*aff17+8*app21*ajj15*aff15*aqq27*aff17+4*aff14*ajj15*aff15*aqq27*aff17-
        3*app14*azz19*app16*aqq27*aff17-abb4*aaa2*abb6*aff17-4*app21*ajj15*aff15*aff17+4*aff14*ajj15*aff15*aff17-
        3*app14*azz19*app16*aff17-12*aff05*abb3*aff10*abb8*abb9+12*bll16*att19*app10*abb8*abb9-6*aff05*abb3*aff10*agg17*abb9+
        6*bll16*att19*app10*agg17*abb9-abb4*aaa2*abb6*aff13-4*app21*ajj15*aff15*aff13+4*aff14*ajj15*aff15*aff13-
        3*app14*azz19*app16*aff13+abb4*aaa2*abb6-2*app21*ajj15*aff15-4*aff14*ajj15*aff15+3*app14*azz19*app16;
      f3  <-  (6*app21*ajj15*aff15*avv19*abb9^4)/abb8^4-6*aff05*abb3*aff10*avv19*bmm34*bmm35+
        6*app08*att19*app10*avv19*bmm34*bmm35-6*bjj24*azz19*app16*bmm34*bmm35+abb4*aaa2*abb6*avv19*ass28*aff17-
        8*app21*ajj15*aff15*avv19*ass28*aff17-4*aff14*ajj15*aff15*avv19*ass28*aff17+3*bjj26*azz19*app16*avv19*ass28*aff17+
        12*app08*att19*app10*ass28*aff17-12*bjj08*btt24*bjj11*ass28*aff17+abb4*aaa2*abb6*avv19*aff17+
        4*app21*ajj15*aff15*avv19*aff17-4*aff14*ajj15*aff15*avv19*aff17+3*bjj26*azz19*app16*avv19*aff17+
        12*app08*att19*app10*aff17-12*bjj08*btt24*bjj11*aff17+12*aff05*abb3*aff10*avv19*abb8*abb9-
        12*app08*att19*app10*avv19*abb8*abb9+9*aff14*ajj15*aff15*abb8*abb9-24*bjj26*azz19*app16*abb8*abb9+
        12*bjj24*azz19*app16*abb8*abb9+15*bjj18*cdd24*bjj21*abb8*abb9+6*aff05*abb3*aff10*avv19*ass23*abb9-
        6*app08*att19*app10*avv19*ass23*abb9-9*aff14*ajj15*aff15*ass23*abb9+24*bjj26*azz19*app16*ass23*abb9+
        6*bjj24*azz19*app16*ass23*abb9-15*bjj18*cdd24*bjj21*ass23*abb9+abb4*aaa2*abb6*avv19*aff13+
        4*app21*ajj15*aff15*avv19*aff13-4*aff14*ajj15*aff15*avv19*aff13+3*bjj26*azz19*app16*avv19*aff13+
        12*app08*att19*app10*aff13-12*bjj08*btt24*bjj11*aff13-abb4*aaa2*abb6*avv19+2*app21*ajj15*aff15*avv19+
        4*aff14*ajj15*aff15*avv19-3*bjj26*azz19*app16*avv19-12*app08*att19*app10+12*bjj08*btt24*bjj11-(8*aqq03*abb3)/aqq05+
        (56*exp1^((-4*t)-4*p)*att19)/aqq05^2-(96*exp1^((-6*t)-6*p)*btt24)/aqq05^3+(48*exp1^((-8*t)-8*p)*aaa2^8)/aqq05^4;
      g3  <-  (-(6*aff05*abb3*aff10*abb9^4)/abb8^4)+2*abb4*aaa2*abb6*bkk33*bkk34-2*aff14*ajj15*aff15*bkk33*bkk34+
        (8*aff05*abb3*aff10*aff17)/aff13-4*aff05*abb3*aff10*aff17-4*abb4*aaa2*abb6*abb8*abb9+4*aff14*ajj15*aff15*abb8*abb9-
        2*abb4*aaa2*abb6*agg17*abb9+2*aff14*ajj15*aff15*agg17*abb9-4*aff05*abb3*aff10*aff13-2*aff05*abb3*aff10;
      h3  <-  (6*aff05*abb3*aff10*avv19*abb9^4)/abb8^4-2*abb4*aaa2*abb6*avv19*bmm34*bmm35+2*aff14*ajj15*aff15*avv19*bmm34*bmm35-
        4*bll16*att19*app10*bmm34*bmm35-8*aff05*abb3*aff10*avv19*ass28*aff17+3*aff14*ajj15*aff15*ass28*aff17-
        3*app14*azz19*app16*ass28*aff17+4*aff05*abb3*aff10*avv19*aff17+3*aff14*ajj15*aff15*aff17-3*app14*azz19*app16*aff17+
        4*abb4*aaa2*abb6*avv19*abb8*abb9-4*aff14*ajj15*aff15*avv19*abb8*abb9+8*bll16*att19*app10*abb8*abb9+
        2*abb4*aaa2*abb6*avv19*ass23*abb9-2*aff14*ajj15*aff15*avv19*ass23*abb9+4*bll16*att19*app10*ass23*abb9+
        4*aff05*abb3*aff10*avv19*aff13+3*aff14*ajj15*aff15*aff13-3*app14*azz19*app16*aff13+2*aff05*abb3*aff10*avv19-
        3*aff14*ajj15*aff15+3*app14*azz19*app16;
      i3  <-  (-(6*aff05*abb3*aff10*bss29*abb9^4)/abb8^4)+2*abb4*aaa2*abb6*bss29*bmm34*bmm35-
        2*aff14*ajj15*aff15*bss29*bmm34*bmm35+2*aff05*abb3*aff10*bss26*bmm34*bmm35+8*app08*att19*app10*bss25*bmm34*bmm35+
        8*aff05*abb3*aff10*bss29*ass28*aff17-abb4*aaa2*abb6*bss26*ass28*aff17+aff14*ajj15*aff15*bss26*ass28*aff17-
        6*aff14*ajj15*aff15*bss25*ass28*aff17+6*bjj26*azz19*app16*bss25*ass28*aff17+4*app08*att19*app10*ass28*aff17-
        8*bjj08*btt24*bjj11*ass28*aff17-4*aff05*abb3*aff10*bss29*aff17-abb4*aaa2*abb6*bss26*aff17+
        aff14*ajj15*aff15*bss26*aff17-6*aff14*ajj15*aff15*bss25*aff17+6*bjj26*azz19*app16*bss25*aff17+
        4*app08*att19*app10*aff17-8*bjj08*btt24*bjj11*aff17-4*abb4*aaa2*abb6*bss29*abb8*abb9+
        4*aff14*ajj15*aff15*bss29*abb8*abb9-4*aff05*abb3*aff10*bss26*abb8*abb9-16*app08*att19*app10*bss25*abb8*abb9+
        6*aff14*ajj15*aff15*abb8*abb9-21*bjj26*azz19*app16*abb8*abb9+15*bjj18*cdd24*bjj21*abb8*abb9-
        2*abb4*aaa2*abb6*bss29*ass23*abb9+2*aff14*ajj15*aff15*bss29*ass23*abb9-2*aff05*abb3*aff10*bss26*ass23*abb9-
        8*app08*att19*app10*bss25*ass23*abb9-6*aff14*ajj15*aff15*ass23*abb9+21*bjj26*azz19*app16*ass23*abb9-
        15*bjj18*cdd24*bjj21*ass23*abb9-4*aff05*abb3*aff10*bss29*aff13-abb4*aaa2*abb6*bss26*aff13+
        aff14*ajj15*aff15*bss26*aff13-6*aff14*ajj15*aff15*bss25*aff13+6*bjj26*azz19*app16*bss25*aff13+
        4*app08*att19*app10*aff13-8*bjj08*btt24*bjj11*aff13-2*aff05*abb3*aff10*bss29+abb4*aaa2*abb6*bss26-
        aff14*ajj15*aff15*bss26+6*aff14*ajj15*aff15*bss25-6*bjj26*azz19*app16*bss25-4*app08*att19*app10+
        8*bjj08*btt24*bjj11-(8*aqq03*abb3)/aqq05+(56*exp1^((-4*t)-4*p)*att19)/aqq05^2-
        (96*exp1^((-6*t)-6*p)*btt24)/aqq05^3+(48*exp1^((-8*t)-8*p)*aaa2^8)/aqq05^4;
      j3  <-  (-(6*abb4*aaa2*abb6*abb9^4)/abb8^4)+(8*abb4*aaa2*abb6*bzz9)/bzz7-4*abb4*aaa2*abb6*bzz9-
        4*abb4*aaa2*abb6*bzz7-2*abb4*aaa2*abb6;
      k3  <-  (6*abb4*aaa2*bdd14*bdd15*abb9^4)/abb8^4-(2*aff14*ajj15*bdd08*abb9^3)/abb8^3-
        (8*abb4*aaa2*bdd14*bdd15*aff17)/aff13+4*abb4*aaa2*bdd14*bdd15*aff17+4*aff14*ajj15*bdd08*abb8*abb9+
        (2*aff14*ajj15*bdd08*abb9)/abb8+4*abb4*aaa2*bdd14*bdd15*aff13+2*abb4*aaa2*bdd14*bdd15;
      cll08 <- 1/bdd07^5;
      cll16 <- aii11+bhh13;
      cll17 <- cll16^2;
      cll18 <- 2*aff14*ajj15*bdd08-3*app14*azz19*cll08;
      cll24 <- aii11+bhh13-aff14*ajj15*bdd08;
      l3  <-  (-(6*abb4*aaa2*bdd14*cll17*abb9^4)/abb8^4)+2*abb4*aaa2*bdd14*cll24*bmm34*bmm35+
        4*aff14*ajj15*bdd08*cll16*bmm34*bmm35+8*abb4*aaa2*bdd14*cll17*ass28*aff17+cll18*ass28*aff17-
        4*abb4*aaa2*bdd14*cll17*aff17+cll18*aff17-4*abb4*aaa2*bdd14*cll24*abb8*abb9-8*aff14*ajj15*bdd08*cll16*abb8*abb9-
        2*abb4*aaa2*bdd14*cll24*ass23*abb9-4*aff14*ajj15*bdd08*cll16*ass23*abb9-4*abb4*aaa2*bdd14*cll17*aff13+
        cll18*aff13-2*abb4*aaa2*bdd14*cll17-2*aff14*ajj15*bdd08+3*app14*azz19*cll08;
      cmm12 <- -3*app14*azz19*cbb09;
      cmm16 <- 2*aff14*ajj15*aii17+cmm12;
      cmm23 <- aii11+ayy14+aff14*ajj15*aii17+cmm12;
      cmm28 <- (-4*aff14*ajj15*aii17)+18*app14*azz19*cbb09-(15*exp1^(aaa1+6*abb5)*aaa2^7)/agg15^7;
      m3  <-  (6*abb4*aaa2*abb6*ccc24*abb9^4)/abb8^4-6*aff14*ajj15*aii17*ayy24*bmm34*bmm35-
        6*abb4*aaa2*abb6*ayy16*ayy17*bmm34*bmm35-8*abb4*aaa2*abb6*ccc24*ass28*aff17+abb4*aaa2*abb6*cmm23*ass28*aff17+
        3*aff14*ajj15*aii17*ayy17*ass28*aff17-3*cmm16*ayy16*ass28*aff17+4*abb4*aaa2*abb6*ccc24*aff17+
        abb4*aaa2*abb6*cmm23*aff17+3*aff14*ajj15*aii17*ayy17*aff17-3*cmm16*ayy16*aff17+12*aff14*ajj15*aii17*ayy24*abb8*abb9+
        12*abb4*aaa2*abb6*ayy16*ayy17*abb8*abb9-cmm28*abb8*abb9+6*aff14*ajj15*aii17*ayy24*ass23*abb9+
        6*abb4*aaa2*abb6*ayy16*ayy17*ass23*abb9+cmm28*ass23*abb9+4*abb4*aaa2*abb6*ccc24*aff13+abb4*aaa2*abb6*cmm23*aff13+
        3*aff14*ajj15*aii17*ayy17*aff13-3*cmm16*ayy16*aff13+2*abb4*aaa2*abb6*ccc24-abb4*aaa2*abb6*cmm23-
        3*aff14*ajj15*aii17*ayy17+3*cmm16*ayy16-(8*abb1*abb3)/aff04+(56*exp1^((-4*t)-4*p)*aaa2^4)/aff04^2-
        (96*exp1^((-6*t)-6*p)*aaa2^6)/aff04^3+(48*exp1^((-8*t)-8*p)*aaa2^8)/aff04^4;
      cnn3 <- abb8^2;
      cnn5 <- abb9^2;
      n3  <-  (-(6*abb9^4)/abb8^4)+(8*cnn5)/cnn3-4*cnn5-4*cnn3-2;
      coo7 <- abb8^2;
      coo9 <- abb9^2;
      o3  <-  (6*bgg4*abb9^4)/abb8^4-(8*bgg4*coo9)/coo7+4*bgg4*coo9+4*bgg4*coo7+2*bgg4;
      cpp06 <- -aaa2/(exp1^t*bdd07);
      cpp08 <- (cpp06+aii11)^2;
      cpp12 <- aii11+cpp06-(exp1^(aaa1+aff08)*aaa2^3)/bdd07^3;
      p3  <-  (-(6*cpp08*abb9^4)/abb8^4)+(2*cpp12*abb9^3)/abb8^3+(8*cpp08*aff17)/aff13-4*cpp08*aff17-
        4*cpp12*abb8*abb9-(2*cpp12*abb9)/abb8-4*cpp08*aff13-2*cpp08;
      cqq12 <- -aff14*ajj15*bdd08;
      cqq19 <- bhh14+bhh13;
      cqq20 <- cqq19^3;
      cqq21 <- bhh14+bhh13+aff14*ajj15*bdd08-3*app14*azz19*cll08;
      cqq25 <- bhh14+bhh13+cqq12;
      cqq28 <- 1/aff13;
      q3  <-  (6*cqq20*abb9^4)/abb8^4-(6*cqq19*cqq25*abb9^3)/abb8^3-8*cqq20*cqq28*aff17+cqq21*cqq28*aff17+
        4*cqq20*aff17+cqq21*aff17+12*cqq19*cqq25*abb8*abb9+(6*cqq19*cqq25*abb9)/abb8+4*cqq20*aff13+
        cqq21*aff13+2*cqq20-ann05*ann06+abb4*aaa2*bdd14+cqq12+3*app14*azz19*cll08;
      crr18 <- aii11+aoo09+aff14*ajj15*aii17-3*app14*azz19*cbb09;
      crr19 <- bii11^4;
      crr21 <- bii15^2;
      crr25 <- aii11+aoo09-3*aff14*ajj15*aii17+15*app14*azz19*cbb09-(15*exp1^(aaa1+6*abb5)*aaa2^7)/agg15^7;
      crr28 <- bii11^2;
      r3  <-  (-(6*crr19*abb9^4)/abb8^4)+(12*crr28*bii15*abb9^3)/abb8^3-3*crr21*ass28*aff17+8*crr19*ass28*aff17-
        4*bii11*crr18*ass28*aff17-3*crr21*aff17-4*crr19*aff17-4*bii11*crr18*aff17-24*crr28*bii15*abb8*abb9-
        crr25*abb8*abb9-12*crr28*bii15*ass23*abb9+crr25*ass23*abb9-3*crr21*aff13-4*crr19*aff13-4*bii11*crr18*aff13+
        3*crr21-2*crr19+4*bii11*crr18-(8*abb1*abb3)/aff04+(56*exp1^((-4*t)-4*p)*aaa2^4)/aff04^2-
        (96*exp1^((-6*t)-6*p)*aaa2^6)/aff04^3+(48*exp1^((-8*t)-8*p)*aaa2^8)/aff04^4;
      
      ## .... end auto code
      
      L4 <- cbind(j2,k2,l2,m2,n2,o2,p2,q2,r2,s2,t2,u2,v2,w2,x2,y2,z2,a3,
                  b3,c3,d3,e3,f3,g3,h3,i3,j3,k3,l3,m3,n3,o3,p3,q3,r3)
      
      G4 <- cbind(family$linfo[[1]]$d4link(mu), family$linfo[[2]]$d4link(tau), 
                  family$linfo[[3]]$d4link(eps), family$linfo[[4]]$d4link(phi))
    }
    
    if (deriv) {
      I2 <- family$tri$i2
      I3 <- family$tri$i3
      I4 <- family$tri$i4
      
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- gamlss.etamu(L1,L2,L3,L4,IG1,G2,G3,G4,I2,I3,I4,deriv-1)
      
      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,de$l1,de$l2,I2,l3=de$l3,i3=I3,l4=de$l4,i4=I4,
                       d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D) 
      
    } else ret <- list()
    ret$l <- l; ret
  } ## end ll
  
  initialize <- expression({
    ## idea is to regress g(y) on model matrix for mean, and then 
    ## to regress the corresponding log absolute residuals on 
    ## the model matrix for log(sigma) - may be called in both
    ## gam.fit5 and initial.spg... note that appropriate E scaling
    ## for full calculation may be inappropriate for initialization 
    ## which is basically penalizing something different here.
    ## best we can do here is to use E only as a regularizer.
    n <- rep(1, nobs)
    ## should E be used unscaled or not?..
    use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {
      jj <- attr(x,"lpi")
      start <- rep(0,ncol(x))
      yt1 <- y
      x1 <- x[ , jj[[1]], drop=FALSE]
      e1 <- E[ , jj[[1]], drop=FALSE] ## square root of total penalty
      #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
      # 1) Ridge regression for the location parameter
      if (use.unscaled) {
        qrx <- qr(rbind(x1, e1))
        x1 <- rbind(x1, e1)
        startMu <- qr.coef(qr(x1), c(yt1,rep(0,nrow(E))))
        startMu[ !is.finite(startMu) ] <- 0       
      } else { startMu <- pen.reg(x1, e1, yt1) }
      start[jj[[1]]] <- startMu
      
      # 2) Ridge regression using log absolute residuals
      lres1 <- log( abs(y-family$linfo[[1]]$linkinv(x[,jj[[1]],drop=FALSE]%*%startMu)) )
      x1 <-  x[,jj[[2]],drop=FALSE]; e1 <- E[,jj[[2]],drop=FALSE]
      #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
      if (use.unscaled) {
        x1 <- rbind(x1,e1)
        startTau <- qr.coef(qr(x1),c(lres1,rep(0,nrow(E))))   
        startTau[!is.finite(startTau)] <- 0
      } else { startTau <- pen.reg(x1,e1,lres1) }
      start[jj[[2]]] <- startTau
      
      # 3) Skewness and kurtosis as for Gaussian density: skewness set to zero (identity link)
      # and log-kurtosis set to zero.
      start[jj[[3]]] <- 0
      start[jj[[4]]] <- c(family$linfo[[4]]$linkfun(0), rep(0, length(jj[[4]])-1))}
  }) ## initialize 
  
  #   postproc <- expression({  ####### XXX ??? #######
  #   })
  
  rd <- function(mu, wt, scale) { 
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }
    muE <- mu[ , 1]
    sigE <- exp(mu[ , 2])
    epsE <- mu[ , 3]
    delE <- exp(mu[ , 4])
    n <- length(muE)
    
    .r <- muE + (delE * sigE) * sinh((1/delE) * asinh(qnorm(runif(n))) + (epsE/delE))
    
    return( .r )
  }
  
  qf <- function(p, mu, wt, scale) {
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }
    muE <- mu[ , 1]
    sigE <- exp(mu[ , 2])
    epsE <- mu[ , 3]
    delE <- exp(mu[ , 4])
    
    q <- muE + (delE * sigE) * sinh((1/delE) * asinh(qnorm(p)) + (epsE/delE))
    
    return( q)
  }
  
  pf <- function(q, mu, wt, scale) {
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }
    muE <- mu[ , 1]
    sigE <- exp(mu[ , 2])
    epsE <- mu[ , 3]
    delE <- exp(mu[ , 4])
    
    #q <- muE + (delE * sigE) * sinh((1/delE) * asinh(qnorm(p)) + (epsE/delE))
    
    p <- pnorm( sinh((asinh( (q-muE)/(delE * sigE) )  - epsE/delE) * delE) )
    
    return( p )
  }
  
  
  structure(list(family="shash",ll=ll, link=paste(link), nlp=npar,
                 tri = trind.generator(npar), ## symmetric indices for accessing derivative arrays
                 initialize=initialize,
                 #postproc=postproc,
                 residuals=residuals,
                 rd=rd,
                 qf=qf,
                 pf=pf,
                 #predict=predict,
                 linfo = stats, ## link information list
                 d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
                 ls=1, ## signals that ls not needed here
                 available.derivs = 2 ## can use full Newton here
  ),class = c("general.family","extended.family","family"))
} ## shash

