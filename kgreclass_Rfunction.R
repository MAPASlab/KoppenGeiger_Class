################################################################################
##########################  "kg_reclass" function  #############################
################################################################################

# R functions based on Beck et al. (2018) "Koppengeiger" function on Matlab.
# "Temp" and "Prec" arguments represent temperature and precipitation regimenes, 
# and should be provided as three dimentional arrays with the third dimension 
# representing time (12 months).
# Temp data should be provided in degrees Celsius, and Prep data in mm/months.

# "type" argument has to be "class" (for results indicating the numeric indentifier
# of the specific class 1-30) or "broadclass" (for results indicating the numeric 
# identifier of the broad class (1-5).

kg_reclass <- function(Temp, Prec, type = c("broadclass", "class")) {
  type <- match.arg(type)
  
  if (identical(dim(Temp), dim(Prec))==F){
    stop("Data matrices have not the same dimensions")
  }
  
  if (any(Prec<0, na.rm = T)==T) {
    stop("Precipitation data can not include negative values") #Added na.rm
  }
  
  Temp[is.na(Temp)] <- -99 #Added
  Prec[is.na(Prec)] <- -99 #Added
  
  T_ONDJFM <- apply(Temp[, , c(1,2,3,10,11,12)], c(1,2), mean) 
  T_AMJJAS <- apply(Temp[, , c(4, 5, 6, 7, 8, 9)], c(1,2), mean)
  tmp <- T_AMJJAS>T_ONDJFM
  SUM_SEL <- array(as.logical(0), dim(Temp))
  SUM_SEL[,, c(4, 5, 6, 7, 8, 9)] <- rep(tmp, 6)
  SUM_SEL[,, c(10, 11, 12, 1, 2, 3)] <- rep(!tmp,6)
  rm(tmp)
  
  Pw <- apply(Prec*!SUM_SEL, c(1,2), sum) 
  Ps <- apply(Prec*SUM_SEL, c(1,2), sum) 
  
  Pdry <- apply (Prec, c(1,2), min) 
  
  tmp <- SUM_SEL
  tmp[tmp==0] = NA
  Psdry <- apply(Prec*tmp, c(1,2), min, na.rm = TRUE) 
  Pswet <- apply(Prec*tmp, c(1,2), max, na.rm = TRUE) 
  
  tmp <- !SUM_SEL
  tmp[tmp==0] = NA
  Pwdry <- apply(Prec*tmp, c(1,2), min, na.rm = TRUE) 
  Pwwet <- apply(Prec*tmp, c(1,2), max, na.rm = TRUE) 
  
  MAT <- apply(Temp, c(1,2), mean)
  MAP <- apply(Prec, c(1,2), sum) 
  Tmon10 <- apply(Temp > 10, c(1,2), sum) 
  Thot <- apply(Temp, c(1,2), max)
  Tcold <- apply(Temp, c(1,2), min)
  
  Pthresh <- 2*MAT+14 #where temp = -99, this is -184
  Pthresh[Pw>Ps*2.333] <- 2*MAT[Pw>Ps*2.333]       
  Pthresh[Ps>Pw*2.333] <- 2*MAT[Ps>Pw*2.333]+28 
  
  B <- MAP < 10*Pthresh 
  BW <- B & MAP < 5*Pthresh 
  BWh <- BW & MAT >= 18
  BWk <- BW & MAT < 18
  BS <- B & MAP >= 5*Pthresh
  BSh <- BS & MAT >= 18
  BSk <- BS & MAT < 18
  
  A <- Tcold >= 18 & !B 
  Af <- A & Pdry >= 60
  Am <- A & !Af & Pdry >= 100-MAP/25
  Aw <- A & !Af & Pdry < 100-MAP/25
  
  C <- Thot > 10 & Tcold > 0 & Tcold < 18 & !B 
  Cs <- C & Psdry<40 & Psdry<Pwwet/3
  Cw <- C & Pwdry<Pswet/10
  overlap <- Cs & Cw
  Cs[overlap & Ps>Pw] <- 0
  Cw[overlap & Ps<=Pw] <- 0
  Csa <- Cs & Thot >= 22
  Csb <- Cs & !Csa & Tmon10 >= 4
  Csc <- Cs & !Csa & !Csb & Tmon10>=1 & Tmon10<4
  Cwa <- Cw & Thot >= 22
  Cwb <- Cw & !Cwa & Tmon10 >= 4
  Cwc <- Cw & !Cwa & !Cwb & Tmon10>=1 & Tmon10<4
  Cf <- C & !Cs & !Cw
  Cfa <- Cf & Thot >= 22
  Cfb <- Cf & !Cfa & Tmon10 >= 4
  Cfc <- Cf & !Cfa & !Cfb & Tmon10>=1 & Tmon10<4
  
  D <- Thot>10 & Tcold<=0 & !B   
  Ds <- D & Psdry<40 & Psdry<Pwwet/3
  Dw <- D & Pwdry<Pswet/10
  overlap <- Ds & Dw
  Ds[overlap & Ps>Pw] <- 0
  Dw[overlap & Ps<=Pw] <- 0
  Dsa <- Ds & Thot>=22
  Dsb <- Ds & !Dsa & Tmon10 >= 4
  Dsd <- Ds & !Dsa & !Dsb & Tcold < (-38) 
  Dsc <- Ds & !Dsa & !Dsb & !Dsd
  
  Dwa <- Dw & Thot>=22
  Dwb <- Dw & !Dwa & Tmon10 >= 4
  Dwd <- Dw & !Dwa & !Dwb & Tcold < (-38)
  Dwc <- Dw & !Dwa & !Dwb & !Dwd
  Df <- D & !Ds & !Dw
  Dfa <- Df & Thot>=22
  Dfb <- Df & !Dfa & Tmon10 >= 4
  Dfd <- Df & !Dfa & !Dfb & Tcold < (-38)
  Dfc <- Df & !Dfa & !Dfb & !Dfd
  
  E <- Thot <= 10 & Thot > (-90) & !B    #Added & Thot > (-90)
  ET <- E & Thot>0
  EF <- E & Thot<=0 & Thot> (-90) # Added & Thot> (-90)
  
  
  if (type == "class") {
    Class <- list(Af, Am, Aw, BWh, BWk, BSh, BSk, Csa, Csb, Csc, Cwa, Cwb,
                  Cwc, Cfa, Cfb, Cfc, Dsa, Dsb, Dsc, Dsd, Dwa, Dwb, Dwc, Dwd, Dfa,
                  Dfb, Dfc, Dfd, ET, EF)
    Class_cont <- array(0, dim(Temp[,,1]))
    for (i in 1:30){
      Class_cont[Class[[i]]==1] <- i
    }
    Class_cont[Class_cont == 0] <- NA #Added
    return(Class_cont)
  } 
  if (type == "broadclass") {
    Broadclass <- list(A, B, C, D, E)
    BroadClass_cont <- array(0, dim(Temp[,,1]))
    for (i in 1:5){
      BroadClass_cont[Broadclass[[i]]==1] <- i
    }
    BroadClass_cont[BroadClass_cont == 0] <- NA #Added 
    return(BroadClass_cont)
  }   
}
