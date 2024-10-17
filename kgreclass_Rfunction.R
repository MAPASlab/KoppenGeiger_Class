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

  if (!identical(dim(Temp), dim(Prec))) {
    stop("Data matrices do not have the same dimensions")
  }

  if (any(Prec < 0, na.rm = TRUE)) {
    stop("Precipitation data cannot include negative values") #Added na.rm
  }
  
  Temp[is.na(Temp)] <- -99 #Added
  Prec[is.na(Prec)] <- -99 #Added

  temp_by_month <- lapply(seq_len(dim(Temp)[3]), function(i) Temp[, , i])
  prec_by_month <- lapply(seq_len(dim(Prec)[3]), function(i) Prec[, , i])

  T_ONDJFM <- rowMeans(Temp[, , c(1, 2, 3, 10, 11, 12)], dims = 2)
  T_AMJJAS <- rowMeans(Temp[, , c(4, 5, 6, 7, 8, 9)], dims = 2)
  tmp <- T_AMJJAS>T_ONDJFM
  SUM_SEL <- array(FALSE, dim(Temp))
  SUM_SEL[,, c(4, 5, 6, 7, 8, 9)] <- rep(tmp, 6)
  SUM_SEL[,, c(10, 11, 12, 1, 2, 3)] <- rep(!tmp,6)

  Pw <- rowSums(Prec * !SUM_SEL, dims = 2)
  Ps <- rowSums(Prec * SUM_SEL, dims = 2)

  Pdry <- do.call(pmin, prec_by_month)

  tmp <- SUM_SEL
  tmp[tmp == 0] <- NA
  prec_masked <- Prec * tmp
  prec_masked_by_month <- lapply(
    seq_len(dim(prec_masked)[3]),
    function(i) prec_masked[, , i]
  )
  Psdry <- do.call(pmin, c(prec_masked_by_month, na.rm = TRUE))
  Pswet <- do.call(pmax, c(prec_masked_by_month, na.rm = TRUE))

  tmp <- !SUM_SEL
  tmp[tmp == 0] <- NA
  prec_masked <- Prec * tmp
  prec_masked_by_month <- lapply(
    seq_len(dim(prec_masked)[3]),
    function(i) prec_masked[, , i]
  )
  Pwdry <- do.call(pmin, c(prec_masked_by_month, na.rm = TRUE))
  Pwwet <- do.call(pmax, c(prec_masked_by_month, na.rm = TRUE))

  MAT <- rowMeans(Temp, dims = 2)
  MAP <- rowSums(Prec, dims = 2)
  Tmon10 <- rowSums(Temp > 10, dims = 2)
  Thot <- do.call(pmax, temp_by_month)
  Tcold <- do.call(pmin, temp_by_month)

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
  } else if (type == "broadclass") {
    Class <- list(A, B, C, D, E)
  }

  Class_cont <- array(NA_integer_, dim(Temp[, , 1]))
  for (i in seq_along(Class)) {
    Class_cont[Class[[i]]] <- i
  }

  Class_cont
}
