## code for fast REML computation. key feature is that first and 
## second derivatives come at no increase in leading order 
## computational cost, relative to evaluation! 
## (c) Simon N. Wood, 2010-2019

singleStrans <- function(S,rank=NULL,ldet=FALSE) {
## transform single penalty matrix to partial identity using Cholesky
## t(D)%*%S%*%D transforms to a partial identity (rank rank)
## if ldet==TRUE then the generalized determinant of S is also returned
  Ri <- R <- suppressWarnings(chol(S,pivot=TRUE))
  k <- ncol(S)
  if (is.null(rank)) {
    rank <- Rrank(R)
    normS <- norm(S)
    while (rank>1&&max(abs(crossprod(R[1:(rank-1),rank:k,drop=FALSE])-S[rank:k,rank:k,drop=FALSE]))
           <.Machine$double.eps^.75*normS) rank <- rank - 1
  }
  piv <- attr(R,"pivot")
  if (rank<k) {
    R[(rank+1):k,] <- 0
    diag(R)[(rank+1):k] <- 1
  }
  ## compute the log determinant if required... 
  ldet <- if (ldet) 2*sum(log(diag(chol(tcrossprod(R[1:rank,]))))) else 0
  Ri[piv,] <- solve(R)
  R[,piv] <- R
  list(D=Ri,Di=R,rank=rank,ldet=ldet)
} ## singleStrans


iniStrans <- function(S,rank = NULL,trans.ldet=FALSE) {
## let a block penalty matrix be \sum_i \lambda_i S_i, and of rank r.
## This routine finds a reparameterization such that the Si are confined to
## the initial r by r block.
  St <- S[[1]]/norm(S[[1]]); m <- length(S)
  if (m>1) for (i in 2:m) St <- St + S[[i]]/norm(S[[i]])
  R <- suppressWarnings(chol(St,pivot=TRUE))
  p <- nrow(St)
  if (is.null(rank)) {
    rank <- Rrank(R)
    normS <- norm(St) ## direct check that rank not over-estimated
    while (rank>1&&max(abs(crossprod(R[1:(rank-1),rank:p,drop=FALSE])-St[rank:p,rank:p,drop=FALSE]))
           <.Machine$double.eps^.75*normS) rank <- rank - 1
  }
  piv <- attr(R,"pivot")
  ipiv <- piv; ipiv[piv] <- 1:p
  if (rank==p) { ## nothing to do - penalty is full rank
    return(list(S=S,T=diag(p),Ti=diag(p),rank=p,trans.ldet=0))
  }
  ind <- (rank+1):p
  R[ind,ind] <- diag(length(ind))
  # S1 <- S # save original for checking
  for (i in 1:m) {
    S[[i]] <- forwardsolve(t(R),t(forwardsolve(t(R),S[[i]][piv,piv])))[1:rank,1:rank]
    S[[i]] <- (S[[i]] + t(S[[i]]))*.5
  }
  T <- backsolve(R,diag(p))[ipiv,]
  #range((t(T)%*%S1[[1]]%*%T)[1:rank,1:rank]-S[[1]])
  #range(T%*%R[,ipiv]-diag(p))
  ## So T will map original S to transformed S, and R[,ipiv] is inverse transform.
  ## now get the correction to be added to the log|sum_i S_i| to get log generalized
  ## determinant in orginal parameterization...
  tldet <- if (trans.ldet) 2*sum(log(diag(chol(tcrossprod(R[1:rank,]))))) else 0
  list(S=S,T=T,Ti=R[,ipiv],rank=rank,trans.ldet=tldet)
} ## iniStrans


Sl.setup <- function(G,cholesky=FALSE,no.repara=FALSE,sparse=FALSE,keepS=FALSE) {
## Sets up a list representing a block diagonal penalty matrix.
## from the object produced by `gam.setup'.
## Uses only pivoted Cholesky if cholesky==TRUE.
## flags all blocks as not to be reparameterized if no.repara=TRUE,
## in which case initial reparameterization is only used to
## compute log generalized determinant, but not otherwise.
## Return object is a list, Sl, with an element for each block.
## For block, b, Sl[[b]] is a list with the following elements
## * repara - should re-parameterization be applied to model matrix etc 
##            usually false if  non-linear in coefs
## * start, stop: start:stop indexes the parameters of this block
## * S a list of penalty matrices for the block (dim = stop-start+1)
##   - If length(S)==1 then this will be an identity penalty.
##   - Otherwise it is a multiple penalty, and an rS list of square
##     root penalty matrices will be added. S (if repara) and rS (always) 
##     will be projected into range space of total penalty matrix.
##     If cholesky==TRUE then rS contains the projected S if !repara and otherwise
##     rS not returned (since S contains same thing).
## * rS sqrt penalty matrices if it's a multiple penalty and cholesky==FALSE.
##      projected penalty if cholesky==TRUE and !repara. NULL otherwise. 
## * D a reparameterization matrix for the block
##   - Applies to cols/params from start:stop.
##   - If numeric then X[,start:stop]%*%diag(D) is repara X[,start:stop],
##     b.orig = D*b.repara
##   - If matrix then X[,start:stop]%*%D is repara X[,start:stop],
##     b.orig = D%*%b.repara
## * Di is inverse of D, but is only supplied if D is not orthogonal, or
##   diagonal.
## The penalties in Sl are in the same order as those in G
## Also returns attribute "E" a square root of the well scaled total
## penalty, suitable for rank deficiency testing, and attribute "lambda"
## the corresponding smoothing parameters.
## keepS==TRUE causes original penalties to be stored in S0 (after any splitting
## of multiple penalties into singletons, but before reparameterization).
  ##if (!is.null(G$H)) stop("paraPen min sp not supported")
  Sl <- list()
  b <- 1 ## block counter
  if (G$n.paraPen) { ## Have to proccess paraPen stuff first
    off <- unique(G$off[1:G$n.paraPen]) ## unique offset lists relating to paraPen
    for (i in 1:length(off)) { ## loop over blocks
      ind <- (1:G$n.paraPen)[G$off[1:G$n.paraPen]%in%off[i]] ## terms in same block
      if (length(ind)>1) { ## additive block
        nr <- 0;for (k in 1:length(ind)) nr <- max(nr,nrow(G$S[[ind[k]]])) ## get block size
        ## now fill Sl[[b]]$S, padding out any penalties that are not "full size"
        Sl[[b]] <- list()
        Sl[[b]]$S <- list()
        St <- matrix(0,nr,nr) ## accumulate a total matrix for rank determination 
        for (k in 1:length(ind)) { ## work through all penalties for this block
          nk <- nrow(G$S[[ind[k]]])
          if (nr>nk) { ## have to pad out this one
            Sl[[b]]$S[[k]] <- matrix(0,nr,nr)
            Sl[[b]]$S[[k]][1:nk,1:nk] <- G$S[[ind[k]]]
          } else Sl[[b]]$S[[k]] <- G$S[[ind[[k]]]]
          St <- St + Sl[[b]]$S[[k]] 
        }
        Sl[[b]]$start <- off[ind[1]]
        Sl[[b]]$stop <- Sl[[b]]$start + nr - 1
	Sl[[b]]$lambda <- rep(1,length(ind)) ## dummy at this stage
	Sl[[b]]$repara <- FALSE
      } else { ## singleton
        Sl[[b]] <- list(start=off[ind], stop=off[ind]+nrow(G$S[[ind]])-1,
                        rank=G$rank[ind],S=list(G$S[[ind]]))
        Sl[[b]]$S <- list(G$S[[ind]])
	Sl[[b]]$lambda <- 1 ## dummy at this stage
	Sl[[b]]$repara <- !no.repara ## allow repara unless completely turned off
      } ## finished singleton
      Sl[[b]]$linear <- TRUE
      b <- b + 1 
    } ## finished this block
  } ## finished paraPen

  ## now work through the smooths....
  if (length(G$smooth)) for (i in 1:length(G$smooth)) {

    if (!is.null(G$smooth[[i]]$fixed)&&G$smooth[[i]]$fixed) m <- 0 else
    m <- length(G$smooth[[i]]$S)

    if (m>0) {
      Sl[[b]] <- list()
      Sl[[b]]$nl.reg <- if (is.null(G$smooth[[i]]$nl.reg)||G$smooth[[i]]$nl.reg<=0) NULL else G$smooth[[i]]$nl.reg ## any fixed regularizer for non-linear coefs
      Sl[[b]]$start <- G$smooth[[i]]$first.para
      Sl[[b]]$stop <- G$smooth[[i]]$last.para
      Sl[[b]]$srank <- G$smooth[[i]]$rank ## original rank of penalty matrix
      ## if the smooth has a g.index field it indicates non-linear params,
      ## in which case re-parameterization will usually break the model.
      ## Or global supression of reparameterization may be requested...
      # Sl[[b]]$repara <- if (is.null(G$smooth[[i]]$g.index) && !no.repara) TRUE else FALSE
      Sl[[b]]$repara <- G$smooth[[i]]$repara && !no.repara  
      if (!is.null(G$smooth[[i]]$updateS)) { ## then this block is nonlinear in smoothing parameters
        Sl[[b]]$repara <-FALSE
	Sl[[b]]$linear <- FALSE
	labs <- c("inisp","updateS","AS","AdS","ldS","St","n.sp","nlinfo")
	Sl[[b]][labs] <- G$smooth[[i]][labs] ## copy the non-linear interface functions
	Sl[[b]]$lambda <- rep(0,Sl[[b]]$n.sp) ## dummy
      } else Sl[[b]]$linear <- TRUE
    }

    if (m==0) {} else ## fixed block
    if (m==1) { ## singleton
     
        Sl[[b]]$rank <- G$smooth[[i]]$rank  
        Sl[[b]]$S <- G$smooth[[i]]$S
	Sl[[b]]$lambda <- 1
        b <- b + 1
     
    } else { ## additive block...
      ## first test whether block can *easily* be split up into singletons
      ## easily here means no overlap in penalties 
      Sl[[b]]$S <- G$smooth[[i]]$S
      Sl[[b]]$lambda <- rep(1,m)
      if (Sl[[b]]$linear) {
        nb <- nrow(Sl[[b]]$S[[1]])      
        sbdiag <- sbStart <- sbStop <- rep(NA,m)
        ut <- upper.tri(Sl[[b]]$S[[1]],diag=FALSE) 
        ## overlap testing requires the block ranges 
        for (j in 1:m) { ## get block range for each S[[j]]
          sbdiag[j] <- sum(abs(Sl[[b]]$S[[j]][ut]))==0 ## is penalty diagonal??
          ir <- range((1:nb)[rowSums(abs(Sl[[b]]$S[[j]]))>0])
          sbStart[j] <- ir[1];sbStop[j] <- ir[2] ## individual ranges
        } 
        split.ok <- TRUE
        for (j in 1:m) { ## test for overlap
          itot <- rep(FALSE,nb)
          if (all(sbdiag)) { ## it's all diagonal - can allow interleaving
            for (k in 1:m) if (j!=k) itot[diag(Sl[[b]]$S[[k]])!=0] <- TRUE
            if (sum(itot[diag(Sl[[b]]$S[[j]])!=0])>0) { ## no good, some overlap detected
              split.ok <- FALSE; break
            }
          } else { ## not diagonal - really need  overlapping blocks
            for (k in 1:m) if (j!=k) itot[sbStart[k]:sbStop[k]] <- TRUE
            if (sum(itot[sbStart[j]:sbStop[j]])>0) { ## no good, some overlap detected
              split.ok <- FALSE; break
            }
          }
        }
      } else split.ok <- FALSE	
      if (split.ok) { ## can split this block into m separate singleton blocks
        for (j in 1:m) {
          Sl[[b]] <- list()
          ind <- sbStart[j]:sbStop[j]
          Sl[[b]]$S <- list(G$smooth[[i]]$S[[j]][ind,ind,drop=FALSE])
          Sl[[b]]$start <- G$smooth[[i]]$first.para + sbStart[j]-1
          Sl[[b]]$stop <- G$smooth[[i]]$first.para + sbStop[j]-1
          Sl[[b]]$rank <- G$smooth[[i]]$rank[j]
	  Sl[[b]]$lambda <- 1 ## dummy here
          Sl[[b]]$repara <- !no.repara ## signals ok to linearly reparameterize, now check this is really ok...
          if (!is.null(G$smooth[[i]]$g.index)) { ## then some parameters are non-linear - can't re-param
            if (any(G$smooth[[i]]$g.index[ind])) Sl[[b]]$repara <- FALSE
          }
	  Sl[[b]]$linear <- TRUE ## linear in smoothing params - assumed always true for multi-blocks
          b <- b + 1
        }
      } else { ## not possible to split
        Sl[[b]]$S <- G$smooth[[i]]$S  
        b <- b + 1 ## next block!!
      } ## additive block finished
    } ## additive block finished
  }

  ## At this stage Sl contains the penalties, identified as singletons or 
  ## multiple S blocks. Now the blocks need re-parameterization applied.
  ## Singletons need to be transformed to identity penalties, while 
  ## multiples need to be projected into total penalty range space. 

  if (length(Sl)==0) return(Sl) ## nothing to do

  np <- ncol(G$X)

  ## In dense case E is a dense matrix and its block get filled in. In sparse
  ## case E is a list of lists of rows (i), cols (j) and elements (x) defining each
  ## block. So E$i[[b]],E$j[[b]],E$x[[b]] defines block b...
  
  S <- E <- if (sparse) list(i=list(),j=list(),x=list()) else matrix(0,np,np) ## well scaled penalty & square root penalty
  lambda <- rep(0,0)

  ## NOTE: computing transforms for repara=FALSE blocks and then not using them
  ##       looks wasteful - remove this??

  for (b in 1:length(Sl)) { ## once more into the blocks, dear friends...
    if (keepS) Sl[[b]]$S0 <- Sl[[b]]$S ## keep untransformed version
    if (!Sl[[b]]$linear) { ## nonlinear term
      ## never re-parameterized, but need contribution to penalty square root, E
      Sl[[b]] <- Sl[[b]]$updateS(Sl[[b]]$lambda,Sl[[b]])
      lambda <- c(lambda,Sl[[b]]$lambda)
      tmp <- Sl[[b]]$St(Sl[[b]],2)
      if (sparse) {
        ## E0 <- as(Sl[[b]]$St(Sl[[b]],1)$E,"dgTMatrix") deprecated
	E0 <- as(as(as(tmp$E, "dMatrix"), "generalMatrix"), "TsparseMatrix")
	E$i[[b]] <- E0@i + Sl[[b]]$start
	E$j[[b]] <- E0@j + Sl[[b]]$start;E$x[[b]] <- E0@x
	E0 <- as(as(as(tmp$St, "dMatrix"), "generalMatrix"), "TsparseMatrix")
	S$i[[b]] <- E0@i + Sl[[b]]$start
	S$j[[b]] <- E0@j + Sl[[b]]$start;S$x[[b]] <- E0@x
      } else { ## dense
        ind <- Sl[[b]]$start:Sl[[b]]$stop
        E[ind,ind] <- tmp$E ## add to the square root
	S[ind,ind] <- tmp$St
      }	
    } else if (length(Sl[[b]]$S)==1) { ## then we have a singleton
      if (sum(abs(Sl[[b]]$S[[1]][upper.tri(Sl[[b]]$S[[1]],diag=FALSE)]))==0) { ## S diagonal
        ## Reparameterize so that S has 1's or zero's on diagonal
        ## In new parameterization smooth specific model matrix is X%*%diag(D)
        ## ind indexes penalized parameters from this smooth's set.

        D <- diag(Sl[[b]]$S[[1]])
        ind <- D > 0 ## index penalized elements
	Sl[[b]]$rank <- sum(ind)
	Sl[[b]]$ldet <- if (cholesky) sum(log(D[ind])) else 0
        D[ind] <- 1/sqrt(D[ind]);D[!ind] <- 1 ## X' = X%*%diag(D) 
        Sl[[b]]$D <- D; Sl[[b]]$ind <- ind
      } else { ## S is not diagonal
        if (cholesky&&is.null(Sl[[b]]$nl.reg)) { ## use Cholesky based reparameterization
          tr <- singleStrans(Sl[[b]]$S[[1]],Sl[[b]]$rank,ldet=!Sl[[b]]$repara)
	  ind <- rep(FALSE,ncol(tr$D))
	  ind[1:tr$rank] <- TRUE
	  Sl[[b]]$D <- tr$D
	  Sl[[b]]$Di <- tr$Di
	  Sl[[b]]$rank <- tr$rank
	  Sl[[b]]$ldet=tr$ldet
        } else { ## use eigen based re-parameterization
          es <- eigen(Sl[[b]]$S[[1]],symmetric=TRUE)
          U <- es$vectors;D <- es$values
          if (is.null(Sl[[b]]$rank)) { ## need to estimate rank
            Sl[[b]]$rank <- sum(D>.Machine$double.eps^.8*max(D))
          }
	  ## non-cholesky stabilized method code ignores the log det
	  ## of transform as it cancels between the two log det terms...
	  Sl[[b]]$ldet <- 0 ## sum(log(D[1:Sl[[b]]$rank]))
	  ind <- rep(FALSE,length(D))
          ind[1:Sl[[b]]$rank] <- TRUE ## index penalized elements
          if (is.null(Sl[[b]]$nl.reg)) {
            D[ind] <- 1/sqrt(D[ind]);D[!ind] <- 1
            Sl[[b]]$D <- t(D*t(U)) ## D <- U%*%diag(D)
            Sl[[b]]$Di <- t(U)/D
	  } else {
            Sl[[b]]$repara <- FALSE ## should be FALSE anyway to get here
	    Sl[[b]]$ev <- D ## store the penalty eigenvalues to allow penalized ldet computation
	    Sl[[b]]$U <- U ## store the eigen-vectors to allow computation of regularized penalty square root
          }
	}  
        ## so if X is smooth model matrix X%*%D is re-parameterized form
	## and t(D)%*%Sl[[b]]$S[[1]]%*%D is the reparameterized penalty
	## -- a partial identity matrix.
        ## Di is the inverse of D and crossprod(Di[1:rank,]) is the original
	## penalty matrix
        Sl[[b]]$ind <- ind
      }
      ## add penalty into S and square root into E  
      if (Sl[[b]]$repara) { ## then it is just the identity
        ind <- (Sl[[b]]$start:Sl[[b]]$stop)[Sl[[b]]$ind]
	if (sparse) {
	  E$j[[b]] <- E$i[[b]] <- ind; E$x[[b]] <- ind*0 + 1
	  S$j[[b]] <- S$i[[b]] <- ind; S$x[[b]] <- ind*0 + 1
	} else { ## dense
          diag(E)[ind] <- 1; diag(S)[ind] <- 1
	}  
        lambda <- c(lambda,1) ## record corresponding lambda
      } else { ## need scaled root penalty in *original* parameterization
        #D <- Sl[[b]]$Di[1:Sl[[b]]$rank,]
	D <- if (is.null(Sl[[b]]$nl.reg)) Sl[[b]]$Di[1:Sl[[b]]$rank,] else
	        sqrt(Sl[[b]]$ev+Sl[[b]]$nl.reg)*t(Sl[[b]]$U)
        D.norm <- norm(D); D <- D/D.norm
        if (sparse) {
	  ## D <- as(D,"dgTMatrix") deprecated
	  D <- as(as(as(D, "dMatrix"), "generalMatrix"), "TsparseMatrix")
	  E$i[[b]] <- D@i + Sl[[b]]$start
	  E$j[[b]] <- D@j + Sl[[b]]$start; E$x[[b]] <- D@x
        } else {
	  indc <- Sl[[b]]$start:(Sl[[b]]$start+ncol(D)-1)
          indr <- Sl[[b]]$start:(Sl[[b]]$start+nrow(D)-1); E[indr,indc] <- D
	}
	## and then the total penalty matrix ...
        D <- if (is.null(Sl[[b]]$nl.reg)) crossprod(Sl[[b]]$Di[1:Sl[[b]]$rank,]) else
	        Sl[[b]]$U %*% ((Sl[[b]]$ev+Sl[[b]]$nl.reg)*t(Sl[[b]]$U))
        D <- D/D.norm^2
	if (sparse) {
	  ## D <- as(D,"dgTMatrix") deprecated
	  D <- as(as(as(D, "dMatrix"), "generalMatrix"), "TsparseMatrix")
	  S$i[[b]] <- D@i + Sl[[b]]$start
	  S$j[[b]] <- D@j + Sl[[b]]$start; S$x[[b]] <- D@x
        } else {
	  indc <- Sl[[b]]$start:(Sl[[b]]$start+ncol(D)-1)
          indr <- Sl[[b]]$start:(Sl[[b]]$start+nrow(D)-1); S[indr,indc] <- D
	}

        lambda <- c(lambda,1/D.norm^2)
      }
    } else { ## multiple S block
      ## must be in range space of total penalty...
      Sl[[b]]$ind <- rep(FALSE,ncol(Sl[[b]]$S[[1]]))
      if (cholesky) {
        tr <- iniStrans(Sl[[b]]$S,Sl[[b]]$rank,trans.ldet=!Sl[[b]]$repara)
	if (!Sl[[b]]$repara) Sl[[b]]$rS <- list() 
	for (i in 1:length(tr$S)) {
          if (Sl[[b]]$repara) Sl[[b]]$S[[i]] <- tr$S[[i]] else
	  Sl[[b]]$rS[[i]] <- tr$S[[i]] ## only need to store here if !repara
        }
	ind <- 1:tr$rank
	Sl[[b]]$rank <- tr$rank
	Sl[[b]]$D <- tr$T
	Sl[[b]]$Di <- tr$Ti
	Sl[[b]]$ldet = tr$trans.ldet
      } else {
        Sl[[b]]$ldet = 0 ## this is pure orthogonal transform
        Sl[[b]]$rS <- list() ## needed for adaptive re-parameterization
        St <- Sl[[b]]$S[[1]]
        for (j in 2:length(Sl[[b]]$S)) St <- St + Sl[[b]]$S[[j]] ## scaled total penalty
        es <- eigen(St,symmetric=TRUE);U <- es$vectors; D <- es$values
        Sl[[b]]$D <- U
        if (is.null(Sl[[b]]$rank)) { ## need to estimate rank
            Sl[[b]]$rank <- sum(D>.Machine$double.eps^.8*max(D))
        }
        ind <- 1:Sl[[b]]$rank
        for (j in 1:length(Sl[[b]]$S)) { ## project penalties into range space of total penalty
          bob <- t(U[,ind])%*%Sl[[b]]$S[[j]]%*%U[,ind]
          bob <- (t(bob) + bob)/2 ## avoid over-zealous chol sym check
          if (Sl[[b]]$repara) { ## otherwise want St and E in original parameterization
            Sl[[b]]$S[[j]] <- bob
          } 
          Sl[[b]]$rS[[j]] <- mroot(bob,Sl[[b]]$rank)
        }
      }	
      #Sl[[b]]$ind <- rep(FALSE,ncol(Sl[[b]]$S[[1]]))
      Sl[[b]]$ind[ind] <- TRUE ## index penalized within sub-range
     
      ## now compute well scaled sqrt
      S.norm <- norm(Sl[[b]]$S[[1]])
      St <- Sl[[b]]$S[[1]]/S.norm
      lambda <- c(lambda,1/S.norm)
      for (j in 2:length(Sl[[b]]$S)) { 
        S.norm <- norm(Sl[[b]]$S[[j]])
        St <- St + Sl[[b]]$S[[j]]/S.norm
        lambda <- c(lambda,1/S.norm)
      } 
      St <- (t(St) + St)/2  ## avoid over-zealous chol sym check
      Sr <- t(mroot(St,Sl[[b]]$rank))
      if (sparse) {
        ## St <- as(St,"dgTMatrix") - deprecated
	Sr <- as(as(as(Sr, "dMatrix"), "generalMatrix"), "TsparseMatrix")
        E$i[[b]] <- Sr@i + Sl[[b]]$start
	E$j[[b]] <- Sr@j + Sl[[b]]$start; E$x[[b]] <- Sr@x
	St <- as(as(as(St, "dMatrix"), "generalMatrix"), "TsparseMatrix")
        S$i[[b]] <- St@i + Sl[[b]]$start
	S$j[[b]] <- St@j + Sl[[b]]$start; S$x[[b]] <- St@x
      } else { ## dense
        indc <- Sl[[b]]$start:(Sl[[b]]$start+ncol(St)-1)
        indr <- Sl[[b]]$start:(Sl[[b]]$start+nrow(St)-1)
        E[indr,indc] <- Sr; S[indr,indc] <- St
      }	
    }
  } ## re-para finished
  if (sparse) {
    E <- sparseMatrix(i=unlist(E$i),j=unlist(E$j),x=unlist(E$x),dims=c(np,np))
    S <- sparseMatrix(i=unlist(S$i),j=unlist(S$j),x=unlist(S$x),dims=c(np,np))
  }
  attr(Sl,"E") <- E ## E'E = scaled total penalty
  attr(Sl,"S") <- S ## scaled total penalty
  attr(Sl,"lambda") <- lambda ## smoothing parameters corresponding to E
  attr(Sl,"cholesky") <- cholesky ## store whether this is Cholesky based or not
  Sl ## the penalty list
} ## end of Sl.setup

Sl.Sb <- function(Sl,rho,beta) {
## computes S %*% beta where S is total penalty matrix defined by Sl and rho,
## the log smoothing parameters. Assumes initial re-parameterization has taken
## place, so single penalties are multiples of identity and uses S for
## multi-S blocks. Logic is identical to Sl.addS.
  k <- 1
  a <- beta * 0
  if (length(Sl)>0) for (b in 1:length(Sl)) {
    ind <- (Sl[[b]]$start:Sl[[b]]$stop)[Sl[[b]]$ind] 
    if (length(Sl[[b]]$S)==1) { ## singleton - multiple of identity
      a[ind] <- a[ind] + beta[ind] * exp(rho[k])
      k <- k + 1
    } else { ## multi-S block
      for (j in 1:length(Sl[[b]]$S)) {
        a[ind] <- a[ind] + exp(rho[k]) * (Sl[[b]]$S[[j]] %*% beta[ind])
        k <- k + 1
      }
    }
  }
  a
} ## Sl.Sb

Sl.rSb <- function(Sl,rho,beta) {
## Computes vector 'a' containing all terms rS %*% beta stacked end to end.
## sum of squares of 'a' this is bSb, but 'a' is linear in beta
## Assumes initial re-parameterization has taken
## place, so single penalties are multiples of identity and uses S for
## multi-S blocks. Logic is identical to Sl.addS.
  k <- 1 ## sp counter
  kk <- 0 ## total length of returned vector
  if (length(Sl)>0) for (b in 1:length(Sl)) {
    ind <- (Sl[[b]]$start:Sl[[b]]$stop)[Sl[[b]]$ind]
    kk <- kk + length(Sl[[b]]$S)*length(ind)
  }
  a <- rep(0,kk)
  kk <- 0
  if (length(Sl)>0) for (b in 1:length(Sl)) {
    ind <- (Sl[[b]]$start:Sl[[b]]$stop)[Sl[[b]]$ind] 
    if (length(Sl[[b]]$S)==1) { ## singleton - multiple of identity
      a[kk + 1:length(ind)] <- beta[ind] * exp(rho[k]/2)
      k <- k + 1
      kk <- kk + length(ind)
    } else { ## multi-S block
      for (j in 1:length(Sl[[b]]$S)) {
        a[kk + 1:length(ind)] <- exp(rho[k]/2) * (beta[ind] %*% Sl[[b]]$rS[[j]])
        k <- k + 1
	kk <- kk + length(ind)
      }
    }
  }
  a
} ## Sl.rSb


Sl.inirep <- function(Sl,X,l=0,r=0,nt=1) {
## Re-parameterize X using initial Sl reparameterization info.
## l,r = -2,-1,0,1,2. O is do not apply, negative to apply inverse transform
##       Di (t(D) if no Di indicating D orthogonal),
##       positive for transform D, 1 for transform, 2 for its transpose.
## Aim is for simpler and cleaner than Sl.initial.repara
  if (length(Sl)==0 && !l && !r) return(X) ## nothing to do
  if (is.matrix(X)) {
    for (b in 1:length(Sl)) if (Sl[[b]]$repara) {
      nuDi <- is.null(Sl[[b]]$Di)
      ind <- Sl[[b]]$start:Sl[[b]]$stop
      if (l) X[ind,] <- if (l == 1) Sl[[b]]$D%*%X[ind,,drop=FALSE] else if (l == 2) t(Sl[[b]]$D)%*%X[ind,,drop=FALSE] else
             if (l == -1) { if (nuDi) t(Sl[[b]]$D)%*%X[ind,,drop=FALSE] else Sl[[b]]$Di%*%X[ind,,drop=FALSE] } else
			  { if (nuDi) Sl[[b]]$D%*%X[ind,,drop=FALSE] else t(Sl[[b]]$Di)%*%X[ind,,drop=FALSE]} 
      if (r) X[,ind] <- if (l == 1) X[,ind,drop=FALSE]%*%Sl[[b]]$D else if (l == 2) X[,ind,drop=FALSE]%*%t(Sl[[b]]$D) else
             if (l == -1) { if (nuDi) X[,ind,drop=FALSE]%*%t(Sl[[b]]$D) else X[,ind,drop=FALSE]%*%Sl[[b]]$Di} else
			{ if (nuDi) X[,ind,drop=FALSE]%*%Sl[[b]]$D else X[,ind,drop=FALSE]%*%t(Sl[[b]]$Di)}			
    }
  } else { ## it's a vector
    for (b in 1:length(Sl)) if (Sl[[b]]$repara) {
      ind <- Sl[[b]]$start:Sl[[b]]$stop; nuDi <- is.null(Sl[[b]]$Di)
      if (l) X[ind] <- if (l == 1) Sl[[b]]$D%*%X[ind] else if (l == 2) t(Sl[[b]]$D)%*%X[ind] else
                        if (l == -1) { if (nuDi) t(Sl[[b]]$D)%*%X[ind] else Sl[[b]]$Di%*%X[ind] } else
			{ if (nuDi) Sl[[b]]$D%*%X[ind] else t(Sl[[b]]$Di)%*%X[ind] }
      if (r) X[ind] <- if (l == 1) X[ind]%*%Sl[[b]]$D else if (l == 2) X[ind]%*%t(Sl[[b]]$D) else
                        if (l == -1) { if (nuDi) X[ind]%*%t(Sl[[b]]$D) else X[ind]%*%Sl[[b]]$Di } else
			{ if (nuDi) X[ind]%*%Sl[[b]]$D else X[ind]%*%t(Sl[[b]]$Di)}			
    }
  }
  X
} ## Sl.inirep

Sl.initial.repara <- function(Sl,X,inverse=FALSE,both.sides=TRUE,cov=TRUE,nt=1) {
## Routine to apply initial Sl re-parameterization to model matrix X,
## or, if inverse==TRUE, to apply inverse re-para to parameter vector 
## or cov matrix. If inverse is TRUE and both.sides=FALSE then 
## re-para only applied to rhs, as appropriate for a choleski factor.
## If both.sides==FALSE, X is a vector and inverse==FALSE then X is
## taken as a coefficient vector (so re-para is inverse of that for model
## matrix...)
  if (length(Sl)==0) return(X) ## nothing to do
  if (inverse) { ## apply inverse re-para
    if (is.matrix(X)) { 
      if (cov) { ## then it's a covariance matrix
        for (b in 1:length(Sl)) if (Sl[[b]]$repara) { 
          ind <- Sl[[b]]$start:Sl[[b]]$stop
          if (is.matrix(Sl[[b]]$D)) { 
            if (both.sides) X[ind,] <- if (nt==1) Sl[[b]]$D%*%X[ind,,drop=FALSE] else 
                                       pmmult(Sl[[b]]$D,X[ind,,drop=FALSE],FALSE,FALSE,nt=nt)
            X[,ind] <- if (nt==1) X[,ind,drop=FALSE]%*%t(Sl[[b]]$D) else
                       pmmult(X[,ind,drop=FALSE],Sl[[b]]$D,FALSE,TRUE,nt=nt)
          } else { ## Diagonal D
            X[,ind] <- t(Sl[[b]]$D * t(X[,ind,drop=FALSE]))
            if (both.sides) X[ind,] <- Sl[[b]]$D * X[ind,,drop=FALSE]
          } 
        } 
      } else { ## regular matrix: need to use Di
         for (b in 1:length(Sl)) if (Sl[[b]]$repara) { 
          ind <- Sl[[b]]$start:Sl[[b]]$stop
          if (is.matrix(Sl[[b]]$D)) { 
            Di <- if(is.null(Sl[[b]]$Di)) t(Sl[[b]]$D) else Sl[[b]]$Di
            if (both.sides) X[ind,] <- if (nt==1) t(Di)%*%X[ind,,drop=FALSE] else
                            pmmult(Di,X[ind,,drop=FALSE],TRUE,FALSE,nt=nt)
            X[,ind] <- if (nt==1) X[,ind,drop=FALSE]%*%Di else
                       pmmult(X[,ind,drop=FALSE],Di,FALSE,FALSE,nt=nt)
          } else { ## Diagonal D
            Di <- 1/Sl[[b]]$D
            X[,ind] <- t(Di * t(X[,ind,drop=FALSE]))
            if (both.sides) X[ind,] <- Di * X[ind,,drop=FALSE]
          } 
        } 
      }
    } else { ## it's a parameter vector
      for (b in 1:length(Sl)) if (Sl[[b]]$repara) { 
        ind <- Sl[[b]]$start:Sl[[b]]$stop
        if (is.matrix(Sl[[b]]$D)) X[ind] <- Sl[[b]]$D%*%X[ind] else 
        X[ind] <- Sl[[b]]$D*X[ind] 
      }
    }
  } else for (b in 1:length(Sl)) if (Sl[[b]]$repara) { ## model matrix re-para
    ind <- Sl[[b]]$start:Sl[[b]]$stop
    if (is.matrix(X)) { 
      if (is.matrix(Sl[[b]]$D)) { 
        if (both.sides)  X[ind,] <- if (nt==1) t(Sl[[b]]$D)%*%X[ind,,drop=FALSE] else 
                         pmmult(Sl[[b]]$D,X[ind,,drop=FALSE],TRUE,FALSE,nt=nt)
        X[,ind] <- if (nt==1) X[,ind,drop=FALSE]%*%Sl[[b]]$D else
                   pmmult(X[,ind,drop=FALSE],Sl[[b]]$D,FALSE,FALSE,nt=nt)
      } else { 
        if (both.sides) X[ind,] <- Sl[[b]]$D * X[ind,,drop=FALSE]
        X[,ind] <- t(Sl[[b]]$D*t(X[,ind,drop=FALSE])) ## X[,ind]%*%diag(Sl[[b]]$D)
      }
    } else {
      if (both.sides) { ## signalling vector to be treated like model matrix X... 
        if (is.matrix(Sl[[b]]$D)) X[ind] <- t(Sl[[b]]$D)%*%X[ind] else 
        X[ind] <- Sl[[b]]$D*X[ind]
      } else { ## both.sides == FALSE is just a signal that X is a parameter vector
        if (is.matrix(Sl[[b]]$D)) X[ind] <-
	  if (is.null(Sl[[b]]$Di)) t(Sl[[b]]$D)%*%X[ind] else Sl[[b]]$Di%*%X[ind] else
        X[ind] <- X[ind]/Sl[[b]]$D
      }
    }
  }
  X
} ## end Sl.initial.repara




ldetSblock <- function(rS,rho,deriv=2,root=FALSE,nt=1) {
## finds derivatives wrt rho of log|S| where 
## S = sum_i tcrossprod(rS[[i]]*exp(rho[i]))
## when S is full rank +ve def and no 
## reparameterization is required....
  lam <- exp(rho)
  S <- pcrossprod(rS[[1]],trans=TRUE,nt=nt)*lam[1]
      ##tcrossprod(rS[[1]])*lam[1] ## parallel
  p <- ncol(S)
  m <- length(rS)
  if (m > 1) for (i in 2:m) S <- S + pcrossprod(rS[[i]],trans=TRUE,nt=nt)*lam[i]
  ## S <- S + tcrossprod(rS[[i]])*lam[i] ## parallel
  if (!root) E <- S
  d <- diag(S);d[d<=0] <- 1;d <- sqrt(d)
  S <- t(S/d)/d ## diagonally pre-condition
  R <- if (nt>1) pchol(S,nt) else suppressWarnings(chol(S,pivot=TRUE))
  piv <- attr(R,"pivot")
  r <- attr(R,"rank")
  if (r<p) R[(r+1):p,(r+1):p] <- 0 ## fix chol bug
  if (root) {
    rp <- piv;rp[rp] <- 1:p ## reverse pivot
    E <- t(t(R[,rp])*d)
  } 
  if (r<p) { ## rank deficiency
    R <- R[1:r,1:r]
    piv <- piv[1:r]
  }
  RrS <- list();dS1 <- rep(0,m);dS2 <- matrix(0,m,m)
  ## use dlog|S|/drhoi = lam_i tr(S^{-1}S_i) = tr(R^{-T}rS[[i]]rS[[i]]R^{-1} etc...
  for (i in 1:m) {
    RrS[[i]] <- pforwardsolve(R,rS[[i]][piv,]/d[piv],nt=nt) ## note R transposed internally - unlike forwardsolve!!
    dS1[i] <- sum(RrS[[i]]^2)*lam[i] ## dlog|S|/drhoi
    if (deriv==2) { 
      RrS[[i]] <- pcrossprod(RrS[[i]],trans=TRUE,nt=nt)
      #tcrossprod(RrS[[i]]) ## parallel
      for (j in 1:i) {
        dS2[i,j] <- dS2[j,i] <- -sum(RrS[[i]]*RrS[[j]])*lam[i]*lam[j]
      }
      dS2[i,i] <- dS2[i,i] + dS1[i]
    }
  }
  list(det = 2*sum(log(diag(R))+log(d[piv])),det1=dS1,det2=dS2,E=E)
} ## ldetSblock

ginv1 <- function(a) {
## Cholesky generalized inverse and log det function
  eps <- norm(a)*.Machine$double.eps^.8
  da <- diag(a+eps)^-.5
  R <- suppressWarnings(chol(t(da*t(da*a)),pivot=TRUE))
  piv <- attr(R,"pivot")
  r <- attr(R,"rank")
  r <- Rrank(R)
  b <- a*0
  E <- R[1:r,]*0
  E[,piv] <- R[1:r,]
  E <- t(t(E)/da)
  piv <- piv[1:r]
  b[piv,piv] <- chol2inv(R[1:r,1:r])
  b <- t(t(da*b)*da)
  ldet <- sum(log(diag(R)[1:r])) - sum(log(da))
  list(inv = b, ldet = ldet,E=E)
} ## ginv1

ldetSt <- function(S,lam,deriv=0,repara=TRUE) {
## completely Cholesky based computation of the determinant of
## an additive block. Assumes initial reparameterization using
## cholesky approach in Sl.setup.
## Formally, St, the sum of S[[i]]*lam[i] must be full rank.
## The S[[i]] are not all full rank (or this is pointless)
## Returns a transformed St which can be used for stable
## determinant evaluation by QR or Cholesky *without pivoting*
## --- pivoting the Cholesky will mess up the very structure
## that ensures stable computation. Also returns the matrix
## mapping original St to its transformed version, and the
## corresponding inverse transform. The returned S is such that
## the sum of S[[i]]*lam[i] is *exactly* the transformed St.
## Set repara to FALSE to get determinant computations without
## stabilizing transforms.
## NOTE: If T is the transformation matrix, log|T| is ommited from
##       returned log|S| as if will cancel with log|X'X+S| in
##       LAML. Can be re-instated by commenting in below. 
  dominant.set <- function(nos,nD) {
  ## on entry 'nD' is the index of the current set and 'nos' the
  ## norms for these. On exit 'D' contains the indices (from nD)
  ## of the dominant terms, and nD the remainder...
    j <- which(nos>max(nos)*1e-5)
    D <- nD[j]; nD <- nD[-j]
    return(list(D=D,nD=nD))
  } ## dominant.set

  Rldet0 <- 0
  m <- length(S)
  nD <- 1:m ## index of terms not yet dealt with
  a <- 1 ## starting row/col
  p <- ncol(S[[1]]) ## dimension
  T <- Ti <- diag(p) 
  while (repara&&length(nD)>1 && a <= p) {
    ## get indices of dominant terms and remainder
    nos <- rep(0,length(nD)) ## for dominance determining norms 
    j <- 1
    for (i in nD) { nos[j] <- norm(S[[i]][a:p,a:p,drop=FALSE])*lam[i]; j <- j + 1}
    ds <- dominant.set(nos,nD) 
    nD <- ds$nD ## not dominant set
    D <- ds$D   ## dominant set
    ## Form the dominant term and its pivoted Cholesky
    k <- p-a+1 ## current block dimension
    Sd <- matrix(0,k,k)
    for (i in D) Sd <- Sd + lam[i] * S[[i]][a:p,a:p,drop=FALSE]
    R <- suppressWarnings(chol(Sd,pivot=TRUE))
    rank <- min(Rrank(R),attr(R,"rank"))
    piv <- attr(R,"pivot")
    ipiv <- piv; ipiv[piv] <- 1:k
    Sp <- Sd[piv,piv,drop=FALSE]; normS <- norm(Sp)
    ## more expensive refinement of rank... 
    while (rank>1&&max(abs(crossprod(R[1:(rank-1),rank:k,drop=FALSE])-Sp[rank:k,rank:k,drop=FALSE]))<.Machine$double.eps^.75*normS) rank <- rank - 1

    if (rank < k) {
      ind <- (rank+1):k
      R[ind,ind] <- diag(length(ind)) ## augment to full rank k
    }
    ## Apply transform to components of D, suppressing definite machine zeroes
    for (i in D) {
      S[[i]][,a:p] <- t(forwardsolve(t(R),t(S[[i]][,a:p,drop=FALSE][,piv,drop=FALSE]))) ## SRi
      S[[i]][a:p,] <- forwardsolve(t(R),S[[i]][a:p,,drop=FALSE][piv,,drop=FALSE]) ## Ri'S
      if (rank < k) S[[i]][a:p,a:p][ind,] <- S[[i]][a:p,a:p][,ind] <- 0   
    }
    ## Apply transform to components of nD
    for (i in nD) {
      S[[i]][,a:p] <- t(forwardsolve(t(R),t(S[[i]][,a:p,drop=FALSE][,piv,drop=FALSE]))) ## SRi
      S[[i]][a:p,] <- forwardsolve(t(R),S[[i]][a:p,,drop=FALSE][piv,,drop=FALSE]) ## Ri'S
    }
    ## Update the total transform matrix, its log determinant and inverse...
    Rldet0 <- Rldet0 + sum(log(diag(R)))
    ## Accumulate T such that |sum_i lam_i*S_i| = |T' sum_i lam_i * St_i T|
    ## St_i being transformed versions...
    Ti[,a:p] <- t(forwardsolve(t(R),t(Ti[,a:p,drop=FALSE][,piv,drop=FALSE]))) ## this is inverse
    T[a:p,] <- R %*% T[a:p,,drop=FALSE][piv,,drop=FALSE]
    a <- a + rank 
  } ## finished transforming
  ## compute the log determinant
  St <- matrix(0,p,p)
  for (i in 1:m) St <- St + lam[i]*S[[i]]
  if (repara) {
    E <- R1 <- chol(St)
    Rldet <- sum(log(diag(R1))) # + Rldet ## note: no log|T| - cancels in REML
  } else { ## use stabilized generalized inverse
    gi <- ginv1(St)
    E <- gi$E
    Rldet <- gi$ldet #+ Rldet ## note: no log|T| - cancels in REML
  }
  det1 <- det2 <- NULL
  if (deriv>0) {
    R1 <- if (repara) chol2inv(R1) else gi$inv
    det1 <- lam*0
    for (i in 1:m) det1[i] <- sum(R1*S[[i]]*lam[i])
  }
  if (deriv>1) {
    SiS <- list()
    det2 <- matrix(0,m,m)
    for (i in 1:m) {
      SiS[[i]] <- R1 %*% S[[i]]
      for (j in 1:i) det2[i,j] <- det2[j,i] <- -sum(SiS[[i]]*t(SiS[[j]]))*lam[i]*lam[j]
      det2[i,i] <- det2[i,i] + det1[i]
    }  
  }
  list(det=2*Rldet,det0=2*Rldet0,T=T,S=S,Ti=Ti,det1=det1,det2 = det2,St=St,E=E,kappa=kappa(St))
} ## ldetSt


ldetS <- function(Sl,rho,fixed,np,root=FALSE,Stot=FALSE,repara=TRUE,nt=1,deriv=2,sparse=FALSE) {
## Get log generalized determinant of S stored blockwise in an Sl list.
## If repara=TRUE multi-term blocks will be re-parameterized using gam.reparam, and
## a re-parameterization object supplied in the returned object.
## rho contains log smoothing parameters, fixed is an array indicating whether they 
## are fixed (or free). np is the number of coefficients. root indicates
## whether or not to return E, and sparse whether or not it should be sparse. 
## Returns: Sl, with modified rS terms, if needed and rho added to each block
##          rp, a re-parameterization list
##          E a total penalty square root such that E'E = S_tot (if root==TRUE)
##          S the total penalty matrix if Stot==TRUE.
##          ldetS,ldetS1,ldetS2 the value, grad vec and Hessian
  n.deriv <- sum(!fixed)
  k.deriv <- k.sp <- k.rp <- 1
  ldS <- 0
  d1.ldS <- rep(0,n.deriv)
  d2.ldS <- matrix(0,n.deriv,n.deriv)
  cholesky <- attr(Sl,"cholesky")
  rp <- list() ## reparameterization list

  ## In dense case E is a dense matrix and its block get filled in. In sparse
  ## case E is a list of lists of rows (i), cols (j) and elements (x) defining each
  ## block. So E$i[[b]],E$j[[b]],E$x[[b]] defines block b...

  if (root) { E <- if (sparse) list(i=list(),j=list(),x=list()) else matrix(0,np,np) } else E <- NULL
  if (Stot) { S <- if (sparse) list(i=list(),j=list(),x=list()) else matrix(0,np,np) } else S <- NULL

  if (length(Sl)>0) for (b in 1:length(Sl)) { ## work through blocks
    
    if (!Sl[[b]]$linear) { ## non-linear block
      ind <- k.sp + 1:Sl[[b]]$n.sp - 1 ## smoothing param index
      Sl[[b]] <- Sl[[b]]$updateS(rho[ind],Sl[[b]]) ## update the block with current params 
      nldS <- Sl[[b]]$ldS(Sl[[b]],deriv) ## get the log determinant and any derivatives
      ldS <- ldS + nldS$ldS
      nldS$ldS1 <- nldS$ldS1[!fixed[ind]] ## discard fixed param derivatives
      
      nderiv <- length(nldS$ldS1)
      if (nderiv) d1.ldS[k.deriv+1:nderiv-1] <- nldS$ldS1
      if (deriv>1) {
        nldS$ldS2 <- nldS$ldS2[!fixed[ind],!fixed[ind]]
	d2.ldS[k.deriv+1:nderiv-1,k.deriv+1:nderiv-1] <- nldS$ldS2
      }
      k.deriv <- k.deriv + nderiv
      k.sp <- k.sp + Sl[[b]]$n.sp
      Sl[[b]]$lambda <- rho[ind] ## not really used in non-linear interface
      if (root||Stot) {
        tmp <- if (root&&Stot) 2 else if (root) 1 else 0
        tmp <- Sl[[b]]$St(Sl[[b]],tmp)
      }	
      if (root) if (sparse) {
        ## E0 <- as(Sl[[b]]$St(Sl[[b]],1)$E,"dgTMatrix") - deprecated
	E0 <- as(as(as(tmp$E, "dMatrix"), "generalMatrix"), "TsparseMatrix")
	E$i[[b]] <- E0@i + Sl[[b]]$start
	E$j[[b]] <- E0@j + Sl[[b]]$start
	E$x[[b]] <- E0@x
      } else {
        ind <- Sl[[b]]$start:Sl[[b]]$stop
	E[ind,ind] <- tmp$E
      }
      if (Stot) if (sparse) {
        ## E0 <- as(Sl[[b]]$St(Sl[[b]],1)$E,"dgTMatrix") - deprecated
	E0 <- as(as(as(tmp$S, "dMatrix"), "generalMatrix"), "TsparseMatrix")
	S$i[[b]] <- E0@i + Sl[[b]]$start
	S$j[[b]] <- E0@j + Sl[[b]]$start
	S$x[[b]] <- E0@x
      } else {
        ind <- Sl[[b]]$start:Sl[[b]]$stop
	S[ind,ind] <- tmp$S
      }
      
    } else if (length(Sl[[b]]$S)==1) { ## linear singleton
      ldS <- ldS + Sl[[b]]$ldet ## initial repara log det correction for this block
      ldS <- ldS + if (is.null(Sl[[b]]$nl.reg)) rho[k.sp] * Sl[[b]]$rank else
                   sum(log(exp(rho[k.sp])*Sl[[b]]$ev + Sl[[b]]$nl.reg))
      if (!fixed[k.sp]) {
        d1.ldS[k.deriv] <- if (is.null(Sl[[b]]$nl.reg)) Sl[[b]]$rank else
	        sum(exp(rho[k.sp])*Sl[[b]]$ev/(exp(rho[k.sp])*Sl[[b]]$ev + Sl[[b]]$nl.reg))
	k.deriv <- k.deriv + 1
      } 
      if (root) { 
        ## insert diagonal from block start to end
        if (Sl[[b]]$repara) {
	  ind <- (Sl[[b]]$start:Sl[[b]]$stop)[Sl[[b]]$ind]
	  if (sparse) {
	    E$i[[b]] <- E$j[[b]] <- ind; E$x[[b]] <- ind*0 + exp(rho[k.sp]*.5)
	  } else {
            ## diag(E)[ind] <- exp(rho[k.sp]*.5) ## sqrt smoothing param
	    .Call(C_wdiag,E,ind,rep(exp(rho[k.sp]*.5),length(ind)))
	  }  
        } else { ## root has to be in original parameterization...
          if (sparse) {
	    ## dgTMatrix is triplet form, which makes combining easier...
	    D <- if (is.null(Sl[[b]]$nl.reg)) as(as(as(Sl[[b]]$Di[1:Sl[[b]]$rank,]* exp(rho[k.sp]*.5), "dMatrix"),
	         "generalMatrix"), "TsparseMatrix") else as(as(as(sqrt(Sl[[b]]$ev*exp(rho[k.sp])+Sl[[b]]$nl.reg)*
		 t(Sl[[b]]$U) , "dMatrix"), "generalMatrix"), "TsparseMatrix")	 
	    E$i[[b]] <- D@i + Sl[[b]]$start
	    E$j[[b]] <- D@j + Sl[[b]]$start
	    E$x[[b]] <- D@x
	  } else {
            D <- if (is.null(Sl[[b]]$nl.reg)) Sl[[b]]$Di[1:Sl[[b]]$rank,]* exp(rho[k.sp]*.5) else
	         sqrt(Sl[[b]]$ev*exp(rho[k.sp])+Sl[[b]]$nl.reg)*t(Sl[[b]]$U)
            indc <- Sl[[b]]$start:(Sl[[b]]$start+ncol(D)-1)
            indr <- Sl[[b]]$start:(Sl[[b]]$start+nrow(D)-1)
            E[indr,indc] <- D 
	  }  
        }
      } ## if (root)
      if (Stot) { 
        ## insert diagonal from block start to end
        if (Sl[[b]]$repara) {
	  ind <- (Sl[[b]]$start:Sl[[b]]$stop)[Sl[[b]]$ind]
	  if (sparse) {
	    S$i[[b]] <- S$j[[b]] <- ind; S$x[[b]] <- ind*0 + exp(rho[k.sp])
	  } else {
            ## diag(S)[ind] <- exp(rho[k.sp]) ## smoothing param
	    .Call(C_wdiag,S,ind,rep(exp(rho[k.sp]),length(ind)))
	  }  
        } else { ## root has to be in original parameterization...
          if (sparse) {
	    ## dgTMatrix is triplet form, which makes combining easier...
	    D <- if (is.null(Sl[[b]]$nl.reg)) as(as(as(crossprod(Sl[[b]]$Di[1:Sl[[b]]$rank,]) * exp(rho[k.sp]), "dMatrix"),
	         "generalMatrix"), "TsparseMatrix") else as(as(as(Sl[[b]]$U%*%((Sl[[b]]$ev*exp(rho[k.sp])+Sl[[b]]$nl.reg)*
		 t(Sl[[b]]$U)), "dMatrix"), "generalMatrix"), "TsparseMatrix")	 
	    S$i[[b]] <- D@i + Sl[[b]]$start
	    S$j[[b]] <- D@j + Sl[[b]]$start
	    S$x[[b]] <- D@x
	  } else {
            D <- if (is.null(Sl[[b]]$nl.reg)) crossprod(Sl[[b]]$Di[1:Sl[[b]]$rank,]) * exp(rho[k.sp]) else
	         Sl[[b]]$U%*%((Sl[[b]]$ev*exp(rho[k.sp])+Sl[[b]]$nl.reg)*t(Sl[[b]]$U))
            indc <- Sl[[b]]$start:(Sl[[b]]$start+ncol(D)-1)
            indr <- Sl[[b]]$start:(Sl[[b]]$start+nrow(D)-1)
            S[indr,indc] <- D 
	  }  
        }
      }
      Sl[[b]]$lambda <- exp(rho[k.sp])
      k.sp <- k.sp + 1 
    } else { ## linear multi-S block
      ldS <- ldS + Sl[[b]]$ldet ## initial repara log det correction for this block
      m <- length(Sl[[b]]$S) ## number of components for this block
      ind <- k.sp:(k.sp+m-1) ## index for smoothing parameters
      ## call gam.reparam to deal with this block
      ## in a stable way...
      if (cholesky) {
	grp <-  if (Sl[[b]]$repara) ldetSt(Sl[[b]]$S,lam=exp(rho[ind]),deriv=deriv,repara) else
	        ldetSt(Sl[[b]]$rS,lam=exp(rho[ind]),deriv=deriv,repara)  
      } else {
        grp <- if (repara) gam.reparam(Sl[[b]]$rS,lsp=rho[ind],deriv=deriv) else 
             ldetSblock(Sl[[b]]$rS,rho[ind],deriv=deriv,root=root,nt=nt)
	grp$St <- if (repara) grp$S else grp$E     
	grp$det0 <- 0     
      }
      Sl[[b]]$lambda <- exp(rho[ind])
      ## If stabilizing repara not applied to coefs, then need to correct log det
      ## for transform, as it won't cancel, hence grp$det0...
      ldS <- ldS + grp$det + if (Sl[[b]]$repara) 0 else grp$det0 
      ## next deal with the derivatives...
      grp$det1 <- grp$det1[!fixed[ind]] ## discard derivatives for fixed components
      grp$det2 <- if (deriv>1) grp$det2[!fixed[ind],!fixed[ind]] else 0 ##NULL
      nd <- length(grp$det1)
      if (nd>0) { ## then not all sp's are fixed
        dind <- k.deriv:(k.deriv+nd-1)
        d1.ldS[dind] <- grp$det1
        d2.ldS[dind,dind] <- grp$det2
        k.deriv <- k.deriv + nd
      }
      ## now store the reparameterization information   
      if (repara) {
        ## note that Ti is equivalent to Qs...
        rp[[k.rp]] <- if (cholesky) list(block =b,ind = (Sl[[b]]$start:Sl[[b]]$stop)[Sl[[b]]$ind],T=grp$T,
	              Ti=grp$Ti,repara=Sl[[b]]$repara) else list(block =b,ind = (Sl[[b]]$start:Sl[[b]]$stop)[Sl[[b]]$ind],
	              Qs = grp$Qs,repara=Sl[[b]]$repara) 
        k.rp <- k.rp + 1
        for (i in 1:m) {
          Sl[[b]]$Srp[[i]] <- if (cholesky) Sl[[b]]$lambda[i]*grp$S[[i]] else
	                      Sl[[b]]$lambda[i]*grp$rS[[i]]%*%t(grp$rS[[i]])
        }
      }
      k.sp <- k.sp + m
      if (Sl[[b]]$repara) {
        if (root) { ## unpack the square root E'E
          if (sparse) {
            # E0 <- as(grp$E,"dgTMatrix") deprecated
	    E0 <- as(as(as(grp$E, "dMatrix"), "generalMatrix"), "TsparseMatrix")
	    E$i[[b]] <- E0@i + Sl[[b]]$start
	    E$j[[b]] <- E0@j + Sl[[b]]$start
	    E$x[[b]] <- E0@x
          } else {
	    ic <- Sl[[b]]$start:(Sl[[b]]$start+ncol(grp$E)-1)
            ir <- Sl[[b]]$start:(Sl[[b]]$start+nrow(grp$E)-1)
            E[ir,ic] <- grp$E
	  }  
          #Sl[[b]]$St <- crossprod(grp$E)
	  Sl[[b]]$St <- grp$St
        } else {  
          ## gam.reparam always returns root penalty in E, but 
          ## ldetSblock returns penalty in E if root==FALSE
	  #if (cholesky) Sl[[b]]$St <- grp$St else
          #Sl[[b]]$St <- if (repara) crossprod(grp$E) else grp$E
	  Sl[[b]]$St <- grp$St
        }
	if (Stot) { ## unpack the square root E'E
          if (sparse) {
	    E0 <- as(as(as(Sl[[b]]$St, "dMatrix"), "generalMatrix"), "TsparseMatrix")
	    S$i[[b]] <- E0@i + Sl[[b]]$start
	    S$j[[b]] <- E0@j + Sl[[b]]$start
	    S$x[[b]] <- E0@x
          } else {
	    ic <- Sl[[b]]$start:(Sl[[b]]$start+ncol(grp$St)-1)
            ir <- Sl[[b]]$start:(Sl[[b]]$start+nrow(grp$St)-1)
            S[ir,ic] <- grp$St
	  }  
        }
      } else { ## square root block and St need to be in original parameterization...
        Sl[[b]]$St <- Sl[[b]]$lambda[1]*Sl[[b]]$S[[1]]
        for (i in 2:m) {
          Sl[[b]]$St <- Sl[[b]]$St + Sl[[b]]$lambda[i]*Sl[[b]]$S[[i]]
        }
        if (Stot) {
	  if (sparse) {
	    Eb <- as(as(as(Sl[[b]]$St, "dMatrix"), "generalMatrix"), "TsparseMatrix")
	    S$i[[b]] <- Eb@i + Sl[[b]]$start
	    S$j[[b]] <- Eb@j + Sl[[b]]$start
	    S$x[[b]] <- Eb@x
	  } else {
            indc <- Sl[[b]]$start:(Sl[[b]]$start+ncol(Sl[[b]]$St)-1)
            indr <- Sl[[b]]$start:(Sl[[b]]$start+nrow(Sl[[b]]$St)-1)
            S[indr,indc] <- Sl[[b]]$St
	  }  
        }
	if (root) {
	  Eb <- t(mroot(Sl[[b]]$St,Sl[[b]]$rank))
	  if (sparse) {
	    # Eb <- as(Eb,"dgTMatrix") - deprecated
	    Eb <- as(as(as(Eb, "dMatrix"), "generalMatrix"), "TsparseMatrix")
	    E$i[[b]] <- Eb@i + Sl[[b]]$start
	    E$j[[b]] <- Eb@j + Sl[[b]]$start
	    E$x[[b]] <- Eb@x
	  } else {
            indc <- Sl[[b]]$start:(Sl[[b]]$start+ncol(Eb)-1)
            indr <- Sl[[b]]$start:(Sl[[b]]$start+nrow(Eb)-1)
            E[indr,indc] <- Eb
	  }  
        } 
      }   
    } ## end of multi-S block branch
  } ## end of block loop
  if (root) E <- if (sparse) sparseMatrix(i=unlist(E$i),j=unlist(E$j),x=unlist(E$x),dims=c(np,np)) else
                             E[rowSums(abs(E))!=0,,drop=FALSE] ## drop zero rows.
  if (Stot&&sparse) S <- sparseMatrix(i=unlist(S$i),j=unlist(S$j),x=unlist(S$x),dims=c(np,np)) 			     
  list(ldetS=ldS,ldet1=d1.ldS,ldet2=d2.ldS,Sl=Sl,rp=rp,E=E,S=S)
} ## end ldetS


Sl.addS <- function(Sl,A,rho) {
## Routine to add total penalty to matrix A. Sl is smooth penalty
## list from Sl.setup, so initial reparameterizations have taken place,
## and should have already been applied to A using Sl.initial.repara
  k <- 1
  A <- A*1 ## force a copy to be made so that A not modified in calling env!!
  if (length(Sl)>0) for (b in 1:length(Sl)) {
    ind <- (Sl[[b]]$start:Sl[[b]]$stop)[Sl[[b]]$ind] 
    if (length(Sl[[b]]$S)==1) { ## singleton
      B <- exp(rho[k]);diag <- -1
      dummy <- .Call(C_mgcv_madi,A,B,ind,diag)
      ## diag(A)[ind] <-  diag(A)[ind] + exp(rho[k]) ## penalty is identity times sp
      k <- k + 1
    } else {
      for (j in 1:length(Sl[[b]]$S)) {
        B <- exp(rho[k]) * Sl[[b]]$S[[j]]; diag <- 0
        .Call(C_mgcv_madi,A,B,ind,diag)
        ## A[ind,ind] <- A[ind,ind] + exp(rho[k]) * Sl[[b]]$S[[j]]
        k <- k + 1
      }
    }
  }
  A
} ## Sl.addS

Sl.addS0 <- function(Sl,A,rho) {
## Routine to add total penalty to matrix A. Sl is smooth penalty
## list from Sl.setup, so initial reparameterizations have taken place,
## and should have already been applied to A using Sl.initial.repara
## inefficient prototype of Sl.addS 
  k <- 1
  for (b in 1:length(Sl)) {
    ind <- (Sl[[b]]$start:Sl[[b]]$stop)[Sl[[b]]$ind] 
    if (length(Sl[[b]]$S)==1) { ## singleton
      diag(A)[ind] <-  diag(A)[ind] + exp(rho[k]) ## penalty is identity times sp
      k <- k + 1
    } else {
      for (j in 1:length(Sl[[b]]$S)) {
        A[ind,ind] <- A[ind,ind] + exp(rho[k]) * Sl[[b]]$S[[j]]
        k <- k + 1
      }
    }
  }
  A
} ## Sl.addS0

Sl.repa <- function(rp,X,l=0,r=0) {
## Re-parameterize X using rp reparameterization info.
## l,r = -2,-1,0,1,2. O is do not apply, negative to apply inverse transform Di,
##       positive for transform D, 1 for transform, 2 for its transpose.
## If b' is transformed and b  orginal. Di b' = b and b' = D b.
## If T present D=T and Di = Ti. Otherwise D = t(Qs) and Di = Qs.
## Aim is for simpler and cleaner than Sl.repara
  nr <- length(rp);if (nr==0) return(X)
  for (i in 1:nr) if (rp[[i]]$repara) {
    if (l) {
      T <- if (is.null(rp[[i]]$Qs)) { if (l<0) { if (l==-2) t(rp[[i]]$Ti) else rp[[i]]$Ti
           } else { if (l==2) t(rp[[i]]$T) else rp[[i]]$T }} else { if (l<0) { if (l==-2) t(rp[[i]]$Qs) else rp[[i]]$Qs
           } else { if (l==2) rp[[i]]$Qs else t(rp[[i]]$Qs) }}
     if (is.matrix(X)) X[rp[[i]]$ind,] <- T %*% X[rp[[i]]$ind,] else X[rp[[i]]$ind] <- T %*% X[rp[[i]]$ind]	   	   
    }
    if (r) {
      T <- if (is.null(rp[[i]]$Qs)) { if (r<0) { if (r==-2) t(rp[[i]]$Ti) else rp[[i]]$Ti
           } else { if (r==2) t(rp[[i]]$T) else rp[[i]]$T }} else {if (r<0) { if (r==-2) t(rp[[i]]$Qs) else rp[[i]]$Qs
           } else { if (r==2) rp[[i]]$Qs else t(rp[[i]]$Qs) }}
      if (is.matrix(X)) X[,rp[[i]]$ind] <- X[,rp[[i]]$ind]%*%T else X[rp[[i]]$ind] <- t(X[rp[[i]]$ind])%*%T	   	   
    }
  }
  X
} ## Sl.repa

Sl.repara <- function(rp,X,inverse=FALSE,both.sides=TRUE) {
## Apply re-parameterization from ldetS to X, blockwise.
## If X is a matrix it is assumed to be a model matrix
## whereas if X is a vector it is assumed to be a parameter vector.
## If inverse==TRUE applies the inverse of the re-para to
## parameter vector X or cov matrix X...
## beta_trans = Ti beta_original T is inverse Ti
  nr <- length(rp);if (nr==0) return(X)
  if (inverse) {
    if (is.matrix(X)) { ## X is a cov matrix
      for (i in 1:nr) if (rp[[i]]$repara) {
        if (both.sides) X[rp[[i]]$ind,]  <- if (is.null(rp[[i]]$Qs)) rp[[i]]$Ti %*% X[rp[[i]]$ind,,drop=FALSE] else
                     rp[[i]]$Qs %*% X[rp[[i]]$ind,,drop=FALSE]
        X[,rp[[i]]$ind]  <-  if (is.null(rp[[i]]$Qs)) X[,rp[[i]]$ind,drop=FALSE] %*% t(rp[[i]]$Ti) else
                     X[,rp[[i]]$ind,drop=FALSE] %*% t(rp[[i]]$Qs)
      }
    } else { ## X is a vector
     for (i in 1:nr) if (rp[[i]]$repara) X[rp[[i]]$ind]  <- if (is.null(rp[[i]]$Qs))
                     rp[[i]]$Ti %*% X[rp[[i]]$ind] else rp[[i]]$Qs %*% X[rp[[i]]$ind]
    }
  } else { ## apply re-para to X
    if (is.matrix(X)) {
      for (i in 1:nr) if (rp[[i]]$repara) X[,rp[[i]]$ind]  <- if (is.null(rp[[i]]$Qs))
                      X[,rp[[i]]$ind]%*%rp[[i]]$Ti else X[,rp[[i]]$ind]%*%rp[[i]]$Qs
    } else {
      for (i in 1:nr) if (rp[[i]]$repara) X[rp[[i]]$ind]  <- if (is.null(rp[[i]]$Qs))
      rp[[i]]$T %*% X[rp[[i]]$ind] else t(rp[[i]]$Qs) %*% X[rp[[i]]$ind]
    }
  }
  X
} ## end Sl.repara

Sl.mult <- function(Sl,A,k = 0,full=TRUE) {
## Sl contains the blocks of block diagonal penalty S.
## If k<=0 this routine forms S%*%A.
## If k>0 then the routine forms S_k%*%A, zero padded 
## if full==TRUE, but in smallest number of rows form otherwise.
  nb <- length(Sl) ## number of blocks
  if (nb==0) return(A*0)
  Amat <- is.matrix(A) 
  if (k<=0) { ## apply whole penalty
    B <- A*0
    for (b in 1:nb) { ## block loop
      ind <- Sl[[b]]$start:Sl[[b]]$stop ## index of coeffs for this bock
      if (!Sl[[b]]$linear) { ## non-linear block
        if (Amat)  B[ind,] <- t(Sl[[b]]$AS(t(A[ind,]),Sl[[b]])) else
	           B[ind] <- drop(Sl[[b]]$AS(A[ind],Sl[[b]]))
      } else if (length(Sl[[b]]$S)==1) { ## singleton
        if (Sl[[b]]$repara) {
          ind <- ind[Sl[[b]]$ind]
          if (Amat) B[ind,] <- Sl[[b]]$lambda*A[ind,] else
	            B[ind] <- Sl[[b]]$lambda*A[ind]
        } else { ## original penalty has to be applied 
          if (Amat) {
	    B[ind,] <- Sl[[b]]$lambda*Sl[[b]]$S[[1]] %*% A[ind,] + if (is.null(Sl[[b]]$nl.reg)) 0 else Sl[[b]]$nl.reg * A[ind,]
	  } else {
            B[ind] <- Sl[[b]]$lambda*Sl[[b]]$S[[1]] %*% A[ind] + if (is.null(Sl[[b]]$nl.reg)) 0 else Sl[[b]]$nl.reg * A[ind]
	  }  
        }
      } else { ## multi-S block
	if (Sl[[b]]$repara) ind <- ind[Sl[[b]]$ind]
        if (Amat) B[ind,] <- Sl[[b]]$St %*% A[ind,] else
                  B[ind] <- Sl[[b]]$St %*% A[ind]
      }
    } ## end of block loop
    A <- B
  } else { ## single penalty matrix selected
    j <- 0 ## S counter
    for (b in 1:nb) { ## block loop
      for (i in 1:length(Sl[[b]]$S)) { ## S loop within blocks
        j <- j + 1
        if (j==k) { ## found block
          if (!Sl[[b]]$linear || length(Sl[[b]]$S)==1) { ## singleton
            ind <- Sl[[b]]$start:Sl[[b]]$stop
            if (Sl[[b]]$repara) {
              ind <- ind[Sl[[b]]$ind]
              if (full) { ## return zero answer with all zeroes in place
                B <- A*0
                if (Amat) B[ind,] <- Sl[[b]]$lambda*A[ind,] else  
                          B[ind] <- Sl[[b]]$lambda*A[ind]
                A <- B
              } else { ## strip zero rows from answer
                A <- if (Amat) Sl[[b]]$lambda*A[ind,] else as.numeric(Sl[[b]]$lambda*A[ind])
              } 
            } else if (Sl[[b]]$linear) { ## not reparameterized version, but linear 
              if (full) {
                B <- A*0
                if (Amat) {
		  B[ind,] <- Sl[[b]]$lambda*Sl[[b]]$S[[1]]%*%A[ind,] + if (is.null(Sl[[b]]$nl.reg)) 0 else Sl[[b]]$nl.reg * A[ind,] 
		} else {  
                  B[ind] <- Sl[[b]]$lambda*Sl[[b]]$S[[1]]%*%A[ind] + if (is.null(Sl[[b]]$nl.reg)) 0 else Sl[[b]]$nl.reg * A[ind]
	        }		  
                A <- B
              } else {
                 A <- if (Amat) {
		   Sl[[b]]$lambda*Sl[[b]]$S[[1]]%*%A[ind,] + if (is.null(Sl[[b]]$nl.reg)) 0 else Sl[[b]]$nl.reg * A[ind,]
		 } else {
                   as.numeric(Sl[[b]]$lambda*Sl[[b]]$S[[1]]%*%A[ind]) + if (is.null(Sl[[b]]$nl.reg)) 0 else Sl[[b]]$nl.reg * A[ind]
		 }  
              }
            } else { ## non-linear
              if (full) {
                B <- A*0
		if (Amat) B[ind,] <- t(Sl[[b]]$AS(t(A[ind,]),Sl)) else
	                  B[ind] <- drop(Sl[[b]]$AS(A[ind],Sl))
                A <- B
              } else {
                A <- if (Amat) t(Sl[[b]]$AS(t(A[ind,]),Sl)) else drop(Sl[[b]]$AS(A[ind],Sl))
              }
            } ## non-linear
          } else { ## multi-S block
            ind <- if (Sl[[b]]$repara) (Sl[[b]]$start:Sl[[b]]$stop)[Sl[[b]]$ind] else Sl[[b]]$start:Sl[[b]]$stop
            if (full) { ## return zero answer with all zeroes in place
              B <- A*0
              if (Amat) {
                B[ind,] <- if (is.null(Sl[[b]]$Srp)||!Sl[[b]]$repara) Sl[[b]]$lambda[i]*(Sl[[b]]$S[[i]]%*%A[ind,]) else 
                                                     Sl[[b]]$Srp[[i]]%*%A[ind,]
              } else {
                B[ind] <- if (is.null(Sl[[b]]$Srp)||!Sl[[b]]$repara) Sl[[b]]$lambda[i]*(Sl[[b]]$S[[i]]%*%A[ind]) else 
                                                    Sl[[b]]$Srp[[i]]%*%A[ind]
              }
              A <- B
            } else { ## strip zero rows from answer
              if (is.null(Sl[[b]]$Srp)||!Sl[[b]]$repara) {
                A <- if (Amat)  Sl[[b]]$lambda[i]*(Sl[[b]]$S[[i]]%*%A[ind,]) else  
                                Sl[[b]]$lambda[i]*as.numeric(Sl[[b]]$S[[i]]%*%A[ind])
              } else {
                A <- if (Amat) Sl[[b]]$Srp[[i]]%*%A[ind,] else as.numeric(Sl[[b]]$Srp[[i]]%*%A[ind])
              }
            }
          }
          break
        } 
      } ## end of S loop
      if (j==k) break
    } ## end of block loop
  }
  A
} ## end Sl.mult

Sl.termMult <- function(Sl,A,full=FALSE,nt=1) {
## returns a list containing the derivative of the penalty w.r.t
## each smoothing/variance parameter multiplied by A. For linear terms 
## this is just the product of each element S of Sl
## with A. If full==TRUE then the results include the zero rows
## otherwise these are stripped out, but in that case each element 
## of the return object contains an "ind" attribute, indicating 
## which rows of the full matrix it relates to.
  Amat <- !is.null(dim(A)) #is.matrix(A) ## allow sparse matrices also
  SA <- list()
  k <- 0 ## component counter
  nb <- length(Sl) ## number of blocks
  nli <- list()
  if (nb>0) for (b in 1:nb) { ## block loop
      if (!Sl[[b]]$linear) { ## non-linear term
        ind <- Sl[[b]]$start:Sl[[b]]$stop
	kk <- length(nli)+1
	nli[[kk]] <- c(k+1,b) ## record where nl term starts in SA 
        for (i in 1:Sl[[b]]$n.sp) { ## loop over its paramaters
	  k <- k + 1  
          if (full) {
	    B <- A*0
            if (Amat) B[ind,] <- t(Sl[[b]]$AdS(t(A[ind,,drop=FALSE]),Sl[[b]],i)) else
	              B[ind] <- drop(Sl[[b]]$AdS(A[ind],Sl[[b]],i))
            SA[[k]] <- B
	  } else {
	    SA[[k]] <- if (Amat) t(Sl[[b]]$AdS(t(A[ind,,drop=FALSE]),Sl[[b]],i)) else
	                         drop(Sl[[b]]$AdS(A[ind],Sl[[b]],i))
            attr(SA[[k]],"ind") <- ind
	  }
        }

      } else if (length(Sl[[b]]$S)==1) { ## singleton
      k <- k + 1
      ind <- Sl[[b]]$start:Sl[[b]]$stop
      if (Sl[[b]]$repara) {
        ind <- ind[Sl[[b]]$ind]
        if (full) { ## return zero answer with all zeroes in place
          B <- A*0
          if (Amat) B[ind,] <- Sl[[b]]$lambda*A[ind,,drop=FALSE] else  
                    B[ind] <- Sl[[b]]$lambda*A[ind]
          SA[[k]] <- B
        } else { ## strip zero rows from answer
          SA[[k]] <- if (Amat)  Sl[[b]]$lambda*A[ind,,drop=FALSE] else
                                drop(Sl[[b]]$lambda*A[ind])
          attr(SA[[k]],"ind") <- ind
        }
      } else {
        if (full) { ## return zero answer with all zeroes in place
          B <- A*0
          if (Amat) B[ind,] <- Sl[[b]]$lambda*Sl[[b]]$S[[1]]%*%A[ind,,drop=FALSE] else  
                    B[ind] <- Sl[[b]]$lambda*Sl[[b]]$S[[1]]%*%A[ind]
          SA[[k]] <- B
        } else { ## strip zero rows from answer
          SA[[k]] <- if (Amat)  Sl[[b]]$lambda*Sl[[b]]$S[[1]]%*%A[ind,,drop=FALSE] else
                                drop(Sl[[b]]$lambda*Sl[[b]]$S[[1]]%*%A[ind])
          attr(SA[[k]],"ind") <- ind
        }
      }
    } else { ## multi-S block
      ind <- if (Sl[[b]]$repara) (Sl[[b]]$start:Sl[[b]]$stop)[Sl[[b]]$ind] else Sl[[b]]$start:Sl[[b]]$stop
      for (i in 1:length(Sl[[b]]$S)) { ## work through S terms
        k <- k + 1
        if (full) { ## return answer with all zeroes in place
          B <- A*0
          if (is.null(Sl[[b]]$Srp)||!Sl[[b]]$repara) {
            if (Amat) { 
              B[ind,] <- if (nt==1)  Sl[[b]]$lambda[i]*(Sl[[b]]$S[[i]]%*%A[ind,,drop=FALSE]) else 
                         Sl[[b]]$lambda[i]*pmmult(Sl[[b]]$S[[i]],A[ind,,drop=FALSE],nt=nt) 
            } else B[ind] <-  Sl[[b]]$lambda[i]*(Sl[[b]]$S[[i]]%*%A[ind])
          } else {
            if (Amat) { 
              B[ind,] <- if (nt==1) Sl[[b]]$Srp[[i]]%*%A[ind,,drop=FALSE] else 
                       pmmult(Sl[[b]]$Srp[[i]],A[ind,,drop=FALSE],nt=nt) 
            } else B[ind] <- Sl[[b]]$Srp[[i]]%*%A[ind]
          }
          SA[[k]] <- B
        } else { ## strip zero rows from answer
          if (is.null(Sl[[b]]$Srp)||!Sl[[b]]$repara) {
            if (Amat) {
              SA[[k]] <- if (nt==1)  Sl[[b]]$lambda[i]*(Sl[[b]]$S[[i]]%*%A[ind,,drop=FALSE]) else
                         Sl[[b]]$lambda[i]*pmmult(Sl[[b]]$S[[i]],A[ind,,drop=FALSE],nt=nt)
            } else SA[[k]] <-  Sl[[b]]$lambda[i]*as.numeric(Sl[[b]]$S[[i]]%*%A[ind])
          } else {
            if (Amat) {
              SA[[k]] <- if (nt==1) Sl[[b]]$Srp[[i]]%*%A[ind,,drop=FALSE] else
                       pmmult(Sl[[b]]$Srp[[i]],A[ind,,drop=FALSE],nt=nt)
            } else SA[[k]] <- as.numeric(Sl[[b]]$Srp[[i]]%*%A[ind])
          }
          attr(SA[[k]],"ind") <- ind
        }
      } ## end of S loop for block b
    }
  } ## end block loop
  nl <- matrix(0,2,length(SA)) ## 0 indicates linear term
  if (length(nli)) { ## record which element of Sl each non-linear SA element relates to 
    for (i in 1:length(nli)) nl[,nli[[i]][1]+1:Sl[[nli[[i]][2]]]$n.sp-1] <- nli[[i]] 
  }
  attr(SA,"nli") <- nl ## record starts of non-linear terms in SA
  SA
} ## end Sl.termMult

d.detXXS <- function(Sl,PP,nt=1,deriv=2,SPP=FALSE) {
## function to obtain derivatives of log |X'X+S| given unpivoted PP' where 
## P is inverse of R from the QR of the augmented model matrix.
## Note that d1[k] is also the EDF suppressed by the kth smoothing penalty,
## for linear smoothing parameters.
  spp <- SPP
  SPP <- Sl.termMult(Sl,PP,full=FALSE,nt=nt) ## SPP[[k]] is S_k PP' where S_k is derivative in nl case
  ## if PP is sparse over block covered by penalty, implying that it is from the ISA and not a
  ## full inverse, then second order terms like tr(PPdS_jPPdS_k) will not be correct - hence will
  ## be omitted...
  if (inherits(PP,"Matrix")) { 
    sparse <- sapply(SPP,function(x,PP) {ii <- attr(x,"ind");length(ii)^2!=length(PP[ii,ii]@x)},PP=PP)
  } else sparse <- rep(FALSE,length(SPP))
  nd <- length(SPP)
  d1 <- rep(0,nd);d2 <- matrix(0,nd,nd)
  if (deriv==2) nli <- attr(SPP,"nli")
  if (nd>0) for (i in 1:nd) { 
    indi <- attr(SPP[[i]],"ind")
    d1[i] <- sum(diag(SPP[[i]][,indi,drop=FALSE]))
    if (deriv==2) {
      b <- nli[2,i] ## non-linear Sl[[b]] b==0 => linear
      for (j in i:nd) {
        indj <- attr(SPP[[j]],"ind")
	## following is correct if only one term is sparse (involves incomplete portion of inverse Hessian PP),
	## but not both. Use as approximation if i and j are for same block, but zero otherwise
        xx <- if (sparse[i]&&sparse[j]&&!all.equal(indj,indi)) 0 else -sum(t(SPP[[i]][,indj,drop=FALSE])*SPP[[j]][,indi,drop=FALSE])  
        if (b && b == nli[2,j]) { ## non-linear second derivative needed
	  k0 <- nli[1,i] - 1 ## where sps start in overall sp vector - 1
	  ind <- Sl[[b]]$start:Sl[[b]]$stop
	  xx <- xx + sum(diag(Sl[[b]]$AdS(t(PP[ind,,drop=FALSE]),Sl[[b]],j-k0,i-k0)[ind,]))
	}
	d2[i,j] <- d2[j,i] <- xx 
      }
      if (nli[2,i]==0) d2[i,i] <- d2[i,i] + d1[i] 
    }  
  }
  if (!spp) SPP <- NULL else attr(SPP,"sparse") <- sparse
  list(d1=d1,d2=d2,SPP=SPP)
} ## end d.detXXS

Sl.ift <- function(Sl,R,X,y,beta,piv,rp) {
## function to obtain derviatives of \hat \beta by implicit differentiation
## and to use these directly to evaluate derivs of b'Sb and the RSS.
## piv and rp are the pivots and inverse pivots from the qr that produced R.
## rssj and bSbj only contain the terms that will not cancel in rssj + bSbj
  beta <- beta[rp] ## unpivot
  Sb <- Sl.mult(Sl,beta,k = 0)          ## unpivoted
  Skb <- Sl.termMult(Sl,beta,full=TRUE) ## unpivoted
  rsd <- (X%*%beta - y)
  #Xrsd <- t(X)%*%rsd ## X'Xbeta - X'y
  nd <- length(Skb)
  np <- length(beta)
  db <- matrix(0,np,nd)
  rss1 <- bSb1 <- rep(0,nd)  
  
  for (i in 1:nd) { ## compute the first derivatives
    db[,i] <- -backsolve(R,forwardsolve(t(R),Skb[[i]][piv]))[rp] ## d beta/ d rho_i
    ## rss1[i] <- 0* 2 * sum(db[,i]*Xrsd)                      ## d rss / d rho_i
    bSb1[i] <- sum(beta*Skb[[i]]) ## + 2 * sum(db[,i]*Sb)   ## d b'Sb / d_rho_i
  }
  XX.db <- t(X)%*%(X%*%db)
  S.db <- Sl.mult(Sl,db,k=0)

  rss2 <- bSb2 <- matrix(0,nd,nd)
  for (k in 1:nd) { ## second derivative loop 
    for (j in k:nd) {
      ## d2b <- (k==j)*db[,k] - backsolve(R,forwardsolve(t(R),Sk.db[[j]][piv,k]+Sk.db[[k]][piv,j]))[rp]
      rss2[j,k] <- rss2[k,j] <- 2 * sum(db[,j]*XX.db[,k]) ## + 2 * sum(d2b*Xrsd) 
      bSb2[j,k] <- bSb2[k,j] <-  (k==j)*sum(beta*Skb[[k]])  + 2*(sum(db[,k]*(Skb[[j]]+S.db[,j])) + 
                                 sum(db[,j]*Skb[[k]])) ## + 2 * (sum(d2b*Sb)   
                               
    }
  }
  list(bSb=sum(beta*Sb),bSb1=bSb1,bSb2=bSb2,d1b=db,rss =sum(rsd^2),rss1=rss1,rss2=rss2)
} ## end Sl.ift

Sl.iftChol <- function(Sl,XX,R,d,beta,piv,nt=1) {
## function to obtain derviatives of \hat \beta by implicit differentiation
## and to use these directly to evaluate derivs of b'Sb and the RSS.
## piv contains the pivots from the chol that produced R.
## rssj and bSbj only contain the terms that will not cancel in rssj + bSbj
## Col j of S1b is dS/drho_j beta
  Sb <- Sl.mult(Sl,beta,k = 0)          ## unpivoted
  Skb <- Sl.termMult(Sl,beta,full=TRUE) ## unpivoted dS/drho beta
  nd <- length(Skb) ## number of derivatives
  cd <- FALSE
  if (cd) { ## check derivatives
    k <- 0
    eps <- 1e-6
    eeps <- exp(eps)
    fdSk <- list()
    for (b in 1:length(Sl)) {
      if (Sl[[b]]$linear) {
        ind <- 1:length(Sl[[b]]$S)
	for (i in ind) {
	  Sl[[b]]$lambda[i] <- Sl[[b]]$lambda[i]*eeps
	  Sb1 <- Sl.mult(Sl,beta,k = 0)
	  Sl[[b]]$lambda[i] <- Sl[[b]]$lambda[i]/eeps
          k <- k + 1
          fdSk[[k]] <- (Sb1-Sb)/eps
	}
      } else {
        ind <- 1:Sl[[b]]$n.sp
	theta <- Sl[[b]]$lambda
	for (i in ind) {
	  theta[i] <- theta[i] + eps
	  Sl[[b]] <- Sl[[b]]$updateS(theta,Sl[[b]])
	  theta[i] <- theta[i] - eps
	  Sb1 <- Sl.mult(Sl,beta,k = 0)
	  k <- k + 1
          fdSk[[k]] <- (Sb1-Sb)/eps
	}
	Sl[[b]] <- Sl[[b]]$updateS(theta,Sl[[b]])
      }
    }
    plot(Skb[[1]],fdSk[[1]])
  } ## if cd deriv check
  np <- length(beta)
  db <- matrix(0,np,nd)
  rss1 <- bSb1 <- rep(0,nd)
  if (nd==0) return(list(bSb=0,bSb1=rep(0,nd),bSb2=rep(0,nd),d1b=db,rss1=rss1,rss2=rss1))

  ## alternative all in one code - matches loop results, but
  ## timing close to identical - modified for parallel exec
  D <- matrix(unlist(Skb),length(beta),nd) ## cols of this are dS/drho_j beta
  bSb1 <- colSums(beta*D)
  if (inherits(R,c("dgCMatrix","dtCMatrix"))) { ## sparse R
    Dd <- Diagonal(length(d),1/d)
    db[piv,] <- -as(solve(R,solve(t(R),(Dd %*% D)[piv,])),"matrix")
    db <- as(Dd %*% db,"matrix")
  } else { ## dense R
    D1 <- .Call(C_mgcv_Rpforwardsolve,R,D[piv,]/d[piv],nt) ## note R transposed internally unlike forwardsolve
    db[piv,] <- -.Call(C_mgcv_Rpbacksolve,R,D1,nt)/d[piv]
  }
  S.db <- Sl.mult(Sl,db,k=0)
  nli <- attr(Skb,"nli") ## index of non-linear terms npi[1,] is first sp of term, npi[2,] is block index 0s for linear 
  bSb2 <- diag(x=colSums(beta*D)*(nli[1,]==0),nrow=nd)  ## linear only - only applies if npi[1,]==0
  bSb2 <- bSb2 + 2 * (.Call(C_mgcv_pmmult2,db,D+S.db,1,0,nt) + .Call(C_mgcv_pmmult2,D,db,1,0,nt))
  if (any(nli[1,]!=0)) { ## non-linear 2nd deriv part needed
    kk <- which(nli[2,]!=0)
    for (i in kk) {
      b <- nli[2,i];k0 <- nli[1,i]-1
      ind <- Sl[[b]]$start:Sl[[b]]$stop
      for (j in kk) if (nli[2,j]==b && i<=j) {
        x <- bSb2[i,j] + sum(Sl[[b]]$AdS(beta[ind],Sl[[b]],i-k0,j-k0)*beta[ind])
	bSb2[i,j] <- bSb2[j,i] <- x
      }
    }  
  }
  
  if (is.null(XX)) return(list(bSb2=bSb2,bSb1=bSb1,db=db,S1b=D)) ## return early
  
  ## XX.db <- XX%*%db
  XX.db <- .Call(C_mgcv_pmmult2,XX,db,0,0,nt)

  rss2 <- 2 * .Call(C_mgcv_pmmult2,db,XX.db,1,0,nt)
 
  list(bSb=sum(beta*Sb),bSb1=bSb1,bSb2=bSb2,
       d1b=db ,rss1=rss1,rss2=rss2,S1b=D)
} ## end Sl.iftChol

Sl.ncv <- function(y,Xd,k,ks,ts,dt,v,qc,nei,Sl,XX,w,f,rho,nt=c(1,1),L=NULL,rho0=0,drop=NULL,tol=0,nthreads=1) {
## given X'WX in XX and f=X'Wy solves the penalized least squares problem
## with penalty defined by Sl and rho, and evaluates an NCV Newton step, the NCV 
## gradient and the estimated coefs beta. rsd contains cross validated weighted residuals for first
## element in prediction neighbourhood on exit. 
  rho <- if (is.null(L)) rho + rho0 else L%*%rho + rho0
  if (length(rho)<length(rho0)) rho <- rho0
  sp <- exp(rho)
  XXp <- Sl.addS(Sl,XX,rho) ## penalized Hessian

  d <- diag(XXp);ind <- d<=0
  d[ind] <- 1;d[!ind] <- sqrt(d[!ind]) ## diagonal pre-conditioner
  ## get preconditioned Cholesky factor
  R <- if (nt[2]>1) pchol(t(XXp/d)/d,nt[2]) else suppressWarnings(chol(t(XXp/d)/d,pivot=TRUE))
  r <- min(attr(R,"rank"),Rrank(R))
  p <- ncol(XXp)
  piv <- attr(R,"pivot") #;rp[rp] <- 1:p
  if (r<p) { ## drop rank deficient terms...
    R <- R[1:r,1:r]
    piv <- piv[1:r]
  }
  beta <- rep(0,p)
  beta[piv] <- backsolve(R,(forwardsolve(t(R),f[piv]/d[piv])))/d[piv]
  betar <- Sl.initial.repara(Sl,beta,inverse=TRUE,both.sides=FALSE,cov=FALSE)
  mu <- Xbd(Xd,betar,k,ks,ts,dt,v,qc,drop)
  rsd <- y-mu
  G <- matrix(0,p,p)
  if (nt[2]>1) {
    P <- pbsi(R,nt=nt[2],copy=TRUE) ## invert R 
    G[piv,piv] <-  pRRt(P,nt[2]) ## PP'
  } else G[piv,piv] <- chol2inv(R)
  G0 <- t(G/d)/d ## inverse of penalized (scaled) Hessian

  
  G <- Sl.initial.repara(Sl,G0,inverse=TRUE,both.sides=TRUE,cov=TRUE,nt=nthreads) 
  

  NCV <- 0; nsp <- length(rho)
  NCV1 <- numeric(nsp); NCV2 <- numeric(nsp*nsp)
  m <- as.integer(unlist(lapply(Xd,nrow)));p <- as.integer(unlist(lapply(Xd,ncol)))
  beta1 <- numeric(length(beta)*nsp)
  ## profiling indicates that by far the dominant cost here is the diagXVXt calls in the C code...
  tt <- system.time(.Call(C_CNCV, NCV, NCV1, NCV2, G, rsd, betar, beta1, w, sp, as.double(unlist(Xd)), k-1L, as.integer(ks-1L), m, p,
        as.integer(ts-1L), as.integer(dt), as.double(unlist(v)), as.integer(qc), as.integer(nthreads), nei, Sl))
  NCV2 <- matrix(NCV2,nsp,nsp)
  beta1 <- matrix(beta1,length(beta),nsp)
  ## NOTE: should modify to pass back dbeta/drho from CNCV
  if (ncol(beta1)>0) for (i in 1:ncol(beta1)) beta1[,i] <- ## d beta / d rho matrix
        Sl.initial.repara(Sl,as.numeric(beta1[,i]),inverse=FALSE,both.sides=FALSE,cov=FALSE,nt=nthreads) 
 
  db <- FALSE
  if (db) { ## NOTE: above .Call now modifies rsd - need a forced copy and save to use here! 
    tt1 <- system.time(dA <- diagXVXd(Xd,G,k,ks,ts,dt,v,qc,drop=NULL,nthreads=1,lt=NULL,rt=NULL)*w)
    rcv <- rsd/(1-dA)
    NCV0 <- sum(rcv^2*w)
    mu1 <- Xbd(Xd,beta1,k,ks,ts,dt,v,qc)
    b1 <-  matrix(0,length(y),length(sp))
    NCV10 <- rep(0,length(sp))
    for (i in 1:length(sp)) {
      ii <- (Sl[[i]]$start:Sl[[i]]$stop)#[Sl[[i]]$ind]
      A1 <- diagXVXd(Xd,-sp[i]*G[,ii]%*%Sl[[i]]$S0[[1]]%*%G[ii,],k,ks,ts,dt,v,qc,drop=NULL,nthreads=1,lt=NULL,rt=NULL)*w
      b1[,i] <- A1/(1-dA)^2*rsd-mu1[,i]/(1-dA) ## C code for first term suspicious
      NCV10[i] <- 2*sum(w*b1[,i]*rcv)
    }
  }
  ## DEBUG end

  if (!is.null(L)) {
    NCV1 <- t(L) %*% NCV1
    NCV2 <- t(L) %*% NCV2 %*% L
  }
  uconv.ind <- (abs(NCV1) > tol)|(abs(diag(NCV2))>tol)
  robj <- list(beta=beta,grad=NCV1,db=beta1,PP=G0,R=R,piv=piv,rank=r,
               hess=NCV2,NCV=NCV,rsd=rsd[1:length(nei$md)])
  if (length(NCV1)>0 && sum(uconv.ind)>0) {
    if (sum(uconv.ind)!=ncol(NCV2)) { 
      NCV1 <- NCV1[uconv.ind]
      NCV2 <- NCV2[uconv.ind,uconv.ind]
    }

    er <- eigen(NCV2,symmetric=TRUE)
    er$values <- abs(er$values)
    me <- max(er$values)*.Machine$double.eps^.5
    er$values[er$values<me] <- me
    step <- rep(0,length(uconv.ind))
    step[uconv.ind] <- -er$vectors %*% ((t(er$vectors) %*% NCV1)/er$values)

    ## limit the step length...
    ms <- max(abs(step))
    if (ms>4) step <- 4*step/ms
  } else step <- 0
  robj$step <- step
  robj
} ## Sl.ncv

Sl.fitChol <- function(Sl,XX,f,rho,yy=0,L=NULL,rho0=0,log.phi=0,phi.fixed=TRUE,
                       nobs=0,Mp=0,nt=c(1,1),tol=0,gamma=1) {
## given X'WX in XX and f=X'Wy solves the penalized least squares problem
## with penalty defined by Sl and rho, and evaluates a REML Newton step, the REML 
## gradient and the estimated coefs beta. If phi.fixed=FALSE then we need 
## yy = y'Wy in order to get derivatives w.r.t. phi.
## NOTE: with an optimized BLAS nt==1 is likely to be much faster than
##       nt > 1
  tol <- as.numeric(tol)
  rho <- if (is.null(L)) rho + rho0 else L%*%rho + rho0
  if (length(rho)<length(rho0)) rho <- rho0 ## ncol(L)==0 or length(rho)==0
  ## get log|S|_+ without stability transform... 
  fixed <- rep(FALSE,length(rho))
  ldS <- ldetS(Sl,rho,fixed,np=ncol(XX),root=FALSE,repara=FALSE,nt=nt[1])
  
  ## now the Cholesky factor of the penalized Hessian... 
  XXp <- Sl.addS(Sl,XX,rho) ## penalized Hessian

  d <- diag(XXp);ind <- d<=0
  d[ind] <- 1;d[!ind] <- sqrt(d[!ind]) ## diagonal pre-conditioner
  ## get preconditioned Cholesky factor
  R <- if (nt[2]>1) pchol(t(XXp/d)/d,nt[2]) else suppressWarnings(chol(t(XXp/d)/d,pivot=TRUE))
  r <- min(attr(R,"rank"),Rrank(R))
  p <- ncol(XXp)
  piv <- attr(R,"pivot") #;rp[rp] <- 1:p
  if (r<p) { ## drop rank deficient terms...
    R <- R[1:r,1:r]
    piv <- piv[1:r]
  }
  beta <- rep(0,p)
  beta[piv] <- backsolve(R,(forwardsolve(t(R),f[piv]/d[piv])))/d[piv]

  ## get component derivatives based on IFT (noting that ldS$Sl has s.p.s updated to current)
  dift <- Sl.iftChol(ldS$Sl,XX,R,d,beta,piv,nt=nt[1])
 
  ## now the derivatives of log|X'X+S|
  PP <- matrix(0,p,p)
  if (nt[2]>1) {
    P <- pbsi(R,nt=nt[2],copy=TRUE) ## invert R 
    PP[piv,piv] <-  pRRt(P,nt[2]) ## PP'
  } else PP[piv,piv] <- chol2inv(R)
  PP <- t(PP/d)/d
  ldetXXS <- 2*sum(log(diag(R))+log(d[piv])) ## log|X'X+S|
  dXXS <- d.detXXS(ldS$Sl,PP,nt=nt[1]) ## derivs of log|X'X+S|

  phi <- exp(log.phi)  

  nrho <- length(rho)

  reml1 <- if (nrho==0) rep(0,0) else (dXXS$d1[!fixed] - ldS$ldet1 + 
            (dift$rss1[!fixed] + dift$bSb1[!fixed])/(phi*gamma))/2

  reml2 <- if (nrho==0) matrix(0,0,0) else (dXXS$d2[!fixed,!fixed] - ldS$ldet2 +  
           (dift$rss2[!fixed,!fixed] + dift$bSb2[!fixed,!fixed])/(phi*gamma))/2 

  if (nrho==0&&phi.fixed) { ## need to return now - nothing else to do
    return(list(beta=beta,grad=0,step=0,db=dift$d1b,PP=PP,R=R,piv=piv,rank=r, hess=reml2,ldetS=ldS$ldetS,ldetXXS=ldetXXS))
  }  

  if (!phi.fixed) {
    n <- length(reml1)
    rss.bSb <- yy - sum(beta*f) ## use identity ||y-Xb|| + b'Sb = y'y - b'X'y (b is minimizer)
    reml1[n+1] <- (-rss.bSb/(phi*gamma) + nobs/gamma - Mp)/2
    d <- c(-(dift$rss1[!fixed] + dift$bSb1[!fixed]),rss.bSb)/(2*phi*gamma)
    reml2 <- if (n>0) rbind(cbind(reml2,d[1:n]),d) else matrix(d,1,1) 
    if (!is.null(L)) L <- rbind(cbind(L,rep(0,nrow(L))),c(rep(0,ncol(L)),1))
  }

  if (!is.null(L)) {
    reml1 <- t(L)%*%reml1
    reml2 <- t(L)%*%reml2%*%L
  }
  uconv.ind <- (abs(reml1) > tol)|(abs(diag(reml2))>tol)
  hess <- reml2
  grad <- reml1
  if (length(grad)>0&&sum(uconv.ind)>0) {
    if (sum(uconv.ind)!=ncol(reml2)) { 
      reml1 <- reml1[uconv.ind]
      reml2 <- reml2[uconv.ind,uconv.ind]
    }

    er <- eigen(reml2,symmetric=TRUE)
    er$values <- abs(er$values)
    me <- max(er$values)*.Machine$double.eps^.5
    er$values[er$values<me] <- me
    step <- rep(0,length(uconv.ind))
    step[uconv.ind] <- -er$vectors%*%((t(er$vectors)%*%reml1)/er$values)

    ## limit the step length...
    ms <- max(abs(step))
    if (ms>4) step <- 4*step/ms
  } else step <- 0
  ## return the coefficient estimate, the reml grad and the Newton step...
  list(beta=beta,grad=grad,step=step,db=dift$d1b,PP=PP,R=R,piv=piv,rank=r,
       hess=hess,ldetS=ldS$ldetS,ldetXXS=ldetXXS)
} ## Sl.fitChol

Sl.fit <- function(Sl,X,y,rho,fixed,log.phi=0,phi.fixed=TRUE,rss.extra=0,nobs=NULL,Mp=0,nt=1,gamma=1) {
## fits penalized regression model with model matrix X and 
## initialised block diagonal penalty Sl to data in y, given 
## log smoothing parameters rho. 
## Returns coefs, reml score + grad and Hessian.
  np <- ncol(X) ## number of parameters
  n <- nrow(X) 
  phi <- exp(log.phi)
  if (is.null(nobs)) nobs <- n
  ## get log|S|_+ stably...
  ldS <- ldetS(Sl,rho,fixed,np,root=TRUE,nt=nt)
  ## apply resulting stable re-parameterization to X...
  X <- Sl.repara(ldS$rp,X)
  ## get pivoted QR decomp of augmented model matrix (in parallel if nt>1)
  qrx <- if (nt>1) pqr2(rbind(X,ldS$E),nt=nt) else qr(rbind(X,ldS$E),LAPACK=TRUE)
  rp <- qrx$pivot;rp[rp] <- 1:np ## reverse pivot vector
  ## find pivoted \hat beta...
  R <- qr.R(qrx)
  Qty0 <- qr.qty(qrx,c(y,rep(0,nrow(ldS$E))))
  beta <- backsolve(R,Qty0)[1:np]
  rss.bSb <- sum(Qty0[-(1:np)]^2) + rss.extra
  ## get component derivatives based on IFT...
  dift <- Sl.ift(ldS$Sl,R,X,y,beta,qrx$pivot,rp)
  ## and the derivatives of log|X'X+S|...
  P <- pbsi(R,nt=nt,copy=TRUE) ## invert R 
  ## P <- backsolve(R,diag(np))[rp,] ## invert R and row unpivot
  ## crossprod and unpivot (don't unpivot if unpivoting P above)
  PP <- if (nt==1) tcrossprod(P)[rp,rp] else pRRt(P,nt)[rp,rp] ## PP'
  ldetXXS <- 2*sum(log(abs(diag(R)))) ## log|X'X+S|
  dXXS <- d.detXXS(ldS$Sl,PP,nt=nt) ## derivs of log|X'X+S|
  ## all ingredients are now in place to form REML score and 
  ## its derivatives....
  reml <- (rss.bSb/(phi*gamma) + (nobs/gamma-Mp)*log(2*pi*phi) + Mp*log(gamma) +
           ldetXXS - ldS$ldetS)/2
  reml1 <-  (dXXS$d1[!fixed] - ldS$ldet1 + # dift$bSb1[!fixed]/phi)/2 
            (dift$rss1[!fixed] + dift$bSb1[!fixed])/(phi*gamma))/2

  reml2 <- (dXXS$d2[!fixed,!fixed] - ldS$ldet2 + #dift$bSb2[!fixed,!fixed]/phi)/2 
           (dift$rss2[!fixed,!fixed] + dift$bSb2[!fixed,!fixed])/(phi*gamma))/2 
  ## finally add in derivatives w.r.t. log.phi
  if (!phi.fixed) {
    n <- length(reml1)
    reml1[n+1] <- (-rss.bSb/(phi*gamma) + nobs/gamma - Mp)/2
    #d <- c(-(dift$bSb1[!fixed]),rss.bSb)/(2*phi)
    d <- c(-(dift$rss1[!fixed] + dift$bSb1[!fixed]),rss.bSb)/(2*phi*gamma)
    reml2 <- rbind(cbind(reml2,d[1:n]),d)
  } 
  ## following are de-bugging lines for testing derivatives of components...
  #list(reml=ldetXXS,reml1=dXXS$d1,reml2=dXXS$d2)
  #list(reml=ldS$ldetS,reml1=ldS$ldet1,reml2=ldS$ldet2)
  #list(reml=dift$rss,reml1=dift$rss1,reml2=dift$rss2)
  #list(reml=dift$bSb,reml1=dift$bSb1,reml2=dift$bSb2) 
  list(reml=as.numeric(reml),reml1=reml1,reml2=reml2,beta=beta[rp],PP=PP,
       rp=ldS$rp,rss=dift$rss+rss.extra,nobs=nobs,d1b=dift$d1b)
} ## Sl.fit

fast.REML.fit <- function(Sl,X,y,rho,L=NULL,rho.0=NULL,log.phi=0,phi.fixed=TRUE,
                 rss.extra=0,nobs=NULL,Mp=0,conv.tol=.Machine$double.eps^.5,nt=1,gamma=gamma) {
## estimates log smoothing parameters rho, by optimizing fast REML 
## using Newton's method. On input Sl is a block diagonal penalty 
## structure produced by Sl.setup, while X is a model matrix 
## reparameterized to correspond to any re-parameterization 
## used in Sl. Both will have had been modified to drop any 
## structurally un-identifiable coefficients. 
## Note that lower bounds on smoothing parameters are not handled.
  maxNstep <- 5  
  
  if (is.null(nobs)) nobs <- nrow(X)
  np <- ncol(X)
  if (nrow(X) > np) { ## might as well do an initial QR step
    qrx <- if (nt>1) pqr2(X,nt=nt) else qr(X,LAPACK=TRUE)
    rp <- qrx$pivot
    rp[rp] <- 1:np
    X <- qr.R(qrx)[,rp]
    y <- qr.qty(qrx,y)
    rss.extra <- rss.extra + sum(y[-(1:np)]^2)
    y <- y[1:np]
  }

  if (is.null(L)) {
    L <- diag(length(rho))
    if (is.null(rho.0)) rho.0 <- rep(0,nrow(L))
  } else { ## convert intial s.p.s to account for L
    if (is.null(rho.0)) rho.0 <- rep(0,nrow(L))
    rho <- as.numeric(coef(lm(rho ~ L-1+offset(rho.0))))
  }

  fixed <- rep(FALSE,nrow(L))
 
  
  best <- Sl.fit(Sl,X,y,L%*%rho+rho.0,fixed,log.phi,phi.fixed,rss.extra,nobs,Mp,nt=nt,gamma=gamma)
  ## get a typical scale for the reml score... 
  reml.scale <- abs(best$reml) + best$rss/best$nobs
 
  nr <- length(rho.0)
  if (!phi.fixed) { 
    rho <- c(rho,log.phi) ## append log.phi for fitting
    rho.0 <- c(rho.0,0)
    L <- rbind(cbind(L,L[,1]*0),c(L[1,]*0,1))
  }
  grad <- as.numeric(t(L)%*% best$reml1)
  hess <- t(L)%*% best$reml2%*%L
  grad2 <- diag(hess)  

  ## create and index for the unconverged... 
  ## idea in following is only to exclude terms with zero first and second derivative
  ## from optimization, as it is only these that slow things down if included...  
  uconv.ind <- (abs(grad) > reml.scale*conv.tol*.1)|(abs(grad2)>reml.scale*conv.tol*.1)
  ## if all appear converged at this stage, then there is probably something wrong,
  ## but reset anyway to see if situation can be recovered. If we don't do this then
  ## need to abort immediately, otherwise fails trying to eigen a 0 by 0 matrix
  if (sum(uconv.ind)==0) {
    warning("Possible divergence detected in fast.REML.fit",call.=FALSE,immediate.=TRUE)
    uconv.ind <- rep(TRUE,length(grad)) 
  }
  step.failed <- FALSE
  for (iter in 1:200) { ## the Newton loop
    ## Work only with unconverged (much quicker under indefiniteness)
    hess <- (t(L)%*% best$reml2%*%L)[uconv.ind,uconv.ind]
    grad <- as.numeric(t(L)%*%best$reml1)[uconv.ind]
    ## check that Hessian is +ve def. Fix if not. 
    eh <- eigen(hess,symmetric=TRUE)
    ## flip negative eigenvalues to get +ve def...
    ind <- eh$values < 0
    eh$values[ind] <- -eh$values[ind]
    ## avoid indefiniteness by further manipulation...
    thresh <- max(abs(eh$values))*.Machine$double.eps^.5
    ind <- eh$values < thresh
    eh$values[ind] <- thresh 
    ## get the Newton direction, -ve inverse hessian times gradient
    uc.step <- - eh$vectors%*%((t(eh$vectors)%*%grad)/eh$values)
    
    ## now make sure step is not too far...
    ms <- max(abs(uc.step))
    if (ms>maxNstep) uc.step <- maxNstep * uc.step/ms

    step <- rep(0,length(uconv.ind)); ## full step (including converged)
    step[uconv.ind] <- uc.step ## step includes converged
    ## try out the step...
    rho1 <- L%*%(rho + step)+rho.0; if (!phi.fixed) log.phi <- rho1[nr+1]
    trial <- Sl.fit(Sl,X,y,rho1[1:nr],fixed,log.phi,phi.fixed,rss.extra,nobs,Mp,nt=nt,gamma=gamma)
    k <- 0
    not.moved <- 0 ## count number of consecutive steps of essentially no change from best
    while (trial$reml>best$reml) { ## step half until improvement or failure
      ## idea with the following is to count the number of consecutive step halvings
      ## without a numerically significant change from best$reml, since
      ## this is an early indicator of step failure.
      if (trial$reml-best$reml < conv.tol*reml.scale) not.moved <- not.moved + 1 else not.moved <- 0
      if (k==25||sum(rho != rho + step)==0||not.moved>3) {
        step.failed <-  TRUE
	break
      }	
      step <- step/2;k <- k + 1
      rho1 <- L%*%(rho + step)+rho.0; if (!phi.fixed) log.phi <- rho1[nr+1]
      trial <- Sl.fit(Sl,X,y,rho1[1:nr],fixed,log.phi,phi.fixed,rss.extra,nobs,Mp,nt=nt,gamma=gamma)
    }
    if (step.failed) break ## can get no further
    #if ((k==35 && trial$reml>best$reml)||(sum(rho != rho + step)==0)) { ## step has failed
    #  step.failed <- TRUE
    #  break ## can get no further
    #}
    ## At this stage the step has been successful. 
    ## Need to test for convergence...
    converged <- TRUE
    grad <- as.numeric(t(L)%*%trial$reml1)
    hess <- t(L)%*%trial$reml2%*%L;grad2 <- diag(hess)
    ## idea in following is only to exclude terms with zero first and second derivative
    ## from optimization, as it is only these that slow things down if included...    
    uconv.ind <- (abs(grad) > reml.scale*conv.tol*.1)|(abs(grad2)>reml.scale*conv.tol*.1)
    ## now do the convergence testing...
    ## First check gradiants...
    if (sum(abs(grad)>reml.scale*conv.tol)) converged <- FALSE
    ## Now check change in REML values
    if (abs(best$reml-trial$reml)>reml.scale*conv.tol) { 
      if (converged) uconv.ind <- uconv.ind | TRUE ## otherwise can't progress
      converged <- FALSE      
    }
    best <- trial ## trial becomes current best.
    rho <- rho + step ## and new log sp accepted.
    if (converged) break ## ok done, leave loop
    reml.scale <- abs(best$reml) + best$rss/best$nobs ## update for next iterate
  } ## end of Newton loop
  if (iter==200) warning("fast REML optimizer reached iteration limit")
  if (step.failed) best$conv <- "step failed" else
  if (iter==200) best$conv <- "no convergence in 200 iterations" else
  best$conv <- "full convergence"
  best$iter <- iter
  best$outer.info <- list(conv = best$conv, iter = best$iter,grad = grad,hess = hess)
  best$rho <- rho
  best$rho.full <-  as.numeric(L%*%rho+rho.0)
  best ## return the best fit (note that it will need post-processing to be useable)
} ## end fast.REML.fit

ident.test <- function(X,E,nt=1) {
## routine to identify structurally un-identifiable coefficients
## for model with model matrix X and scaled sqrt penalty matrix E
## lambda is smoothing parameter vector corresponding to E, 
## and the routine also suggests starting values for lambda
## based on scaling of X and E. 
## If length(drop)>0 then X[,-drop] is new model matrix
## if beta contains coefs with unidentified dropped, and 
## if beta.full is a vector of zeroes for each original coef
## then beta.full[undrop] <- beta, is the full, zero padded 
## coeff vector, with dropped coefs re-nstated as zeroes. 
  Xnorm <- norm(X,type="F")
  qrx <- if (nt>1) pqr2(rbind(X/Xnorm,E),nt=nt) else qr(rbind(X/Xnorm,E),LAPACK=TRUE) ## pivoted QR
  rank <- Rrank(qr.R(qrx),tol=.Machine$double.eps^.75)
  drop <- qrx$pivot[-(1:rank)] ## index of un-identifiable coefs
  undrop <- 1:ncol(X) 
  if (length(drop)>0) undrop <- undrop[-drop]
  list(drop=drop,undrop=undrop)
} ## ident.test



Sl.drop <- function(Sl,drop,np) {
## routine to drop coefs in drop, from block diagonal penalty
## stored in Sl. np is total number of coeffs
  if (length(drop)==0) return(Sl)
  keep <- rep(TRUE,np)
  keep[drop] <- FALSE  ## logical indexing of retained coefs
  ## let b1 be coefs after dropping and b0 be full vector before.
  ## new.loc gives location in b1 of ith element in b0. If i is 
  ## in drop then new.loc[i] is position of last b0[j] not dropped.
  ## i.e. for i not in drop, b0[i] = b1[new.loc[i]].
  ##      for i in drop, b1[new.loc[i]] = b0[j] where j is largest
  ##      j < i s.t. j not in drop. 
  ## These indices facilitate easy dropping from parts of blocks 
  ## corresponding to coef indices in drop.
  cholesky <- attr(Sl,"cholesky") ## is setup all Cholesky based?
  new.loc <- cumsum(keep) 
  dropped.blocks <- FALSE
  for (b in 1:length(Sl)) {
    ind <- (Sl[[b]]$start:Sl[[b]]$stop)##[Sl[[b]]$ind]
    if (length(Sl[[b]]$S)==1) { ## singleton
      ## need to count terms dropped from penalty,
      ## adjusting rank, ind, start and stop
      bdrop <- ind%in%drop ## which are dropped here?
      npd <- sum(bdrop[Sl[[b]]$ind]) ## number of penalized dropped
      Sl[[b]]$ind <- Sl[[b]]$ind[!bdrop] ## retain not dropped
      Sl[[b]]$rank <- Sl[[b]]$rank - npd ## reduce rank by penalized dropped
      if (Sl[[b]]$rank <=0) dropped.blocks <- TRUE 
      Sl[[b]]$start <- new.loc[Sl[[b]]$start]
      Sl[[b]]$stop <- new.loc[Sl[[b]]$stop]
    } else { ## multi-S 
      bdrop <- ind%in%drop ## which are dropped here?
      keep <- !bdrop[Sl[[b]]$ind] ## index of what to keep in the penalties
      npd <- sum(!keep) ## number of penalized dropped
      Sl[[b]]$ind <- Sl[[b]]$ind[!bdrop] ## retain not dropped
      Sl[[b]]$rank <- Sl[[b]]$rank - npd ## reduce rank by penalized dropped
      if (Sl[[b]]$rank <=0) dropped.blocks <- TRUE 
      ## need to drop rows and cols from S and and rows from rS        
      for (i in 1:length(Sl[[b]]$S)) {
        if (length(Sl[[b]]$rS)) Sl[[b]]$rS[[i]] <- if (cholesky) Sl[[b]]$rS[[i]][keep,keep] else Sl[[b]]$rS[[i]][keep,]
        Sl[[b]]$S[[i]] <- Sl[[b]]$S[[i]][keep,keep] 
      }
      Sl[[b]]$start <- new.loc[Sl[[b]]$start]
      Sl[[b]]$stop <- new.loc[Sl[[b]]$stop]
    }
  }
  if (dropped.blocks) {
    new.b <- 1
    for (i in 1:length(Sl)) {
      if (Sl[[b]]$rank>0) {
        Sl[[new.b]] <- Sl[[b]]
        new.b <- new.b + 1
      }
    }
  }
  attr(Sl,"drop") <- drop
  Sl
} ## Sl.drop

Sl.Xprep <- function(Sl,X,nt=1) { 
## Sl is block diag object from Sl.setup, X is a model matrix
## this routine applies preliminary Sl transformations to X
## tests for structural identifibility problems and drops
## un-identifiable parameters.
  X <- Sl.initial.repara(Sl,X,inverse=FALSE,both.sides=FALSE,cov=FALSE,nt=nt) ## apply re-para used in Sl to X
  id <- ident.test(X,attr(Sl,"E"),nt=nt) ## deal with structural identifiability
  ## id contains drop, undrop, lambda
  if (length(id$drop)>0) { ## then there is something to do here 
    Sl <- Sl.drop(Sl,id$drop,ncol(X)) ## drop unidentifiable from Sl
    X <- X[,-id$drop] ## drop unidentifiable from X 
  }
  rank <- 0
  for (b in 1:length(Sl)) rank <- rank + Sl[[b]]$rank ## the total penalty rank
  ## Also add Mp, the total null space dimension to return list.  
  list(X=X,Sl=Sl,undrop=id$undrop,rank=rank,Mp=ncol(X)-rank) 
} ## end Sl.Xprep


Sl.postproc <- function(Sl,fit,undrop,X0,cov=FALSE,scale = -1,L,nt=1) {
## reverse the various fitting re-parameterizations.
## X0 is the orginal model matrix before any re-parameterization
## or parameter dropping. Sl is also the original *before parameter 
## dropping*
  np <- ncol(X0)
  beta <- rep(0,np)
  beta[undrop] <- Sl.repara(fit$rp,fit$beta,inverse=TRUE)
  beta <- Sl.initial.repara(Sl,beta,inverse=TRUE,both.sides=TRUE,cov=TRUE,nt=nt)
 
  if (cov) { 
    d1b <- matrix(0,np,ncol(fit$d1b))
    ## following construction a bit ugly due to Sl.repara assumptions...
    d1b[undrop,] <- t(Sl.repara(fit$rp,t(fit$d1b),inverse=TRUE,both.sides=FALSE))
    for (i in 1:ncol(d1b)) d1b[,i] <- 
        Sl.initial.repara(Sl,as.numeric(d1b[,i]),inverse=TRUE,both.sides=TRUE,cov=TRUE,nt=nt) ## d beta / d rho matrix
    PP <- matrix(0,np,np)
    PP[undrop,undrop] <-  Sl.repara(fit$rp,fit$PP,inverse=TRUE)
    PP <- Sl.initial.repara(Sl,PP,inverse=TRUE,both.sides=TRUE,cov=TRUE,nt=nt)
    #XPP <- crossprod(t(X0),PP)*X0
    #hat <- rowSums(XPP);edf <- colSums(XPP)
    XPP <- crossprod(t(X0),PP)
    hat <- rowSums(XPP*X0)
    F <- crossprod(XPP,X0)
    edf <- diag(F)
    edf1 <- 2*edf - rowSums(t(F)*F) 
    ## edf <- rowSums(PP*crossprod(X0)) ## diag(PP%*%(t(X0)%*%X0))
    if (scale<=0) { 
      scale <- fit$rss/(fit$nobs - sum(edf))
    }
    Vp <- PP * scale ## cov matrix
    ## sp uncertainty correction... 
    if (!is.null(L)) d1b <- d1b%*%L
    M <- ncol(d1b) 
    ev <- eigen(fit$outer.info$hess,symmetric=TRUE)
    ind <- ev$values <= 0
    ev$values[ind] <- 0;ev$values[!ind] <- 1/sqrt(ev$values[!ind])
    rV <- (ev$values*t(ev$vectors))[,1:M]
    V.sp <- crossprod(rV)
    attr(V.sp,"L") <- L
    Vc <- crossprod(rV%*%t(d1b))
    Vc <- Vp + Vc  ## Bayesian cov matrix with sp uncertainty
    edf2 <- rowSums(Vc*crossprod(X0))/scale

    ##bias <- as.numeric(beta-F%*%beta) ## estimate of smoothing bias in beta
    return(list(beta=beta,Vp=Vp,Vc=Vc,Ve=F%*%Vp,V.sp=V.sp,edf=edf,edf1=edf1,edf2=edf2,hat=hat,F=F))
  } else return(list(beta=beta))
} ## Sl.postproc


## USEAGE SEQUENCE:
## 1. Use gam.setup to setup gam object, G, say, as usual
## 2. Call Sl.setup which uses info in G$smooth and G$paraPen
##    to set up a block diagonal penalty structure, Sl, say
## 3. At this stage reweight the model matrix in G if needed 
##    (e.g. in IRLS) to get X
## 4. um <- Sl.Xprep(Sl,X) to deal with identifiability and re-para.
## 5. initial smoothing parameters from initial.sp(X,G$S,G$off),
##    initial phi from, say variance of y over 10??
## 6. fit <- fast.REML.fit(um$Sl,um$X,G$y,rho,L=G$L,rho.0=G$lsp0,
##                         log.phi=log.phi,phi.fixed=FALSE/TRUE,Mp=um$Mp)
## 7. res <- Sl.postproc(Sl,fit,um$undrop,X,cov=TRUE), to get postproc 
##    stuff 
## Notice: that only steps 3-7 are needed in an IRLS loop and  cov=TRUE
## is only needed after convergence of such a loop. 
## Warning: min.sp not handled by this approach.

