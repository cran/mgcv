## (c) Simon N. Wood (2013-2022) distributed under GPL2
## Code for the gamlss families.
## idea is that there are standard functions converting
## derivatives w.r.t. mu to derivatives w.r.t. eta, given 
## given the links and derivatives. 
## Then there are standard routines to take the family 
## specific derivatives and the model matrices, and convert 
## these to the required gradient, hessian, etc...

## Storage convections:
## 1. A single model matrix is stored, along with a single param vector.
##    an index object associates columns with particular gamlss predictors.
## 2. Distribution specific derivatives are stored in d1l-d4l.  

## Need to somehow record block starts... 
## idea is that if n blocks are stored using loops with the 
## given l >= k >= j >= i structure then the block for index
## i,j,k,l starts at i4[i,j,k,l]*n+1, given symmetry over the indices. 

trind.generator <- function(K=2,ifunc=FALSE,reverse=!ifunc) {
## Generates index arrays for 'upper triangular' storage up to order 4
## Suppose you fill an array using code like...
## m = 1
## for (i in 1:K) for (j in i:K) for (k in j:K) for (l in k:K) {
##   a[,m] <- something; m <- m+1 }
## ... and do this because actually the same 'something' would 
## be stored for any permutation of the indices i,j,k,l.
## Clearly in storage we have the restriction l>=k>=j>=i,
## but for access we want no restriction on the indices.
## i4[i,j,k,l] produces the appropriate m for unrestricted 
## indices. i3 and i2 do the same for 3d and 2d arrays.
## if ifunc==TRUE then rather than index arrays, index functions
## are returned, so e.g.i4(i,j,k,l) is equivalent to above.
## Index functions require less storage for high K.
## ixr will extract the unique elements from an x dimensional
## upper triangular array in the correct order.
## Note that there are K*(K+1)/2 unique entries in a 2d array
## (K*(K+3)+2)*K/6 in a 2d array and (6+K*11+K^2*6+K^3)*K/24
## in a 4d array. 
  m.start <- 1
  if (ifunc) { ## return index functions
   eval(parse(text= paste("i2 <- function(i,j) {\n",
        "  if (i>j) {ii <- i;i <- j;j <- ii}\n",
        "  (i-1)*(",2*K+2,"-i)/2 +j-i+1\n}",sep="")))

   eval(parse(text=paste("i3 <- function(i,j,k) {\n",
        "  if (i>j||j>k) { \n    ii <- sort(c(i,j,k))\n",
	"    i <- ii[1];j <- ii[2];k <- ii[3]\n  }\n",
        "  (i-1)*(",3*K*(K+1),"+(i-2)*(i-",3*(K+1),"))/6 + (j-i)*(",2*K+3,"-i-j)/2+k-j+1 \n}",
	sep="")))

   eval(parse(text=paste("i4 <- function(i,j,k,l) {\n",
        "  if (i>j||j>k||k>l) { \n    ii <- sort(c(i,j,k,l))\n",
	"    i <- ii[1];j <- ii[2];k <- ii[3];l <- ii[4]\n  }\n",
        "  i1 <- i-1;i2 <- i-2; i1i2 <- i1*i2/2\n",
	"  l-k+1 + (k-j)*(",2*K+3,"-j-k)/2 +",
	"  (j-i)*(3*(",K+1,"-i)^2+3*(",K+1,"-i) + (j-i-1)*(j+2*i-",3*K+5,"))/6 +\n",
	"  (i1*",K^3+3*K^2+2*K,"+i1i2*(",K+1,"*(2*i-3) - ",3*K^2+6*K+2,"-i1i2))/6\n}",
	sep="")))
  } else { ## return index arrays
    i4 <- array(0,dim=c(K,K,K,K))  
    m <- m.start
    for (i in 1:K) for (j in i:K) for (k in j:K) for (l in k:K) {
      i4[i,j,k,l] <- i4[i,j,l,k] <- i4[i,k,l,j] <- i4[i,k,j,l] <- i4[i,l,j,k] <- 
      i4[i,l,k,j] <- i4[j,i,k,l] <- i4[j,i,l,k] <- i4[j,k,l,i] <- i4[j,k,i,l] <-
      i4[j,l,i,k] <- i4[j,l,k,i] <- i4[k,j,i,l] <- i4[k,j,l,i] <- i4[k,i,l,j] <-
      i4[k,i,j,l] <- i4[k,l,j,i] <- i4[k,l,i,j] <- i4[l,j,k,i] <- i4[l,j,i,k] <-
      i4[l,k,i,j] <- i4[l,k,j,i] <- i4[l,i,j,k] <- i4[l,i,k,j] <- m
      m <- m + 1
    }

    i3 <- array(0,dim=c(K,K,K))
    m <- m.start
    for (j in 1:K) for (k in j:K) for (l in k:K) {
      i3[j,k,l] <- i3[j,l,k] <- i3[k,l,j] <- i3[k,j,l] <- i3[l,j,k] <- 
      i3[l,k,j] <- m
      m <- m + 1
    }

    i2 <- array(0,dim=c(K,K))
    m <- m.start
    for (k in 1:K) for (l in k:K) {
      i2[k,l] <- i2[l,k] <- m
      m <- m + 1
    }
  }  
  ## now create the reverse indices...
  if (reverse) {
    m <- m.start
    maxi4 <- if (ifunc) i4(K,K,K,K) else i4[K,K,K,K]
    i4r <- rep(0,maxi4) ## extracts the unique elements from a symmetric array in packing order.
    for (i in 1:K) for (j in i:K) for (k in j:K) for (l in k:K) {
      i4r[m] <- l + (k-1)*K + (j-1)*K^2 + (i-1)*K^3
      m <- m + 1
    }
    m <- m.start
    maxi3 <- if (ifunc) i3(K,K,K) else i3[K,K,K]
    i3r <- rep(0,maxi3) ## extracts the unique elements from a symmetric array in packing order.
    for (j in 1:K) for (k in j:K) for (l in k:K) {
      i3r[m] <- l + (k-1)*K + (j-1)*K^2
      m <- m + 1
    }
    m <- m.start
    maxi2 <- if (ifunc) i2(K,K) else i2[K,K]
    i2r <- rep(0,maxi2) ## extracts the unique elements from a symmetric array in packing order.
    for (k in 1:K) for (l in k:K) {
      i2r[m] <- l + (k-1)*K
      m <- m + 1
    }
  } else i2r <- i3r <- i4r <- NULL  
  list(i2=i2,i3=i3,i4=i4,i2r=i2r,i3r=i3r,i4r=i4r)
} ## trind.generator

gamlss.ncv <- function(X,y,wt,nei,beta,family,llf,H=NULL,Hi=NULL,R=NULL,offset=NULL,dH=NULL,db=NULL,deriv=FALSE,nt=1) {
## computes the neighbourhood cross validation score and its derivative for a
## gamlss model. llf is what was returned by family$ll when ncv info requested.
## If derivs not required then ll must be called with deriv >=1, otherwise deriv >= 3.
## To enable NCV for a gamlss family:
## 1. the 'll' function must be modified to have an 'ncv' argument. When this is TRUE and
##    deriv!=0 then ll should return l1, l2 and l3 the derivatives of the log likelihood
##    w.r.t. the linear preditors (typically returned from gamlss.mueta).
## 2. the 'll' function must have an eta argument allowing the linear predictors to be
##    supplied directly, rather than being computed from X and beta. 
## 3. The family must contain an 'ncv' wrapper function, which simply calls this function.
## ... gaulss provides an example.
  jj <- attr(X,"lpi") ## extract linear predictor index
  nlp <- length(jj); n <- nrow(X)
  if (deriv>0) {
    nsp <- ncol(db)
    deta <- matrix(0,n*nlp,nsp)
    ind <- 1:n
    for (i in 1:nlp) {
      deta[ind,] <- X[,jj[[i]],drop=FALSE] %*% db[jj[[i]],,drop=FALSE]
      ind <- ind + n
    }  
  } else deta <- 0.0
  ## debug section
  eta <- matrix(0,n,nlp)
  for (i in 1:nlp) eta[,i] <- X[,jj[[i]],drop=FALSE] %*% beta[jj[[i]]]
  ## end debug
  nm <- length(nei$d)
  eta.cv <- matrix(0,nm,nlp)
  
  deta.cv <- if (deriv>0) matrix(0,nm*nlp,nsp) else if (deriv<0) matrix(0,length(nei$ma),length(beta))  else 0.0
  if (is.null(R)) {
    cg.iter <- .Call(C_ncvls,X,jj,H,Hi,dH,llf$l1,llf$l2,llf$l3,nei$d-1,nei$md,nei$ma,nei$a-1,beta,eta.cv,deta.cv,
                   deta,db,deriv)		   
  } else {
    cg.iter <- .Call(C_Rncvls,X,jj,R,dH,llf$l1,llf$l2,llf$l3,nei$d-1,nei$md,nei$ma,nei$a-1,beta,eta.cv,deta.cv,
                   deta,db,deriv,.Machine$double.eps,nt)
  }
  if (!is.null(offset)) {
    for (i in 1:ncol(eta.cv)) if (i <= length(offset)&&!is.null(offset[[i]])) eta.cv[,i] <- eta.cv[,i] + offset[[i]][nei$d]
  }
  ## ll must be set up to return l1..l3 as derivs w.r.t. linear predictors if ncv=TRUE
  ncv1 <- NULL
  gamma <- llf$gamma;qapprox <- family$qapprox
  dev <- if (gamma!=1||qapprox) -sum(family$ll(y,X,beta,wt,family,offset,deriv=0,db,eta=eta,ncv=TRUE)$l0[nei$d]) else 0 
  if (qapprox) { ## quadratic approximate version
    ncv <-  dev - gamma*sum(llf$l1[nei$d,]*(eta.cv-eta[nei$d,]))
    k <- 0
    for (i in 1:nlp) for (j in i:nlp) {
      k <- k  + 1
      ncv <- ncv - 0.5*gamma*(1+(i!=j))*sum(llf$l2[nei$d,k]*(eta.cv[,i]-eta[nei$d,i])*(eta.cv[,j]-eta[nei$d,j])) ## symmetric term
    }
    if (deriv>0) {
      #ncv1 <- -colSums(as.numeric(llf$l1[nei$d,])*(deta.cv*gamma+(1-gamma)*deta))
      rowk <- 1:nm
      ncv1 <- -as.numeric(llf$l1[nei$d,1])*(deta.cv[rowk,]*gamma+(1-gamma)*deta[rowk,])
      if (nlp>1) for (j in 2:nlp) {
        rowk <- rowk + nm 
        ncv1 <- ncv1 - as.numeric(llf$l1[nei$d,j])*(deta.cv[rowk,]*gamma+(1-gamma)*deta[rowk,])
      }
      kk <- 0;jj <- 0
      for (j in 1:nlp) for (k in j:nlp) {
        kk <- kk  + 1
	rowk <- 1:nm+(k-1)*nm
	rowj <- 1:nm+(j-1)*nm
	#ncv1 <- ncv1 - colSums(llf$l2[nei$d,kk]*((deta[1:nm+(k-1)*nm,] + deta.cv[1:nm+(k-1)*nm,])*(eta.cv[,j]-eta[nei$d,j]) + 
	#                   (eta.cv[,k]-eta[nei$d,k])*(deta.cv[1:nm+(j-1)*nm,] - deta[nei$d+(j-1)*n,])))*gamma*.5
        ncv1 <- ncv1 - llf$l2[nei$d,kk]*((deta[rowk,] + deta.cv[rowk,])*(eta.cv[,j]-eta[nei$d,j]) +
	                       (eta.cv[,k]-eta[nei$d,k])*(deta.cv[rowj,] - deta[nei$d+(j-1)*n,]))*gamma*.5
        #if (j!=k) ncv1 <- ncv1 - colSums(llf$l2[nei$d,kk]*((deta[1:nm+(j-1)*nm,] + deta.cv[1:nm+(j-1)*nm,])*(eta.cv[,k]-eta[nei$d,k]) + 
	#                   (eta.cv[,j]-eta[nei$d,j])*(deta.cv[1:nm+(k-1)*nm,] - deta[nei$d+(k-1)*n,])))*gamma*.5		  
        if (j!=k) ncv1 <- ncv1 - llf$l2[nei$d,kk]*((deta[rowj,] + deta.cv[rowj,])*(eta.cv[,k]-eta[nei$d,k]) +
	          (eta.cv[,j]-eta[nei$d,j])*(deta.cv[rowk,] - deta[nei$d+(k-1)*n,]))*gamma*.5
        for (l in k:nlp) {
          jj <- jj + 1
	  #ncv1 <- ncv1 - (1+(j!=k)) * gamma*.5 * colSums(
   	  #        llf$l3[nei$d,jj]*deta[nei$d+(l-1)*n,]*(eta.cv[,k]-eta[nei$d,k])*(eta.cv[,j]-eta[nei$d,j]))
	  ncv1 <- ncv1 - (1+(j!=k)) * gamma*.5 * llf$l3[nei$d,jj]*deta[nei$d+(l-1)*n,]*(eta.cv[,k]-eta[nei$d,k])*(eta.cv[,j]-eta[nei$d,j])
	  #if (l!=k) ncv1 <- ncv1 - (1+(l!=j&&j!=k)) * gamma * .5 * colSums(
	  #        llf$l3[nei$d,jj]*deta[nei$d+(k-1)*n,]*(eta.cv[,l]-eta[nei$d,l])*(eta.cv[,j]-eta[nei$d,j]))
	  if (l!=k) ncv1 <- ncv1 - (1+(l!=j&&j!=k)) * gamma * .5 * 
	            llf$l3[nei$d,jj]*deta[nei$d+(k-1)*n,]*(eta.cv[,l]-eta[nei$d,l])*(eta.cv[,j]-eta[nei$d,j])
	  #if (l!=j) ncv1 <- ncv1 - (1+(l!=k&&j!=k)) * gamma * .5 * colSums(
	  #        llf$l3[nei$d,jj]*deta[nei$d+(j-1)*n,]*(eta.cv[,k]-eta[nei$d,k])*(eta.cv[,l]-eta[nei$d,l]))
          if (l!=j) ncv1 <- ncv1 - (1+(l!=k&&j!=k)) * gamma * .5 * 
	          llf$l3[nei$d,jj]*deta[nei$d+(j-1)*n,]*(eta.cv[,k]-eta[nei$d,k])*(eta.cv[,l]-eta[nei$d,l])
        }
      }
    } 
  } else { ## exact
    offi <- offset
    if (!is.null(offset)) for (i in 1:length(offset)) if (!is.null(offset[[i]])) offi[[i]] <- offset[[i]][nei$d]
    ll <- family$ll(y[nei$d],X[nei$d,],beta,wt[nei$d],family,offi,deriv=1,db,eta=eta.cv,ncv=TRUE)
    ncv <- -ll$l
    ncv <- gamma*ncv - (gamma-1)*dev
    if (deriv>0) {
      dev1 <- ncv1 <- matrix(0,nm,nsp) #rep(0,nsp)
      ind <- 1:nm; iin <- 1:n
      for (i in 1:nlp) {
        #ncv1 <- ncv1 - colSums(ll$l1[,i]*deta.cv[ind,])
        ncv1 <- ncv1 - ll$l1[,i]*deta.cv[ind,]
        #if (gamma!=1) dev1 <- dev1 - colSums((llf$l1[,i]*deta[iin,])[nei$d,,drop=FALSE])
	if (gamma!=1) dev1 <- dev1 - (llf$l1[,i]*deta[iin,])[nei$d,,drop=FALSE]
        ind <- ind + nm; iin <- iin + n
      }
      ncv1 <- gamma*ncv1 - (gamma-1)*dev1
    } 
  }
  if (deriv>0) {
    Vg <- crossprod(ncv1)
    ncv1 <- colSums(ncv1)
  } else Vg <- NULL
  attr(ncv,"eta.cv") <- eta.cv
  if (deriv!=0) attr(ncv,"deta.cv") <- deta.cv ## actually the perturbations if deriv<0
  return(list(NCV=ncv,NCV1=ncv1,error=cg.iter,Vg=Vg))
} ## gamlss.ncv

gamlss.etamu <- function(l1,l2,l3=NULL,l4=NULL,ig1,g2,g3=NULL,g4=NULL,i2,i3=NULL,i4=NULL,deriv=0) {
## lj is the array of jth order derivatives of l
## gj[,k] contains the jth derivatives for the link of the kth lp
## ig1 is one over first deriv of link
## kth parameter. This routine transforms derivatives 
## w.r.t. the parameters (mu_1..mu_K) to derivatives
## w.r.t. the linear predictors (eta_1.. eta_K)
## i2, i3 and i4 are the upper triangular indexing arrays
## e.g. l4[,i4[i,j,l,m]] contains the partial w.r.t.  
## params indexed by i,j,l,m with no restriction on
## the index values except that they are in 1..K
## lj may have an attribute "remap". If so then columns containing
## only zeros are not actually stored. So, for example, if
## remap is the attribute for l4 and k = remap[i4[i,j,l,m]] then
## the derivative is zero if k==0 and l4[,k] otherwise. This
## allows for situations in which K is quite high, but there
## are many zero derivatives. Note, however that a zero col
## in l4 does not always imply a zero column in d4 (deriv w.r.t.
## eta), since the latter often involves lower order derivatives
## l3, l2 etc. The same goes for l3 and l2.
## Returned arrays have remap attributes if input arrays do -
## they will generally be different to the input versions. 

  ordf <- function(i,j,l=NULL,m=NULL) {
  ## helper function to work out derivative orders
   if (is.null(l)) { ## 2d
     ord <- rep(1,2)
     if (i==j) {ord[1] <- ord[1] + 1; ord[2] <- 0 }
   } else if (is.null(m)) { ## 3 d 
      ord <- rep(1,3)
      if (i==j) {ord[1] <- ord[1] + 1; ord[2] <- 0 }
      if (i==l) {ord[1] <- ord[1] + 1; ord[3] <- 0 }
      if (ord[2]) {
        if (j==l) {ord[2] <- ord[2] + 1; ord[3] <- 0 }
      }
    } else { ## 4 d
      ord <- rep(1,4)
      if (i==j) {ord[1] <- ord[1] + 1; ord[2] <- 0 }
      if (i==l) {ord[1] <- ord[1] + 1; ord[3] <- 0 }
      if (i==m) {ord[1] <- ord[1] + 1; ord[4] <- 0 }
      if (ord[2]) {
        if (j==l) {ord[2] <- ord[2] + 1; ord[3] <- 0 }
        if (j==m) {ord[2] <- ord[2] + 1; ord[4] <- 0 }
      }
      if (ord[3]&&l==m) { ord[3] <- ord[3] + 1; ord[4] <- 0 }
    }
    ord
  } ## ordf

  K <- ncol(l1) ## number of parameters of distribution
  l1map <- attr(l1,"remap") ## lmap[i] is col of l1 storing ith deriv or zero if deriv zero
  d1 <- l1 ## not "remap" matches l1
  if (is.null(l1map)) { ## all derivs stored explicitly in l1
    for (i in 1:K) { ## first derivative loop
      d1[,i] <- l1[,i]*ig1[,i]
    }
  } else { ## some derivative are zero and not stored in l1
    for (ii in 1:K) { ## first derivative loop
      i <- l1map[ii] ## actual column in which iith deriv stored, or 0 if deriv zero.
      if (i>0) d1[,i] <- l1[,i]*ig1[,i]
    }
  }

  ifunc <- !is.array(i2) ## are index functions provided in place of index arrays?
 
  l2map <- attr(l2,"remap")
  g2zero <- colMeans(abs(g2))==0
  if (!is.null(l2map)||!is.null(l1map)) { ## l1 and or l2 are supplied with missing zero cols
    if (is.null(l1map)) l1map <- 1:K
    K2 <- ((K+1)*K)/2 ## number of second derivatives in total
    d2map <- 1:K2
    if (is.null(l2map)) l2map <- 1:K2 else { ## need to do a dummy run to establish which elements of d2 are non-zero
      k <- 0
      for (i in 1:K) for (j in i:K) {
        ord <- ordf(i,j); k <- k+1
        mo <- max(ord)
        if (mo==2) { ## pure 2nd derivative transform
          if (l2map[k]==0 && (l1map[i]==0||g2zero[i])) d2map[k] <- 0 ## d2[,k] zero as components zero.
        } else { ## all first derivative
          if (l2map[k]==0) d2map[k] <- 0 
        }
      }
      K2 <- sum(d2map!=0)
      d2map[d2map!=0] <- 1:K2
    }
    ## Now know d2map and l1map and l2map both exist. Do transforms...
    d2 <- matrix(0,nrow(l2),K2)
    for (i in 1:K) for (j in i:K) {
      ord <- ordf(i,j);k <- k+1
      mo <- max(ord)
      if (d2map[k]>0) { ## non-zero term
        if (mo==2) { ## pure 2nd derivative transform
          a <- if (l2map[k]>0) l2[,l2map[k]] else 0
	  b <- if (l1map[i]>0) l1[,l1map[i]]*g2[,i]*ig1[,i] else 0
          d2[,d2map[k]] <- (a - b)*ig1[,i]^2   
        } else { ## all first derivative
          d2[,d2map[k]] <- l2[,l2map[k]]*ig1[,i]*ig1[,j]
        }
      }
    }
    attr(d2,"remap") <- d2map
  } else { ## l1 and or l2 are supplied complete
    k <- 0; d2 <- l2
    for (i in 1:K) for (j in i:K) {
      ## obtain the order of differentiation associated 
      ## with the i,j derivatives...
      ord <- ordf(i,j);k <- k+1
      ## l2[,k] is derivative to transform
      mo <- max(ord)
      if (mo==2) { ## pure 2nd derivative transform
        d2[,k] <- (l2[,k] - l1[,i]*g2[,i]*ig1[,i])*ig1[,i]^2   
      } else { ## all first derivative
        d2[,k] <- l2[,k]*ig1[,i]*ig1[,j]
      }
    }
  } ## 2nd order transform done

  l3map <- attr(l3,"remap")
  if (deriv>0) g3zero <- colMeans(abs(g3))==0
  if (deriv>0&& (!is.null(l3map)||!is.null(l2map)||!is.null(l1map))) { ## 3rd order required, but some lj have dropped zero cols
    if (is.null(l1map)) l1map <- 1:K
    if (is.null(l2map)) { K2 <- ((K+1)*K)/2; l2map <- 1:K2 }  
    K3 <- (K*(K+3)+2)*K/6
    d3map <- 1:K3
    if (is.null(l3map)) l3map <- 1:K3 else { ## dummy run to work out d3map
      k <- 0
      for (i in 1:K) for (j in i:K) for (l in j:K) {
        k <- k + 1
	ord <- ordf(i,j,l)
        ii <- c(i,j,l)
        ## l3[,k] is derivative to transform
        mo <- max(ord)
        if (mo==3) { ## pure 3rd derivative transform
          mind <- if (ifunc) i2(i,i) else i2[i,i]
	  if (l3map[k]==0&&(l2map[mind]==0||g2zero[i])&&(l1map[i]==0||(g3zero[i]&&g2zero[i]))) d3map[k] <- 0      
        } else if (mo==1) { ## all first derivative
          if (l3map[k]==0) d3map[k] <- 0
        } else { ## 2,1 deriv
          k1 <- ii[ord==1] ## index of order 1 deriv
          k2 <- ii[ord==2] ## index of order 2 part
          mind <- if (ifunc) i2(k2,k1) else i2[k2,k1]
          if (l3map[k]==0&&(l2map[mind]==0||g2zero[k2])) d3map[k] <- 0
        }
      }
      K3 <- sum(d3map!=0)
      d3map[d3map!=0] <- 1:K3
    }
    ## now create and fill in non-zero cols of d3... 
    d3 <- matrix(0,nrow(l3),K3)
    k <- 0
    for (i in 1:K) for (j in i:K) for (l in j:K) {
      ## obtain the order of differentiation associated 
      ## with the i,j,l derivatives...
      k <- k+1
      if (d3map[k]>0) {
        ord <- ordf(i,j,l)
        ii <- c(i,j,l)
        ## l3[,k] is derivative to transform
        mo <- max(ord)
        if (mo==3) { ## pure 3rd derivative transform
          mind <- if (ifunc) i2(i,i) else i2[i,i]
	  aa <- if (l3map[k]>0) l3[,l3map[k]] else 0
	  bb <- if (l2map[mind]>0) -3*l2[,l2map[mind]]*g2[,i]*ig1[,i] else 0
	  cc <- if (l1map[i]>0) l1[,l1map[i]]*(3*g2[,i]^2*ig1[,i]^2 - g3[,i]*ig1[,i]) else 0
          d3[,d3map[k]] <- (aa + bb + cc)*ig1[,i]^3         
        } else if (mo==1) { ## all first derivative
          d3[,d3map[k]] <- l3[,l3map[k]]*ig1[,i]*ig1[,j]*ig1[,l]
        } else { ## 2,1 deriv
          k1 <- ii[ord==1] ## index of order 1 deriv
          k2 <- ii[ord==2] ## index of order 2 part
          mind <- if (ifunc) i2(k2,k1) else i2[k2,k1]
          aa <- if (l3map[k]>0) l3[,l3map[k]] else 0
	  bb <- if (l2map[mind]>0) -l2[,l2map[mind]]*g2[,k2]*ig1[,k2] else 0
          d3[,d3map[k]] <- (aa+bb)*ig1[,k1]*ig1[,k2]^2
        } 
      } ## d3map[k]>0
    } ## loop
    attr(d3,"remap") <- d3map
  } else { ## l3, l2 and l1 all supplied complete without zero cols dropped
    k <- 0
    d3 <- l3
    if (deriv>0) for (i in 1:K) for (j in i:K) for (l in j:K) {
      ## obtain the order of differentiation associated 
      ## with the i,j,l derivatives...
      ord <- ordf(i,j,l);k <- k+1
      ii <- c(i,j,l)
      ## l3[,k] is derivative to transform
      mo <- max(ord)
      if (mo==3) { ## pure 3rd derivative transform
        mind <- if (ifunc) i2(i,i) else i2[i,i]
        d3[,k] <- (l3[,k] - 3*l2[,mind]*g2[,i]*ig1[,i] +
                l1[,i]*(3*g2[,i]^2*ig1[,i]^2 - g3[,i]*ig1[,i]))*ig1[,i]^3          
      } else if (mo==1) { ## all first derivative
        d3[,k] <- l3[,k]*ig1[,i]*ig1[,j]*ig1[,l]
      } else { ## 2,1 deriv
        k1 <- ii[ord==1] ## index of order 1 deriv
        k2 <- ii[ord==2] ## index of order 2 part
        mind <- if (ifunc) i2(k2,k1) else i2[k2,k1]
        d3[,k] <- (l3[,k] - l2[,mind]*g2[,k2]*ig1[,k2])*
                ig1[,k1]*ig1[,k2]^2
      } 
    }
  }  ## 3rd order transform done

  l4map <- attr(l4,"remap")
  if (deriv>2&& (!is.null(l4map)||!is.null(l3map)||!is.null(l2map)||!is.null(l1map))) { ## 4th order required, but some lj have dropped zero cols
    g4zero <- colMeans(abs(g4))==0
    if (is.null(l1map)) l1map <- 1:K
    if (is.null(l2map)) { K2 <- ((K+1)*K)/2; l2map <- 1:K2 }  
    if (is.null(l3map)) { K3 <- (K*(K+3)+2)*K/6;l3map <- 1:K3}
    K4 <- (6+K*11+K^2*6+K^3)*K/24
    d4map <- 1:K4
    if (is.null(l4map)) l4map <- 1:K4 else { ## dummy run to create d4map
      k <- 0
      for (i in 1:K) for (j in i:K) for (l in j:K) for (m in l:K) { 
        ## obtain the order of differentiation associated 
        ## with the i,j,l & m derivatives...
        ord <- ordf(i,j,l,m);k <- k+1
        ii <- c(i,j,l,m)
        ## l4[,k] is derivative to transform
        mo <- max(ord)
        if (mo==4) { ## pure 4th derivative transform
          mi2 <- if (ifunc) i2(i,i) else i2[i,i]
          mi3 <- if (ifunc) i3(i,i,i) else i3[i,i,i]
          if (l4map[k]==0&&(l3map[mi3]==0||g2zero[i])&&(l2map[mi2]==0||(g2zero[i]&&g3zero[i]))&&(l1map[i]==0||(g4zero[i]&&g2zero[i]))) d4map[k] <- 0
        } else if (mo==1) { ## all first derivative
          if (l4map[k]==0) d4map[k] <- 0
        } else if (mo==3) { ## 3,1 deriv
          k1 <- ii[ord==1] ## index of order 1 deriv
          k3 <- ii[ord==3] ## index of order 3 part
          mi2 <- if (ifunc) i2(k3,k1) else i2[k3,k1]
          mi3 <- if (ifunc) i3(k3,k3,k1) else i3[k3,k3,k1]
	  if (l4map[k]==0&&(l3map[mi3]==0||g2zero[k3])&&(l2map[mi2]==0||(g2zero[k3]&&g3zero[k3]))) d4map[k] <- 0
        } else { 
          if (sum(ord==2)==2) { ## 2,2
            k2a <- (ii[ord==2])[1];k2b <- (ii[ord==2])[2]
	    mi2 <- if (ifunc) i2(k2a,k2b) else i2[k2a,k2b]
	    mi3 <- if (ifunc) i3(k2a,k2b,k2b) else i3[k2a,k2b,k2b]
	    mi3a <- if (ifunc) i3(k2a,k2a,k2b) else i3[k2a,k2a,k2b]
	    if (l4map[k]==0&&(l3map[mi3]==0||g2zero[k2a])&&(l3map[mi3a]==0||g2zero[k2b])&&(l2map[mi2]==0||g2zero[k2a]||g2zero[k2b])) d4map[k] <- 0
          } else { ## 2,1,1
            k2 <- ii[ord==2] ## index of order 2 derivative
            k1a <- (ii[ord==1])[1];k1b <- (ii[ord==1])[2]
	    mi3 <- if (ifunc) i3(k2,k1a,k1b) else i3[k2,k1a,k1b]
	    if (l4map[k]==0&&(l3map[mi3]==0||g2zero[k2])) d4map[k] <- 0
          }
        }
      } ## loop	
    }
    K4 <- sum(d4map!=0)
    d4map[d4map!=0] <- 1:K4
    d4 <- matrix(0,nrow(l4),K4)
    k <- 0
    for (i in 1:K) for (j in i:K) for (l in j:K) for (m in l:K) { ## fill in d4
      ## obtain the order of differentiation associated 
      ## with the i,j,l & m derivatives...
      k <- k+1
      if (d4map[k]>0) {
        ord <- ordf(i,j,l,m);
        ii <- c(i,j,l,m)
        ## l4[,k] is derivative to transform
        mo <- max(ord)
        if (mo==4) { ## pure 4th derivative transform
          mi2 <- if (ifunc) i2(i,i) else i2[i,i]
          mi3 <- if (ifunc) i3(i,i,i) else i3[i,i,i]
          aa <- if (l4map[k]>0) l4[,l4map[k]] else 0
	  bb <- if (l3map[mi3]>0) -6*l3[,l3map[mi3]]*g2[,i]*ig1[,i] else 0
	  cc <- if (l2map[mi2]>0) l2[,l2map[mi2]]*(15*g2[,i]^2*ig1[,i]^2 - 4*g3[,i]*ig1[,i]) else 0
	  dd <- if (l1map[i]>0) -l1[,l1map[i]]*(15*g2[,i]^3*ig1[,i]^3 - 10*g2[,i]*g3[,i]*ig1[,i]^2 + g4[,i]*ig1[,i]) else 0 
          d4[,d4map[k]] <- (aa+bb+cc+dd)*ig1[,i]^4    
        } else if (mo==1) { ## all first derivative
          d4[,d4map[k]] <- l4[,l4map[k]]*ig1[,i]*ig1[,j]*ig1[,l]*ig1[,m]
        } else if (mo==3) { ## 3,1 deriv
          k1 <- ii[ord==1] ## index of order 1 deriv
          k3 <- ii[ord==3] ## index of order 3 part
          mi2 <- if (ifunc) i2(k3,k1) else i2[k3,k1]
          mi3 <- if (ifunc) i3(k3,k3,k1) else i3[k3,k3,k1]
	  aa <- if (l4map[k]>0) l4[,l4map[k]] else 0
	  bb <- if (l3map[mi3]>0) -3*l3[,l3map[mi3]]*g2[,k3]*ig1[,k3] else 0
	  cc <- if (l2map[mi2]>0) l2[,l2map[mi2]]*(3*g2[,k3]^2*ig1[,k3]^2 - g3[,k3]*ig1[,k3]) else 0
          d4[,d4map[k]] <- (aa+bb+cc)*ig1[,k1]*ig1[,k3]^3
        } else { 
          if (sum(ord==2)==2) { ## 2,2
            k2a <- (ii[ord==2])[1];k2b <- (ii[ord==2])[2]
	    mi2 <- if (ifunc) i2(k2a,k2b) else i2[k2a,k2b]
	    mi3 <- if (ifunc) i3(k2a,k2b,k2b) else i3[k2a,k2b,k2b]
	    mi3a <- if (ifunc) i3(k2a,k2a,k2b) else i3[k2a,k2a,k2b]
	    aa <- if (l4map[k]>0) l4[,l4map[k]] else 0
	    bb <- if (l3map[mi3]>0) -l3[,l3map[mi3]]*g2[,k2a]*ig1[,k2a] else 0
	    cc <- if (l3map[mi3a]>0) -l3[,l3map[mi3a]]*g2[,k2b]*ig1[,k2b] else 0
	    dd <- if (l2map[mi2]>0) l2[,l2map[mi2]]*g2[,k2a]*g2[,k2b]*ig1[,k2a]*ig1[,k2b] else 0
            d4[,d4map[k]] <- (aa+bb+cc+dd)*ig1[,k2a]^2*ig1[,k2b]^2
          } else { ## 2,1,1
            k2 <- ii[ord==2] ## index of order 2 derivative
            k1a <- (ii[ord==1])[1];k1b <- (ii[ord==1])[2]
	    mi3 <- if (ifunc) i3(k2,k1a,k1b) else i3[k2,k1a,k1b]
	    aa <- if (l4map[k]>0) l4[,l4map[k]] else 0
	    bb <- if (l3map[mi3]>0) -l3[,l3map[mi3]]*g2[,k2]*ig1[,k2] else 0
            d4[,d4map[k]] <- (aa+bb)*ig1[,k1a]*ig1[,k1b]*ig1[,k2]^2
          }
        }
      } ## if d4map[k]>0	
    } ## loop  
    attr(d4,"remap") <- d4map
  } else { ## l1-l4 are all supplied complete with no dropping 
    k <- 0
    d4 <- l4
    if (deriv>2) for (i in 1:K) for (j in i:K) for (l in j:K) for (m in l:K) {
      ## obtain the order of differentiation associated 
      ## with the i,j,l & m derivatives...
      ord <- ordf(i,j,l,m);k <- k+1
      ii <- c(i,j,l,m)
      ## l4[,k] is derivative to transform
      mo <- max(ord)
      if (mo==4) { ## pure 4th derivative transform
      mi2 <- if (ifunc) i2(i,i) else i2[i,i]
      mi3 <- if (ifunc) i3(i,i,i) else i3[i,i,i]
      d4[,k] <-  (l4[,k] - 6*l3[,mi3]*g2[,i]*ig1[,i] + 
        l2[,mi2]*(15*g2[,i]^2*ig1[,i]^2 - 4*g3[,i]*ig1[,i]) - 
        l1[,i]*(15*g2[,i]^3*ig1[,i]^3 - 10*g2[,i]*g3[,i]*ig1[,i]^2 
         + g4[,i]*ig1[,i]))*ig1[,i]^4    
      } else if (mo==1) { ## all first derivative
        d4[,k] <- l4[,k]*ig1[,i]*ig1[,j]*ig1[,l]*ig1[,m]
      } else if (mo==3) { ## 3,1 deriv
        k1 <- ii[ord==1] ## index of order 1 deriv
        k3 <- ii[ord==3] ## index of order 3 part
        mi2 <- if (ifunc) i2(k3,k1) else i2[k3,k1]
        mi3 <- if (ifunc) i3(k3,k3,k1) else i3[k3,k3,k1]
        d4[,k] <- (l4[,k] - 3*l3[,mi3]*g2[,k3]*ig1[,k3] +
        l2[,mi2]*(3*g2[,k3]^2*ig1[,k3]^2 - g3[,k3]*ig1[,k3])         
        )*ig1[,k1]*ig1[,k3]^3
      } else { 
        if (sum(ord==2)==2) { ## 2,2
          k2a <- (ii[ord==2])[1];k2b <- (ii[ord==2])[2]
	  mi2 <- if (ifunc) i2(k2a,k2b) else i2[k2a,k2b]
	  mi3 <- if (ifunc) i3(k2a,k2b,k2b) else i3[k2a,k2b,k2b]
	  mi3a <- if (ifunc) i3(k2a,k2a,k2b) else i3[k2a,k2a,k2b]
          d4[,k] <- (l4[,k] - l3[,mi3]*g2[,k2a]*ig1[,k2a]
          -l3[,mi3a]*g2[,k2b]*ig1[,k2b] + 
           l2[,mi2]*g2[,k2a]*g2[,k2b]*ig1[,k2a]*ig1[,k2b]
          )*ig1[,k2a]^2*ig1[,k2b]^2
        } else { ## 2,1,1
          k2 <- ii[ord==2] ## index of order 2 derivative
          k1a <- (ii[ord==1])[1];k1b <- (ii[ord==1])[2]
	  mi3 <- if (ifunc) i3(k2,k1a,k1b) else i3[k2,k1a,k1b]
          d4[,k] <- (l4[,k] - l3[,mi3]*g2[,k2]*ig1[,k2] 
                   )*ig1[,k1a]*ig1[,k1b]*ig1[,k2]^2
        }
      }
    }  
  } ## 4th order transform done

  list(l1=d1,l2=d2,l3=d3,l4=d4)
} # gamlss.etamu


gamlss.gH <- function(X,jj,l1,l2,i2,l3=0,i3=0,l4=0,i4=0,d1b=0,d2b=0,deriv=0,fh=NULL,D=NULL,sandwich=FALSE) {
## X[,jj[[i]]] is the ith model matrix.
## lj contains jth derivatives of the likelihood for each datum,
## columns are w.r.t. different combinations of parameters.
## ij is the symmetric array indexer for the jth order derivs...
## e.g. l4[,i4[i,j,l,m]] contains derivatives with 
## respect to parameters indexed by i,j,l,m
## d1b and d2b are first and second derivatives of beta w.r.t. sps.
## fh is a factorization of the penalized hessian, while D contains the corresponding
##    Diagonal pre-conditioning weights.
## deriv: 0 - just grad and Hess
##        1 - tr(Hp^{-1} dH/drho_j) vector - Hp^{-1} must be supplied in fh
##        2 - first deriv of Hess
##        3 - everything.
  K <- length(jj)
  if (is.list(X)) {
    discrete <- TRUE
    p <- X$p;n <- nrow(X$kd)
    sparse <- !inherits(X$Xd[[1]],"matrix")
  } else {
    sparse <- discrete <- FALSE
    p <- ncol(X);n <- nrow(X)
  }  
  trHid2H <- d1H <- d2H <- NULL ## defaults
  ifunc <- !is.array(i2) ## are index functions provided in place of index arrays?

  ## the gradient...
  l1map <- attr(l1,"remap")
  lb <- rep(0,p)
  if (is.null(l1map)) { ## all cols of l1 supplied
    for (i in 1:K) { ## first derivative loop
      lb[jj[[i]]] <- lb[jj[[i]]] + if (discrete) XWyd(X$Xd,rep(1,n),l1[,i],X$kd,X$ks,X$ts,X$dt,X$v,X$qc,X$drop,lt=X$lpid[[i]]) else
                     colSums(l1[,i]*X[,jj[[i]],drop=FALSE]) ## !
    }
  } else { ## only non-zero cols of l1 supplied (only supplied for completeness - unclear any sensible model could ever want this)
    for (i in 1:K) if (l1map[i]>0) { ## first derivative loop
      lb[jj[[i]]] <- lb[jj[[i]]] + if (discrete) XWyd(X$Xd,rep(1,n),l1[,l1map[i]],X$kd,X$ks,X$ts,X$dt,X$v,X$qc,X$drop,lt=X$lpid[[i]]) else
                     colSums(l1[,l1map[i]]*X[,jj[[i]],drop=FALSE]) ## !
    }
  }
  
  ## the Hessian...
  lbb <- if (sparse) Matrix(0,p,p) else matrix(0,p,p)
  if (sandwich) { ## reset l2 so that Hessian becomes 'filling' for sandwich estimate
    if (deriv>0) warning("sandwich requested with higher derivatives")
    if (!is.null(l1map)) stop("sandwich requested with structurally zero first derivatives - can't be sensible")
    k <- 0;
    for (i in 1:K) for (j in i:K) { k <- k + 1;l2[,k] <- l1[,i]*l1[,j] }
    attr(l2,"remap") <- NULL ## l2 has to be full now.
  }

  l2map <- attr(l2,"remap")
  if (is.null(l2map)) { ## l2 supplied with all columns
    for (i in 1:K) for (j in i:K) {
      ## A <- t(X[,jj[[i]],drop=FALSE])%*%(l2[,i2[i,j]]*X[,jj[[j]],drop=FALSE])
      mi2 <- if (ifunc) i2(i,j) else i2[i,j] 
      A <- if (discrete) XWXd(X$Xd,w=l2[,mi2],k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,nthreads=1,drop=X$drop,lt=X$lpid[[i]],rt=X$lpid[[j]]) else
             crossprod(X[,jj[[i]],drop=FALSE],l2[,mi2]*X[,jj[[j]],drop=FALSE])
      lbb[jj[[i]],jj[[j]]] <- lbb[jj[[i]],jj[[j]]] + A 
      if (j>i) lbb[jj[[j]],jj[[i]]] <- lbb[jj[[j]],jj[[i]]] + t(A) 
    } 
  } else { ## l2 supplied with zero columns dropped
     for (i in 1:K) for (j in i:K) {
      mi2 <- if (ifunc) i2(i,j) else i2[i,j]
      if (l2map[mi2]>0) {
        A <- if (discrete) XWXd(X$Xd,w=l2[,l2map[mi2]],k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,nthreads=1,drop=X$drop,lt=X$lpid[[i]],rt=X$lpid[[j]]) else
             crossprod(X[,jj[[i]],drop=FALSE],l2[,l2map[mi2]]*X[,jj[[j]],drop=FALSE])
        lbb[jj[[i]],jj[[j]]] <- lbb[jj[[i]],jj[[j]]] + A 
        if (j>i) lbb[jj[[j]],jj[[i]]] <- lbb[jj[[j]],jj[[i]]] + t(A)
      }	
    }
  }


  if (deriv>0) {
    ## the first derivative of the Hessian, using d1b
    ## the first derivates of the coefficients wrt the sps
    m <- ncol(d1b) ## number of smoothing parameters
    ## stack the derivatives of the various linear predictors on top
    ## of each other...
    d1eta <- matrix(0,n*K,m)
    ind <- 1:n
    for (i in 1:K) { 
      d1eta[ind,] <- if (discrete) Xbd(X$Xd,d1b,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[i]]) else
                                   X[,jj[[i]],drop=FALSE]%*%d1b[jj[[i]],]
      ind <- ind + n
    }
    l3map <- attr(l3,"remap")
    null3map <- is.null(l3map)
  }

  if (deriv==1) { 
   ## assuming fh contains the inverse penalized Hessian, Hp, forms tr(Hp^{-1}dH/drho_j) for each j
   g.index <- attr(d1b,"g.index") ## possible index indicating log parameterization
   if (!is.null(g.index)) { ## then several transform related quantities are required 
     beta <- attr(d1b,"beta") ##  regression coefficients 
     d1g <- d1b; d1g[g.index,] <- d1g[g.index,]/beta[g.index] ## derivartive w.r.t. working parameters
     hess.diag <- attr(d1b,"hess.diag") ## should diagonal correction terms be included?
   }
   d1H <- rep(0,m)
   if (discrete) {
     ## lpi <- attr(X,"lpi") ## this line was in original code for this discrete section, and lpi replaced jj below - mistake, I think
     for (i in 1:K) for (j in i:K) { ## lp block loop
       for (l in 1:m) { ## sp loop
         v <- rep(0,n);ind <- 1:n
         for (q in 1:K) { ## diagonal accumulation loop
	   mi3 <- if (ifunc) i3(i,j,q) else i3[i,j,q]
	   if (!null3map) mi3 <- l3map[mi3]
           if (mi3>0) v <- v + l3[,mi3] * d1eta[ind,l]
           ind <- ind + n
         }
	 XVX <- XWXd(X$Xd,w=v,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,nthreads=1,drop=X$drop,lt=X$lpid[[i]],rt=X$lpid[[j]])
	 if (!is.null(g.index)) { ## non-linear correction terms required
           gi <- g.index[jj[[i]]];gj <- g.index[jj[[j]]]
	   if (any(gi)) XVX[gi,] <- beta[jj[[i]]][gi]*XVX[gi,]
	   if (any(gj)) XVX[,gj] <- t(beta[jj[[j]]][gj]*t(XVX[,gj]))
	   if (any(gi)) {
	     XWX <- beta[jj[[i]]][gi]*d1g[jj[[i]],l][gi]*lbb[jj[[i]],jj[[j]]][gi,]
	     if (any(gj)) XWX[,gj] <- t(beta[jj[[j]]][gj]*t(XWX[,gj]))
	     XVX[gi,] <- XVX[gi,] + XWX   
	   }  
	   if (any(gj)) {
	     XWX <- t(beta[jj[[j]]][gj]*d1g[jj[[j]],l][gj]*t(lbb[jj[[i]],jj[[j]]][,gj]))
	     if (any(gi)) XWX[gi,] <- beta[jj[[i]]][gi]*XWX[gi,]
	     XVX[,gj] <- XVX[,gj] + XWX
	     if (i==j&&hess.diag) { ## add diagonal corrections
	       dd <- beta[jj[[i]]][gi]*(lbb[jj[[i]][gi],] %*% d1b[,l] + lb[jj[[i]]][gi]*d1g[jj[[i]],l][gi])
	       XVX[gi,gj] <- XVX[gi,gj] + diag(drop(dd),nrow=sum(gi))
             }
           }
         } ## end of non-linear corrections
	 mult <- if (i==j) 1 else 2
	 d1H[l] <- d1H[l] + mult * sum(XVX * fh[jj[[i]],jj[[j]]]) ## accumulate tr(Hp^{-1}dH/drho_l)
       }
     }  
   } else for (i in 1:K) for (j in i:K) { ## lp block loop
      Hpi <- fh[jj[[i]],jj[[j]]] ## correct component of inverse Hessian
      d1hc <- rep(0,m)
      if (!is.null(g.index)) { ## correct for non-linearity
        gi <- g.index[jj[[i]]];gj <- g.index[jj[[j]]]
        for (l in 1:m) { ## s.p. loop
          dcor <- 0
	  if (any(gi)) {
	    XWX <- beta[jj[[i]]][gi]*d1g[jj[[i]],l][gi]*lbb[jj[[i]],jj[[j]]][gi,]
	    if (any(gj)) XWX[,gj] <- t(beta[jj[[j]]][gj]*t(XWX[,gj]))
	    dcor <- dcor + sum(XWX * Hpi[gi,])
	  }  
	  if (any(gj)) {
	    XWX <- t(beta[jj[[j]]][gj]*d1g[jj[[j]],l][gj]*t(lbb[jj[[i]],jj[[j]]][,gj]))
	    if (any(gi)) XWX[gi,] <- beta[jj[[i]]][gi]*XWX[gi,]
	    dcor <- dcor + sum(XWX * Hpi[,gj])
	    if (i==j&&hess.diag) { ## diagonal correction
               dd <- beta[jj[[i]]][gi]*(lbb[jj[[i]][gi],] %*% d1b[,l] + lb[jj[[i]]][gi]*d1g[jj[[i]],l][gi])
	       dcor <- dcor + sum(dd*diag(Hpi)[gi])
            }
	  }
	  d1hc[l] <- dcor
	} ## s.p. loop end
        if (any(gi)) Hpi[gi,] <- Hpi[gi,]*beta[jj[[i]]][gi]
        if (any(gj)) Hpi[,gj] <- t(t(Hpi[,gj])*beta[jj[[j]]][gj]) ## was jj[i] -- wrong
      } ## end of non-linearity correction
      a <- rowSums((X[,jj[[i]]] %*% Hpi) * X[,jj[[j]]])
      for (l in 1:m) { ## sp loop
        v <- rep(0,n);ind <- 1:n
        for (q in 1:K) { ## diagonal accumulation loop
	  mi3 <- if (ifunc) i3(i,j,q) else i3[i,j,q]
	  if (!null3map) mi3 <- l3map[mi3]
          if (mi3>0) v <- v + l3[,mi3] * d1eta[ind,l]
          ind <- ind + n
        }
	mult <- if (i==j) 1 else 2
	d1H[l] <- d1H[l] + mult * (sum(a*v) + d1hc[l]) ## accumulate tr(Hp^{-1}dH/drho_l)
      }
    }
  } ## if deriv==1

  if (deriv>1) {
    if (discrete) stop("er... no discrete methods for higher derivatives")
    d1H <- list()
    for (l in 1:m) {
      d1H[[l]] <- matrix(0,p,p)
      for (i in 1:K) for (j in i:K) {
        v <- rep(0,n);ind <- 1:n
        for (q in 1:K) {
	  mi3 <- if (ifunc) i3(i,j,q) else i3[i,j,q] 
          if (!null3map) mi3 <- l3map[mi3]
          if (mi3>0) v <- v + l3[,mi3] * d1eta[ind,l]
          ind <- ind + n
        }
        ## d1H[[l]][jj[[j]],jj[[i]]] <- 
        ## A <- t(X[,jj[[i]],drop=FALSE])%*%(v*X[,jj[[j]],drop=FALSE])
        A <- crossprod(X[,jj[[i]],drop=FALSE],v*X[,jj[[j]],drop=FALSE])
        d1H[[l]][jj[[i]],jj[[j]]] <- d1H[[l]][jj[[i]],jj[[j]]] + A
        if (j>i) d1H[[l]][jj[[j]],jj[[i]]] <- d1H[[l]][jj[[j]],jj[[i]]] + t(A)
      } 
    }
  } ## if deriv>1

  if (deriv>2) {
    l4map <- attr(l4,"remap")
    null4map <- is.null(l4map)
    ## need tr(Hp^{-1} d^2H/drho_k drho_j)
    ## First form the expanded model matrix...
    VX <- Xe <- matrix(0,K*n,ncol(X))
    ind <- 1:n 
    for (i in 1:K) { 
      Xe[ind,jj[[i]]] <- X[,jj[[i]]]
      ind <- ind + n
    }
    ## Now form Hp^{-1} Xe'...
    if (is.list(fh)) { ## then the supplied factor is an eigen-decomposition
       d <- fh$values;d[d>0] <- 1/d[d>0];d[d<=0] <- 0
       Xe <- t(D*((fh$vectors%*%(d*t(fh$vectors)))%*%(D*t(Xe))))
    } else { ## the supplied factor is a choleski factor
       ipiv <- piv <- attr(fh,"pivot");ipiv[piv] <- 1:p
       Xe <- t(D*(backsolve(fh,forwardsolve(t(fh),(D*t(Xe))[piv,]))[ipiv,]))
    }
    ## now compute the required trace terms
    d2eta <- matrix(0,n*K,ncol(d2b))
    ind <- 1:n
    for (i in 1:K) { 
      d2eta[ind,] <- X[,jj[[i]],drop=FALSE]%*%d2b[jj[[i]],]
      ind <- ind + n
    }
    trHid2H <- rep(0,ncol(d2b))
    kk <- 0 ## counter for second derivatives
    for (k in 1:m) for (l in k:m) { ## looping over smoothing parameters...
      kk <- kk + 1
      for (i in 1:K) for (j in 1:K) {
        v <- rep(0,n);ind <- 1:n
        for (q in 1:K) { ## accumulate the diagonal matrix for X_i'diag(v)X_j
	  mi3 <- if (ifunc) i3(i,j,q) else i3[i,j,q]
	  if (!null3map) mi3 <- l3map[mi3]
          if (mi3>0) v <- v + d2eta[ind,kk]*l3[,mi3]
          ins <- 1:n
          for (s in 1:K) {
	    mi4 <- if (ifunc) i4(i,j,q,s) else i4[i,j,q,s]
	    if (!null4map) mi4 <- l4map[mi4]
            if (mi4>0) v <- v + d1eta[ind,k]*d1eta[ins,l]*l4[,mi4]
            ins <- ins + n
          }
          ind <- ind + n
        }
        if (i==j) {
          rind <- 1:n + (i-1)*n
          VX[rind,jj[[i]]] <- v * X[,jj[[i]]]
        } else {
          rind1 <- 1:n + (i-1)*n
          rind2 <- 1:n + (j-1)*n
          VX[rind2,jj[[i]]] <- v * X[,jj[[i]]]
          VX[rind1,jj[[j]]] <- v * X[,jj[[j]]]
        }
      }
      trHid2H[kk] <- sum(Xe*VX)
    }
  } ## if deriv>2

  list(lb=lb,lbb=lbb,d1H=d1H,d2H=d2H,trHid2H=trHid2H)
} ## end of gamlss.gH 



gaulss <- function(link=list("identity","logb"),b=0.01) {
## Extended family for Gaussian location scale model...
## so mu is mu1 and tau=1/sig is mu2
## tau = 1/(b + exp(eta)) eta = log(1/tau - b)
## 1. get derivatives wrt mu, tau
## 2. get required link derivatives and tri indices.
## 3. transform derivs to derivs wrt eta (gamlss.etamu).
## 4. get the grad and Hessian etc for the model
##    via a call to gamlss.gH  
## the first derivatives of the log likelihood w.r.t
## the first and second parameters...
  
  ## first deal with links and their derivatives...
  if (length(link)!=2) stop("gaulss requires 2 links specified as character strings")
  okLinks <- list(c("inverse", "log", "identity","sqrt"),"logb")
  stats <- list()
  if (link[[1]] %in% okLinks[[1]]) stats[[1]] <- make.link(link[[1]]) else 
  stop(link[[1]]," link not available for mu parameter of gaulss")
  fam <- structure(list(link=link[[1]],canonical="none",linkfun=stats[[1]]$linkfun,
           mu.eta=stats[[1]]$mu.eta),
           class="family")
  fam <- fix.family.link(fam)
  stats[[1]]$d2link <- fam$d2link
  stats[[1]]$d3link <- fam$d3link
  stats[[1]]$d4link <- fam$d4link
  if (link[[2]] %in% okLinks[[2]]) { ## creating the logb link
    stats[[2]] <- list()
    stats[[2]]$valideta <- function(eta) TRUE 
    stats[[2]]$link = link[[2]]
    stats[[2]]$linkfun <- eval(parse(text=paste("function(mu) log(1/mu -",b,")")))
    stats[[2]]$linkinv <- eval(parse(text=paste("function(eta) 1/(exp(eta) +",b,")")))
    stats[[2]]$mu.eta <- eval(parse(text=
                         paste("function(eta) { ee <- exp(eta); -ee/(ee +",b,")^2 }")))
    stats[[2]]$d2link <-  eval(parse(text=
    paste("function(mu) { mub <- pmax(1 - mu *",b,",.Machine$double.eps);(2*mub-1)/(mub*mu)^2}" )))
    stats[[2]]$d3link <-  eval(parse(text=
    paste("function(mu) { mub <-  pmax(1 - mu *",b,",.Machine$double.eps);((1-mub)*mub*6-2)/(mub*mu)^3}" )))
    stats[[2]]$d4link <-  eval(parse(text=
    paste("function(mu) { mub <- pmax(1 - mu *",b,",.Machine$double.eps);(((24*mub-36)*mub+24)*mub-6)/(mub*mu)^4}")))
  } else stop(link[[2]]," link not available for precision parameter of gaulss")
  
  residuals <- function(object,type=c("deviance","pearson","response")) { 
      type <- match.arg(type)
      rsd <- object$y-object$fitted[,1]
      if (type=="response") return(rsd) else
      return((rsd*object$fitted[,2])) ## (y-mu)/sigma 
    } ## gaulss residuals
    
  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    ## in principle the following seems reasonable, but because no
    ## price is paid for the high null variance, it leads to silly
    ## % deviance explained...
    #er <- fitNull(G$y,G$family,G$w,G$offset,nlp=length(attr(G$X,"lpi")),tol=1e-7)
    #object$null.deviance <- sum(((object$y-er$mu[,1])*er$mu[,2])^2*G$w)
    object$null.deviance <- sum(((object$y-mean(object$y))*object$fitted[,2])^2)
  })

  ncv <- function(X,y,wt,nei,beta,family,llf,H=NULL,Hi=NULL,R=NULL,offset=NULL,dH=NULL,db=NULL,deriv=FALSE,nt=1) {
    gamlss.ncv(X,y,wt,nei,beta,family,llf,H=H,Hi=Hi,R=R,offset=offset,dH=dH,db=db,deriv=deriv,nt=nt)
  } ## ncv  

  ll <- function(y,X,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL,eta=NULL,ncv=FALSE,sandwich=FALSE) {
  ## function defining the gamlss Gaussian model log lik. 
  ## N(mu,sigma^2) parameterized in terms of mu and log(sigma)
  ## deriv: 0 - eval
  ##        1 - grad and Hess
  ##        2 - diagonal of first deriv of Hess
  ##        3 - first deriv of Hess
  ##        4 - everything.
    if (!is.null(offset)) offset[[3]] <- 0
    discrete <- is.list(X)
    jj <- attr(X,"lpi") ## extract linear predictor index
    if (is.null(eta)) {
      eta <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[1]]) else X[,jj[[1]],drop=FALSE]%*%coef[jj[[1]]]
      if (!is.null(offset[[1]])) eta <- eta + offset[[1]]
      eta1 <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[2]]) else X[,jj[[2]],drop=FALSE]%*%coef[jj[[2]]]
      if (!is.null(offset[[2]])) eta1 <- eta1 + offset[[2]]
    } else { ## eta supplied directly
      eta1 <- eta[,2]
      eta <- eta[,1]
    }
    
    mu <- family$linfo[[1]]$linkinv(eta)
    tau <-  family$linfo[[2]]$linkinv(eta1) ## tau = 1/sig here
    
    n <- length(y)
    l1 <- matrix(0,n,2)
    ymu <- y-mu;ymu2 <- ymu^2;tau2 <- tau^2
 
    l0 <- -.5 * ymu2 * tau2 - .5 * log(2*pi) + log(tau)
    l <- sum(l0)

    if (deriv>0) {

      l1[,1] <- tau2*ymu
      l1[,2] <- 1/tau - tau*ymu2  

      ## the second derivatives
      ## order mm,ms,ss
      l2 <- cbind(-tau2,2*l1[,1]/tau,-ymu2 - 1/tau2)
     
      ## need some link derivatives for derivative transform
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),family$linfo[[2]]$d2link(tau))
    }
 
    l3 <- l4 <- g3 <- g4 <- 0 ## defaults

    if (deriv>1) {
      ## the third derivatives
      ## order mmm,mms,mss,sss
      if (TRUE) {
        l3 <- cbind(0,-2*tau,2*ymu,2/tau^3)
      } else { ## test infrastructure for dropping zero columns
        l3 <- cbind(-2*tau,2*ymu,2/tau^3)
	attr(l3,"remap") <- c(0,1:3)
      }
      g3 <- cbind(family$linfo[[1]]$d3link(mu),family$linfo[[2]]$d3link(tau))
    }

    if (deriv>3) {
      ## the fourth derivatives
      ## order mmmm,mmms,mmss,msss,ssss
      if (TRUE) {
        l4 <- cbind(0,0,-2,0,-6/tau2^2)
      } else { ## illustrates/tests 0 col dropping
        l4 <- cbind(-2,-6/tau2^2)
	attr(l4,"remap") <- c(0,0,1,0,2)
      }
      g4 <- cbind(family$linfo[[1]]$d4link(mu),family$linfo[[2]]$d4link(tau))
    }
    if (deriv) {
      i2 <- family$tri$i2; i3 <- family$tri$i3
      i4 <- family$tri$i4
   
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)

      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                      d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D,sandwich=sandwich)
      if (ncv) {
        ret$l1 <- de$l1; ret$l2 = de$l2; ret$l3 = de$l3
      }
    } else ret <- list()
    ret$l <- l; ret$l0 <- l0; ret
  } ## end ll gaulss

  sandwich <- function(y,X,coef,wt,family,offset=NULL) {
  ## compute filling for sandwich estimate of cov matrix
    ll(y,X,coef,wt,family,offset=NULL,deriv=1,sandwich=TRUE)$lbb
  }
  
  initialize <- expression({ ## init gaulss
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
	if (!is.null(offset)) offset[[3]] <- 0
        yt1 <- if (family$link[[1]]=="identity") y else 
               family$linfo[[1]]$linkfun(abs(y)+max(y)*1e-7)
	if (!is.null(offset[[1]])) yt1 <- yt1 - offset[[1]]
        if (is.list(x)) { ## discrete case
	  start <- rep(0,max(unlist(jj)))
	  R <- suppressWarnings(mchol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,
	              v=x$v,qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[1]])+crossprod(E[,jj[[1]]])))
	  Xty <- XWyd(x$Xd,rep(1,length(y)),yt1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[1]])
          piv <- attr(R,"pivot")
	  rrank <- attr(R,"rank") 
	  startji <- rep(0,ncol(R))
	  if (rrank<ncol(R)) {
            R <- R[1:rrank,1:rrank]
	    piv <- piv[1:rrank]
          }
          startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
	  startji[!is.finite(startji)] <- 0
	  start[jj[[1]]] <- startji
	  eta1 <- Xbd(x$Xd,start,k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,qc=x$qc,drop=x$drop,lt=x$lpid[[1]])
	  lres1 <- log(abs(y-family$linfo[[1]]$linkinv(eta1)))
	  if (!is.null(offset[[2]])) lres1 <- lres1 - offset[[2]]
	  R <- suppressWarnings(mchol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,
	       v=x$v,qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[2]])+crossprod(E[,jj[[2]]])))
	  Xty <- XWyd(x$Xd,rep(1,length(y)),lres1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[2]])
	  piv <- attr(R,"pivot")
	  rrank <- attr(R,"rank")
	  startji <- rep(0,ncol(R))
          if (rrank<ncol(R)) {
            R <- R[1:rrank,1:rrank]
	    piv <- piv[1:rrank]
          }
          startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
          start[jj[[2]]] <- startji
        } else { ## regular case
	  start <- rep(0,ncol(x))
	  x1 <- x[,jj[[1]],drop=FALSE]
          e1 <- E[,jj[[1]],drop=FALSE] ## square root of total penalty
          #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
          if (use.unscaled) {
            qrx <- qr(rbind(x1,e1))
            x1 <- rbind(x1,e1)
            startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
            startji[!is.finite(startji)] <- 0       
          } else startji <- pen.reg(x1,e1,yt1)
          start[jj[[1]]] <- startji
          lres1 <- log(abs(y-family$linfo[[1]]$linkinv(x[,jj[[1]],drop=FALSE]%*%start[jj[[1]]])))
	  if (!is.null(offset[[2]])) lres1 <- lres1 - offset[[2]]
          x1 <-  x[,jj[[2]],drop=FALSE];e1 <- E[,jj[[2]],drop=FALSE]
          #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
          if (use.unscaled) {
            x1 <- rbind(x1,e1)
            startji <- qr.coef(qr(x1),c(lres1,rep(0,nrow(E))))   
            startji[!is.finite(startji)] <- 0
          } else startji <- pen.reg(x1,e1,lres1)
          start[jj[[2]]] <- startji
	}  
      }
  }) ## initialize gaulss

  rd <- function(mu,wt,scale) {
  ## simulate responses 
    return( rnorm(nrow(mu), mu[ , 1], sqrt(scale/wt)/mu[ , 2]) )
  } ## gaulss rd


  structure(list(family="gaulss",ll=ll,link=paste(link),ncv=ncv,nlp=2,sandwich=sandwich,
    tri = trind.generator(2), ## symmetric indices for accessing derivative arrays
    initialize=initialize,postproc=postproc,residuals=residuals,
    linfo = stats,rd=rd, ## link information list
    d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
    ls=1, ## signals that ls not needed here
    available.derivs = 2, ## can use full Newton here
    discrete.ok = TRUE
    ),class = c("general.family","extended.family","family"))
} ## end gaulss


multinom <- function(K=1) {
## general family for multinomial logistic regression model...
## accepts no links as parameterization directly in terms of 
## linear predictor. 
## 1. get derivatives wrt mu, tau
## 2. get required link derivatives and tri indices.
## 3. transform derivs to derivs wrt eta (gamlss.etamu).
## 4. get the grad and Hessian etc for the model
##    via a call to gamlss.gH  
## the first derivatives of the log likelihood w.r.t
## the first and second parameters...
 
  if (K<1) stop("number of categories must be at least 2") 
   stats <- list()
  
  for (i in 1:K) {
    stats[[i]] <- make.link("identity")
    fam <- structure(list(link="identity",canonical="none",linkfun=stats[[i]]$linkfun,
           mu.eta=stats[[i]]$mu.eta),
           class="family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  } 

  residuals <- function(object,type=c("deviance")) {
  ## Deviance residuals where sign depends on whether classification correct (+ve)
  ## or not (-ve)...
      type <- match.arg(type)
      ## get category probabilities...
      p <- object$family$predict(object$family,eta=object$fitted.values)[[1]] ## changed from linear predictor for consistency
      ## now get most probable category for each observation
      pc <- apply(p,1,function(x) which(max(x)==x)[1])-1 
      n <- length(pc)
      ## +ve sign if class correct, -ve otherwise
      sgn <- rep(-1,n); sgn[pc==object$y] <- 1
      ## now get the deviance...
      sgn*sqrt(-2*log(pmax(.Machine$double.eps,p[1:n + object$y*n]))) 
  } ## multinom residuals

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL) {
  ## optional function to give predicted values - idea is that 
  ## predict.gam(...,type="response") will use this, and that
  ## either eta will be provided, or {X, beta, off, Vb}. family$data
  ## contains any family specific extra information. 
  ## if se = FALSE returns one item list containing matrix otherwise 
  ## list of two matrices "fit" and "se.fit"... 

    if (is.null(eta)) {
      discrete <- is.list(X)
      lpi <- attr(X,"lpi") 
      if (is.null(lpi)) {
        lpi <- list(1:ncol(X))
      } 
      K <- length(lpi) ## number of linear predictors
      nobs <- if (discrete) nrow(X$kd) else nrow(X)
      eta <- matrix(0,nobs,K)
      if (se) { 
        ve <- matrix(0,nobs,K) ## variance of eta
        ce <- matrix(0,nobs,K*(K-1)/2) ## covariance of eta_i eta_j
      }
      ii <- 0
      for (i in 1:K) {
        if (discrete) {
	  eta[,i] <- Xbd(X$Xd,beta,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[i]]) 
        } else {
          Xi <- X[,lpi[[i]],drop=FALSE]
          eta[,i] <- Xi%*%beta[lpi[[i]]] ## ith linear predictor
        }
        if (!is.null(off[[i]])) eta[,i] <- eta[,i] + off[[i]]
        if (se) { ## variance and covariances for kth l.p.
	  
          ve[,i] <- if (discrete) diagXVXd(X$Xd,Vb,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,nthreads=1,
	                   lt=X$lpid[[i]],rt=X$lpid[[i]]) else drop(pmax(0,rowSums((Xi%*%Vb[lpi[[i]],lpi[[i]]])*Xi)))
          ## ii <- 0 BUGGY location!
          if (i<K) for (j in (i+1):K) {
            ii <- ii + 1
            ce[,ii] <- if (discrete) diagXVXd(X$Xd,Vb,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,nthreads=1,
	                   lt=X$lpid[[i]],rt=X$lpid[[j]]) else drop(pmax(0,rowSums((Xi%*%Vb[lpi[[i]],lpi[[j]]])*X[,lpi[[j]]])))
          }
        }
      }
    } else { 
      se <- FALSE
    }
    gamma <- cbind(1,exp(eta))
    beta <- rowSums(gamma)
    gamma <- gamma/beta ## category probabilities
    vp <- gamma*0
    if (se) { ## need to loop to find se of probabilities...
      for (j in 1:(K+1)) {
        ## get dp_j/deta_k...
        if (j==1) dp <- -gamma[,-1,drop=FALSE]/beta else { 
          dp <- -gamma[,j]*gamma[,-1,drop=FALSE]
          dp[,j-1] <- gamma[,j]*(1-gamma[,j]) 
        }
        ## now compute variance... 
        vp[,j] <- rowSums(dp^2*ve)
        ii <- 0
        for (i in 1:K) if (i<K) for (k in (i+1):K) {
          ii <- ii + 1
          vp[,j] <- vp[,j] + 2 * dp[,i]*dp[,k]*ce[,ii] 
        }
        vp[,j] <- sqrt(pmax(0,vp[,j])) ## transform to se
      }
      return(list(fit=gamma,se.fit=vp))
    } ## if se
    list(fit=gamma)
  } ## multinom predict

  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    multinom <- list()
    object$y <- round(object$y) 
    multinom$nj <- tabulate(object$y+1) ## count each class membership
    multinom$n <- sum(multinom$nj) ## total number
    multinom$K <- length(multinom$nj)-1 ## number of linear predictors
    ## compute exp(eta) for categories 1..K
    multinom$gamma <- c(1,solve(diag(multinom$n/multinom$nj[-1],multinom$K)-
                        matrix(1,multinom$K,multinom$K),rep(1,multinom$K)))
    multinom$gamma <- log(multinom$gamma/sum(multinom$gamma))
    object$null.deviance <- -2*sum(multinom$gamma[object$y+1])
  })
  
  ncv <- function(X,y,wt,nei,beta,family,llf,H=NULL,Hi=NULL,R=NULL,offset=NULL,dH=NULL,db=NULL,deriv=FALSE,nt=1) {
    gamlss.ncv(X,y,wt,nei,beta,family,llf,H=H,Hi=Hi,R=R,offset=offset,dH=dH,db=db,deriv=deriv,nt=nt)
  } ## ncv  

  ll <- function(y,X,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL,eta=NULL,ncv=FALSE,sandwich=FALSE) {
  ## Function defining the logistic multimomial model log lik. 
  ## Assumption is that coding runs from 0..K, with 0 class having no l.p.
  ## ... this matches binary log reg case... 
  ## deriv: 0 - eval
  ##        1 - grad and Hess
  ##        2 - diagonal of first deriv of Hess
  ##        3 - first deriv of Hess
  ##        4 - everything.
    n <- length(y)
    jj <- attr(X,"lpi") ## extract linear predictor index
    if (is.null(eta)) {
      discrete <- is.list(X)
      ##return.l <- FALSE
      K <- length(jj) ## number of linear predictors 
      eta <- matrix(1,n,K+1) ## linear predictor matrix (dummy 1's in first column)
      if (is.null(offset)) offset <- list()
      offset[[K+1]] <- 0
      for (i in 1:K) if (is.null(offset[[i]])) offset[[i]] <- 0
      for (i in 1:K) eta[,i+1] <- offset[[i]] + if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,
                                  drop=X$drop,lt=X$lpid[[i]]) else X[,jj[[i]],drop=FALSE]%*%coef[jj[[i]]]
    } else { l2 <- 0;K <- ncol(eta);eta <- cbind(1,eta)} ##; return.l <- TRUE}
 
    if (K!=family$nlp) stop("number of linear predictors doesn't match")
    y <- round(y) ## just in case
    if (min(y)<0||max(y)>K) stop("response not in 0 to number of predictors + 1")
    
    ee <- exp(eta[,-1,drop=FALSE])
    beta <- 1 + rowSums(ee); alpha <- log(beta)
    
    l0 <- eta[1:n+y*n] - alpha ## log likelihood
    l <- sum(l0)    

    l1 <- matrix(0,n,K) ## first deriv matrix
 
    if (deriv>0) {
      for (i in 1:K) l1[,i] <- ee[,i]/beta ## alpha1
    
      ## the second derivatives...
    
      l2 <- matrix(0,n,K*(K+1)/2)
      ii <- 0; b2 <- beta^2
      for (i in 1:K) for (j in i:K) {
        ii <- ii + 1 ## column index
        l2[,ii] <- if (i==j) -l1[,i] + ee[,i]^2/b2 else (ee[,i]*ee[,j])/b2
      }

      ## finish first derivatives...
      for (i in 1:K) l1[,i] <- as.numeric(y==i) - l1[,i] 

    } ## if (deriv>0)
 
    l3 <- l4 <- 0 ## defaults
    tri <- family$tri ## indices to facilitate access to earlier results
    
    if (deriv>1) { ## the third derivatives...
      l3 <- matrix(0,n,(K*(K+3)+2)*K/6)
      ii <- 0; b3 <- b2 * beta
      for (i in 1:K) for (j in i:K) for (k in j:K) {
        ii <- ii + 1 ## column index
        if (i==j&&j==k) { ## all same
           l3[,ii] <- l2[,tri$i2[i,i]] + 2*ee[,i]^2/b2 - 2*ee[,i]^3/b3
        } else if (i!=j&&j!=k&i!=k) { ## all different
           l3[,ii] <- -2*(ee[,i]*ee[,j]*ee[,k])/b3
        } else { ## two same one different
           kk <- if (i==j) k else j ## get indices for differing pair
           l3[,ii] <- l2[,tri$i2[i,kk]] - 2*(ee[,i]*ee[,j]*ee[,k])/b3
        }
      }
    } ## if (deriv>1)

    if (deriv>3) { ## the fourth derivatives...
      l4 <- matrix(0,n,(6+K*11+K^2*6+K^3)*K/24)
      ii <- 0; b4 <- b3 * beta
      for (i in 1:K) for (j in i:K) for (k in j:K) for (l in k:K) {
        ii <- ii + 1 ## column index
        uni <- unique(c(i,j,k,l));
        nun <- length(uni) ## number of unique indices
        if (nun==1) { ## all equal
          l4[,ii] <- l3[,tri$i3[i,i,i]] + 4*ee[,i]^2/b2 - 10*ee[,i]^3/b3 + 6*ee[,i]^4/b4
        } else if (nun==4) { ## all unequal
          l4[,ii] <- 6*ee[,i]*ee[,j]*ee[,k]*ee[,l]/b4
        } else if (nun==3) { ## 2 same 2 different
          l4[,ii] <- l3[,tri$i3[uni[1],uni[2],uni[3]]]  +6*ee[,i]*ee[,j]*ee[,k]*ee[,l]/b4
        } else if (sum(uni[1]==c(i,j,k,l))==2) { ## 2 unique (2 of each)
          l4[,ii] <- l3[,tri$i3[uni[1],uni[2],uni[2]]] - 2 * ee[,uni[1]]^2*ee[,uni[2]]/b3 + 
                     6*ee[,i]*ee[,j]*ee[,k]*ee[,l]/b4
        } else { ## 3 of one 1 of the other
          if (sum(uni[1]==c(i,j,k,l))==1) uni <- uni[2:1] ## first index is triple repeat index
          l4[,ii] <- l3[,tri$i3[uni[1],uni[1],uni[2]]] - 4 * ee[,uni[1]]^2*ee[,uni[2]]/b3 + 
                     6*ee[,i]*ee[,j]*ee[,k]*ee[,l]/b4
        }
      }
    } ## if deriv>3

    ##if (return.l) return(list(l=l0,l1=l1,l2=l2,l3=l3,l4=l4)) ## for testing...

    if (deriv) {
      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,l1,l2,tri$i2,l3=l3,i3=tri$i3,l4=l4,i4=tri$i4,
                      d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D,sandwich=sandwich)
      if (ncv) { ret$l1=l1; ret$l2=l2; ret$l3=l3 }		      
    } else ret <- list()
    ret$l <- l; ret
  } ## end ll multinom

  sandwich <- function(y,X,coef,wt,family,offset=NULL) {
  ## compute filling for sandwich estimate of cov matrix
    ll(y,X,coef,wt,family,offset=NULL,deriv=1,sandwich=TRUE)$lbb
  }

  rd <- function(mu,wt,scale) {
    ## simulate data given fitted linear predictor matrix in mu 
    p <- exp(cbind(0,mu))
    p <- p/rowSums(p)
    cp <- t(apply(p,1,cumsum))
    apply(cp,1,function(x) min(which(x>runif(1))))-1    
  } ## rd

  initialize <- expression({ ## for multinom
  ## Binarize each category and lm on 6*y-3 by category.
 
      n <- rep(1, nobs)
      ## should E be used unscaled or not?..
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      if (is.null(start)) {
        jj <- attr(x,"lpi")
	if (is.list(x)) { ## discrete case
          start <- rep(0,max(unlist(jj)))
          for (k in 1:length(jj)) { ## loop over the linear predictors
            yt1 <- 6*as.numeric(y==k)-3
	    R <- suppressWarnings(mchol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,
	         v=x$v,qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[k]])+crossprod(E[,jj[[k]]])))
	    Xty <- XWyd(x$Xd,rep(1,length(y)),yt1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[k]])
	    piv <- attr(R,"pivot")
	    rrank <- attr(R,"rank")
	    startji <- rep(0,ncol(R))
            if (rrank<ncol(R)) {
              R <- R[1:rrank,1:rrank]
	      piv <- piv[1:rrank]
            }
            startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
	    startji[!is.finite(startji)] <- 0
	    start[jj[[k]]] <- startji
            
          } ## lp loop
        } else { ## regular case
          start <- rep(0,ncol(x))
          for (k in 1:length(jj)) { ## loop over the linear predictors      
            yt1 <- 6*as.numeric(y==k)-3
            x1 <- x[,jj[[k]],drop=FALSE]
            e1 <- E[,jj[[k]],drop=FALSE] ## square root of total penalty
            if (use.unscaled) {
              qrx <- qr(rbind(x1,e1))
              x1 <- rbind(x1,e1)
              startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
              startji[!is.finite(startji)] <- 0       
            } else startji <- pen.reg(x1,e1,yt1)
            start[jj[[k]]] <- startji ## copy coefficients back into overall start coef vector
          } ## lp loop
        }
      }		
  }) ## initialize multinom

  structure(list(family="multinom",ll=ll,link=NULL,#paste(link),
    nlp=round(K),rd=rd,ncv=ncv,sandwich=sandwich,
    tri = trind.generator(K), ## symmetric indices for accessing derivative arrays
    initialize=initialize,postproc=postproc,residuals=residuals,predict=predict,
    linfo = stats, ## link information list
    d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
    ls=1, ## signals that ls not needed here
    available.derivs = 2, ## can use full Newton here
    discrete.ok = TRUE
    ),class = c("general.family","extended.family","family"))
} ## end multinom



pen.reg <- function(x,e,y) {
## get coefficients of penalized regression of y on matrix x
## where e is a square root penalty. Idea is to use e mainly for 
## regularization, so that edf is close to rank of x.  
  if (sum(abs(e))==0) { ## no penalization - easy
    b <- qr.coef(qr(x),y);b[!is.finite(b)] <- 0
    return(b)
  }
  ## need to adjust degree of penalization, so best to QR
  ## the x matrix up front...
  qrx <- qr(x,LAPACK=TRUE)
  R <- qr.R(qrx)
  r <- ncol(R)
  rr <- Rrank(R) ## rank of R/X
  R[,qrx$pivot] <- R ## unpivot
  Qy <- qr.qty(qrx,y)[1:ncol(R)]  
  ## now we want estimates with penalty weight low enough 
  ## EDF is k * rr where k is somewhere in e.g. (.7,.9)
  k <- .01 * norm(R)/norm(e)
  qrr <- qr(rbind(R,e*k));
  edf <- sum(qr.Q(qrr)[1:r,]^2)
  ## compute rank of e less rank of the space penalized by e not in
  ## range space of x, this is how much penalty can in principle change
  ## edf by. Needed for corner cases where e.g. penalty is imposed for
  ## identifiability reasons and then only penalizes null space of x...
  re <- min(sum(colSums(abs(e))!=0),nrow(e)) - Rrank(qr.R(qrr)) + rr
  while (edf > rr-.1*re) { ## increase penalization
    k <- k*10
    qrr <- qr(rbind(R,e*k));
    edf <- sum(qr.Q(qrr)[1:r,]^2)
  } 
  while (edf<.85*rr) { ## reduce penalization (was .7)
    k <- k/5 ## was 20! 
    qrr <- qr(rbind(R,e*k));
    edf <- sum(qr.Q(qrr)[1:r,]^2)
  } 
  b <- qr.coef(qrr,c(Qy,rep(0,nrow(e))));b[!is.finite(b)] <- 0
  b
} ## pen.reg

## code for zero inflated Poisson models


#log1ex <- function(x) {
## evaluate log(1+exp(x)) accurately and avoiding overflow
#  y <- x
#  big <- -log(.Machine$double.eps)+5 ## exp(big) overwhelms 1
#  ind <- x > big
#  y[ind] <- x[ind] ## log(1+exp(x)) = x to machine precision
#  ## compute number below which log(1+exp(x)) = exp(x) to
#  ## machine precision... 
#  small <- log(sqrt(.Machine$double.eps))
#  ind1 <- x < small
#  y[ind1] <- exp(x[ind1])
#  ind <- !ind&!ind1 ## the moderate size elements 
#  y[ind] <- log(1+exp(x[ind]))
#  y 
#}

#logist <- function(x) {
## overflow proof logistic
#  ind <- x > 0; y <- x
#  y[ind] <- 1/(exp(-x[ind])+1) 
#  ex <- exp(x[!ind])
#  y[!ind] <- ex/(1+ex)
#  y
#}

l1ee <- function(x) {
## log(1-exp(-exp(x)))...
  ind <- x < log(.Machine$double.eps)/3
  ex <- exp(x);exi <- ex[ind]
  l <- log(1-exp(-ex))
  l[ind] <- log(exi-exi^2/2+exi^3/6)
  ind <- x < -log(.Machine$double.xmax)
  l[ind] <- x[ind]
  l
}

lee1 <- function(x) {
## log(exp(exp(x))-1)...
  ind <- x < log(.Machine$double.eps)/3
  ex <- exp(x);exi <- ex[ind]
  l <- log(exp(ex)-1)
  l[ind] <- log(exi+exi^2/2+exi^3/6)
  ind <- x < -log(.Machine$double.xmax)
  l[ind] <- x[ind]
  ind <- x > log(log(.Machine$double.xmax))
  l[ind] <- ex[ind]
  l
}

ldg <- function(g,deriv=4) {
  alpha <- function(g) {
    ind <- g > log(.Machine$double.eps)/3
    eg <- exp(g)
    g[ind] <- eg[ind]/(1-exp(-eg[ind]))
    g[!ind] <- 1+eg[!ind]/2 + eg[!ind]^2/12
    g
  }
  ind <- g < log(.Machine$double.eps)/3
  ghi <- log(log(.Machine$double.xmax)) + 1 
  ## ... above ghi alpha(g) is simply exp(g) 
  ii <- g>ghi 
  a <- alpha(g)
  eg <- exp(g)
  l2 <- a*(a-eg-1)
  egi <- eg[ind]
  ## in the lower tail alpha = 1 + b, where b = eg/2 + eg^2/12
  ## so l'' = alpha*(b-eg)...
  b <- egi*(1+egi/6)/2
  l2[ind] <- a[ind]*(b-egi)
  l2[ii] <- -exp(g[ii])
  l3 <- l4 <- NULL
  ## in a similar vein l3 can be robustified...  
  if (deriv>1) {
    l3 <- a*(a*(-2*a + 3*(eg+1)) - 3*eg - eg^2 - 1) 
    l3[ind] <- a[ind]*(-b-2*b^2+3*b*egi-egi^2)
    l3[ii] <- -exp(g[ii])
  }
  ## finally l4, which requires a similar approach...
  if (deriv>2) {
    l4 <- a*(6*a^3 - 12*(eg+1)*a^2+4*eg*a+7*(eg+1)^2*a-(4+3*eg)*eg -(eg+1)^3)
    l4[ind] <- a[ind]*(6*b*(3+3*b+b^2) - 12*egi*(1+2*b+b^2) - 12*b*(2-b) + 4*egi*(1+b)+
                     7*(egi^2+2*egi+b*egi^2+2*b*egi+b)-(4+3*egi)*egi-egi*(3+3*egi+egi^2))
   
    l4[ii] <- -exp(g[ii])
  }
  l1=-a
  ghi <- log(.Machine$double.xmax)/5
  ii <- g > ghi
  if (sum(ii)) {
    l1[ii] <- l2[ii] <- l3[ii] <- l4[ii] <- -exp(ghi)
  }
  list(l1=l1,l2=l2,l3=l3,l4=l4)
} ## ldg

lde <- function(eta,deriv=4) {
  ## llog lik derivs w.r.t. eta
  ind <- eta < log(.Machine$double.eps)/3
  ii <- eta > log(.Machine$double.xmax)
  l1 <- et <- exp(eta);eti <- et[ind]
  l1[!ind] <- et[!ind]/(exp(et[!ind])-1)
  b <- -eti*(1+eti/6)/2
  l1[ind] <- 1+b
  l1[ii] <- 0
  ## l2 ...   
  l2 <- l1*((1-et)-l1)
  l2[ind] <- -b*(1+eti+b) - eti
  l2[ii] <- 0
  l3 <- l4 <- NULL
  ## l3 ...
  if (deriv>1) {
    ii <- eta > log(.Machine$double.xmax)/2
    l3 <- l1*((1-et)^2-et - 3*(1-et)*l1 + 2*l1^2)
    l3[ind] <- l1[ind]*(-3*eti+eti^2 -3*(-eti+b-eti*b) + 2*b*(2+b))
    l3[ii] <- 0
  }
  ## l4 ...
  if (deriv>2) {
    ii <- eta > log(.Machine$double.xmax)/3
    l4 <- l1*((3*et-4)*et + 4*et*l1 + (1-et)^3 - 7*(1-et)^2*l1 + 12*(1-et)*l1^2 
            - 6*l1^3)
    l4[ii] <- 0
    l4[ind] <- l1[ind]*(4*l1[ind]*eti - eti^3 - b -7*b*eti^2 - eti^2 - 5*eti -
                 10*b*eti - 12*eti*b^2 - 6*b^2 - 6*b^3)
  }
  list(l1=l1,l2=l2,l3=l3,l4=l4)
} ## lde


zipll <- function(y,g,eta,deriv=0) {
## function to evaluate zero inflated Poisson log likelihood
## and its derivatives w.r.t. g/gamma and eta where 
## 1-p = exp(-exp(eta)) and lambda = exp(gamma), for each datum in vector y.
## p is probability of potential presence. lambda is Poisson mean
## given potential presence. 
## deriv: 0 - eval
##        1 - grad (l,p) and Hess (ll,lp,pp)
##        2 - third derivs lll,llp,lpp,ppp
##        4 - 4th derivs. llll,lllp,llpp,lppp,pppp

   l1 <- El2 <- l2 <- l3 <- l4 <- NULL
   zind <- y == 0 ## the index of the zeroes
   ## yz <- y[zind];
   yp <- y[!zind]
   l <- et <- exp(eta)
   l[zind] <- -et[zind] # -exp(eta[ind])
   l[!zind] <- l1ee(eta[!zind]) + yp*g[!zind] - lee1(g[!zind]) - lgamma(yp+1)
   p <- 1-exp(-et) ## probablity of non-zero

   if (deriv>0) { ## get first and second derivs...
     n <- length(y)
     l1 <- matrix(0,n,2)
     le <- lde(eta,deriv) ## derivs of ll wrt eta     
     lg <- ldg(g,deriv) ## derivs of ll wrt gamma
     l1[!zind,1] <- yp + lg$l1[!zind]  ## l_gamma, y>0
     l1[zind,2] <- l[zind] ## l_eta, y==0
     l1[!zind,2] <- le$l1[!zind]  ## l_eta, y>0    

     El2 <- l2 <- matrix(0,n,3)
     ## order gg, ge, ee... 
     l2[!zind,1] <- lg$l2[!zind]   ## l_gg, y>0
     l2[!zind,3] <- le$l2[!zind]   ## l_ee, y>0
     l2[zind,3]  <- l[zind] ## l_ee, y=0
     El2[,1] <- p*lg$l2             ## E(l_gg)
     El2[,3] <- -(1-p)*et + p*le$l2 ## E(l_ee)
   }
   if (deriv>1) {
      ## the third derivatives
      ## order ggg,gge,gee,eee
      l3 <- matrix(0,n,4) 
      l3[!zind,1] <- lg$l3[!zind]   ## l_ggg, y>0
      l3[!zind,4] <- le$l3[!zind]   ## l_eee, y>0
      l3[zind,4]  <- l[zind]        ## l_eee, y=0
   }
   if (deriv>3) {
      ## the fourth derivatives
      ## order gggg,ggge,ggee,geee,eeee
      l4 <- matrix(0,n,5) 
      l4[!zind,1] <- lg$l4[!zind]   ## l_gggg, y>0
      l4[!zind,5] <- le$l4[!zind]   ## l_eeee, y>0
      l4[zind,5]  <- l[zind]        ## l_eeee, y=0
   }
   list(l=l,l1=l1,l2=l2,l3=l3,l4=l4,El2=El2)
} ## zipll


ziplss <-  function(link=list("identity","identity")) {
## Extended family for Zero Inflated Poisson fitted as gamlss 
## type model.
## mu1 is Poisson mean, while mu2 is zero inflation parameter.
  ## first deal with links and their derivatives...
  if (length(link)!=2) stop("ziplss requires 2 links specified as character strings")
  okLinks <- list(c("identity"),c("identity"))
  stats <- list()
  param.names <- c("Poisson mean","binary probability")
  for (i in 1:2) {
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- make.link(link[[i]]) else 
    stop(link[[i]]," link not available for ",param.names[i]," parameter of ziplss")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
           mu.eta=stats[[i]]$mu.eta),
           class="family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  } 

  residuals <- function(object,type=c("deviance","response")) {
      ls <- function(y) {
        ## compute saturated likelihood for ziplss model 
        l <- y;l[y<2] <- 0
        ind <- y > 1 & y < 18
        ## lambda maximizing likelihood for y = 2 to 17 
        glo <- c(1.593624,2.821439,3.920690,4.965114,5.984901,6.993576,
                 7.997309,8.998888,9.999546,10.999816,11.999926,12.999971,
                 13.999988,14.999995,15.999998,16.999999)
        g <- y ## maximizing lambda essentially y above this
        g[ind] <- glo[y[ind]-1]
        ind <- y > 1
        l[ind] <- zipll(y[ind],log(g[ind]),g[ind]*0+1e10,deriv=0)$l
        l
      } ## ls

      type <- match.arg(type)
      p <-  1 - exp(-exp(object$fitted[,2]));
      lam <- exp(object$fitted[,1])
      ind <- lam > .Machine$double.eps^.5
      ## compute E(y)
      Ey <- p ## very small lambda causes conditional expectation to be 1
      Ey[ind] <- p[ind]*lam[ind]/(1-exp(-lam[ind])) 
      rsd <- object$y - Ey ## raw residuals
      if (type=="response") return(rsd)
      else { ## compute deviance residuals
        sgn <- sign(rsd)
        ind <- object$y == 0
        rsd <- pmax(0,2*(ls(object$y) - zipll(object$y,object$fitted[,1],object$fitted[,2],deriv=0)$l))
        rsd <- sqrt(rsd)*sgn
      }
      rsd
  } ## ziplss residuals

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL) {
  ## optional function to give predicted values - idea is that 
  ## predict.gam(...,type="response") will use this, and that
  ## either eta will be provided, or {X, beta, off, Vb}. family$data
  ## contains any family specific extra information. 
  ## if se = FALSE returns one item list containing matrix otherwise 
  ## list of two matrices "fit" and "se.fit"... 

    if (is.null(eta)) {
      if (is.null(off)) off <- list(0,0)
      off[[3]] <- 0
      for (i in 1:2) if (is.null(off[[i]])) off[[i]] <- 0
      lpi <- attr(X,"lpi") 
      X1 <- X[,lpi[[1]],drop=FALSE]
      X2 <- X[,lpi[[2]],drop=FALSE]
      gamma <- drop(X1%*%beta[lpi[[1]]] + off[[1]]) ## linear predictor for poisson parameter 
      eta <- drop(X2%*%beta[lpi[[2]]] + off[[2]])  ## linear predictor for presence parameter 
      if (se) {
        v.g <- drop(pmax(0,rowSums((X1%*%Vb[lpi[[1]],lpi[[1]]])*X1))) ## var of gamma
        v.e <- drop(pmax(0,rowSums((X1%*%Vb[lpi[[1]],lpi[[1]]])*X1))) ## var of eta
        v.eg <- drop(pmax(0,rowSums((X1%*%Vb[lpi[[1]],lpi[[2]]])*X2))) ## cov of eta, gamma
      }
    } else { 
      se <- FALSE
      gamma <- eta[,1]
      eta <- eta[,2]
    }
    et <- exp(eta)
    mu <- p <- 1 - exp(-et)
    fv <- lambda <- exp(gamma)  
    ind <- gamma < log(.Machine$double.eps)/2
    mu[!ind] <- lambda[!ind]/(1-exp(-lambda[!ind]))
    mu[ind] <- 1
    fv <- list(p*mu)    ## E(y)    
    if (!se) return(fv) else {
      df.de <- p  
      ind <- eta < log(.Machine$double.xmax)/2
      df.de[!ind] <- 0
      df.de[ind] <- exp(-et[ind])*et[ind]
      df.de <- df.de * mu
      df.dg <- ((lambda + 1)*mu - mu^2)*p
      fv[[2]] <- sqrt(df.dg^2*v.g + df.de^2*v.e + 2 * df.de * df.dg * v.eg)   
      names(fv) <- c("fit","se.fit")
      return(fv)
    }
  } ## ziplss predict


  rd <- function(mu,wt,scale) {
    ## simulate data given fitted latent variable in mu 
      rzip <- function(gamma,eta) { ## generate ziP deviates according to model and lp gamma
        y <- gamma; n <- length(y)
        lambda <- exp(gamma)
        p <- 1- exp(-exp(eta)) ## prob present
        ind <- p > runif(n)
        y[!ind] <- 0
        np <- sum(ind)
        ## generate from zero truncated Poisson, given presence...
	u <- runif(np,dpois(0,lambda[ind]),1)
	## qpois can produce infinite answers at low lambda
	## if u too close to one, and it can get very close at low
	## lambda!
	one.eps <- 1 - .Machine$double.eps^.75
	u[u>one.eps] <- one.eps 
        y[ind] <- qpois(u,lambda[ind])
        y
      } 
      rzip(mu[,1],mu[,2])
  } ## rd

  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    ## null model really has two parameters... probably need to newton iterate
    ls <- function(y) {
        ## compute saturated likelihood for ziplss model 
        l <- y;l[y<2] <- 0
        ind <- y > 1 & y < 18
        ## lambda maximizing likelihood for y = 2 to 17 
        glo <- c(1.593624,2.821439,3.920690,4.965114,5.984901,6.993576,
                 7.997309,8.998888,9.999546,10.999816,11.999926,12.999971,
                 13.999988,14.999995,15.999998,16.999999)
        g <- y ## maximizing lambda essentially y above this
        g[ind] <- glo[y[ind]-1]
        ind <- y > 1
        l[ind] <- zipll(y[ind],log(g[ind]),g[ind]*0+1e10,deriv=0)$l
        l
    } ## ls

    fp <- function(p,y) {
    ## compute zero related part of log likelihood
       eps <- .Machine$double.eps^.5 
       l1p <- if (p>eps) log(1-p) else -p - p^2/2 
       l1p*sum(y==0) + log(p)*sum(y>0)
    } ## fp

    flam <- function(lam,y) {
      ## compute >0 part of log likelihood
      y <- y[y>0]
      sum(y*log(lam) - log(exp(lam)-1) - lgamma(y+1))
    } ## flam

    ## optimize zero repated part of likelihood w.r.t. p...
    lnull <- optimize(fp,interval=c(1e-60,1-1e-10),y=object$y,maximum=TRUE)$objective
    ## optimize >0 part for lambda...
    my <- mean(object$y[object$y>0])
    lnull <- lnull + optimize(flam,interval=c(my/2,my*2),y=object$y,maximum=TRUE)$objective
    object$null.deviance <- 2*(sum(ls(object$y)) - lnull)
   
  }) ## postproc
  
  ncv <- function(X,y,wt,nei,beta,family,llf,H=NULL,Hi=NULL,R=NULL,offset=NULL,dH=NULL,db=NULL,deriv=FALSE,nt=1) {
    gamlss.ncv(X,y,wt,nei,beta,family,llf,H=H,Hi=Hi,R=R,offset=offset,dH=dH,db=db,deriv=deriv,nt=nt)
  } ## ncv  


  ll <- function(y,X,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL,eta=NULL,ncv=FALSE,sandwich=FALSE) {
  ## function defining the gamlss ZIP model log lik. 
  ## First l.p. defines Poisson mean, given presence (lambda)
  ## Second l.p. defines probability of presence (p)
  ## deriv: 0 - eval
  ##        1 - grad and Hess
  ##        2 - diagonal of first deriv of Hess
  ##        3 - first deriv of Hess
  ##        4 - everything.
    if (is.null(offset)) offset <- list(0,0) else offset[[3]] <- 0
    for (i in 1:2) if (is.null(offset[[i]])) offset[[i]] <- 0
    jj <- attr(X,"lpi") ## extract linear predictor index
    if (is.null(eta)) {
      eta <- X[,jj[[1]],drop=FALSE]%*%coef[jj[[1]]] + offset[[1]]
      eta1 <- X[,jj[[2]],drop=FALSE]%*%coef[jj[[2]]] +offset[[2]]
    } else { ## eta supplied
      eta1 <- eta[,2]
      eta <- eta[,1]
    }
    lambda <- family$linfo[[1]]$linkinv(eta)
    p <-  family$linfo[[2]]$linkinv(eta1) 
    
    ##n <- length(y)
    ## l1 <- matrix(0,n,2)
    zl <- zipll(y,lambda,p,deriv)

    if (deriv>0) {
      ## need some link derivatives for derivative transform
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(lambda),family$linfo[[2]]$d2link(p))
    }
 
    ## l3 <- l4 <- 
    g3 <- g4 <- 0 ## defaults

    if (deriv>1) {
      ## the third derivatives
      ## order lll,llp,lpp,ppp
      g3 <- cbind(family$linfo[[1]]$d3link(lambda),family$linfo[[2]]$d3link(p))
    }

    if (deriv>3) {
      ## the fourth derivatives
      ## order llll,lllp,llpp,lppp,pppp
      g4 <- cbind(family$linfo[[1]]$d4link(lambda),family$linfo[[2]]$d4link(p))
    }
    if (deriv) {
      i2 <- family$tri$i2; i3 <- family$tri$i3
      i4 <- family$tri$i4
   
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- gamlss.etamu(zl$l1,zl$l2,zl$l3,zl$l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)

      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                      d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D,sandwich=sandwich)
      if (ncv) {
        ret$l1 <- de$l1; ret$l2 = de$l2; ret$l3 = de$l3
      }		      
    } else ret <- list()
    ret$l <- sum(zl$l); ret
  } ## end ll for ZIP

  sandwich <- function(y,X,coef,wt,family,offset=NULL) {
  ## compute filling for sandwich estimate of cov matrix
    ll(y,X,coef,wt,family,offset=NULL,deriv=1,sandwich=TRUE)$lbb
  }

  initialize <- expression({ ## for ZIP
  ## Idea is to regress binarized y on model matrix for p. 
  ## Then downweight any y=0 with p<0.5 and regress g(y) on 
  ## the model matrix for lambda - don't drop as this may
  ## induce rank deficiency in model matrix! 
  ## May be called in both gam.fit5 and initial.spg... 
  ## note that appropriate E scaling
  ## for full calculation may be inappropriate for initialization 
  ## which is basically penalizing something different here.
  ## best we can do here is to use E only as a regularizer.
      n <- rep(1, nobs)
      if (all.equal(y,round(y))!=TRUE) {
          stop("Non-integer response variables are not allowed with ziplss ")
      }
      if ((min(y)==0&&max(y)==1)) stop("Using ziplss for binary data makes no sense")
      ## should E be used unscaled or not?..
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      if (is.null(start)) {
        jj <- attr(x,"lpi")
        start <- rep(0,ncol(x))
        x1 <- x[,jj[[2]],drop=FALSE]
        e1 <- E[,jj[[2]],drop=FALSE] ## square root of total penalty
        yt1 <- as.numeric(as.logical(y)) ## binarized response
        if (use.unscaled) {
          qrx <- qr(rbind(x1,e1))
          x1 <- rbind(x1,e1)
          startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
          startji[!is.finite(startji)] <- 0       
        } else startji <- pen.reg(x1,e1,yt1)
        start[jj[[2]]] <- startji
        p <- drop(x1[1:nobs,,drop=FALSE] %*% startji) ## probability of presence
        ind <- y==0 & p < 0.5 ## downweight these for estimating lambda
        w <- rep(1,nobs); w[ind] <- .1

        ## note assumption that working scale is log...
        yt1 <-  family$linfo[[1]]$linkfun(log(abs(y)+(y==0)*.2))

        yt1 <- yt1*w        

        x1 <-  w*x[,jj[[1]],drop=FALSE];e1 <- E[,jj[[1]],drop=FALSE]
        if (use.unscaled) {
          x1 <- rbind(x1,e1)
          startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))   
          startji[!is.finite(startji)] <- 0
        } else startji <- pen.reg(x1,e1,yt1)
        start[jj[[1]]] <- startji
      }
  }) ## initialize ziplss

  structure(list(family="ziplss",ll=ll,link=paste(link),nlp=2,ncv=ncv,sandwich=sandwich,
    tri = trind.generator(2), ## symmetric indices for accessing derivative arrays
    initialize=initialize,postproc=postproc,residuals=residuals,rd=rd,predict=predict,
    linfo = stats, ## link information list
    d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
    ls=1, ## signals that ls not needed here
    available.derivs = 2 ## can use full Newton here
    ),class = c("general.family","extended.family","family"))
} ## ziplss





gevlss <- function(link=list("identity","identity","logit")) {
## General family for GEV location scale model...
## so mu is mu1, rho = log(sigma) is mu2 and xi is mu3
## 1. get derivatives wrt mu, rho and xi.
## 2. get required link derivatives and tri indices.
## 3. transform derivs to derivs wrt eta (gamlss.etamu).
## 4. get the grad and Hessian etc for the model
##    via a call to gamlss.gH  
  
  ## first deal with links and their derivatives...
  if (length(link)!=3) stop("gevlss requires 3 links specified as character strings")
  okLinks <- list(c("log", "identity"),"identity",c("identity","logit"))
  stats <- list()
  for (i in 1:3) {
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- make.link(link[[i]]) else 
    stop(link[[i]]," link not available for gevlss")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
           mu.eta=stats[[i]]$mu.eta),
           class="family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }
  
  if (link[[3]]=="logit") { ## shifted logit link to confine xi to (-1,.5)
    ## Smith '85 Biometrika shows that -1 limit needed for MLE consistency
    ## but would need -0.5 for normality...
    stats[[3]]$linkfun <- function(mu) binomial()$linkfun((mu + 1)/1.5)
    stats[[3]]$mu.eta <- function(eta) binomial()$mu.eta(eta)*1.5
    stats[[3]]$linkinv <- function(eta) 1.5* binomial()$linkinv(eta) - 1
    stats[[3]]$d2link <- function(mu) { mu <- (mu+ 1)/1.5; (1/(1 - mu)^2 - 1/mu^2)/1.5^2}
    stats[[3]]$d3link <- function(mu) { mu <- (mu+ 1)/1.5; (2/(1 - mu)^3 + 2/mu^3)/1.5^3}
    stats[[3]]$d4link <- function(mu) {  mu <- (mu+ 1)/1.5; (6/(1-mu)^4 - 6/mu^4)/1.5^4}
  }

  residuals <- function(object,type=c("deviance","pearson","response")) {
      mu <- object$fitted[,1]
      rho <- object$fitted[,2]
      xi <- object$fitted[,3]
      y <- object$y
      fv <- mu + exp(rho)*(gamma(1-xi)-1)/xi 
      eps <- 1e-7; xi[xi>=0&xi<eps] <- eps; xi[xi<0&xi>-eps] <- -eps
      type <- match.arg(type)
      if (type=="deviance") {
        rsd <- (xi+1)/xi * log(1+(y-mu)*exp(-rho)*xi) + (1+(y-mu)*exp(-rho)*xi)^(-1/xi) +
	       (1+xi)*log(1+xi) - (1 + xi) ## saturated part
        rsd <- sqrt(pmax(0,rsd))*sign(y-fv)
      } else if (type=="pearson") {
        sd <- exp(rho)/xi*sqrt(pmax(0,gamma(1-2*xi)-gamma(1-xi)^2))
        rsd <- (y-fv)/sd
      } else {
        rsd <- y-fv
      }
      rsd
    } ## gevlss residuals
    
  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    ## It's difficult to define a sensible version of this that ensures
    ## that the data fall in the support of the null model, whilst being
    ## somehow equivalent to the full fit
    object$null.deviance <- NA
    
  })
  
  ncv <- function(X,y,wt,nei,beta,family,llf,H=NULL,Hi=NULL,R=NULL,offset=NULL,dH=NULL,db=NULL,deriv=FALSE,nt=1) {
    gamlss.ncv(X,y,wt,nei,beta,family,llf,H=H,Hi=Hi,R=R,offset=offset,dH=dH,db=db,deriv=deriv,nt=nt)
  } ## ncv  

  ll <- function(y,X,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL,eta=NULL,ncv=FALSE,sandwich=FALSE) {
  ## function defining the gamlss GEV model log lik. 
  ## deriv: 0 - eval
  ##        1 - grad and Hess
  ##        2 - diagonal of first deriv of Hess
  ##        3 - first deriv of Hess
  ##        4 - everything.
    if (!is.null(offset)) offset[[4]] <- 0
    discrete <- is.list(X)
    jj <- attr(X,"lpi") ## extract linear predictor index
    if (is.null(eta)) {
      eta <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[1]]) else X[,jj[[1]],drop=FALSE]%*%coef[jj[[1]]]
      if (!is.null(offset[[1]])) eta <- eta + offset[[1]] 
      etar <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[2]]) else X[,jj[[2]],drop=FALSE]%*%coef[jj[[2]]] ## log sigma
      if (!is.null(offset[[2]])) etar <- etar + offset[[2]]
      etax <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[3]]) else X[,jj[[3]],drop=FALSE]%*%coef[jj[[3]]] ## shape parameter
      if (!is.null(offset[[3]])) etax <- etax + offset[[3]]
    } else { ## eta supplied
      etar <- eta[,2]
      etax <- eta[,3]
      eta <- eta[,1]
    }
    mu <- family$linfo[[1]]$linkinv(eta) ## mean
    rho <- family$linfo[[2]]$linkinv(etar) ## log sigma
    xi <- family$linfo[[3]]$linkinv(etax) ## shape parameter
    
    ## Avoid xi == 0 - using a separate branch for xi==0 requires
    ## seperate treatment of derivative w.r.t. xi, and statistically
    ## brings us nothing. 
    eps <- 1e-7 
 
    xi[xi>=0&xi<eps] <- eps
    xi[xi<0&xi>-eps] <- -eps

    n <- length(y)
    l1 <- matrix(0,n,3)
    
    ## note that the derivative code is largely auto-generated, and
    ## auto-simplified. Optimized Maxima derivs exported as Maxima
    ## code, translated to R code in R, then processed in R to
    ## remove redundant auxiliary variables and their definitions.
    ## Modifications of auto code (but not the consequent substitutions) are
    ## flagged '## added'. Code post auto and non-auto modification has
    ## been tested against raw translated code. 

    exp1 <- exp(1); ## facilitates lazy auto-translation
    ymu <- y - mu
    aa0 <- (xi*ymu)/exp1^rho # added
    ind <- which(aa0 <= -1) ## added
    if (FALSE&&length(ind)>0) { ## all added
      xii <- xi[ind] ## this idea is really not a good one - messes up derivatives when triggered
      erho <- exp1^rho[ind]
      eps1 <- 1-.Machine$double.eps^.25
      ymu[ind] <- -erho/xii*eps1
      aa0[ind] <- -eps1
    }
    log.aa1 <- log1p(aa0) ## added
    aa1 <- aa0 + 1 # (xi*(y-mu))/exp1^rho+1;
    aa2 <- 1/xi;
    l0  <- (-aa2*(1+xi)*log.aa1)-1/aa1^aa2-rho;
    l <- sum(l0)

    if (deriv>0) {
      ## first derivatives m, r, x...
      bb1 <- 1/exp1^rho;
      bb2 <- bb1*xi*ymu+1;
     
      l1[,1]  <-  (bb1*(xi+1))/bb2-bb1*bb2^((-1/xi)-1);
      cc2 <- ymu;
      cc0 <- bb1*xi*cc2 ## added
      
      log.cc3 <- log1p(cc0) ## added
      cc3 <- cc0 + 1 ##bb1*xi*cc2+1;
      l1[,2]  <-  (-bb1*cc2*cc3^((-1/xi)-1))+(bb1*(xi+1)*cc2)/cc3-1;
      dd3 <- xi+1;
      dd6 <- 1/cc3;
      dd7 <- log.cc3;
      dd8 <- 1/xi^2;
      l1[,3]  <-  (-(dd8*dd7-bb1*aa2*cc2*dd6)/cc3^aa2)+dd8*dd3*dd7-
                 aa2*dd7-bb1*aa2*dd3*cc2*dd6;

      ## the second derivatives mm mr mx rr rx xx
    
      l2 <- matrix(0,n,6)
      ee1 <- 1/exp1^(2*rho);
      ee3 <- -1/xi;
      l2[,1]  <-  ee1*(ee3-1)*xi*aa1^(ee3-2)+(ee1*xi*(xi+1))/aa1^2;
      ff7 <- ee3-1;
      l2[,2]  <-  bb1*cc3^ff7+ee1*ff7*xi*cc2*cc3^(ee3-2)-(bb1*dd3)/cc3+
                  (ee1*xi*dd3*cc2)/cc3^2;
      gg7 <- -aa2;
      l2[,3]  <-  (-bb1*cc3^(gg7-1)*(log.cc3/xi^2-bb1*aa2*cc2*dd6))+
                  ee1*cc2*cc3^(gg7-2)+bb1*dd6-(ee1*(xi+1)*cc2)/cc3^2;
      hh4 <- cc2^2;
      l2[,4]  <-  bb1*cc2*cc3^ff7+ee1*ff7*xi*hh4*cc3^(ee3-2)-
              (bb1*dd3*cc2)/cc3+(ee1*xi*dd3*hh4)/cc3^2;
      l2[,5]  <-  (-bb1*cc2*cc3^(gg7-1)*(log.cc3/xi^2-bb1*aa2*cc2*dd6))+
                  ee1*hh4*cc3^(gg7-2)+bb1*cc2*dd6-(ee1*(xi+1)*hh4)/cc3^2;
      jj08 <- 1/cc3^2;
      jj12 <- 1/xi^3;
      jj13 <- 1/cc3^aa2;
      l2[,6]  <-  (-jj13*(dd8*dd7-bb1*aa2*cc2*dd6)^2)-jj13*(ee1*aa2*hh4*jj08+
                   2*bb1*dd8*cc2*dd6-2*jj12*dd7)-2*jj12*dd3*dd7+2*dd8*dd7+
		   2*bb1*dd8*dd3*cc2*dd6-2*bb1*aa2*cc2*dd6+ee1*aa2*dd3*hh4*jj08;

      ## need some link derivatives for derivative transform
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),family$linfo[[2]]$mu.eta(etar),
                   family$linfo[[3]]$mu.eta(etax))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),family$linfo[[2]]$d2link(rho),
                  family$linfo[[3]]$d2link(xi))
    }
 
    l3 <- l4 <- g3 <- g4 <- 0 ## defaults

    if (deriv>1) {
      ## the third derivatives
      ## order mmm mmr mmx mrr mrx mxx rrr rrx rxx xxx
      l3 <- matrix(0,n,10) 
      kk1 <- 1/exp1^(3*rho);
      kk2 <- xi^2;
      l3[,1]  <-  (2*kk1*kk2*(xi+1))/aa1^3-kk1*(ee3-2)*(ee3-1)*kk2*aa1^(ee3-3);
      ll5 <- (xi*cc2)/exp1^rho+1;
      ll8 <- ee3-2;
      l3[,2]  <-  (-2*ee1*ff7*xi*ll5^ll8)-kk1*ll8*ff7*kk2*cc2*ll5^(ee3-3)-
                  (2*ee1*xi*dd3)/ll5^2+(2*kk1*kk2*dd3*cc2)/ll5^3;
      mm10 <- cc3^(gg7-3);
      mm11 <- gg7-2;
      mm12 <- cc3^mm11;
      l3[,3]  <-  ee1*(gg7-1)*xi*mm12*(log.cc3/xi^2-(bb1*aa2*cc2)/cc3)-ee1*mm12-
                  kk1*mm11*xi*cc2*mm10+kk1*cc2*mm10+ee1*dd3*jj08+ee1*xi*jj08-
		  (2*kk1*xi*dd3*cc2)/cc3^3;
      l3[,4]  <-  (-bb1*cc3^ff7)-3*ee1*ff7*xi*cc2*cc3^ll8-kk1*ll8*ff7*kk2*hh4*
                  cc3^(ee3-3)+(bb1*dd3)/cc3-(3*ee1*xi*dd3*cc2)/cc3^2+
		  (2*kk1*kk2*dd3*hh4)/cc3^3;
      oo10 <- gg7-1;
      oo13 <- log.cc3/xi^2;
      l3[,5]  <-  bb1*cc3^oo10*(bb1*oo10*cc2*dd6+oo13)+ee1*oo10*xi*cc2*mm12*
                  (bb1*mm11*cc2*dd6+oo13)+ee1*aa2*cc2*mm12+ee1*oo10*cc2*mm12-
		  bb1*dd6+2*ee1*dd3*cc2*jj08+ee1*xi*cc2*jj08-
		  (2*xi*dd3*cc2^2)/(exp1^(3*rho)*cc3^3);
      pp07 <- (-1/xi)-1;
      pp08 <- cc3^pp07;
      l3[,6]  <-  (-bb1*pp08*(bb1*pp07*cc2*dd6+dd8*dd7)^2)-bb1*pp08*
                  ((-ee1*pp07*hh4*jj08)+2*bb1*dd8*cc2*dd6-(2*dd7)/xi^3)-
		  2*ee1*cc2*jj08+(2*(xi+1)*hh4)/(exp1^(3*rho)*cc3^3);
      qq05 <- cc2^3;
      l3[,7]  <-  (-bb1*cc2*cc3^ff7)-3*ee1*ff7*xi*hh4*cc3^ll8-
                  kk1*ll8*ff7*kk2*qq05*cc3^(ee3-3)+(bb1*dd3*cc2)/cc3-
		  (3*ee1*xi*dd3*hh4)/cc3^2+(2*kk1*kk2*dd3*qq05)/cc3^3;
      rr17 <- log.cc3/xi^2-bb1*aa2*cc2*dd6;
      l3[,8]  <-  bb1*cc2*cc3^oo10*rr17+ee1*oo10*xi*hh4*mm12*rr17-2*ee1*hh4*mm12-
                  kk1*mm11*xi*qq05*mm10+kk1*qq05*mm10-bb1*cc2*dd6+
		  2*ee1*dd3*hh4*jj08+ee1*xi*hh4*jj08-(2*kk1*xi*dd3*qq05)/cc3^3;
      l3[,9]  <-  (-bb1*cc2*pp08*(bb1*pp07*cc2*dd6+dd8*dd7)^2)-bb1*cc2*pp08*
                  ((-ee1*pp07*hh4*jj08)+2*bb1*dd8*cc2*dd6-(2*dd7)/xi^3)-
		  2*ee1*hh4*jj08+(2*(xi+1)*cc2^3)/(exp1^(3*rho)*cc3^3);
      tt08 <- 1/cc3^3;
      tt16 <- 1/xi^4;
      tt18 <- dd8*dd7-bb1*aa2*cc2*dd6;
      l3[,10]  <-  (-jj13*tt18^3)-3*jj13*(ee1*aa2*hh4*jj08+2*bb1*dd8*cc2*dd6-2*jj12*dd7)*
              tt18-jj13*((-2*kk1*aa2*qq05*tt08)-3*ee1*dd8*hh4*jj08-6*bb1*jj12*cc2*dd6+
	      6*tt16*dd7)+6*tt16*dd3*dd7-6*jj12*dd7-6*bb1*jj12*dd3*cc2*dd6+
	      6*bb1*dd8*cc2*dd6-3*ee1*dd8*dd3*hh4*jj08+3*ee1*aa2*hh4*jj08-
	      2*kk1*aa2*dd3*qq05*tt08;
 
      g3 <- cbind(family$linfo[[1]]$d3link(mu),family$linfo[[2]]$d3link(rho),
                  family$linfo[[3]]$d3link(xi))
    }

    if (deriv>3) {
      ## the fourth derivatives
      ## mmmm mmmr mmmx mmrr mmrx mmxx mrrr mrrx mrxx mxxx
      ## rrrr rrrx rrxx rxxx xxxx
      l4 <- matrix(0,n,15)
      uu1 <- 1/exp1^(4*rho);
      uu2 <- xi^3;
      l4[,1]  <-  uu1*(ee3-3)*(ee3-2)*(ee3-1)*uu2*aa1^(ee3-4)+(6*uu1*uu2*(xi+1))/aa1^4;
      vv09 <- ee3-3;
      l4[,2]  <-  3*kk1*ll8*ff7*kk2*ll5^vv09+uu1*vv09*ll8*ff7*uu2*cc2*ll5^(ee3-4)-
                 (6*kk1*kk2*dd3)/ll5^3+(6*uu1*uu2*dd3*cc2)/ll5^4;
      ww11 <- gg7-3;
      ww12 <- cc3^(gg7-4);
      ww15 <- cc3^ww11;
      l4[,3]  <- (-kk1*mm11*oo10*kk2*ww15*(log.cc3/kk2-(bb1*aa2*cc2)/cc3))+
                 2*kk1*mm11*xi*ww15-kk1*ww15+uu1*ww11*mm11*kk2*cc2*ww12-
		 uu1*oo10*xi*cc2*ww12-uu1*ww11*xi*cc2*ww12+2*kk1*kk2*tt08+
		 4*kk1*xi*dd3*tt08-(6*uu1*kk2*dd3*cc2)/cc3^4;
      l4[,4]  <-  4*ee1*ff7*xi*ll5^ll8+5*kk1*ll8*ff7*kk2*cc2*ll5^vv09+
                  uu1*vv09*ll8*ff7*uu2*hh4*ll5^(ee3-4)+(4*ee1*xi*dd3)/ll5^2-
		  (10*kk1*kk2*dd3*cc2)/ll5^3+(6*uu1*uu2*dd3*hh4)/ll5^4;
      yy18 <- log.cc3/kk2;
      l4[,5]  <-  (-2*ee1*oo10*xi*mm12*(bb1*mm11*cc2*dd6+yy18))-
                   kk1*mm11*oo10*kk2*cc2*ww15*(bb1*ww11*cc2*dd6+yy18)-
		   2*ee1*aa2*mm12-2*ee1*oo10*mm12-2*kk1*mm11*oo10*xi*cc2*ww15-
		   kk1*oo10*cc2*ww15-kk1*mm11*cc2*ww15-2*ee1*dd3*jj08-
		   2*ee1*xi*jj08+2*kk1*kk2*cc2*tt08+8*kk1*xi*dd3*cc2*tt08-
		   (6*kk2*dd3*cc2^2)/(exp1^(4*rho)*cc3^4);
      l4[,6]  <-  ee1*oo10*xi*mm12*tt18^2-2*ee1*mm12*tt18-2*kk1*mm11*xi*cc2*ww15*tt18+
                  2*kk1*cc2*ww15*tt18+ee1*oo10*xi*mm12*(ee1*aa2*hh4*jj08+2*bb1*
		  dd8*cc2*dd6-(2*dd7)/xi^3)+4*kk1*cc2*ww15+2*uu1*ww11*xi*hh4*ww12-
		  4*uu1*hh4*ww12+2*ee1*jj08-4*kk1*dd3*cc2*tt08-4*kk1*xi*cc2*tt08+
		  (6*uu1*xi*dd3*hh4)/cc3^4;
      l4[,7]  <-  bb1*cc3^ff7+7*ee1*ff7*xi*cc2*cc3^ll8+6*kk1*ll8*ff7*kk2*hh4*cc3^vv09+
                  uu1*vv09*ll8*ff7*uu2*qq05*cc3^(ee3-4)-(bb1*dd3)/cc3+
		  (7*ee1*xi*dd3*cc2)/cc3^2-(12*kk1*kk2*dd3*hh4)/cc3^3+
		  (6*uu1*uu2*dd3*qq05)/cc3^4;
      l4[,8]  <-  (-bb1*cc3^oo10*(bb1*oo10*cc2*dd6+yy18))-3*ee1*oo10*xi*cc2*mm12*
                  (bb1*mm11*cc2*dd6+yy18)-kk1*mm11*oo10*kk2*hh4*ww15*
		  (bb1*ww11*cc2*dd6+yy18)-3*ee1*aa2*cc2*mm12-3*ee1*oo10*cc2*mm12-
		  2*kk1*mm11*oo10*xi*hh4*ww15-kk1*oo10*hh4*ww15-kk1*mm11*hh4*ww15+
		  bb1*dd6-4*ee1*dd3*cc2*jj08-3*ee1*xi*cc2*jj08+2*kk1*kk2*hh4*tt08+
		  10*kk1*xi*dd3*hh4*tt08-(6*kk2*dd3*cc2^3)/(exp1^(4*rho)*cc3^4); 
      ad17 <- 2*bb1*dd8*cc2*dd6;
      ad19 <- -(2*dd7)/xi^3;
      ad20 <- cc3^oo10;
      ad21 <- dd8*dd7;
      ad22 <- ad21+bb1*mm11*cc2*dd6;
      l4[,9] <-  bb1*ad20*(bb1*oo10*cc2*dd6+ad21)^2+ee1*oo10*xi*cc2*mm12*ad22^2+
                  2*ee1*aa2*cc2*mm12*ad22+2*ee1*oo10*cc2*mm12*ad22+
		  bb1*ad20*((-ee1*oo10*hh4*jj08)+ad17+ad19)+ee1*oo10*xi*cc2*mm12*
		  ((-ee1*mm11*hh4*jj08)+ad17+ad19)+4*ee1*cc2*jj08-6*kk1*dd3*hh4*tt08-
		  4*kk1*xi*hh4*tt08+(6*xi*dd3*cc2^3)/(exp1^(4*rho)*cc3^4);
      ae16 <- dd8*dd7+bb1*pp07*cc2*dd6;
      l4[,10]  <- (-bb1*pp08*ae16^3)-3*bb1*pp08*((-ee1*pp07*hh4*jj08)+
                  2*bb1*dd8*cc2*dd6-2*jj12*dd7)*ae16-bb1*pp08*(2*kk1*pp07*qq05*tt08-
		  3*ee1*dd8*hh4*jj08-6*bb1*jj12*cc2*dd6+(6*dd7)/xi^4)+6*kk1*hh4*tt08-
		  (6*(xi+1)*qq05)/(exp1^(4*rho)*cc3^4);
      af05 <- cc2^4;
      l4[,11]  <-  bb1*cc2*cc3^ff7+7*ee1*ff7*xi*hh4*cc3^ll8+6*kk1*ll8*ff7*kk2*qq05*
                   cc3^vv09+uu1*vv09*ll8*ff7*uu2*af05*cc3^(ee3-4)-(bb1*dd3*cc2)/cc3+
		   (7*ee1*xi*dd3*hh4)/cc3^2-(12*kk1*kk2*dd3*qq05)/cc3^3+
		   (6*uu1*uu2*dd3*af05)/cc3^4;
      ag23 <- log.cc3/kk2-bb1*aa2*cc2*dd6;
      l4[,12]  <- (-bb1*cc2*cc3^oo10*ag23)-3*ee1*oo10*xi*hh4*mm12*ag23-
                  kk1*mm11*oo10*kk2*qq05*ww15*ag23+4*ee1*hh4*mm12+
		  5*kk1*mm11*xi*qq05*ww15-4*kk1*qq05*ww15+uu1*ww11*mm11*kk2*af05*ww12-
		  uu1*oo10*xi*af05*ww12-uu1*ww11*xi*af05*ww12+bb1*cc2*dd6-
		  4*ee1*dd3*hh4*jj08-3*ee1*xi*hh4*jj08+2*kk1*kk2*qq05*tt08+
		  10*kk1*xi*dd3*qq05*tt08-(6*uu1*kk2*dd3*af05)/cc3^4;
      ah24 <- (-(2*dd7)/xi^3)+2*bb1*dd8*cc2*dd6+ee1*aa2*hh4*jj08;
      ah27 <- tt18^2;
      l4[,13]  <- bb1*cc2*ad20*ah27+ee1*oo10*xi*hh4*mm12*ah27-4*ee1*hh4*mm12*tt18-
                  2*kk1*mm11*xi*qq05*ww15*tt18+2*kk1*qq05*ww15*tt18+bb1*cc2*ad20*ah24+
		  ee1*oo10*xi*hh4*mm12*ah24+6*kk1*qq05*ww15+2*uu1*ww11*xi*af05*ww12-
		  4*uu1*af05*ww12+4*ee1*hh4*jj08-6*kk1*dd3*qq05*tt08-
		  4*kk1*xi*qq05*tt08+(6*uu1*xi*dd3*af05)/cc3^4;
      l4[,14]  <- (-bb1*cc2*pp08*ae16^3)-3*bb1*cc2*pp08*((-ee1*pp07*hh4*jj08)+
                  2*bb1*dd8*cc2*dd6-2*jj12*dd7)*ae16-bb1*cc2*pp08*(2*kk1*pp07*qq05*
		  tt08-3*ee1*dd8*hh4*jj08-6*bb1*jj12*cc2*dd6+(6*dd7)/xi^4)+
		  6*kk1*qq05*tt08-(6*(xi+1)*cc2^4)/(exp1^(4*rho)*cc3^4);
      aj08 <- 1/cc3^4;
      aj20 <- 1/xi^5;
      aj23 <- (-2*jj12*dd7)+2*bb1*dd8*cc2*dd6+ee1*aa2*hh4*jj08;
      l4[,15]  <- (-jj13*tt18^4)-6*jj13*aj23*tt18^2-3*jj13*aj23^2-
                  4*jj13*((-2*kk1*aa2*qq05*tt08)-3*ee1*dd8*hh4*jj08-6*bb1*jj12*cc2*dd6+
		  6*tt16*dd7)*tt18-jj13*(6*uu1*aa2*af05*aj08+8*kk1*dd8*qq05*tt08+
		  12*ee1*jj12*hh4*jj08+24*bb1*tt16*cc2*dd6-24*aj20*dd7)-
		  24*aj20*dd3*dd7+24*tt16*dd7+24*bb1*tt16*dd3*cc2*dd6-
		  24*bb1*jj12*cc2*dd6+12*ee1*jj12*dd3*hh4*jj08-12*ee1*dd8*hh4*jj08+
		  8*kk1*dd8*dd3*qq05*tt08-8*kk1*aa2*qq05*tt08+6*uu1*aa2*dd3*af05*aj08;
   
      g4 <- cbind(family$linfo[[1]]$d4link(mu),family$linfo[[2]]$d4link(rho),
                  family$linfo[[3]]$d4link(xi))
    }
    if (deriv) {
      i2 <- family$tri$i2; i3 <- family$tri$i3
      i4 <- family$tri$i4
   
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)

      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                      d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D,sandwich=sandwich)
      if (ncv) {
        ret$l1 <- de$l1; ret$l2 = de$l2; ret$l3 = de$l3
      }		      
    } else ret <- list()
    ret$l <- l; ret$l0 <- l0; ret
  } ## end ll gevlss

  sandwich <- function(y,X,coef,wt,family,offset=NULL) {
  ## compute filling for sandwich estimate of cov matrix
    ll(y,X,coef,wt,family,offset=NULL,deriv=1,sandwich=TRUE)$lbb
  }

  initialize <- expression({
  ## start out with xi close to zero. If xi==0 then
  ## mean is mu + sigma*gamma and var is sigma^2*pi^2/6
  ## where sigma = exp(rho) and gamma is Euler's constant.  
  ## idea is to regress g(y) on model matrix for mean, and then 
  ## to regress the corresponding log absolute residuals on 
  ## the model matrix for log(sigma) - may be called in both
  ## gam.fit5 and initial.spg... note that appropriate E scaling
  ## for full calculation may be inappropriate for initialization 
  ## which is basically penalizing something different here.
  ## best we can do here is to use E only as a regularizer.
      .euler <- 0.5772156649015328606065121 ## Euler's constant
      n <- rep(1, nobs)
      ## should E be used unscaled or not?..
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      if (is.null(start)) {
        jj <- attr(x,"lpi")
	if (is.list(x)) { ## discrete case
	  ## LP 1...
          start <- rep(0,max(unlist(jj)))
          yt1 <- if (family$link[[1]]=="identity") y else 
                 family$linfo[[1]]$linkfun(abs(y)+max(y)*1e-7)
          R <- suppressWarnings(mchol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,
	         v=x$v,qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[1]])+crossprod(E[,jj[[1]]])))
	  Xty <- XWyd(x$Xd,rep(1,length(y)),yt1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[1]])
          piv <- attr(R,"pivot");rrank <- attr(R,"rank");startji <- rep(0,ncol(R))
          if (rrank<ncol(R)) {
              R <- R[1:rrank,1:rrank]
	      piv <- piv[1:rrank]
          }
          startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
	  startji[!is.finite(startji)] <- 0
	  start[jj[[1]]] <- startji
          ## LP 2...
	  lres1 <- Xbd(x$Xd,start,k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,qc=x$qc,drop=x$drop,lt=x$lpid[[1]])
          lres1 <- log(abs(y-family$linfo[[1]]$linkinv(lres1)))
          R <- suppressWarnings(mchol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,
	         v=x$v,qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[2]])+crossprod(E[,jj[[2]]])))
	  Xty <- XWyd(x$Xd,rep(1,length(y)),lres1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[2]])
          piv <- attr(R,"pivot");rrank <- attr(R,"rank");startji <- rep(0,ncol(R))
          if (rrank<ncol(R)) {
              R <- R[1:rrank,1:rrank]
	      piv <- piv[1:rrank]
          }
          startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
	  startji[!is.finite(startji)] <- 0
	  start[jj[[2]]] <- startji
	  ## LP 3...
	  yt1 <- rep(family$linfo[[3]]$linkfun(1e-3),length(yt1))
	  R <- suppressWarnings(mchol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,
	         v=x$v,qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[3]])+crossprod(E[,jj[[3]]])))
	  Xty <- XWyd(x$Xd,rep(1,length(y)),yt1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[3]])
          piv <- attr(R,"pivot");rrank <- attr(R,"rank");startji <- rep(0,ncol(R))
          if (rrank<ncol(R)) {
              R <- R[1:rrank,1:rrank]
	      piv <- piv[1:rrank]
          }
	  fob <- function(m=1) {
            startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]*m))
	    startji[!is.finite(startji)] <- 0
	    start[jj[[3]]] <- startji
	    list(l=family$ll(y,x,start,weights,family,offset)$l,start=start)
	  } ## fob
	  f0b <- fob(); dm <- .2; mm <- 1; up <- FALSE
	  while (mm>-4.2 && mm <4.2) { ## crude search for improved initial xi
            f1b <- fob(mm+dm)
	    if (is.finite(f1b$l) && f1b$l>f0b$l) {
              up <- TRUE;f0b <- f1b;mm <- mm + dm
            } else if (up) { ## f0b best
              break
            } else if (dm>0) dm <- -dm else break
          }
	  if (!is.finite(f1b$l)) { ## move back a bit from non finite regime
            f1b <- fob(mm-dm)
	    if (is.finite(f1b$l)) f0b <- f1b
          }
	  start <- f0b$start
        } else { ## not discrete case
	  start <- rep(0,ncol(x))
          yt1 <- if (family$link[[1]]=="identity") y else 
                 family$linfo[[1]]$linkfun(abs(y)+max(y)*1e-7)
          x1 <- x[,jj[[1]],drop=FALSE]
          e1 <- E[,jj[[1]],drop=FALSE] ## square root of total penalty
          #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
          if (use.unscaled) {
            qrx <- qr(rbind(x1,e1))
            x1 <- rbind(x1,e1)
            startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
            startji[!is.finite(startji)] <- 0       
          } else startji <- pen.reg(x1,e1,yt1)
          start[jj[[1]]] <- startji
          lres1 <- log(abs(y-family$linfo[[1]]$linkinv(x[,jj[[1]],drop=FALSE]%*%start[jj[[1]]])))
          x1 <-  x[,jj[[2]],drop=FALSE];e1 <- E[,jj[[2]],drop=FALSE]
          #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
          if (use.unscaled) {
            x1 <- rbind(x1,e1)
            startji <- qr.coef(qr(x1),c(lres1,rep(0,nrow(E))))   
            startji[!is.finite(startji)] <- 0
          } else startji <- pen.reg(x1,e1,lres1)
          start[jj[[2]]] <- startji
	  x1 <-  x[,jj[[3]],drop=FALSE]
	  qrx1 <- qr(x1);yt1 <- rep(family$linfo[[3]]$linkfun(1e-3),nrow(x1))
	  fob <- function(m=1) {
	    startji <- qr.coef(qrx1,yt1*m)   
            startji[!is.finite(startji)] <- 0
	    start[jj[[3]]] <- startji
	    list(l=family$ll(y,x,start,weights,family,offset)$l,start=start)
	  } ## fob
	  f0b <- fob(); dm <- .2; mm <- 1; up <- FALSE
	  while (mm>-4.2 && mm <4.2) { ## crude search for improved initial xi
            f1b <- fob(mm+dm)
	    if (is.finite(f1b$l) && f1b$l>f0b$l) {
              up <- TRUE;f0b <- f1b;mm <- mm + dm
            } else if (up) { ## f0b best
              break
            } else if (dm>0) dm <- -dm else break
          }
	  if (!is.finite(f1b$l)) { ## move back a bit from non finite regime
            f1b <- fob(mm-dm)
	    if (is.finite(f1b$l)) f0b <- f1b
          }
	  start <- f0b$start  
        } ## non-discrete initiailization
      }	
  }) ## initialize gevlss

  rd <- function(mu,wt,scale) {
    Fi.gev <- function(z,mu,sigma,xi) {
      ## GEV inverse cdf.
      xi[abs(xi)<1e-8] <- 1e-8 ## approximate xi=0, by small xi
      x <- mu + ((-log(z))^-xi-1)*sigma/xi
    }
    Fi.gev(runif(nrow(mu)),mu[,1],exp(mu[,2]),mu[,3])
  } ## gevlss rd

  structure(list(family="gevlss",ll=ll,link=paste(link),nlp=3,ncv=ncv,sandwich=sandwich,
    tri = trind.generator(3), ## symmetric indices for accessing derivative arrays
    initialize=initialize,postproc=postproc,residuals=residuals,
    linfo = stats, ## link information list
    d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
    ls=1, ## signals that ls not needed here
    rd=rd,
    available.derivs = 2, ## can use full Newton here
    discrete.ok = TRUE,qapprox=TRUE
    ),class = c("general.family","extended.family","family"))
} ## end gevlss


## attempt at a Tweedie gamlss family, suitable for extended FS
## estimation (i.e. avoiding horrendous 4th derivative series
## summations)


tw.null.fit <- function(y,a=1.001,b=1.999) {
## MLE of c(mu,p,phi) for sample of Tweedie data, y.
## Stabilized, step controlled, Newton
  th <- c(0,0,0)
  ld <- ldTweedie(y,exp(th[1]),rho=th[3],theta=th[2],a=a,b=b,all.derivs=TRUE)
  lds <- colSums(ld)
  for (i in 1:50) {
    g <- lds[c(7,4,2)]
    if (sum(abs(g)>1e-9*abs(lds[1]))==0) break
    g[1] <- g[1] * exp(th[1]) ## work on log scale for mu
    H <- matrix(0,3,3) ## mu, th, rh 
    diag(H) <- c(lds[8],lds[5],lds[3])
    H[1,2] <- H[2,1] <- lds[9]
    H[1,3] <- H[3,1] <- lds[10]
    H[2,3] <- H[3,2] <- lds[6]
    H[,1] <- H[,1]*exp(th[1])
    H[1,-1] <- H[1,-1] * exp(th[1])
    eh <- eigen(H,symmetric=TRUE) 
    tol <- max(abs(eh$values))*1e-7
    eh$values[eh$values>-tol] <- -tol
    step <- as.numeric(eh$vectors%*%((t(eh$vectors)%*%g)/eh$values))
    ms <- max(abs(step))
    if (ms>3) step <- step*3/ms
    ok <- FALSE
    while (!ok) {
      th1 <- th - step
      ld1 <- ldTweedie(y,exp(th1[1]),rho=th1[3],theta=th1[2],a=a,b=b,all.derivs=TRUE)
      if (sum(ld1[,1])<lds[1]) step <- step/2 else {
        ok <- TRUE
        th <- th1
        lds <- colSums(ld1)
      }
    }
  }
  p <- if (th[2]>0) (b + a*exp(-th[2]))/(1+exp(-th[2])) else (b*exp(th[2])+a)/(exp(th[2])+1)
  c(exp(th[1]),p,exp(th[3])) # mu, p, sigma
} ## tw.null.fit


twlss <- function(link=list("log","identity","identity"),a=1.01,b=1.99) {
## General family for Tweedie location scale model...
## so mu is mu1, rho = log(sigma) is mu2 and transformed p is mu3
## Need to redo ldTweedie to allow vector p and phi
## -- advantage is that there is no point buffering
## 1. get derivatives wrt mu, rho and p.
## 2. get required link derivatives and tri indices.
## 3. transform derivs to derivs wrt eta (gamlss.etamu).
## 4. get the grad and Hessian etc for the model
##    via a call to gamlss.gH  
  
  ## first deal with links and their derivatives...
  if (length(link)!=3) stop("gevlss requires 3 links specified as character strings")
  okLinks <- list(c("log", "identity", "sqrt"),"identity",c("identity"))
  stats <- list()
  for (i in 1:3) {
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- make.link(link[[i]]) else 
    stop(link[[i]]," link not available for mu parameter of twlss")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
           mu.eta=stats[[i]]$mu.eta),
           class="family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }

  env <- new.env(parent = .GlobalEnv)
  assign(".a",a, envir = env);assign(".b",b, envir = env)

  residuals <- function(object,type=c("deviance","pearson","response")) {
      
      a <- get(".a");b <- get(".b")
      type <- match.arg(type)
      mu <- object$fitted.values[,1]
      p <- object$fitted.values[,2]
      ind <- p > 0;
      ethi <- exp(-p[ind]);ethni <- exp(p[!ind])
      p[ind] <- (b+a*ethi)/(1+ethi)
      p[!ind] <- (b*ethni+a)/(ethni+1)
      phi <- exp(object$fitted.values[,3])
      if (type=="pearson") {
        rsd <- (object$y-mu)/sqrt(phi*mu^p) ## Pearson
      } else if (type=="response") rsd <- object$y-mu else {
        y1 <- object$y + (object$y == 0)
        theta <- (y1^(1 - p) - mu^(1 - p))/(1 - p)
        kappa <- (object$y^(2 - p) - mu^(2 - p))/(2 - p)
        rsd <- sign(object$y-mu)*sqrt(pmax(2 * (object$y * theta - kappa) * object$prior.weights/phi,0))
      }
      return(rsd) ## (y-mu)/sigma 
    } ## twlss residuals
    
  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    ## used for dev explained - really a mean scale concept.
    ## makes no sense to use single scale param here...
    tw.para <- tw.null.fit(object$y) ## paramaters mu, p and phi
    tw.y1 <- object$y + (object$y == 0)
    tw.theta <- (tw.y1^(1 - tw.para[2]) - tw.para[1]^(1 - tw.para[2]))/(1 - tw.para[2])
    tw.kappa <- (object$y^(2 - tw.para[2]) - tw.para[1]^(2 - tw.para[2]))/(2 - tw.para[2])
    object$null.deviance <- sum(pmax(2 * (object$y * tw.theta - tw.kappa) * object$prior.weights/exp(object$fitted.values[,3]),0))
  })

  ll <- function(y,X,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL,sandwich=FALSE) {
  ## function defining the gamlss Tweedie model log lik. 
  ## deriv: 0 - eval
  ##        1 - grad and Hess
  ##        2 - diagonal of first deriv of Hess
  ##        3 - first deriv of Hess
  ##        4 - everything.
  ## This family does not have code for 3 and 4
    if (is.null(offset)) offset <- list(0,0,0) else offset[[4]] <- 0
    for (i in 1:3) if (is.null(offset[[i]])) offset[[i]] <- 0
    jj <- attr(X,"lpi") ## extract linear predictor index
    eta <- X[,jj[[1]],drop=FALSE]%*%coef[jj[[1]]] + offset[[1]]
    mu <- family$linfo[[1]]$linkinv(eta) ## mean
    theta <- X[,jj[[2]],drop=FALSE]%*%coef[jj[[2]]] +offset[[2]] ## transformed p
    rho <- X[,jj[[3]],drop=FALSE]%*%coef[jj[[3]]] +offset[[3]] ## log scale parameter
    a <- get(".a");b <- get(".b")
   
    ld <- ldTweedie(y,mu=mu,p=NA,phi=NA,rho=rho,theta=theta,a=a,b=b,all.derivs=TRUE)
    ## m, t, r ; mm, mt, mr, tt, tr, rr
    l0 <- ld[,1]
    l <- sum(l0)
    l1 <- cbind(ld[,7],ld[,4],ld[,2])
    l2 <- cbind(ld[,8],ld[,9],ld[,10],ld[,5],ld[,6],ld[,3])

    ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),family$linfo[[2]]$mu.eta(theta),
                   family$linfo[[3]]$mu.eta(rho))
    g2 <- cbind(family$linfo[[1]]$d2link(mu),family$linfo[[2]]$d2link(theta),
                  family$linfo[[3]]$d2link(rho))

    n <- length(y)
    
    l3 <- l4 <- g3 <- g4 <- 0 ## defaults


    if (deriv) {
      i2 <- family$tri$i2;#i3 <- i4 <- 0
      i3 <- family$tri$i3;i4 <- family$tri$i4
   
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,0)

      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                      d1b=d1b,d2b=d2b,deriv=0,fh=fh,D=D,sandwich=sandwich) 
    } else ret <- list()
    ret$l <- l; ret$l0 <- l0; ret
  } ## end ll twlss

  sandwich <- function(y,X,coef,wt,family,offset=NULL) {
  ## compute filling for sandwich estimate of cov matrix
    ll(y,X,coef,wt,family,offset=NULL,deriv=1,sandwich=TRUE)$lbb
  }

  initialize <- expression({
   ## idea is to regress g(y) on model matrix for mean, and then 
   ## to regress the corresponding log absolute scaled residuals on 
   ## the model matrix for log(sigma) - may be called in both
   ## gam.fit5 and initial.spg... note that appropriate E scaling
   ## for full calculation may be inappropriate for initialization 
   ## which is basically penalizing something different here.
   ## best we can do here is to use E only as a regularizer.
   ## initial theta params are zero for p = 1.5
      n <- rep(1, nobs)
      ## should E be used unscaled or not?..
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      if (is.null(start)) {
        jj <- attr(x,"lpi")
        start <- rep(0,ncol(x))
        yt1 <- if (family$link[[1]]=="identity") y else 
               family$linfo[[1]]$linkfun(abs(y)+max(y)*1e-7)
        x1 <- x[,jj[[1]],drop=FALSE]
        e1 <- E[,jj[[1]],drop=FALSE] ## square root of total penalty
        #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
        if (use.unscaled) {
          qrx <- qr(rbind(x1,e1))
          x1 <- rbind(x1,e1)
          startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
          startji[!is.finite(startji)] <- 0       
        } else startji <- pen.reg(x1,e1,yt1)
        ## now the scale parameter
        start[jj[[1]]] <- startji
	mu1 <- family$linfo[[1]]$linkinv(x[,jj[[1]],drop=FALSE]%*%start[jj[[1]]])
        lres1 <- log(abs((y-mu1)/mu1^1.5))
        x1 <-  x[,jj[[3]],drop=FALSE];e1 <- E[,jj[[3]],drop=FALSE]
        #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
        if (use.unscaled) {
          x1 <- rbind(x1,e1)
          startji <- qr.coef(qr(x1),c(lres1,rep(0,nrow(E))))   
          startji[!is.finite(startji)] <- 0
        } else startji <- pen.reg(x1,e1,lres1)
        start[jj[[3]]] <- startji
      } ## is.null(start)
  }) ## initialize twlss

  environment(ll) <- environment(residuals) <- env

  structure(list(family="twlss",ll=ll,link=paste(link),nlp=3,sandwich=sandwich,
    tri = trind.generator(3), ## symmetric indices for accessing derivative arrays
    initialize=initialize,postproc=postproc,residuals=residuals,
    linfo = stats, ## link information list
    d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
    ls=1, ## signals that ls not needed here
    available.derivs = 0 ## no higher derivs
    ),class = c("general.family","extended.family","family"))
} ## end twlss


lb.linkfun <- function(mu,b=-7) {
## lower bound link function - see gammals for related routines. 
  eta <- mub <- mu-b
  ii <- mub < .Machine$double.eps
  if (any(ii)) eta[ii] <- log(.Machine$double.eps)
  jj <- mub > -log(.Machine$double.eps)
  if (any(jj)) eta[jj] <- mub[jj]
  jj <- !jj & !ii
  if (any(jj)) eta[jj] <- log(exp(mub[jj])-1)
  eta
} ## lb.linkfun


gammals <- function(link=list("identity","log"),b=-7) {
## General family for gamma location scale model...
## parameterization is in terms of log mean and log scale.
## so log(mu) is mu1, scale = log(sigma) is mu2
## 1. get derivatives wrt mu, rho and xi.
## 2. get required link derivatives and tri indices.
## 3. transform derivs to derivs wrt eta (gamlss.etamu).
## 4. get the grad and Hessian etc for the model
##    via a call to gamlss.gH  
  
  ## first deal with links and their derivatives...
  if (length(link)!=2) stop("gammals requires 2 links specified as character strings")
  okLinks <- list(c("identity"),c("identity","log"))
  stats <- list()
  for (i in 1:2) {
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- make.link(link[[i]]) else 
    stop(link[[i]]," link not available for gammals")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
           mu.eta=stats[[i]]$mu.eta),
           class="family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }
  if (link[[2]]=="log") { ## g^{-1}(eta) = b + log(1+exp(eta)) link
    stats[[2]]$valideta <- function(eta) TRUE
    stats[[2]]$linkfun <- eval(parse(text=paste("function(mu,b=",b,") {\n eta <- mub <- mu-b;\n",
      "ii <- mub < .Machine$double.eps;\n if (any(ii)) eta[ii] <- log(.Machine$double.eps);\n",
      "jj <- mub > -log(.Machine$double.eps);if (any(jj)) eta[jj] <- mub[jj];\n",
      "jj <- !jj & !ii;if (any(jj)) eta[jj] <- log(exp(mub[jj])-1);eta }")))

    stats[[2]]$mu.eta <- eval(parse(text=paste("function(eta,b=",b,") {\n",
      "ii <- eta < 0;eta <- exp(-eta*sign(eta))\n",
      "if (any(ii)) { ei <- eta[ii];eta[ii] <- ei/(1+ei)}\n",
      "ii <- !ii;if (any(ii)) eta[ii] <- 1/(1+eta[ii])\n",
      "eta }\n")))

    stats[[2]]$linkinv <- eval(parse(text=paste("function(eta,b=",b,") {\n",
      "mu <- eta;ii <- eta > -log(.Machine$double.eps)\n",
      "if (any(ii)) mu[ii] <- b + eta[ii]\n",
      "ii <- !ii;if (any(ii)) mu[ii] <- b + log(1 + exp(eta[ii]))\n",
      "mu }\n")))

    stats[[2]]$d2link <- eval(parse(text=paste("function(mu,b=",b,") {\n",
      "eta <- lb.linkfun(mu,b=b); ii <- eta > 0\n",
      "eta <- exp(-eta*sign(eta))\n",
      "if (any(ii)) { ei <- eta[ii];eta[ii] <- -(ei^2 + ei) }\n",
      "ii <- !ii;if (any(ii)) { ei <- eta[ii];eta[ii] <- -(1+ei)/ei^2 }\n",
      "eta }\n")))

    stats[[2]]$d3link <- eval(parse(text=paste("function(mu,b=",b,") {\n",
      "eta <- lb.linkfun(mu,b=b);ii <- eta > 0\n",
      "eta <- exp(-eta*sign(eta))\n",
      "if (any(ii)) { ei <- eta[ii]; eta[ii] <- (2*ei^2+ei)*(ei+1) }\n",
      "ii <- !ii;if (any(ii)) { ei <- eta[ii]; eta[ii] <- (2+ei)*(1+ei)/ei^3 }\n",
      "eta }\n")))
    
    stats[[2]]$d4link <- eval(parse(text=paste("function(mu,b=",b,") {\n",
      "eta <- lb.linkfun(mu,b=b);ii <- eta > 0\n",
      "eta <- exp(-eta*sign(eta))\n",
      "if (any(ii)) { ei <- eta[ii];eta[ii] <- -(6*ei^3+6*ei^2+ei)*(ei+1) }\n",
      "ii <- !ii;if (any(ii)) { ei <- eta[ii]; eta[ii] <- -(6+6*ei+ei^2)*(1+ei)/ei^4 }\n",
      "eta }\n")))
  }

  residuals <- function(object,type=c("deviance","pearson","response")) {
      mu <- object$fitted.values[,1]
      rho <- object$fitted.values[,2]
      y <- object$y
      type <- match.arg(type)
      if (type=="deviance") {
        rsd <- 2*((y-mu)/mu-log(y/mu))*exp(-rho)
        rsd <- sqrt(pmax(0,rsd))*sign(y-mu)
      } else if (type=="pearson") {
        rsd <- (y-mu)/(exp(rho*.5)*mu)
      } else {
        rsd <- y-mu
      }
      rsd
    } ## gammls residuals
    
  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    object$fitted.values[,1] <- exp(object$fitted.values[,1])
    .my <- mean(object$y)
    object$null.deviance <- sum(((object$y-.my)/.my-log(object$y/.my))*exp(-object$fitted.values[,2]))*2
  })

  ncv <- function(X,y,wt,nei,beta,family,llf,H=NULL,Hi=NULL,R=NULL,offset=NULL,dH=NULL,db=NULL,deriv=FALSE,nt=1) {
    gamlss.ncv(X,y,wt,nei,beta,family,llf,H=H,Hi=Hi,R=R,offset=offset,dH=dH,db=db,deriv=deriv,nt=nt)
  } ## ncv  

 
  ll <- function(y,X,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL,eta=NULL,ncv=FALSE,sandwich=FALSE) {
  ## function defining the gamlss gamma model log lik. 
  ## deriv: 0 - eval
  ##        1 - grad and Hess
  ##        2 - diagonal of first deriv of Hess
  ##        3 - first deriv of Hess
  ##        4 - everything.

    if (!is.null(offset)) offset[[3]] <- 0
    discrete <- is.list(X)
    jj <- attr(X,"lpi") ## extract linear predictor index
    if (is.null(eta)) {
      eta <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[1]]) else X[,jj[[1]],drop=FALSE]%*%coef[jj[[1]]]
      if (!is.null(offset[[1]])) eta <- eta + offset[[1]] ## log mu
  
      etat <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[2]]) else X[,jj[[2]],drop=FALSE]%*%coef[jj[[2]]]
      if (!is.null(offset[[2]])) etat <- etat + offset[[2]]
    } else {
      etat <- eta[,2]
      eta <- eta[,1]
    }
    mu <- family$linfo[[1]]$linkinv(eta) ## mean
    th <-  family$linfo[[2]]$linkinv(etat) ## log sigma
 
    eth <- exp(-th) ## 1/exp1^th;
    logy <- log(y);
    ethmu <- exp(-th-mu)
    ethmuy <- ethmu*y
    etlymt <- eth*(logy-mu-th)
    n <- length(y)

    l0  <-  etlymt-logy-ethmuy-lgamma(eth) ## l
    l <- sum(l0)

    if (deriv>0) {
      l1 <- matrix(0,n,2)
     
      l1[,1]  <- ethmuy-eth ## lm
      digeth <- digamma(eth)
      l1[,2]  <- -etlymt+ethmuy+eth*digeth-eth; ## lt 

      ## the second derivatives
    
      l2 <- matrix(0,n,3)
      ## order mm,mt,tt
     
      l2[,1]  <- -ethmuy; ## lmm
      l2[,2]  <-  eth-ethmuy; ## lmt
      eth2 <- eth^2;treth <- trigamma(eth)
      l2[,3]  <- etlymt-ethmuy-treth*eth2-eth*digeth+2*eth; #ltt
      ## need some link derivatives for derivative transform
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),family$linfo[[2]]$mu.eta(etat))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),family$linfo[[2]]$d2link(th))
    }

    l3 <- l4 <- g3 <- g4 <- 0 ## defaults

    if (deriv>1) {
      ## the third derivatives
      ## order mmm,mmt,mtt,ttt
      l3 <- matrix(0,n,4)
      l3[,1] <- ethmuy; ## lmmm
      l3[,2] <- ethmuy; ## lmmt
      l3[,3] <- ethmuy-eth; ## lmtt
      eth3 <- eth2*eth; g3eth <- psigamma(eth,deriv=2)
      l3[,4] <- -etlymt+ethmuy+g3eth*eth3+3*treth*eth2+eth*digeth-3*eth; ## lttt
      g3 <- cbind(family$linfo[[1]]$d3link(mu),family$linfo[[2]]$d3link(th))
    }

    if (deriv>3) {
      ## the fourth derivatives
      ## order mmmm,mmmt,mmtt,mttt,tttt
      l4 <- matrix(0,n,5)
      l4[,1] <- -ethmuy; ## lmmmm
      l4[,2] <- -ethmuy; ## lmmmt
      l4[,3] <- -ethmuy; ## lmmtt
      l4[,4] <-  eth-ethmuy; ## lmttt
      eth4 <- eth3*eth
      l4[,5] <- etlymt-ethmuy-psigamma(eth,deriv=3)*eth4-6*g3eth*eth3-
                7*treth*eth2-eth*digeth+4*eth; ## ltttt 
      g4 <- cbind(family$linfo[[1]]$d4link(mu),family$linfo[[2]]$d4link(th))
    }

    if (deriv) {
      i2 <- family$tri$i2; i3 <- family$tri$i3
      i4 <- family$tri$i4
   
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)

      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                      d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D,sandwich=sandwich)
      if (ncv) {
        ret$l1 <- de$l1; ret$l2 = de$l2; ret$l3 = de$l3
      }		      
    } else ret <- list()
    ret$l <- l; ret$l0 <- l0; ret
  } ## end ll gammals

  sandwich <- function(y,X,coef,wt,family,offset=NULL) {
  ## compute filling for sandwich estimate of cov matrix
    ll(y,X,coef,wt,family,offset=NULL,deriv=1,sandwich=TRUE)$lbb
  }

  initialize <- expression({
  ## regress X[,[jj[[1]]] on log(y) then X[,jj[[2]]] on log abs
  ## raw residuals.
  ## note that appropriate E scaling
  ## for full calculation may be inappropriate for initialization 
  ## which is basically penalizing something different here.
  ## best we can do here is to use E only as a regularizer.
      n <- rep(1, nobs)
      ## should E be used unscaled or not?..
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      if (is.null(start)) {
        jj <- attr(x,"lpi")
	if (!is.null(offset)) offset[[3]] <- 0
        yt1 <- log(y+max(y)*.Machine$double.eps^.75)
	if (!is.null(offset[[1]])) yt1 <- yt1 - offset[[1]]
        if (is.list(x)) { ## discrete case
	  start <- rep(0,max(unlist(jj)))
	  R <- suppressWarnings(mchol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,
	            qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[1]])+crossprod(E[,jj[[1]]])))
	  Xty <- XWyd(x$Xd,rep(1,length(y)),yt1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[1]])
          piv <- attr(R,"pivot")
	  rrank <- attr(R,"rank")
	  startji <- rep(0,ncol(R))
	  if (rrank<ncol(R)) {
            R <- R[1:rrank,1:rrank]
	    piv <- piv[1:rrank]
          }
          startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
	  startji[!is.finite(startji)] <- 0
	  start[jj[[1]]] <- startji
	  eta1 <- Xbd(x$Xd,start,k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,qc=x$qc,drop=x$drop,lt=x$lpid[[1]])
	  lres1 <- family$linfo[[2]]$linkfun(log(abs(y-family$linfo[[1]]$linkinv(eta1))))
	  if (!is.null(offset[[2]])) lres1 <- lres1 - offset[[2]]
	  R <- suppressWarnings(mchol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,
	            qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[2]])+crossprod(E[,jj[[2]]])))
	  Xty <- XWyd(x$Xd,rep(1,length(y)),lres1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[2]])
	  piv <- attr(R,"pivot");startji <- rep(0,ncol(R));rrank <- attr(R,"rank")
	  if (rrank<ncol(R)) {
            R <- R[1:rrank,1:rrank]
	    piv <- piv[1:rrank]
          }
	  startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
          start[jj[[2]]] <- startji
        } else { ## regular case
	  start <- rep(0,ncol(x))
	  x1 <- x[,jj[[1]],drop=FALSE]
          e1 <- E[,jj[[1]],drop=FALSE] ## square root of total penalty
          if (use.unscaled) {
            qrx <- qr(rbind(x1,e1))
            x1 <- rbind(x1,e1)
            startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
            startji[!is.finite(startji)] <- 0       
          } else startji <- pen.reg(x1,e1,yt1)
          start[jj[[1]]] <- startji
          lres1 <- family$linfo[[2]]$linkfun(log(abs(y-family$linfo[[1]]$linkinv(x[,jj[[1]],drop=FALSE]%*%start[jj[[1]]]))))
	  if (!is.null(offset[[2]])) lres1 <- lres1 - offset[[2]]
          x1 <-  x[,jj[[2]],drop=FALSE];e1 <- E[,jj[[2]],drop=FALSE]
          if (use.unscaled) {
            x1 <- rbind(x1,e1)
            startji <- qr.coef(qr(x1),c(lres1,rep(0,nrow(E))))   
            startji[!is.finite(startji)] <- 0
          } else startji <- pen.reg(x1,e1,lres1)
          start[jj[[2]]] <- startji
	}  
      }
  }) ## initialize gammals

  rd <- function(mu,wt,scale) {
    ## simulate responses
    phi <- exp(mu[,2])
    rgamma(nrow(mu),shape=1/phi,scale=mu[,1]*phi)
  } ## rd

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL) {
  ## optional function to give predicted values - idea is that 
  ## predict.gam(...,type="response") will use this, and that
  ## either eta will be provided, or {X, beta, off, Vb}. family$data
  ## contains any family specific extra information. 
  ## if se = FALSE returns one item list containing matrix otherwise 
  ## list of two matrices "fit" and "se.fit"... 

    if (is.null(eta)) {
      discrete <- is.list(X) 
      lpi <- attr(X,"lpi") 
      if (is.null(lpi)) {
        lpi <- list(1:ncol(X))
      }
      nobs <- if (discrete) nrow(X$kd) else nrow(X)
      eta <- matrix(0,nobs,2)
      ve <- matrix(0,nobs,2) ## variance of eta 
      for (i in 1:2) {
        if (discrete) {
	  eta[,i] <- Xbd(X$Xd,beta,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[i]])
        } else {
          Xi <- X[,lpi[[i]],drop=FALSE]
          eta[,i] <- Xi%*%beta[lpi[[i]]] ## ith linear predictor
        } 
        if (!is.null(off[[i]])) eta[,i] <- eta[,i] + off[[i]]
        if (se) ve[,i] <- if (discrete) diagXVXd(X$Xd,Vb,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,nthreads=1,lt=X$lpid[[i]]) else
	                  drop(pmax(0,rowSums((Xi%*%Vb[lpi[[i]],lpi[[i]]])*Xi)))
      }
    } else { 
      se <- FALSE
    }
    gamma <- cbind(exp(eta[,1]),family$linfo[[2]]$linkinv(eta[,2]))
   
    if (se) { ## need to loop to find se of probabilities...
      vp <- gamma
      vp[,1] <- abs(gamma[,1])*sqrt(ve[,1])
      vp[,2] <- abs(family$linfo[[2]]$mu.eta(eta[,2]))*sqrt(ve[,2])
      return(list(fit=gamma,se.fit=vp))
    } ## if se
    list(fit=gamma)
  } ## gammals predict

  structure(list(family="gammals",ll=ll,link=paste(link),nlp=2,ncv=ncv,sandwich=sandwich,
    tri = trind.generator(2), ## symmetric indices for accessing derivative arrays
    initialize=initialize,postproc=postproc,residuals=residuals,
    linfo = stats,rd=rd,predict=predict, ## link information list
    d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
    ls=1, ## signals that ls not needed here
    available.derivs = 2, ## can use full Newton here
    discrete.ok = TRUE
    ),class = c("general.family","extended.family","family"))
} ## end gammals




gumbls <- function(link=list("identity","log"),b=-7) {
## General family for Gumbel location scale model...
## parameterization is in terms of mean and log scale (beta).
## so log(mu) is mu1, scale = log(sigma) is mu2
## 1. get derivatives wrt mu, rho and xi.
## 2. get required link derivatives and tri indices.
## 3. transform derivs to derivs wrt eta (gamlss.etamu).
## 4. get the grad and Hessian etc for the model
##    via a call to gamlss.gH  
  
  ## first deal with links and their derivatives...
  if (length(link)!=2) stop("gumbls requires 2 links specified as character strings")
  okLinks <- list(c("identity"),c("identity","log"))
  stats <- list()
  for (i in 1:2) {
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- make.link(link[[i]]) else 
    stop(link[[i]]," link not available for mu parameter of gammals")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
           mu.eta=stats[[i]]$mu.eta),
           class="family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }
  if (link[[2]]=="log") { ## g^{-1}(eta) = b + log(1+exp(eta)) link
    stats[[2]]$valideta <- function(eta) TRUE
    stats[[2]]$linkfun <- eval(parse(text=paste("function(mu,b=",b,") {\n eta <- mub <- mu-b;\n",
      "ii <- mub < .Machine$double.eps;\n if (any(ii)) eta[ii] <- log(.Machine$double.eps);\n",
      "jj <- mub > -log(.Machine$double.eps);if (any(jj)) eta[jj] <- mub[jj];\n",
      "jj <- !jj & !ii;if (any(jj)) eta[jj] <- log(exp(mub[jj])-1);eta }")))

    stats[[2]]$mu.eta <- eval(parse(text=paste("function(eta,b=",b,") {\n",
      "ii <- eta < 0;eta <- exp(-eta*sign(eta))\n",
      "if (any(ii)) { ei <- eta[ii];eta[ii] <- ei/(1+ei)}\n",
      "ii <- !ii;if (any(ii)) eta[ii] <- 1/(1+eta[ii])\n",
      "eta }\n")))

    stats[[2]]$linkinv <- eval(parse(text=paste("function(eta,b=",b,") {\n",
      "mu <- eta;ii <- eta > -log(.Machine$double.eps)\n",
      "if (any(ii)) mu[ii] <- b + eta[ii]\n",
      "ii <- !ii;if (any(ii)) mu[ii] <- b + log(1 + exp(eta[ii]))\n",
      "mu }\n")))

    stats[[2]]$d2link <- eval(parse(text=paste("function(mu,b=",b,") {\n",
      "eta <- lb.linkfun(mu,b=b); ii <- eta > 0\n",
      "eta <- exp(-eta*sign(eta))\n",
      "if (any(ii)) { ei <- eta[ii];eta[ii] <- -(ei^2 + ei) }\n",
      "ii <- !ii;if (any(ii)) { ei <- eta[ii];eta[ii] <- -(1+ei)/ei^2 }\n",
      "eta }\n")))

    stats[[2]]$d3link <- eval(parse(text=paste("function(mu,b=",b,") {\n",
      "eta <- lb.linkfun(mu,b=b);ii <- eta > 0\n",
      "eta <- exp(-eta*sign(eta))\n",
      "if (any(ii)) { ei <- eta[ii]; eta[ii] <- (2*ei^2+ei)*(ei+1) }\n",
      "ii <- !ii;if (any(ii)) { ei <- eta[ii]; eta[ii] <- (2+ei)*(1+ei)/ei^3 }\n",
      "eta }\n")))
    
    stats[[2]]$d4link <- eval(parse(text=paste("function(mu,b=",b,") {\n",
      "eta <- lb.linkfun(mu,b=b);ii <- eta > 0\n",
      "eta <- exp(-eta*sign(eta))\n",
      "if (any(ii)) { ei <- eta[ii];eta[ii] <- -(6*ei^3+6*ei^2+ei)*(ei+1) }\n",
      "ii <- !ii;if (any(ii)) { ei <- eta[ii]; eta[ii] <- -(6+6*ei+ei^2)*(1+ei)/ei^4 }\n",
      "eta }\n")))
  }

  residuals <- function(object,type=c("deviance","pearson","response")) {
      mean <- object$fitted.values[,1]
      beta <- exp(object$fitted.values[,2])
      .euler <- 0.5772156649015328606065121 ## Euler's constant
      mu <- mean - beta * .euler
      sd <- pi*beta/sqrt(6)
      y <- object$y
      type <- match.arg(type)
      if (type=="deviance") {
        z <- (y-mu)/beta
        rsd <- 2*(z+exp(-z)-1)
        rsd <- sqrt(pmax(0,rsd))*sign(y-mu)
      } else if (type=="pearson") {
        rsd <- (y-mean)/sd
      } else {
        rsd <- y-mean
      }
      rsd
    } ## gumbls residuals
    
  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    ## It's difficult to define a sensible version of this that ensures
    ## that the data fall in the support of the null model, whilst being
    ## somehow equivalent to the full fit
    .euler <- 0.5772156649015328606065121 ## Euler's constant
    ## replace mean with expected value
    object$fitted.values[,1] <- object$fitted.values[,1] + exp(object$fitted.values[,2]) * .euler
    .my <- mean(object$y)
    object$null.deviance <- NA
  })

  ncv <- function(X,y,wt,nei,beta,family,llf,H=NULL,Hi=NULL,R=NULL,offset=NULL,dH=NULL,db=NULL,deriv=FALSE,nt=1) {
    gamlss.ncv(X,y,wt,nei,beta,family,llf,H=H,Hi=Hi,R=R,offset=offset,dH=dH,db=db,deriv=deriv,nt=nt)
  } ## ncv  


  ll <- function(y,X,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL,eta=NULL,ncv=FALSE,sandwich=FALSE) {
  ## function defining the gamlss gamma model log lik. 
  ## deriv: 0 - eval
  ##        1 - grad and Hess
  ##        2 - diagonal of first deriv of Hess
  ##        3 - first deriv of Hess
  ##        4 - everything.

    if (!is.null(offset)) offset[[3]] <- 0
    discrete <- is.list(X)
    jj <- attr(X,"lpi") ## extract linear predictor index
    if (is.null(eta)) {
      eta <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[1]]) else
                           X[,jj[[1]],drop=FALSE]%*%coef[jj[[1]]]
      if (!is.null(offset[[1]])) eta <- eta + offset[[1]] ## mu
    
      etab <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[2]]) else
                            X[,jj[[2]],drop=FALSE]%*%coef[jj[[2]]]
      if (!is.null(offset[[2]])) etab <- etab + offset[[2]]
    } else {
      etab <- eta[,2]
      eta <- eta[,1]
    }
    mu <- family$linfo[[1]]$linkinv(eta) ## mean
    beta <-  family$linfo[[2]]$linkinv(etab) ## log beta

    eb <- exp(-beta)
    z <- (y-mu)*eb
    ez <- exp(-z)
    
    l0 <- -beta - z - ez
    l <- sum(l0)
    
    n <- length(y)
     
    if (deriv>0) {
      lz <- ez - 1 ## e^{-z} - 1
      zm <- -eb   ## -e^{-b}
      zb <- -z
     
      l1 <- matrix(0,n,2)
      l1[,1] <- lz*zm    ## lm
      l1[,2] <- lz*zb - 1 ## lb

      lzz <- -ez; zmb <- eb; zbb <- z
      l2 <- matrix(0,n,3)
      ## order mm,mb,bb
      l2[,1] <- lzz*zm^2
      l2[,2] <- lzz*zm*zb + lz*zmb
      l2[,3] <- lzz*zb^2 + lz*zbb
      ## need some link derivatives for derivative transform
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),family$linfo[[2]]$mu.eta(etab))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),family$linfo[[2]]$d2link(beta))
    }

    l3 <- l4 <- g3 <- g4 <- 0 ## defaults

    if (deriv>1) {
      lzzz <- ez;zbbb <- -z;zmbb <- -eb
      ## the third derivatives
      ## order mmm,mmb,mbb,bbb
      l3 <- matrix(0,n,4)
      l3[,1] <- lzzz*zm^3; ## lmmm
      l3[,2] <- lzzz*zm^2*zb + 2*lzz*zm*zmb;   ## lmmb
      l3[,3] <- lzzz*zb^2*zm + 2*lzz*zb*zmb + lzz*zbb*zm + lz*zmbb; ## lmbb
      l3[,4] <- lzzz*zb^3 + 3*lzz*zb*zbb + lz*zbbb; ## lbbb
      g3 <- cbind(family$linfo[[1]]$d3link(mu),family$linfo[[2]]$d3link(beta))
    }

    if (deriv>3) {
      lzzzz <- -ez;zbbbb <- z;zmbbb <- eb
      ## the fourth derivatives
      ## order mmmm,mmmb,mmbb,mbbb,bbbb
      l4 <- matrix(0,n,5)
      l4[,1] <- lzzzz*zm^4; ## lmmmm
      l4[,2] <- lzzzz*zm^3*zb + 3*lzzz*zm^2*zmb ; ## lmmmb
      l4[,3] <- lzzzz*zm^2*zb^2 + 4*lzzz*zm*zb*zmb + lzzz*zm^2*zbb + 2*lzz*zmb^2 + 2*lzz*zm*zmbb; ## lmmbb
      l4[,4] <- lzzzz*zb^3*zm + 3*lzzz*zb^2*zmb + 3*lzzz*zm*zb*zbb + 3*lzz*zmb*zbb + 3*lzz*zb*zmbb + lzz*zm*zbbb + lz*zmbbb; ## lmbbb
      l4[,5] <- lzzzz*zb^4 + 6*lzzz*zb^2*zbb + 3*lzz*zbb^2 + 4*lzz*zb*zbbb + lz*zbbbb ## lbbbb 
      g4 <- cbind(family$linfo[[1]]$d4link(mu),family$linfo[[2]]$d4link(beta))
    }

    if (deriv) {
      i2 <- family$tri$i2; i3 <- family$tri$i3
      i4 <- family$tri$i4
   
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)

      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                      d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D,sandwich=sandwich)
      if (ncv) {
        ret$l1 <- de$l1; ret$l2 = de$l2; ret$l3 = de$l3
      }		      
    } else ret <- list()
    ret$l <- l;ret$l0 <- l0; ret
  } ## end ll gumbls

  sandwich <- function(y,X,coef,wt,family,offset=NULL) {
  ## compute filling for sandwich estimate of cov matrix
    ll(y,X,coef,wt,family,offset=NULL,deriv=1,sandwich=TRUE)$lbb
  }

  initialize <- expression({
  ## regress X[,[jj[[1]]] on y then X[,jj[[2]]] on
  ## log((y-mu)^2)/2 - 0.25 where mu is fit to y
  ## note that appropriate E scaling
  ## for full calculation may be inappropriate for initialization 
  ## which is basically penalizing something different here.
  ## best we can do here is to use E only as a regularizer.
      n <- rep(1, nobs)
      ## should E be used unscaled or not?..
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      if (is.null(start)) {
        jj <- attr(x,"lpi")
	if (!is.null(offset)) offset[[3]] <- 0
        yt1 <- y
	if (!is.null(offset[[1]])) yt1 <- yt1 - offset[[1]]
        if (is.list(x)) { ## discrete case
	  start <- rep(0,max(unlist(jj)))
	  R <- suppressWarnings(mchol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,
	            qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[1]])+crossprod(E[,jj[[1]]])))
	  Xty <- XWyd(x$Xd,rep(1,length(y)),yt1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[1]])
          piv <- attr(R,"pivot")
	  rrank <- attr(R,"rank")
	  startji <- rep(0,ncol(R))
	  if (rrank<ncol(R)) {
            R <- R[1:rrank,1:rrank]
	    piv <- piv[1:rrank]
          }
          startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
	  startji[!is.finite(startji)] <- 0
	  start[jj[[1]]] <- startji
	  eta1 <- Xbd(x$Xd,start,k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,qc=x$qc,drop=x$drop,lt=x$lpid[[1]])
	  lres1 <- family$linfo[[2]]$linkfun(log((y-family$linfo[[1]]$linkinv(eta1))^2)/2 - .25)
	  if (!is.null(offset[[2]])) lres1 <- lres1 - offset[[2]]
	  #lres1 <- family$linfo[[2]]$linkfun(lres1)
	  R1 <- suppressWarnings(mchol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,
	            qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[2]])+crossprod(E[,jj[[2]]])))
	  Xty <- XWyd(x$Xd,rep(1,length(y)),lres1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[2]])
	  piv1 <- attr(R1,"pivot");startji <- rep(0,ncol(R1));rrank <- attr(R1,"rank")
	  if (rrank<ncol(R1)) {
            R1 <- R1[1:rrank,1:rrank]
	    piv1 <- piv1[1:rrank]
          }
	  startji[piv1] <- backsolve(R1,forwardsolve(t(R1),Xty[piv1]))
          start[jj[[2]]] <- startji
	  ## pass 2 at mean param...
	  eta2 <- Xbd(x$Xd,start,k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,qc=x$qc,drop=x$drop,lt=x$lpid[[2]]) +
	          if (is.null(offset[[2]])) 0 else offset[[2]] 
	  yt1 <- yt1 - .57721*exp(family$linfo[[2]]$linkinv(eta2))
	  Xty <- XWyd(x$Xd,rep(1,length(y)),yt1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[1]])
	  startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
	  startji[!is.finite(startji)] <- 0
	  start[jj[[1]]] <- startji
        } else { ## regular case
	  start <- rep(0,ncol(x))
	  x1 <- x[,jj[[1]],drop=FALSE]
          e1 <- E[,jj[[1]],drop=FALSE] ## square root of total penalty
          if (use.unscaled) {
            qrx <- qr(rbind(x1,e1))
            startji <- qr.coef(qrx,c(yt1,rep(0,nrow(E))))
            startji[!is.finite(startji)] <- 0       
          } else startji <- pen.reg(x1,e1,yt1)
          start[jj[[1]]] <- startji
          lres1 <- family$linfo[[2]]$linkfun(log((y-family$linfo[[1]]$linkinv(x[,jj[[1]],drop=FALSE]%*%start[jj[[1]]]))^2)/2 - .25)
	  if (!is.null(offset[[2]])) lres1 <- lres1 - offset[[2]]
	  #lres1 <- family$linfo[[2]]$linkfun(lres1)
          x2 <-  x[,jj[[2]],drop=FALSE];e2 <- E[,jj[[2]],drop=FALSE]
          if (use.unscaled) {
            x2 <- rbind(x2,e2)
            startji <- qr.coef(qr(x2),c(lres1,rep(0,nrow(E))))   
            startji[!is.finite(startji)] <- 0
          } else startji <- pen.reg(x2,e2,lres1)
          start[jj[[2]]] <- startji
	  ## pass 2 at mean parameter...
	  eta2 <- x[,jj[[2]],drop=FALSE]%*%start[jj[[2]]] + if (is.null(offset[[2]])) 0 else offset[[2]] 
	  yt1 <- yt1 - .57721*exp(family$linfo[[2]]$linkinv(eta2))
	  if (use.unscaled) {
            startji <- qr.coef(qrx,c(yt1,rep(0,nrow(E))))
            startji[!is.finite(startji)] <- 0       
          } else startji <- pen.reg(x1,e1,yt1)
	  start[jj[[1]]] <- startji
	}  
      }
  }) ## initialize gumbls

  rd <- function(mu,wt,scale) {
    ## simulate Gumbel responses from quantile function
    .euler <- 0.5772156649015328606065121 ## Euler's constant
    u <- runif(nrow(mu))
    beta <- exp(mu[,2])
    ## assumes mu[,1] is mean not mu param
    mu[,1] - beta*(.euler +log(-log(u)))
  } ## rd

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL) {
  ## optional function to give predicted values - idea is that 
  ## predict.gam(...,type="response") will use this, and that
  ## either eta will be provided, or {X, beta, off, Vb}. family$data
  ## contains any family specific extra information. 
  ## if se = FALSE returns one item list containing matrix otherwise 
  ## list of two matrices "fit" and "se.fit"... 

    if (is.null(eta)) {
      discrete <- is.list(X) 
      lpi <- attr(X,"lpi") 
      if (is.null(lpi)) {
        lpi <- list(1:ncol(X))
      }
      nobs <- if (discrete) nrow(X$kd) else nrow(X)
      eta <- matrix(0,nobs,2)
      ve <- matrix(0,nobs,2) ## variance of eta 
      for (i in 1:2) {
        if (discrete) {
	  eta[,i] <- Xbd(X$Xd,beta,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[i]])
        } else {
          Xi <- X[,lpi[[i]],drop=FALSE]
          eta[,i] <- Xi%*%beta[lpi[[i]]] ## ith linear predictor
        } 
        if (!is.null(off[[i]])) eta[,i] <- eta[,i] + off[[i]]
        if (se) ve[,i] <- if (discrete) diagXVXd(X$Xd,Vb,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,nthreads=1,lt=X$lpid[[i]]) else
	                  drop(pmax(0,rowSums((Xi%*%Vb[lpi[[i]],lpi[[i]]])*Xi)))
      }
    } else { 
      se <- FALSE
    }
    gamma <- cbind(eta[,1],family$linfo[[2]]$linkinv(eta[,2]))
   
    if (se) { ## need to loop to find se of probabilities...
      vp <- gamma
      vp[,1] <- sqrt(ve[,1])
      vp[,2] <- abs(family$linfo[[2]]$mu.eta(eta[,2]))*sqrt(ve[2])
      return(list(fit=gamma,se.fit=vp))
    } ## if se
    list(fit=gamma)
  } ## gumbls predict

  structure(list(family="gumbls",ll=ll,link=paste(link),nlp=2,ncv=ncv,sandwich=sandwich,
    tri = trind.generator(2), ## symmetric indices for accessing derivative arrays
    initialize=initialize,postproc=postproc,residuals=residuals,
    linfo = stats,rd=rd,predict=predict, ## link information list
    d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
    ls=1, ## signals that ls not needed here
    available.derivs = 2, ## can use full Newton here
    discrete.ok = TRUE
    ),class = c("general.family","extended.family","family"))
} ## end gumbls


## shash - Matteo Fasiolo 2020

shash <- function(link = list("identity", "logeb", "identity", "identity"), b = 1e-2, phiPen = 1e-3) 
{ 
  sech <- function(.x){ 1 / cosh(.x) }

  npar <- 4
  if (length(link) != npar) stop("shash requires 4 links specified as character strings")
  okLinks <- list("identity", "logeb", "identity", "identity")
  stats <- list()
  param.names <- c("mu", "tau", "eps", "phi")
  for (i in c(1, 3, 4)) { # Links for mu, eps and phi
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
  
  # variance <- function(mu) exp(get(".Theta"))  ##### XXX ##### Necessary?
  
  # validmu <- function(mu) all( is.finite(mu) )
  
  residuals <- function(object, type = c("deviance", "response")) {
    
    mu <-  object$fitted[ , 1, drop = TRUE]
    tau <- object$fitted[ , 2, drop = TRUE]
    eps <- object$fitted[ , 3, drop = TRUE]
    phi <- object$fitted[ , 4, drop = TRUE]
    
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
      
      l <- - tau - 0.5*log(2*pi) + log(CC) - 0.5*log1p(z^2) - 0.5*SS^2
      
      # By putting ls to zero we are using only log-likelihood
      ls <- 0
      
      rsd <- pmax(0, 2*(ls - l))
      
      rsd <- sqrt(rsd)*sgn
    }
    rsd
  } ## shash residuals

  ncv <- function(X,y,wt,nei,beta,family,llf,H=NULL,Hi=NULL,R=NULL,offset=NULL,dH=NULL,db=NULL,deriv=FALSE,nt=1) {
    gamlss.ncv(X,y,wt,nei,beta,family,llf,H=H,Hi=Hi,R=R,offset=offset,dH=dH,db=db,deriv=deriv,nt=nt)
  } ## ncv  


  ll <- function(y, X, coef, wt, family, offset = NULL, deriv=0, d1b=0, d2b=0, Hp=NULL, rank=0, fh=NULL, D=NULL,
                 eta=NULL,ncv=FALSE,sandwich=FALSE) {
    ## function defining the shash model log lik. 
    ## deriv: 0 - eval
    ##        1 - grad and Hess
    ##        2 - diagonal of first deriv of Hess
    ##        3 - first deriv of Hess
    ##        4 - everything.
    
    ##### We need some helper functions
    # Calculates \code{log(1+exp(x))} in a numerically stable fashion.
    .log1pexp <- function(x)
    {
      indx <- .bincode(x, c(-Inf, -37, 18, 33.3, Inf), T)
      kk <- which(indx==1)
      if( length(kk) ){  x[kk] <- exp(x[kk])  }
      kk <- which(indx==2)
      if( length(kk) ){  x[kk] <- log1p( exp(x[kk]) ) }
      kk <- which(indx==3)
      if( length(kk) ){  x[kk] <- x[kk] + exp(-x[kk]) }
      return(x)
    }
    
    # Compute sqrt(x^2 + m) when |x| >> 0 and m is reasonably small (e.g. + 1 or - 1)
    .sqrtX2pm <- function(x, m){ 
      x <- abs(x)
      kk <- which( x < 1e8 )
      if( length(kk) ){
        x[kk] <- sqrt(x[kk]^2 + m)
      }
      return(x)
    }
    
    # Compute (a*x^2 + m1) / (x^2 + m2)^2 when |x| >> 0 and m1, m2 are reasonably small (e.g. + 1 or - 1)
    .ax2m1DivX2m2SQ <- function(x, m1, m2, a = 1){
      if(a < 0){ stop("'a' has to be positive")  }
      x <- abs(x)
      kk <- (a * x^2 + m1) < 0
      o <- x * 0
      if( any(kk) ){
        o[kk] <- (a * x[kk]^2 + m1) / (x[kk]^2 + m2)^2
      }
      if( sum(kk) < length(x) ){
        o[!kk] <- ((.sqrtX2pm(sqrt(a)*x[!kk], m1) / .sqrtX2pm(x[!kk], m2)) / .sqrtX2pm(x[!kk], m2))^2
      }
      return(o)
    }
  
    # If offset is not null or a vector of zeros, give an error
    if( !is.null(offset[[1]]) && sum(abs(offset)) )  stop("offset not still available for this family")
    
    jj <- attr(X, "lpi") ## extract linear predictor index
    
    npar <- 4
    n <- length(y)
    if (is.null(eta)) {
      eta <-  drop( X[ , jj[[1]], drop=FALSE] %*% coef[jj[[1]]] )
      eta1 <- drop( X[ , jj[[2]], drop=FALSE] %*% coef[jj[[2]]] )
      eta2 <- drop( X[ , jj[[3]], drop=FALSE] %*% coef[jj[[3]]] )
      eta3 <- drop( X[ , jj[[4]], drop=FALSE] %*% coef[jj[[4]]] )
    } else {
      eta1 <- eta[,2]
      eta2 <- eta[,3]
      eta3 <- eta[,4]
      eta <- eta[,1]
    }
    mu <-  family$linfo[[1]]$linkinv( eta )
    tau <- family$linfo[[2]]$linkinv( eta1 )
    eps <- family$linfo[[3]]$linkinv( eta2 )
    phi <- family$linfo[[4]]$linkinv( eta3 )
    
    sig <- exp( tau )
    del <- exp( phi )
    
    z <- (y - mu) / (sig*del)
    
    dTasMe <- del*asinh(z) - eps
    g <- -dTasMe
    CC <- cosh( dTasMe )
    SS <- sinh( dTasMe )
    
    l0 <-  - tau - 0.5*log(2*pi) + log(CC) - 0.5*.log1pexp(2*log(abs(z))) - 0.5*SS^2 - phiPen*phi^2 
    l <- sum(l0)
    
    if (deriv>0) {
      
      zsd <- z*sig*del
      sSp1 <- .sqrtX2pm(z, 1) # sqrt(z^2+1)
      asinhZ <- asinh(z)
      
      ## First derivatives 
      De <- tanh(g) - 0.5*sinh(2*g)
      Dm <- 1/(del*sig*sSp1)*(del*(De)+z/sSp1)
      Dt <- zsd*Dm - 1
      Dp <- Dt + 1 - del*asinhZ*De - 2*phiPen*phi
      
      L1 <- cbind(Dm,Dt,De,Dp)
      
      ## the second derivatives  
      Dme <- (sech(g)^2 - cosh(2*g)) / (sig*sSp1)
      Dte <- zsd*Dme
      Dmm <- Dme/(sig*sSp1) + z*De/(sig^2*del*sSp1^3) + .ax2m1DivX2m2SQ(z, -1, 1)/(del*sig*del*sig)
      Dmt <- zsd*Dmm - Dm
      Dee <- -2*cosh(g)^2 + sech(g)^2 + 1 
      Dtt <-  zsd*Dmt
      Dep <- Dte - del*asinhZ*Dee
      Dmp <- Dmt + De/(sig*sSp1) - del*asinhZ*Dme
      Dtp <- zsd*Dmp
      Dpp <- Dtp - del*asinhZ*Dep + del*(z/sSp1-asinhZ)*De - 2*phiPen
      
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
      Dmme <- Dmee/(sig*sSp1) + z*Dee/(sig*sig*del*sSp1^3)
      Dmmm <- 2*z*Dme/(sig*sig*del*sSp1^3) + Dmme/(sig*sSp1) + 
        .ax2m1DivX2m2SQ(z, -1, 1, 2)*De/(sig^3*del^2*sSp1) + 
        2*(z/sSp1)*.ax2m1DivX2m2SQ(z, -3, 1)/((sig*del)^3*sSp1)
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
                       d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D,sandwich=sandwich) 
      if (ncv) {
        ret$l1 <- de$l1; ret$l2 = de$l2; ret$l3 = de$l3
      }
    } else ret <- list()
    ret$l <- l;ret$l0 <- l0; ret
  } ## end ll

  sandwich <- function(y,X,coef,wt,family,offset=NULL) {
  ## compute filling for sandwich estimate of cov matrix
    ll(y,X,coef,wt,family,offset=NULL,deriv=1,sandwich=TRUE)$lbb
  }

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
        x1 <- rbind(x1, e1)
        startMu <- qr.coef(qr(x1), c(yt1,rep(0,nrow(E))))
        startMu[ !is.finite(startMu) ] <- 0       
      } else { startMu <- pen.reg(x1, e1, yt1) }
      start[jj[[1]]] <- startMu
      
      # 2) Ridge regression using log absolute residuals
      lres1 <- log( abs(y-drop(family$linfo[[1]]$linkinv(x[,jj[[1]],drop=FALSE]%*%startMu))) )
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
      x1 <-  x[ , jj[[3]],drop=FALSE]
      startji <- qr.coef(qr(x1), c(rep(family$linfo[[3]]$linkfun(0),nrow(x1))))   
      startji[!is.finite(startji)] <- 0
      start[jj[[3]]] <- startji
      
      x1 <-  x[ , jj[[4]],drop=FALSE]
      startji <- qr.coef(qr(x1), c(rep(family$linfo[[4]]$linkfun(0),nrow(x1))))   
      startji[!is.finite(startji)] <- 0
      start[jj[[4]]] <- startji
      
      }
  }) ## initialize 
  
  rd <- function(mu, wt, scale) { 
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }
    muE <- mu[ , 1, drop = TRUE]
    sigE <- exp(mu[ , 2, drop = TRUE])
    epsE <- mu[ , 3, drop = TRUE]
    delE <- exp(mu[ , 4, drop = TRUE])
    n <- length(muE)
    
    .r <- muE + (delE * sigE) * sinh((1/delE) * asinh(qnorm(runif(n))) + (epsE/delE))
    
    return( .r )
  }
  
  qf <- function(p, mu, wt, scale) {
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }
    muE <- mu[ , 1, drop = TRUE]
    sigE <- exp(mu[ , 2, drop = TRUE])
    epsE <- mu[ , 3, drop = TRUE]
    delE <- exp(mu[ , 4, drop = TRUE])
    
    q <- muE + (delE * sigE) * sinh((1/delE) * asinh(qnorm(p)) + (epsE/delE))
    
    return( q)
  }
  
  cdf <- function(q, mu, wt, scale, logp) {
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }
    muE <- mu[ , 1, drop = TRUE]
    sigE <- exp(mu[ , 2, drop = TRUE])
    epsE <- mu[ , 3, drop = TRUE]
    delE <- exp(mu[ , 4, drop = TRUE])
    
    p <- pnorm( sinh((asinh( (q-muE)/(delE * sigE) )  - epsE/delE) * delE), log.p = logp )
    
    return( p )
  }
  
  
  structure(list(family="shash",ll=ll, link=paste(link), nlp=npar,ncv=ncv,sandwich=sandwich,
                 tri = trind.generator(npar), ## symmetric indices for accessing derivative arrays
                 initialize=initialize,
                 #postproc=postproc,
                 residuals=residuals,
                 rd=rd,
                 qf=qf,
                 cdf=cdf,
                 #predict=predict,
                 linfo = stats, ## link information list
                 d2link=1, d3link=1, d4link=1, ## signals to fix.family.link that all done    
                 ls=1, ## signals that ls not needed here
                 available.derivs = 2 ## can use full Newton here
  ),class = c("general.family","extended.family","family"))
} ## shash


