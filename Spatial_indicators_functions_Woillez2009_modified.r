"dg2nm" <-
function(x, y = NA, modproj, mlong, mlat)
{
#===============================================================================
# PROJECTION from decimal degrees to nautical miles
#
# Routine from Geostatistics for Estimating Fish Abundance (GEFA)
# & EU program Fisboat, DG-Fish, STREP n? 502572
# Authors : M.Woillez (Mines-ParisTech), N.Bez (IRD) 
#           and J.Rivoirard (Mines-ParisTech)
# Last update : 01 march 2008
#
# Argument:
# x, y     2 vectors of same length or a list with x and y.
# modproj  Type of projection
#          if modproj = "": x and y are NOT changed.
#          if modproj = "mean": longitudes are modified by the same cosine equal
#             to the cosine of the mean latitude of y.
#          if modproj = 0.3: longitudes are modified by the same given cosine.
#          if modproj = "cosine": each longitude is modified according to 
#             its latitude.
# mlong    mean longitude in DEGREES of the data set to be transformed
# mlat     mean latitude in DEGREES of the data set to be transformed
#
#===============================================================================

  miss <- function(x){
	  length(x) == 1 && is.na(x)
  }
  
	if(is.list(x)) {
		y <- x$y
		x <- x$x
	}
	if(!miss(modproj)) {
		x <- x - mlong
		y <- y - mlat
		if(modproj == "mean") {
			x <- x * 60 * cos((mlat * pi)/180)
			y <- y * 60
		}
		else if(is.numeric(modproj)) {
			x <- x * 60 * modproj
			y <- y * 60
		}
		else if(modproj == "cosine") {
			x <- x * 60 * cos((y * pi)/180)
			y <- y * 60
		}
	}
	list(x = x, y = y)
}

"nm2dg" <-
function(x, y = NA, modproj, mlong, mlat)
{
#===============================================================================
# PROJECTION from nautical miles to decimal degrees
#
# Routine from Geostatistics for Estimating Fish Abundance (GEFA)
# & EU program Fisboat, DG-Fish, STREP n? 502572
# Authors : M.Woillez (Mines-ParisTech), N.Bez (IRD) 
#           and J.Rivoirard (Mines-ParisTech)
# Last update : 01 march 2008 
#
# Argument:
# x, y     2 vectors of same length or a list with x and y.
# modproj  if modproj = "": x and y are NOT changed.
#          if modproj = "mean": longitudes are modified by the same cosine equal
#                         to the cosine of the mean latitude of y.
#          if modproj = 0.3: longitudes are modified by the same given cosine.
#          if modproj = "cosine": each longitude is modified according to 
#                         its latitude.
# mlong    mean longitude in DEGREES of the data set to be transformed
# mlat     mean latitude in DEGREES of the data set to be transformed
#
#===============================================================================

  miss <- function(x){
	  length(x) == 1 && is.na(x)
  }
  
	if(is.list(x)) {
		y <- x$y
		x <- x$x
	}
  if(!miss(modproj)) {
  	if(modproj == "mean") {
	  	y <- y/60
		  x <- x/(60 * cos((mlat * pi)/180))
  	}
	  else if(is.numeric(modproj)) {
		  x <- x/(60 * modproj)
  		y <- y/60
	  }
  	else if(modproj == "cosine") {
	  	y <- y/60
		  x <- x/(60 * cos(((y + mlat) * pi)/180))
  	}
  }
	x <- x + mlong
	y <- y + mlat
	list(x = x, y = y)
}

"infl" <- 
function (x, y, z, dlim = NA, pol = NA, ndisc = NA, 
    visu = F, opt = 0, mode = 0) 
{
#===============================================================================
# AREAS OF INFLUENCE
#
# Routine from Geostatistics for Estimating Fish Abundance (GEFA)
# & EU program Fisboat, DG-Fish, STREP n? 502572
# Authors : M.Woillez (Mines-ParisTech), N.Bez (IRD) 
#           and J.Rivoirard (Mines-ParisTech)
# Last update : 01 march 2008  
#
# Arguments:
# x      The x-coordinate of the first population
#        Can be a vector, a matrix or an xyZ list (-,|,[]).
# y      The y-coordinates of the first population
# z      The regionalised variable of the first population.
#        If missing, the results of 'cgi' will concern the samples only. 
# dlim   Limit distance for valuating a node. Default is half the
#        diagonal of the square with contains all the selected data
#        points. GIVEN IN THE PROJECTED SPACE UNITS. Can be a vector
#        with ONE OR TWO components. 
#        When the areas of influence is a circle (opt==0), its
#        radius is given by dlim=c(1)
#        When the areas of influence is a rectangle (opt==1), its (half-)
#        dimensions are given by dlim=c(1,2) 
# pol    Optional polygon definition
# ndisc  Discretization of the data area when constructing the
#        defaulted grid. 
# visu   When TRUE, the results are demonstrated graphically. This
#        plot is only performed in the projected space. 
# opt    When equal to 0, the influence areas is defined as a
#        circle; if different from 0, the areas corresponds to a
#        rectangle 
# mode   When equal to 0, each grid node contains the value of the
#        closest data. When equal to 1, each grid node contains the
#        rank of the closest data. When equal to 2, each grid node
#        contains the areas of the influence polygon to which it
#        belongs. 
#
#===============================================================================

    miss <- function(x){
            length(x) == 1 && is.na(x)
    }

    extend <- rev(sort(c((2 * dlim)/diff(range(x)), (2 * dlim)/diff(range(y)))))[1]
    nd <- length(x)
    if (length(dlim) == 1) 
        dlim <- rep(dlim, 2)
    else opt <- 1

    # Build the grid #

    xrange <- range(x)
    yrange <- range(y)
    dx <- xrange[2] - xrange[1]
    dy <- yrange[2] - yrange[1]
    xmin <- xrange[1] - extend * dx
    xmax <- xrange[2] + extend * dx
    ymin <- yrange[1] - extend * dy
    ymax <- yrange[2] + extend * dy
    nx <- ndisc
    ny <- ndisc
    dx <- (xmax - xmin) / nx
    dy <- (ymax - ymin) / ny
    ng <- nx * ny
    maille <- dx * dy 
  
    xg <- rep(xmin + dx * seq(0, nx - 1), ny)
    yg <- rep(ymin + dy * seq(0, ny - 1), rep(nx, ny))
  
    if(!miss(pol)){
      sel.pol <- inout(data.frame(xg,yg),data.frame(pol$x,pol$y))
      xg[!sel.pol] <- NA
      yg[!sel.pol] <- NA
    }

    # Affectations #

    argin <- function (x, test = 1.234e+30){
        if (length(x) > 0) 
            x[is.na(x)] <- test
        x
    }
    
    argout <- function (x, test = 1.234e+30){
        if (length(x) > 0) 
            x[x == test] <- NA
        x
    }
    
    x <- argin(x)
    y <- argin(y)
    z <- argin(z)
    xg <- argin(xg)
    yg <- argin(yg)
    dlim <- argin(dlim)

    # Calculate the data extension #

    xxmin <- min(x[x!=argin(NA)])
    xxmax <- max(x[x!=argin(NA)])
    deltax <- xxmax - xxmin
    yymin <- min(y[y!=argin(NA)])
    yymax <- max(y[y!=argin(NA)])
    deltay <- yymax - yymin
    if(dlim[1]==argin(NA)) dlim[1] <- deltax * extend
    if(dlim[2]==argin(NA)) dlim[2] <- deltax * extend
    if(opt==0) dmax <- sqrt(dlim[1] * dlim[2])
    if(opt==1) dmax <- sqrt(dlim[1] * dlim[1] + dlim[2] * dlim[2])

    # Initializations #

    surf <- rep(0,nd)
    zg <- rep(argin(NA),ng)

    # Build matrix of data points and of target points #
 
    xd <- matrix(x,ng,nd,byrow=T)
    yd <- matrix(y,ng,nd,byrow=T)
    rkd <- matrix(1:nd,ng,nd,byrow=T)
    xxg <- matrix(xg,ng,1)
    yyg <- matrix(yg,ng,1)
    rkg <- matrix(1:ng,ng,nd,byrow=F)
    
    # Compute distance matrix #    
    
    dxx <- abs(sweep(xd,1,xxg,"-"))
    dyy <- abs(sweep(yd,1,yyg,"-"))   
    dist <- sqrt(dxx^2 + dyy^2)
    misval.dist <- sqrt((1.234e+30)^2+(1.234e+30)^2)
    
    # Find the closest data point of an active target point according to options #

    valid.dg.f<-function(d,g,dist,dxx,dyy){
        valid <- F
        if((dxx[g,d]<=dlim[1] & dyy[g,d]<=dlim[2])==T) valid <- T
        return(valid) } 

    valid.g.f<-function(g,dist,dxx,dyy){
        sapply(1:nd,valid.dg.f,g,dist,dxx,dyy)} 
    
    valid.f<-function(dist,dxx,dyy){
        t(sapply(1:ng,valid.g.f,dist,dxx,dyy))}

    if(opt==0)
        cond <- dist==apply(dist,1,min) & dist!=misval.dist & dist<=dmax
    if(opt==1){
        dist[!valid.f(dist, dxx, dyy)] <- misval.dist
        cond <- dist==apply(dist,1,min) & dist!=misval.dist
    }
    proxd <- rkd[cond]
    proxg <- rkg[cond]   
    
    # Count nodes attributed to data points #

    nbn <- as.data.frame(table(proxd))

    # Scale the count by the mesh unit to get areas of influence #
    
    surf[as.numeric(levels(nbn$proxd))] <- nbn$Freq * maille

    # Evaluate the influence areas at each target #

    if(mode==0) zg[proxg] <- z[proxd]
    if(mode==1) zg[proxg] <- proxd
    if(mode==2) zg[proxg] <- surf[proxd]       

    zg <- argout(zg)
    surf <- argout(surf)
    x <- argout(x)
    y <- argout(y)
    z <- argout(z)
    
    zg <- matrix(zg, nrow = nx, ncol = ny, byrow = F)
    xg <- xmin + dx * seq(0, (nx - 1))
    yg <- ymin + dy * seq(0, (ny - 1))
    
    miss <- function(x){
          length(x) == 1 && is.na(x)
    }
    
    if (visu) {
        image(xg, yg, zg, col = c("blue", rev(rainbow(abs(diff(range(z))), start=0, end=1/6))))
        symbols(x, y, sqrt(surf), fg = 3, inches = 0.2, add = T)
        if (!miss(pol)) 
            lines(pol)
    }
    surf
}

"abundance" <-
function(z, w)
{
#===============================================================================
# ABUNDANCE
#
# Routine from EU program Fisboat, DG-Fish, STREP n? 502572
# Authors : M.Woillez and J.Rivoirard (Mines-ParisTech)
# Last update : 01 march 2008 
#
# Argument:
# z     variable of interest (i.e. fish density)
# w     appropriate areas of influence set as weighted factors
#
#===============================================================================

    AB <- sum(z*w,na.rm=T)
    AB
}

"cgi" <-
function(x = long, y = lat, z = NA, w = NA, modproj = NA, mlong = NA, 
	mlat = NA, col = 1, plot = F)
{
#===============================================================================
# CENTER OF GRAVITY, INERTIA AND ISOTROPY
#
# Routine from Geostatistics for Estimating Fish Abundance (GEFA)
# & EU program Fisboat, DG-Fish, STREP n? 502572
# Authors : M.Woillez (Mines-ParisTech), N.Bez (IRD) 
#           and J.Rivoirard (Mines-ParisTech)
# Last update : 01 march 2008 
#
# Argument:
#	x	      The x-coordinate (MUST be a vector).
#	y	      The y-coordinates (MUST be a vector).
#	z	      The regionalised variable in 2d (MUST be a vector). 
#         If missing, the results of 'cgi' will concern the samples only.
# w	      Optional. A weight or a area of influence. Set to 1 if missing
#	modproj	Optional. Indicates the type of projection to perform.
# mlong   mean longitude in DEGREES of the data set to be transformed
# mlat    mean latitude in DEGREES of the data set to be transformed
#	        See 'dg2nm' for precisions.
#	col	    Color for representing the axes.
#	plot	  If plot=T the principal axes of the inertia are automatically 
#		      plotted on an ALREADY EXISTING figure.
#
#	The output consists in a list with :
#	xcg, ycg	    the coordinates of the center of gravity of z
#	I	            the value of the inertia of z around its center of gravity
# Imax          the value of the inertia of z according to the first principal 
#               axes of the inertia
# Imin          the value of the inertia of z according to the second principal 
#               axes of the inertia
#	Iso           the value of the isotropy of z
# xaxe1, yaxe1  the coordinates of the first principal axes of the inertia of z
# xaxe2, yaxe2	the coordinates of the second principal axes of the inertia of z
#
#===============================================================================

  miss <- function(x){
	  length(x) == 1 && is.na(x)
  }

	if(miss(z))
		z <- rep(1, length(x))
	if(miss(w))
		w <- rep(1, length(x))
  sel <- !is.na(x * y * z * w)
	x <- x[sel]
	y <- y[sel]
	z <- z[sel]
	w <- w[sel]
	if(length(x[!is.na(x)]) > 0) {
		if(!miss(modproj)) {
		  bid <- dg2nm(x = x, y = y, modproj = modproj, mlong = mlong, mlat = mlat)
			x <- bid$x
			y <- bid$y
		}
    # Center of gravity coordinates
		xg <- sum(x * z * w)/sum(z * w)
		yg <- sum(y * z * w)/sum(z * w)
		
	  # Inertia
    dx <- x - xg
		dy <- y - yg
		d <- sqrt(dx^2 + dy^2)
		inert <- sum(z * w * (d^2))/sum(z * w)
		I <- inert	
	  
    # Weigthed PCA 
    if(!is.na(I)) {
			M11 <- sum(dx^2 * z * w)
			M22 <- sum(dy^2 * z * w)
			M21 <- sum(dx * dy * z * w)
			M12 <- M21
			M <- matrix(c(M11, M12, M21, M22), byrow = T, ncol = 2)
			x1 <- eigen(M)$vectors[1, 1]
			y1 <- eigen(M)$vectors[2, 1]
			x2 <- eigen(M)$vectors[1, 2]
			y2 <- eigen(M)$vectors[2, 2]
			r1 <- eigen(M)$values[1]/(eigen(M)$values[1] + eigen(M)$values[2])
      	
	    # Principal axis coordinates
			e1 <- (y1/x1)^2
			sx1 <- x1/abs(x1)
			sy1 <- y1/abs(y1)
			sx2 <- x2/abs(x2)
			sy2 <- y2/abs(y2)
			xa <- xg + sx1 * sqrt((r1 * inert)/(1 + e1))
			ya <- yg + sy1 * sqrt((r1 * inert)/(1 + (1/e1)))
			xb <- 2 * xg - xa
			yb <- 2 * yg - ya
			xc <- xg + sx2 * sqrt(((1 - r1) * inert)/(1 + (1/e1)))
			yc <- yg + sy2 * sqrt(((1 - r1) * inert)/(1 + e1))
			xd <- 2 * xg - xc
			yd <- 2 * yg - yc
			Imax <- r1*inert 
      Imin <- (1-r1)*inert
      Iso <- sqrt(Imin/Imax)
		}
		else {
			xa <- NA
			ya <- NA
			xb <- NA
			yb <- NA
			xc <- NA
			yc <- NA
			xd <- NA
			yd <- NA
			Imax <- NA
			Imin <- NA
			Iso <- NA
		}
		if(!miss(modproj)) {
			bid <- nm2dg(x = c(xg, xa, xb, xc, xd), y = c(yg, ya, yb, yc, yd), 
        modproj = modproj, mlong = mlong, mlat = mlat)
			res <- list(xcg = bid$x[1], ycg = bid$y[1], I = I, Imax = Imax, 
        Imin = Imin, Iso = Iso, xaxe1 = bid$x[2:3], yaxe1 = bid$y[2:3], 
        xaxe2 = bid$x[4:5],	yaxe2 = bid$y[4:5])
		}
		else res <- list(xcg = xg, ycg = yg, I = I, Imax = Imax, Imin = Imin, 
      Iso = Iso, xaxe1 = c(xa, xb), yaxe1 = c(ya, yb), xaxe2 = c(xc, xd), 
      yaxe2 = c(yc, yd))
		if(plot == T) {
			segments(res$xaxe1[1], res$yaxe1[1], res$xaxe1[2], res$yaxe1[2], col = col)
			segments(res$xaxe2[1], res$yaxe2[1], res$xaxe2[2], res$yaxe2[2], col = col)
		}
	}
	else {
		res <- list(xcg = NA, ycg = NA, I = NA, Imax = NA, 
      Imin = NA, Iso = NA, xaxe1 = NA, yaxe1 = NA, xaxe2 = NA, yaxe2 = NA)
	}
	res
}

"gic" <-
function (x1, y1, z1, w1 = NA, x2, y2, z2, w2 = NA, modproj = NA, mlong = NA, 
	mlat = NA)
{
#===============================================================================
# GLOBAL INDEX OF COLLOCATION
#
# Routine from Geostatistics for Estimating Fish Abundance (GEFA)
# & EU program Fisboat, DG-Fish, STREP n? 502572
# Authors : M.Woillez (Mines-ParisTech), N.Bez (IRD) 
#           and J.Rivoirard (Mines-ParisTech)
# Last update : 01 march 2008  
#
# Arguments:
# x1	    The x-coordinate of the first population
#	        Can be a vector, a matrix or an xyZ list (-,|,[]).
# y1	    The y-coordinates of the first population
# z1	    The regionalised variable of the first population.
#	        If missing, the results of 'cgi' will concern the samples only.
# w1	    Optional. A weight or an area of influence of the first population.
#	        Set to 1 if missing
# x2	    The x-coordinate of the second population
#	        Can be a vector, a matrix or an xyZ list (-,|,[]).
# y2	    The y-coordinates of the second population
# z2	    The regionalised variable of the second population.
#	        If missing, the results of 'cgi' will concern the samples only.
# w2	    Optional. A weight or an area of influence of the second population.
#	        Set to 1 if missing
# modproj	Optional. Indicates the type of projection to perform.
# mlong   mean longitude in DEGREES of the data set to be transformed
# mlat    mean latitude in DEGREES of the data set to be transformed
#	        See 'dg2nm' for precisions.    
#
#===============================================================================

    # Compute the centers of gravity and inertia of the two populations
    popZ1 <- cgi(x1, y1, z1, w1, modproj=modproj, mlong=mlong, mlat=mlat, plot=F)
    popZ2 <- cgi(x2, y2, z2, w2, modproj=modproj, mlong=mlong, mlat=mlat, plot=F)
    
    # Perform the projection
    Z1 <- dg2nm(x=popZ1$xcg, y=popZ1$ycg, modproj=modproj, mlong=mlong, mlat=mlat)
    Z2 <- dg2nm(x=popZ2$xcg, y=popZ2$ycg, modproj=modproj, mlong=mlong, mlat=mlat)   
    
    # Compute the 'GIC' index
    GIC <- (((Z1$x-Z2$x)^2+(Z1$y-Z2$y)^2) / (((Z1$x-Z2$x)^2+(Z1$y-Z2$y)^2) 
        + popZ1$I + popZ2$I))
    if(!is.na(GIC))
        GIC <- 1-GIC
    else GIC <- 1
    GIC
}

"spatialpatches" <-
function (x, y, z, w, Lim.D = 100, A.li = 10, modproj = NA, mlong = NA, mlat = NA)
{
#===============================================================================
# IDENTIFICATION OF SPATIAL PATCHES 
# AND COUNT OF CENTERS OF GRAVITY OF SPATIAL PATCHES
#
# Routine from EU program Fisboat, DG-Fish, STREP n? 502572
# Authors : P.Petitgas (Ifremer), M.Woillez and J.Rivoirard (Mines-ParisTech)
# Last update : 01 march 2008  
#
# Arguments:
# x	      The x-coordinate of the first population
#	        Can be a vector, a matrix or an xyZ list (-,|,[]).
# y	      The y-coordinates of the first population
# z	      The regionalised variable of the first population.
#	        If missing, the results of 'cgi' will concern the samples only.
# w	      Optional. A weight or an area of influence of the first population.
#	        Set to 1 if missing
# Lim.D   Select minimum distance from sample to patch centre: to 
#         identify patches (units are those of coordinates)
# A.li    Visualisation of gravity centres for those patches with 
#         abundance > A.li (in %)
# modproj	Optional. Indicates the type of projection to perform.
# mlong   mean longitude in DEGREES of the data set to be transformed
# mlat    mean latitude in DEGREES of the data set to be transformed
#	        See 'dg2nm' for precisions
#
#===============================================================================
 
    patch.id<-function(xx,yy,zz,ww,dli,ali,modproj,mlong,mlat)
    {
    ############################################################################
    # inputs : xx , yy , ww & zz ranked by decreasing order of zz
    #          dli is a distance limit from the patch gravity centre to define
    #          the patch border
    #
    # the function starts from the richest zz value and considers each sample
    # in decreasing order of zz. It tests whether the current value is spatially
    # close enough to the gravity centre of previously formed patches. if not a
    # new patch is considered. and so on until the last value. Patches of nul
    # values are returned with centres as NA and code 0 and their areas are 
    # summed.
    #
    # outputs : the patch number for each sample,
    #           the gravity center for each patch,
    #           the percent abundance and area for each patch
    ############################################################################

        g<-c(1,rep(0,(length(xx)-1)))
        for (j in 1:(length(xx)-1)) {
            # if (j%%50==0) { cat(j,"\n") }
            xg<-tapply(ww[1:j]*zz[1:j]*xx[1:j],g[1:j],sum,na.rm=T)/
                tapply(ww[1:j]*zz[1:j],g[1:j],sum,na.rm=T)
            yg<-tapply(ww[1:j]*zz[1:j]*yy[1:j],g[1:j],sum,na.rm=T)/
                tapply(ww[1:j]*zz[1:j],g[1:j],sum,na.rm=T)
            d<-sqrt(((xg-x1[j+1]))^2+(yg-y1[j+1])^2)
            o<-order(d)
            if (d[o[1]]<dli)
                g[j+1]<-o[1]
            else
                g[j+1]<-max(g[1:j],na.rm=T)+1
        }
        xg<-tapply(ww*zz*xx,g,sum,na.rm=T)/tapply(ww*zz,g,sum,na.rm=T)
        yg<-tapply(ww*zz*yy,g,sum,na.rm=T)/tapply(ww*zz,g,sum,na.rm=T)
        pb<-tapply(ww*zz,g,sum,na.rm=T)/sum(ww*zz,na.rm=T)
        pa<-tapply(ww,g,sum,na.rm=T)/sum(ww,na.rm=T)
        sel<-(tapply(ww*zz,g,sum,na.rm=T)/sum(ww*zz,na.rm=T))*100>=ali
        cat('total nb of patches : ',max(g,na.rm=T),"\n")
        cat('nb of patches with abundance > ',round(ali,2),'% : ',
            max(unique(g)[sel],na.rm=T),"\n")
        res<-max(unique(g)[sel],na.rm=T)
        cat('percent abundance in these patches : ',round(pb[sel],4),"\n")
        cat('percent area in these patches : ',round(pa[sel],4),"\n")
        n<-sort(unique(g)); selna<-is.na(xg)
        if (sum(selna)>1) {
            n<-c(n[!selna],0)
            pb<-c(pb[!selna],0)
            pa<-c(pa[!selna],sum(pa[selna],na.rm=T))
            xg<-c(xg[!selna],NA)
            yg<-c(yg[!selna],NA)
            g[!is.element(g,n)]<-0
        }
        cg<-nm2dg(xg,yg,modproj=modproj,mlong=mlong,mlat=mlat)
        return(list(n=g,mat=cbind(n=n,xg=cg$x,yg=cg$y,pabun=round(pb*100,2),parea=pa),
            nsp=res))
    }

    # prepare data
    bid<-dg2nm(x,y,modproj=modproj,mlong=mlong,mlat=mlat)
    Xsta <- bid$x
    Ysta <- bid$y
    Zsta <- z
    Wsta <- w/sum(w,na.rm=T)

    # exclude NA values
    SEL <- is.na(Xsta)+is.na(Ysta)+is.na(Zsta)+is.na(Wsta)
    Xsta <- Xsta[!SEL]
    Ysta <- Ysta[!SEL]
    Zsta <- Zsta[!SEL]
    Wsta <- Wsta[!SEL] 
 
    # order data by decreasing order of value z
    zz1<-sort(Zsta,decreasing=T,index.return=T)
    w1<-Wsta[zz1$ix]
    x1<-Xsta[zz1$ix]
    y1<-Ysta[zz1$ix]
    z1<-zz1$x
 
    # identify patches around high values
    SP<-patch.id(x1,y1,z1,w1,Lim.D,A.li,modproj=modproj,mlong=mlong,mlat=mlat)
 
    # graphical visualisation : patch identification
    # shows data points, circles for data values, number of patch &
    # crosses for patch gravity centres with abundance > A.li of total
    o<-sort(zz1$ix,index.return=T)
    xy<-nm2dg(x1,y1,modproj=modproj,mlong=mlong,mlat=mlat)
    symbols(x,y,sqrt(z),xlab=" ",ylab=" ",fg=SP$n[o$ix]+1,inches=0.25,asp=1/cos((mlat*pi)/180))
    coast(add=T)
    text(xy$x, xy$y, paste(SP$n), cex=.7, col=SP$n+1)
    points(SP$mat[,2][SP$mat[,4]>A.li],SP$mat[,3][SP$mat[,4]>A.li],pch=3,cex=2,lwd=2)
    SP$n<-SP$n[o$ix]
    SP
}

"positivearea" <-
function (z,w)
{
#===============================================================================
# POSITIVE AREA
#
# Routine from EU program Fisboat, DG-Fish, STREP n? 502572
# Authors : M.Woillez and J.Rivoirard (Mines-ParisTech)
# Last update : 01 march 2008 
#
# Argument:
# z     variable of interest (i.e. fish density)
# w     appropriate areas of influence set as weighted factors
#
#===============================================================================

    PA <- sum(w[z>0],na.rm=T)
    PA
}

"spreadingarea" <-
function (z, w, plot = F)
{
#===============================================================================
# SPREADING AREA
#
# Routine from EU program Fisboat, DG-Fish, STREP n? 502572
# Authors : M.Woillez and J.Rivoirard (Mines-ParisTech)
# Last update : 01 march 2008 
#
# Arguments:
# z     variable of interest (i.e. fish density)
# w     appropriate areas of influence set as weighted factors
# plot  if TRUE, the curve expressing (Q?Q(T))/Q as a function of T is 
#       plotted with T the cumulated area occupied by the density values,
#       ranked in decreasing order, Q(T) the corresponding cumulated abundance,
#       and Q the overall abundance. The spreading area SA (expressed in square
#       nautical miles) is then simply defined as twice the area below this 
#       curve
#
#===============================================================================

    # extract data
    nb<-length(z)
    
    # sort data in increasing order
    zi <- sort(z,index.return=T)
    z<-zi$x
    w<-w[zi$ix]

    # computation of the spreading area 
    Q <- sum(z*w)
    QT <- c(0,cumsum(z*w))
    SA <- sum((QT[1:nb]+QT[2:(nb+1)])*w)/Q

    # computation of (Q?Q(T))/Q as a function of T
    T <- c(0,cumsum(w))
    T <- T[nb+1] - T
    T <- rev(T)
    QT <- QT[nb+1] - QT
    QT <- rev(QT)
        
    # display
    if(plot)
        plot(T, (Q-QT)/Q, main="Curve (Q?Q(T))/Q", type="o", pch="+")
    
    # outputs
    SA
}

"equivalentarea" <-
function(z, w)
{
#===============================================================================
# EQUIVALENT AREA
#
# Routine from EU program Fisboat, DG-Fish, STREP n? 502572
# Authors : M.Woillez and J.Rivoirard (Mines-ParisTech)
# Last update : 01 march 2008 
#
# Arguments:
# z     variable of interest (i.e. fish density)
# w     appropriate areas of influence set as weighted factors
#
#===============================================================================

    EA <- sum(z*w,na.rm=T)^2 / sum(z^2*w,na.rm=T)
    EA
}

"microstructure" <-
function(x, y, z, w, h0, pol, dlim, ndisc, modproj = NA, mlong = NA, mlat = NA)
{
#===============================================================================
# MICROSTRUCTURE INDEX
#
# Routine from EU program Fisboat, DG-Fish, STREP n? 502572
# Authors : M.Woillez and J.Rivoirard (Mines-ParisTech)
# Last update : 01 march 2008 
#
# Arguments:
# x	      The x-coordinate of the first population
#	        Can be a vector, a matrix or an xyZ list (-,|,[]).
# y	      The y-coordinates of the first population
# z	      The regionalised variable of the first population.
#	        If missing, the results of 'cgi' will concern the samples only.
# w	      Optional. A weight or an area of influence of the first population.
#	        Set to 1 if missing
# h0      mean lag between samples
# pol     polygon (in degree) used to delineate the maximal sampled area
# dlim    limit distance for valuating a node, when migrating points 
#         on a grid for computing the covariogram. 
#         GIVEN IN THE PROJECTED SPACE UNITS
# ndisc   Discretization of the data area when constructing the defaulted grid
# modproj	Optional. Indicates the type of projection to perform.
# mlong   mean longitude in DEGREES of the data set to be transformed
# mlat    mean latitude in DEGREES of the data set to be transformed
#	        See 'dg2nm' for precisions
# plot    If TRUE, plot is performed
#
#===============================================================================

    miss <- function(x){
	    length(x) == 1 && is.na(x)
    }

    sel <- !is.na(x * y * z * w)
	  x <- x[sel]
	  y <- y[sel]
	  z <- z[sel]
	  w <- w[sel]
	  nd <- length(x)
    if (length(dlim) == 1){
      opt <- 0 
      dlim <- rep(dlim, 2)
	  }
	  else opt <- 1

	  # Perform the projection if needed
	  
	  if(!miss(modproj)) {
		  bid <- dg2nm(x = x, y = y, modproj = modproj, mlong = mlong, mlat = mlat)
			x <- bid$x
			y <- bid$y
      if(!miss(pol))
        pol <- dg2nm(x = pol$x, y = pol$y, modproj = modproj, mlong = mlong, mlat = mlat)
		}

    # Build the grid 
    
    extend <- rev(sort(c((2*dlim)/diff(range(x)),(2*dlim)/diff(range(y)))))[1]
    xrange <- range(x)
    yrange <- range(y)
    dx <- xrange[2] - xrange[1]
    dy <- yrange[2] - yrange[1]
    xmin <- xrange[1] - extend * dx
    xmax <- xrange[2] + extend * dx
    ymin <- yrange[1] - extend * dy
    ymax <- yrange[2] + extend * dy
    nx <- ndisc
    ny <- ndisc
    dx <- (xmax - xmin) / nx
    dy <- (ymax - ymin) / ny
    ng <- nx * ny
    maille <- dx * dy 
  
    xg <- rep(xmin + dx * seq(0, nx - 1), ny)
    yg <- rep(ymin + dy * seq(0, ny - 1), rep(nx, ny))
  
    # Activate the grid nodes that are inside the polygon
  
    if(!miss(pol)){
      sel.pol <- inout(data.frame(xg,yg),data.frame(pol$x,pol$y))
      xg[!sel.pol] <- NA
      yg[!sel.pol] <- NA
    }

    # Migration of a variable defined on points to a grid
    # 1-Affectations

    argin <- function (x, test = 1.234e+30){
        if (length(x) > 0) 
            x[is.na(x)] <- test
        x
    }
    
    argout <- function (x, test = 1.234e+30){
        if (length(x) > 0) 
            x[x == test] <- NA
        x
    }
    
    x <- argin(x)
    y <- argin(y)
    z <- argin(z)
    xg <- argin(xg)
    yg <- argin(yg)
    zg <- rep(argin(NA),ng)
    dlim <- argin(dlim)

    # 2-Calculate the data extension #

    xxmin <- min(x[x!=argin(NA)])
    xxmax <- max(x[x!=argin(NA)])
    deltax <- xxmax - xxmin
    yymin <- min(y[y!=argin(NA)])
    yymax <- max(y[y!=argin(NA)])
    deltay <- yymax - yymin
    if(dlim[1]==argin(NA)) dlim[1] <- deltax * extend
    if(dlim[2]==argin(NA)) dlim[2] <- deltax * extend
    if(opt==0) dmax <- sqrt(dlim[1] * dlim[2])
    if(opt==1) dmax <- sqrt(dlim[1] * dlim[1] + dlim[2] * dlim[2])

    # 3-Build matrix of data points and of target points #
 
    xd <- matrix(x,ng,nd,byrow=T)
    yd <- matrix(y,ng,nd,byrow=T)
    rkd <- matrix(1:nd,ng,nd,byrow=T)
    xxg <- matrix(xg,ng,1)
    yyg <- matrix(yg,ng,1)
    rkg <- matrix(1:ng,ng,nd,byrow=F)
    
    # 4-Compute distance matrix #    
    
    dxx <- abs(sweep(xd,1,xxg,"-"))
    dyy <- abs(sweep(yd,1,yyg,"-"))   
    dist <- sqrt(dxx^2 + dyy^2)
    misval.dist <- sqrt((1.234e+30)^2+(1.234e+30)^2)
    
    # 5-Find the closest data point of an active target point according to options #
    
    if(opt==0)
        cond <- dist==apply(dist,1,min) & dist!=misval.dist
    if(opt==1)
        cond <- dxx==apply(dxx,1,min) & dyy==apply(dyy,1,min) & dist!=misval.dist
    proxd <- rkd[cond]
    proxg <- rkg[cond]   
    
    # 6-Migration step
    
    zg[proxg] <- z[proxd]
    zg <- argout(zg)
    zg[!sel.pol] <- NA

    # Compute the variogram map lags
       
    nlag <- ceiling((3*h0/2)/min(dy,dx))   
    if(length(nlag) == 1) 
        nlag <- rep(nlag, 2)
    nrow <- 2 * nlag[1] + 1
    ncol <- 2 * nlag[2] + 1
    nsize <- ncol * nrow
    g <- matrix(0,nrow,ncol)

    for (kx in -nlag[1]:nlag[1])
      for (ky in -nlag[2]:nlag[2]){
        xv1 <- (1:nx)-kx
        yv1 <- (1:ny)-ky
        xv2 <- (1:nx)+kx
        yv2 <- (1:ny)+ky
        xv1 <- xv1[xv1>0 & xv1<nx+1]
        yv1 <- yv1[yv1>0 & yv1<ny+1]
        xv2 <- xv2[xv2>0 & xv2<nx+1]
        yv2 <- yv2[yv2>0 & yv2<ny+1]                               
        g[(((ky)+nlag[2]) * (2*nlag[1]+1) + ((kx)+nlag[1]))+1] <- sum(matrix(zg,nx,ny)[yv1,xv1] * matrix(zg,nx,ny)[yv2,xv2],na.rm=T)      
      }
    g <- g*maille
    h <- sqrt(((rep((-nlag[1]:nlag[1]),nrow)*dx)^2)+((rep((-nlag[2]:nlag[2]),each=ncol)*dy)^2))
    gh0 <- mean(as.numeric(t(g))[h>=h0/2 & h<=3*h0/2])
    g0 <- sum(zg*zg*maille,na.rm=T)
    MI <- (g0-gh0)/g0
    MI
}

