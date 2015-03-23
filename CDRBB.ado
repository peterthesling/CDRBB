/*###################################################################################*\
|#####                                                                           #####|
|#####     STATA code to compute block bootstrap panel unit root tests from      #####|
|#####                                                                           #####|
|##### "Cross-Sectional Dependence Robust Block Bootstrap Panel Unit Root Tests" #####|
|#####                                    by                                     #####|
|#####           Franz C. Palm, Stephan Smeekes and Jean-Pierre Urbain           #####|
|#####                                                                           #####|
|#####    Code by Peter Thesling (e-mail: peter.thesling@gmail.com)              #####|
|#####                                                                           #####|
\*###################################################################################*/

/*###################################################################################*\

To improve:
- make sure that matodd and moremata is installed?
\*###################################################################################*/


program CDRBB
	version 13.1
	syntax varlist [if] [in], [DC(integer 0) B(integer 9999) AUTOBL(integer 0) BK(integer 499) A(real 0.05)]
	
	/*
	This code could be implemented, if I am sure that the standard Stata
	does not know certain commands. It probably should be put right where the 
	command is called. So ask right before, if this works out and then prompt with
	installing certain packages. Matodd and moremata are to my knowledge not 
	installed yet and I actually need them. 
	#delimit ;
	cap which mm_median(X);
	if (_rc!=0)
	{
		di "You will need to install matodd and moremate before proceeding"
		di "type: findit moremata"
		di "type: findit matodd"
		exit
	}
	*/
	
	marksample touse // need this for the use of "if/in" commands
	preserve
	quietly keep if `touse'
	matrix test = (1 \ 1 \ 1)
	loc obs = `r(tmax)'-`r(tmin)'+1
	loc ind =  `r(imax)'-`r(imin)'+1
	loc imin = `r(imin)'
	loc imax = `r(imax)'
	loc balance = "`r(balanced)'"
	mata: mata clear
	putmata temp_y=`varlist' panelid=`r(panelvar)' timevar=`r(timevar)', replace omitmissing
	mata: initialize(temp_y, panelid, timevar,`b',`autobl',`dc',"`balance'", `bk',`a')
end

mata:
void initialize(real matrix temp_y, real matrix panelid, real matrix timevar, real scalar B, real scalar autobl, real scalar dc, string scalar balance, real scalar bk, real scalar a)
{
	tstat = J(3,1,0);
	cr = J(1,3,0); pval = J(3,1,0); rej = J(1,3,0)
	ind = max(panelid)-min(panelid)+1
	rows_temp_y = rows(temp_y);

	// Upcoming is the transformation of the data into a matrix.
	// I distinguish here between balanced and unbalanced panels.
	j=1; 
	if (balance == "strongly balanced")
	{
		obs = max(timevar)-min(timevar)+1
		y = J(obs, ind, 0)
		
		for (i=1;i<=rows_temp_y;i++)
		{
			y[i-(j-1)*obs,j] = temp_y[i]
			if (i-(j-1)*obs==obs) j=j+1
		}
	}
	
	else // unbalanced panel
	{
		imax = max(panelid)
		imin = min(panelid)
		current_time = J(2,imax-imin+1,-1)
		for (i=1;i<=rows_temp_y;i++)
		{
			for (j=imin; j<=imax;j++)
			{
				if (panelid[i]==j)
				{
					if(timevar[i]<current_time[1,j-imin+1] || current_time[1,j-imin+1]==-1)
					{
						current_time[1,j-imin+1] = timevar[i,1]
					}
					else if (timevar[i] >current_time[2,j-imin+1] || current_time[2,j-imin+1]==-1)
					{
						current_time[2,j-imin+1] = timevar[i,1]
					}
				}
			}
		}
		current_time
		tmax = min(current_time[2,.]) // min of all end times
		tmin = max(current_time[1,.]) // max of all beginning times
		y = J(tmax-tmin+1,ind,0)
		
		for (i=1;i<=rows_temp_y;i++)
		{
			if (timevar[i,1]>=tmin && timevar[i,1] <=tmax)
			{
				y[timevar[i]-tmin+1,panelid[i]]=temp_y[i]
			}
		}
	}
	
	T = rows(y)
	N = cols(y)
	bmax = ceil(0.75*T); 
	b0 = ceil(1.75*T^(1/3))
	nb = ceil((T-1)/b0)
	ybv = J(N*T,bk,0)
	indm = J(nb,bk,0)
	
	if (autobl==0) bl = b0*J(3,1,1)
	
	// automatic bootstrap selection
	else if (autobl==-1)
	{
		printf("Since you selected an automatic determination of the block length, this might take several minutes.\n")
		bl = J(3,1,0)
		timer_clear()
		for (i=1;i<=3;i++)
		{
			timer_on(i)
			bl[i]=WS_calibration(i,y,a,b0,bmax,bk,dc,ybv,indm,3)
			printf("Finished calibration for %f/3 tests. \n",i)
			timer_off(i)
			timer(i)
		}
	}
	
	else if (autobl>=1) bl = autobl*J(3,1,1)
	// The upcoming programs are very similar to yours and I only had to slightly adjust them.
	pur_b(y,B,bl,tstat,cr,pval,rej,dc,a)
	output_on_screen(tstat,cr,pval,bl)
}


void pur_b(real matrix y, real scalar B, real matrix bl, real matrix tstat, real matrix cr, real matrix pval, real matrix rej, real scalar dc, real scalar a)
{
	for (i=1;i<=3;i++)
	{
		if (i == 1) tstat[i,1] = pooled(y,dc)
		
		if (i == 2) tstat[i,1] = group_mean(y,dc)
		
		if (i == 3) tstat[i,1] = med(y,dc)

		// Have to ask Stephan about this last number
		tb = pur_bl(i,y,B,bl[i],dc)
		cr[.,i] = tb[ceil(a*B),.]
		pval[i,.] = counts(tb,tstat[i,1])/B
		rej[.,i] = tstat[i,1]:<cr[1,i]
	}
}


real scalar pooled(real matrix y,real scalar dc)
{
	real matrix x, xls, dxs
	T = rows(y)
	x = detr(y, dc)
	xls = vec(x[(1::T-1),.])
	dxs = vec(dif(x))
	t_p = T*luinv(xls'*xls)*xls'*dxs
	return(t_p)
}

real matrix group_mean(real matrix y, real scalar dc)
{
	real matrix x, dx, tvec
	N = cols(y)
	T = rows(y)
	x = J(T-1,N,0)
	tvec = J(N,1,0)
	for (i=1;i<=N;i++)
	{
		x = detr(y[.,i],dc)
		dx = dif(x)
		tvec[i,1] = T*luinv(x[(1::T-1),.]'*x[(1::T-1),.])*x[(1::T-1),.]'*dx
	}
	t_gm = mean(tvec)
	return(t_gm)
}


real matrix med(real matrix y, real scalar dc)
{
	real matrix x, dx, tvec
	N = cols(y)
	T = rows(y)
	tvec = J(N,1,0)
	for (i=1;i<=N;i++)
	{
		x = detr(y[.,i],dc)
		dx = dif(x)
		tvec[i,1] = T*luinv(x[(1::T-1),.]'*x[(1::T-1),.])*x[(1::T-1),.]'*dx
	}
	t_med = mm_median(tvec)
	return(t_med)
}


real matrix detr(real matrix y,real scalar dc)
{
	T = rows(y)
	real matrix m, b, yd
	m = J(T,1,1),(1::T)
	m = m[.,(1::ceil((dc+1)/2))]
	b = luinv(m'm)*(m'*y)
	yd = y - m*b*ceil(dc/2)
	return(yd)
}

real matrix dif(real matrix k)
{
	real matrix dk
	dk = k[(2::rows(k)),.]-k[(1::rows(k)-1),.]
	return(dk)
}

real matrix pur_bl(real scalar test, real matrix y, real scalar B, real scalar bl, real scalar dc)
{
	N = cols(y)
	T = rows(y)
	nb = ceil((T-1)/bl)
	u = u_hat(y,dc)
	tb = J(B,1,0)
	for (i=1;i<=B;i++)
	{
		urn = uniform(nb,1)
		Ii = ceil((J(nb,1,1)-urn)*(T-bl))
		ind = seqa(Ii[1],1,bl)
		for (j=2; j<=nb;j++) ind = ind \ seqa(Ii[j,1],1,bl)
		ind = ind[(1::T-1),.]
		ub = u[ind,.]
		yb = recserar(J(1,N,0) \ ub, y[1,.], J(1,N,1))

		if (test == 1) tb[i,.] = pooled(yb,dc)
		
		else if (test == 2) tb[i,.] = group_mean(yb,dc)
		
		else if (test == 3) tb[i,.] = med(yb,dc)
	}
	tb = sort(tb,1)
	return(tb)
}


real matrix u_hat(real matrix y, real scalar dc)
{
	T = rows(y)
	N = cols(y)
	u = J(T-1,N,0)
	
	for (i=1;i<=N;i++) u[.,i] = u_hat_i(y[.,i],dc)
	
	return(u)
}

real matrix u_hat_i(real matrix y, real scalar dc)
{
	real matrix x, dx, lx
	T = rows(y)
	x = detr(y,dc)
	dx = dif(x)
	lx = x[(1::T-1),.]
	r = luinv(lx'*lx)*lx'*x[(2::T),.] // changed this slightly from Stephans code
	u = demean(x[(2::T),.]-lx*r)
	return(u)
}	

real matrix demean(real matrix u)
{
	real matrix udm
	udm = u - mean(u)*J(rows(u),cols(u),1)
	return(udm)
}


real scalar WS_calibration(real scalar test, real matrix y, real scalar a, real scalar b0, real scalar bmax, real scalar bk, real scalar dc, real matrix ybv, real matrix indm, real scalar nt)
{
	T = rows(y)
	N = cols(y)
	tstat = pur_bl_sel(test,y,bk,b0,dc,ybv,indm)
	tsts = sort(tstat,1)
	crb = tsts[ceil(a*bk),.]
	adrf = J(bmax,1,0)
	rejf = J(bmax,1,0)
	
	for (j=1;j<=bmax;j++)
	{
		cr = J(bk,1,0)
		for (i=1;i<=bk;i++)
		{
			yb = reshape(ybv[.,i],N,T)'
			cr[i,.] = pur_bl(test,yb,1,j,dc)
		}
		crs = sort(cr,1)
		adrf[j,.] = abs(crb-crs[ceil(a*bk)])
		rejf[j,.] = mean(tstat:<crs[ceil(a*bk)])
	}
	// adrf
	optb = minindc(adrf)
	return(optb)
}

real matrix pur_bl_sel(real scalar test, real matrix y, real scalar B, real scalar bl, real scalar dc, real matrix ybv, real matrix indm)
{
	N = cols(y)
	T = rows(y)
	nb = ceil((T-1)/bl)
	u = u_hat(y,dc)
	tb = J(B,1,0)
	for(i=1;i<=B;i++)
	{
		urn = uniform(nb,1)
		Ii = ceil((J(nb,1,1)-urn)*(T-bl))
		indm[.,i] = Ii
		ind = seqa(Ii[1,.],1,bl)
		
		for (j=2; j<=nb; j++) ind = ind \ seqa(Ii[j,.],1,bl)
		
		ind = ind[(1::T-1),.]
		ub = u[ind,.]
		yb = recserar(J(1,N,0)\ub,y[1,.],J(1,N,1))
		ybv[.,i] = vec(yb)
		
		if (test == 1) tb[i,.] = pooled(yb,dc)
		
		else if (test == 2) tb[i,.] = group_mean(yb,dc)
		
		else if (test == 3) tb[i,.] = med(yb,dc)		
	}
	
	tb = sort(tb,1)
	return(tb)
}

void output_on_screen(real matrix tstat, real matrix cr, real matrix pval, real matrix bl)
{	
	printf("\n")
	
	printf("{txt} Test {space 2} {c |} {space 1} Statistic \t {space 2} 5%% c.v. \t p-value {space 1} Block length \n")
	printf("{hline 9}{c +}{hline 60}\n")
	
	for (i=1;i<=3;i++)
	{
		if (i==1) printf("{txt}Pooled \t {c |} {res}")
		if (i==2) printf("{txt}GM \t {c |} {res}")
		if (i==3) printf("{txt}Median \t {c |} {res}")
		printf(" %9.3fc \t ", tstat[i])
		printf(" %9.3fc ", cr[i])
		printf(" %9.3fc ", pval[i])
		printf(" %10.0g ", bl[i])
		printf("\n")
	}
}

// Functions that I had to program myself, but are given in GAUSS already

real matrix seqa(real scalar start, real scalar inc, real scalar n)
{
	output = J(n,1,0)
	for (k=1; k<=n; k=k++) output[k,1] = (k-1)*inc+start
	return(output)
}

// Found code for this online: http://www.spatial-econometrics.com/util/recserar.m
real matrix recserar(real matrix x, real matrix y0, real matrix a)
{	
	n1 = rows(x)
	k1 = cols(x)
	p1 = rows(y0)
	k2 = cols(y0)
	p2 = rows(a)
	k3 = cols(a)
	result=J(n1,k1,0)
	for (j=1;j<=p1;j++)
	{
		result[j,.] = y0[j,.]
	}
	for (j=p1+1;j<=n1;j++){
		result[j,.]=x[j,.]
		for (k=1;k<=p1;k++){
			result[j,.] = result[j,.]+a[k,.]:*result[j-k,.]
		}
	}
	return(result)
}

real matrix minindc(real matrix x)
{
	xr = rows(x)
	xc = cols(x)
	current_index = J(xc,1,0)
	for (j=1;j<=xc;j++)
	{
		current_min = -9999
		for (i=1;i<=xr;i++)
		{
			if (x[i,j]< current_min || current_min == -9999)
			{
				current_min = x[i,j]
				current_index[j,1] = i
			}
		}
	}
	return(current_index)
}

real matrix reshape(real matrix x, real scalar r, real scalar c)
{
	output = J(r,c,0)
	xr = rows(x)
	xc = cols(x)
	k=1
	l=1
	for (i=1;i<=r;i++)
	{
		for (j=1;j<=c;j++)
		{
			output[i,j] = x[k,l]
			k=k+1
			if (k==xr && l<xc)
			{
				l=l+1
			}
			else if(k==xr && l==xc)
			{
				j=c+1 
				i=r+1 // meaning stop (maybe warning message, that some stuff deleted)
			}
		}
	}
	return(output)
}

real matrix counts(real matrix x, real matrix v)
{
	if (rows(v)==1)
	{
		counter = (0);
		for (i=1;i<=rows(x);i++)
		{
			if (x[i,1]<=v[1,1]) counter[1,1]=counter[1,1]+1
		}
	}
	
	else
	{
		counter = J(rows(v),1,0)
		for (i=1;i<=rows(v);i++)
		{
			for (j=1;j<=rows(x);j++)
			{
				if (x[j,1]<=v[i,1]) counter[i,1] = counter[i,1]+1

				else if (x[j,1] <= v[i-1,1] && x[j,1] >= v[i,1]) counter[i,1] = counter[i,1]+1
			}
		}
	}
	return(counter)

}

end
exit
