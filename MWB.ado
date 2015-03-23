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


program CDRBB
	version 13.1
	syntax varlist [if] [in], [BS(string) TEST(string) B(integer 999) LVL(real 0.05) KAPPA(integer 0) GLSDT(integer 0) ADFT(integer 1) PMIN(integer 0) PMAX(integer -1) IC(string) LLB(integer 1) LBWB(real -1) LDWB(real -1) KDWB(string)  GAMMAAWB(real -1) KRS(string) HRS(real 0.1)]
	if ("`bs'" == "") loc bs = "bwb,dwb,awb"
	else local bs "`bs'"
	
	if ("`test'" == "") loc test "av,med"
	else local test "`test'"
	
	if ("`ic'" == "") local ic "aic"
	else local ic "`ic'"
	
	if ("`kdwb'" == "") local k_dwb "trapk"
	else local k_dwb "`k_dwb'"
	
	if ("`krs'" == "") local krs "gaussian"
	else local krs "`krs'"

	/*
	#delimit ;
	cap which integrate(X);
	if (_rc!=0)
	{
		di "You will need to install the integrate package before proceeding"
		di "type: findit integrate"
		exit
	}
	*/
	
	marksample touse // need this for the use of "if/in" commands
	preserve
	quietly keep if `touse'
	loc balance = "`r(balanced)'"
	
	mata: mata clear
	mata: _c = 0.43
	mata: _t = 0
	
	quietly putmata temp_y=`varlist' panelid=`r(panelvar)' timevar=`r(timevar)', replace omitmissing
	mata: initialize(temp_y,panelid,timevar,"`balance'","`test'","`bs'",`b',`lvl',`kappa',`glsdt',`adft',`pmin',`pmax',"`ic'",`llb',`lbwb',`ldwb',"`k_dwb'",`gammaawb',"`krs'",`hrs')
end


mata:

void initialize(real matrix temp_y, real matrix panelid, real matrix timevar, string scalar balance, string scalar test, string scalar bs, real scalar B, real scalar lvl, real scalar kappa, real scalar gls_dt, real scalar adf_t, real scalar pmin, real scalar pmax, string scalar ic, real scalar llb, real scalar l_bwb, real scalar l_dwb, string scalar k_dwb, real scalar gamma_awb, string scalar k_rs, real scalar h_rs)
{	
	// Upcoming is the transformation of the data into a matrix.
	// I distinguish here between balanced and unbalanced panels.
	j=1
	ind = max(panelid)-min(panelid)+1
	rows_temp_y = rows(temp_y)
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
	
	if (pmax  == -1)     pmax  = floor(12*(T/100)^(1/3))
	if (l_bwb == -1)     l_bwb = round(1.75*T^(1/3))
	if (l_dwb == -1)     l_dwb = 1.75*T^(1/3)
	if (gamma_awb == -1) gamma_awb = 0.01^l_dwb
	
	bs_m = J(3,1,"")
	if (strpos(bs,"bwb")!= 0) bs_m[1] = "bwb"
	if (strpos(bs,"dwb")!= 0) bs_m[2] = "dwb"
	if (strpos(bs,"awb")!= 0) bs_m[3] = "awb"
	
	test_m = J(2,1,"")
	if (strpos(test,"av")!= 0 || strpos(test,"mean")) test_m[1]="av"
	if (strpos(test,"med")!= 0) test_m[2]="med"
	
	nt = rows(test_m)
    nb = rows(bs_m)
	ts = J(nt,1,0)
	cv = J(nt,nb,0)
	pv = J(nt,nb,0)
	
	pur(ts,cv,pv,test_m,bs_m,y,B,l_bwb,l_dwb,k_dwb,gamma_awb,(pmin\pmax),llb,ic,kappa,gls_dt,adf_t,k_rs,h_rs,lvl)
	output_on_screen(bs_m,test_m,lvl,ts,cv,pv)
}

void pur(real matrix ts, real matrix cv, real matrix pv, string matrix test, string matrix bs, real matrix y, real scalar B, real scalar bl, real scalar h, string scalar k, real scalar ar, real matrix pm, real scalar llb, string scalar ic, real scalar kap, real scalar g, real scalar cft, string scalar ks, real scalar hs, real scalar a)
{
	T = rows(y)
	N = cols(y)
	dc = kap+1
	c_gls = 7 \ 13.5
	cT = (c_gls[ceil((dc+1)/2)]/T)^g
	p = panellags(ic,y,pm[1],pm[2],1,dc,ks,hs)
	s = sqrtsigma(k,h,T)
	ti = istat(y,p,dc,cT,cft)
	nt = rows(test)
	nb = rows(bs)

	for (i=1;i<=nt;i++)
	{
		if (i==1 && test[1] == "av")      ts[i] = av(ti)
		else if (i==2 && test[2] == "med") ts[i] = med(ti)
		
		for (j=1;j<=nb;j++)
		{
			if (bs[j]!="" && test[i]!="")
			{
				temp1 = cv[i,j]
				temp2 = pv[i,j]
				purt(temp1,temp2,test[i],bs[j],y,B,bl,s,ar,p,pm,llb,ic,dc,cT,cft,ks,hs,a,ts[i])
				cv[i,j] = temp1
				pv[i,j] = temp2 
			}
		}
	}
}


// "Panel" test statistic
real scalar av(real matrix x)
{
	return(mean(x))
}

real scalar med(real matrix x)
{
	return(mm_median(x))
}

real matrix istat(real matrix y, real matrix p, real scalar dc, real scalar c, real scalar cft)
{
	T = rows(y)
	N = cols(y)
	tti = J(N,1,0)
	cfi = J(N,1,0)
	e  = J(T,N,0)
	phi = J(max(p)+2,N,0)
	adf_N(tti,cfi,e,phi,y,p,dc,c)
	ti = cft*tti+(1-cft)*cfi
	return(ti)
}

// N-dimensional vector of ADF statistics
void adf_N(real matrix tt, real matrix tc, real matrix e, real matrix phi, real matrix y, real matrix p, real scalar dc, real scalar c)
{
	T = rows(y)
	N = cols(y)
	yd = detr(y,dc,c)
	for (i=1;i<=N;i++)
	{
		temp1 = tt[i]
		temp2 = tc[i]
		temp3 = e[.,i]
		
		adf_i(temp1,temp2,temp3,phii,yd[.,i],p[i],T)
		phi[(1::rows(phii)),i] = phii
		
		tt[i]   = temp1
		tc[i]   = temp2
		e[.,i]  = temp3
	}
}

void adf_i(real scalar tt, real scalar tc, real matrix e, real matrix b, real matrix y, real scalar p, real scalar T)
{
	dy0 = dif0(y)
	dy = dy0[(2+p::T)]
	mm = lgm0(y,1),lgm0(dy0,p+1)
	m = mm[(2+p::T),(1::1+p)]
	b = luinv(m'*m)*m'*dy
	e = dy0 - mm[.,(1::1+p)]*b
	s2 = e[(2+p::T)]'*e[(2+p::T)]/(T-p-1)
	Vb = s2*luinv(m'*m)
	tt = b[1]/sqrt(Vb[1,1])
	tc = T*b[1]/(1-sum(b)+b[1])
}

real matrix detr(real matrix y, real scalar dc, real scalar cn)
{
	T = rows(y)
	N = cols(y)
	dcc = ceil((dc+1)/2)
	z = J(T,1,1), seqa(1,1,T)
	zc = dif0(z) + cn*(J(1,2,0) \ z[(1::T-1),.])
	yc = dif0(y) + cn*(J(1,N,0) \ y[(1::T-1),.])
	z = z[.,(1::dcc)]
	zc = zc[.,(1::dcc)]
	b = luinv(zc'*zc)*zc'*yc
	x = y-ceil(dc/2)*z*b
	return(x)
}

void purt(real scalar cr, real scalar pv, string scalar test, string scalar bs, real matrix y, real scalar B, real scalar bl, real matrix s, real scalar ar, real matrix p, real matrix pm, real scalar llb, string scalar ic, real scalar dc, real scalar c, real scalar cft, string scalar k, real scalar h, real scalar a, real scalar tst)
{
	N = cols(y)
	tb = J(B,1,0)
	pb = J(B,N,0)
	pur_bl(tb,pb,test,bs,y,B,bl,s,ar,p,pm,llb,ic,dc,c,cft,k,h)
	ts = sort(tb,1)
	cr = ts[ceil(a*B)]
	pv = mean(tb:<tst)
}

void pur_bl(real matrix tb, real matrix pb, string scalar test, string scalar bs, real matrix y, real scalar B, real scalar bl, real matrix s, real scalar ar, real matrix p, real matrix pm, real scalar llb, string scalar ic, real scalar dc, real scalar c, real scalar cft, string scalar k, real scalar h)
{
	N = cols(y)
	T = rows(y)
	nb = ceil(T/round(bl))
	tt = J(N,1,0)
	tc = J(N,1,0)
	e  = J(T,N,0)
	phi = J(max(p)+2,N,0)
	adf_N(tt,tc,e,phi,y,p,dc,1)
	r = 1:+phi[1,.]'
	yd = detr(y,dc,1)
	u = yd - lgm0(yd,1):*r'

	for (i=1;i<=B;i++)
	{
		if (bs == "bwb") yb = bwb(u,e,J(1,N,0),nb,bl,s,ar,phi)
		else if (bs == "dwb") yb = dwb(u,e,J(1,N,0),nb,bl,s,ar,phi)
		else if (bs == "awb") yb = awb(u,e,J(1,N,0),nb,bl,s,ar,phi)
		
		if (llb == 1) pb[i,.] = panellags(ic,yb,pm[1],pm[2],1,dc,k,h)'
		else pb[i,.] = p'
		tib = istat(yb,pb[i,.]',dc,c,cft)
		
		if (test == "av")       tb[i] = av(tib)
		else if (test == "med") tb[i] = med(tib)
	}
}

real matrix bwb(real matrix u, real matrix e, real matrix y1, real matrix nb, real scalar bl, real matrix s, real scalar ar, real matrix ph)
{
	x = uniform(round(nb),1)

	xi = x#J(bl,1,1)
	ub = y1 \ (u:*xi[(1::rows(u))])
	yb = recserar(ub,ub[1,.],J(1,cols(u),1))
	yb = yb[(2::rows(u)+1),.]
	return(yb)
}

real matrix dwb(real matrix u, real matrix e, real matrix y1, real matrix nb, real scalar bl, real matrix s, real scalar ar, real matrix ph)
{
	x = uniform(rows(u),1)
	xi = s'x
	ub = y1 \ (u:*xi)
	yb = recserar(ub,ub[1,.],J(1,cols(u),1))
	yb = yb[(2::rows(u)+1),.]
	return(yb)
}

real matrix awb(real matrix u, real matrix e, real matrix y1, real matrix nb, real scalar bl, real matrix s, real scalar ar, real matrix ph)
{
	x = uniform(rows(u),1)
	xi = recserar(sqrt(1-ar^2)*x,x[1],ar)
	ub = y1 \ (u:*xi)
	yb = recserar(ub,ub[1,.],J(1,cols(u),1))
	yb = yb[(2::rows(u)+1),.]
	return(yb)
}

real matrix panellags(string scalar ic, real matrix y, real scalar pmin, real scalar pmax, real scalar c, real scalar dc, string scalar k, real scalar h)
{
	
	N = cols(y)
	p = J(N,1,0)
	for (i=1;i<=N;i++) p[i] = ic_adf(ic,y[.,i],pmin,pmax,c,dc,k,h)
	return(p)
}

real scalar ic_adf(string scalar cri, real matrix y, real scalar pmin, real scalar pmax, real scalar dc, real scalar c, string scalar k, real scalar h)
{
	r = 0
	n = rows(y)
	d = n-pmax
	ys = J(rows(y),cols(y),0)
	rescale(ys,cri,y,pmax,dc,h,k)
	
	yd = detr(ys,dc,c)

	yl = yd[(1+pmax::n-1)]
	
	real matrix e
	e = adf_t_ic(r,yd,pmin,pmax,n)
	if      (cri == "aic")    ic =    aic(e,pmin,d,yl,r)
	else if (cri == "bic")    ic =    bic(e,pmin,d,yl,r)
	else if (cri == "maic")   ic =   maic(e,pmin,d,yl,r)
	else if (cri == "mbic")   ic =   mbic(e,pmin,d,yl,r)
	
	p_ic = seqa(pmin,1,pmax-pmin+1),(ic*J(pmax-pmin+1,1,1))
	for (p=pmin+1;p<=pmax;p++)
	{
		e = adf_t_ic(r,yd,p,pmax,n)
		if      (cri == "aic")    p_ic[p-pmin+1,2] =    aic(e,p,d,yl,r)
		else if (cri == "bic")    p_ic[p-pmin+1,2] =    bic(e,p,d,yl,r)
		else if (cri == "maic")   p_ic[p-pmin+1,2] =   maic(e,p,d,yl,r)
		else if (cri == "mbic")   p_ic[p-pmin+1,2] =   mbic(e,p,d,yl,r)
	}
	mat = sort(p_ic,2)
	return(mat[1,1])
}

void rescale(real matrix ys, string scalar cri,real matrix y, real scalar pmax, real scalar dc, real matrix h, string scalar k)
{
	if (substr(cri,1,2)=="rs")
	{
		cri = substr(cri,3,strlen(cri)-2)
		detr(yd,y,dc,1)
		adfsc(u,phi,yd,0)
		dif0(yb,yd)
		sh = sqrt(npve(u,h,k))
		ys = recserar(yb:/sh,yb[1]/sh[1],1)
	}
	else ys = y
}

void npve(real matrix vm, real matrix u, real matrix h, string scalar k)
{
	n = rows(u)
	r = seqa(1,1,n)/n

	if (k=="gaussian") vm = (gaussian(spec_minus(r)'/h)*u^2):/sum(gaussian((r-r')/h))
	else if (k=="truncated") vm = (truncated(spec_minus(r)'/h)*u^2):/sum(truncated((r-r')/h))
	else if (k=="bartlett") vm = (bartlett(spec_minus(r)'/h)*u^2):/sum(bartlett((r-r')/h))
	else if (k=="parzen") vm = (parzen(spec_minus(r)'/h)*u^2):/sum(parzen((r-r')/h))
	else if (k=="tukey_hanning") vm = (tukey_hanning(spec_minus(r)'/h)*u^2):/sum(tukey_hanning((r-r')/h))
	else if (k=="quadratic_spectral") vm = (quadratic_spectral(spec_minus(r)'/h)*u^2):/sum(quadratic_spectral((r-r')/h))
	
}

real matrix spec_minus(real matrix x)
{
	// calculates t-t', where t is a column vector
	y = J(rows(x),rows(x),0)
	for (i=1;i<=rows(x);i++)
	{
		for (j=1;j<=rows(x);j++)
		{
			y[i,j] = x[j]-x'[i]
		}
	}
	return(y')
}

// ADF estimation for rescaling; can be used to modify RSIC - for details see Cavaliere et al. (2014)
void adfsc(real matrix e, real matrix b, real matrix y, real scalar p)
{
	n = rows(y)
	dy = dif0(y)
	m = lgm0(y,1),lgm0(dy,p+1)
	m = m[.,(1::1+p)]
	m2 = m[(2+p::n),.]
	b = luinv(m2'*m2)*m2'*dy[(2+p::n)]
	e = dy-m*b
}


real matrix adf_t_ic(real scalar r, real matrix y, real scalar p, real scalar q, real scalar n)
{	
	dy = dif(y)
	m = y[(1::n-1)],lgm0(dy,p+1)
	m = m[(1+q::n-1),(1::1+p)]
	b = luinv(m'*m)*m'*dy[(1+q::n-1)]
	e = dy[(1+q::n-1)]-m*b
	r = b[1]
	return(e)
}


real scalar aic(real matrix e, real matrix p, real matrix d, real matrix yl, real matrix r)
{
	return(ln(det(e'*e/d))+2*p/d)
}	

real scalar bic(real matrix e, real matrix p, real matrix d, real matrix yl, real matrix r)
{
	return(ln(det(e'*e/d))+ln(d)*p/d)
}	

real scalar maic(real matrix e, real matrix p, real matrix d, real matrix yl, real matrix r)
{
	t_p = (yl'*yl)*r^2/(e'*e/d)
	return(ln(det(e'*e/d))+2*(p+t_p)/d)
}

real scalar mbic(real matrix e, real matrix p, real matrix d, real matrix yl, real matrix r)
{
	t_p = (yl'*yl)*r^2/(e'*e/d)
	return(ln(det(e'*e/d))+ln(d)*(p+t_p)/d)
}

// Procedures for calculating the DWB covariance matrix
real matrix sqrtsigma(string scalar k, real matrix h, real scalar n)
{
	t = seqa(1,1,n)
	
	if (k=="gaussian") s2 = gaussian(spec_minus(t)/h)
	else if (k=="trapk") s2 = trapk(spec_minus(t)/h)
	else if (k=="truncated") s2 = truncated(spec_minus(t)/h)
	else if (k=="bartlett") s2 = bartlett(spec_minus(t)/h)
	else if (k=="parzen") s2 = parzen(spec_minus(t)/h)
	else if (k=="tukey_hanning") s2 = tukey_hanning(spec_minus(t)/h)
	else if (k=="quadratic_spectral") s2 = quadratic_spectral(spec_minus(t)/h)
	s = pddc(s2,1e-10)
	return(s)
}

real matrix pddc(real matrix m, real matrix e)
{
	y = pdm(m,e)
	s = cholesky(y)
	return(s)
}

real matrix pdm(real matrix m, real scalar e) // G. Rapuch & T. Roncalli (2001) - DefiniteCorrelation
{
	n = rows(m)
	eigensystem(m,ve,va)
	va = Re(va)
	ve = Re(ve)
	b = va:>e
	f = sum(va-(J(rows(b),cols(b),1)-b)*e)/sum(b:*va)
	
	temp = va
	va = J(n,n,0)
	
	_diag(va,f*(b:*temp)+(J(rows(b),cols(b),1)-b)*e)
	k = ve*va*qrinv(ve)
	s = sqrt(diagonal(k))
	k = k:/s:/s'
	return(Re(k))
}

// "trapk" kernel
real matrix trapk(real matrix t)
{
	external _c
	y = self_conv(t,_c)/self_conv(0,_c)
	return(y)
}

// Procedures needed for "trapk"
real matrix self_conv(real matrix t, real scalar c)
{
	y = J(rows(t), cols(t), 0)
	external _c
	external _t
	_c = c
	for (j=1;j<=cols(t);j++)
	{
		for (i=1;i<=rows(t);i++)
		{
			_t = t[i,j]
			y[i,j] = integrate(&conv_k(),-1,1)
		}
	}
	return(y)
}

real matrix conv_k(real matrix x)
{
	external _c
	external _t
	y = w_c(x,_c):*w_c(x+J(rows(x),cols(x),abs(_t)),_c)
	return(y)
}

real matrix w_c(real matrix x, real scalar c)
{
	y = (x/c):*(x:>=0):*(x:<c) + (x:>=c):*(x:<=1-c) + ((J(rows(x),cols(x),1)-x)/c):*(x:>1-c):*(x:<=1)
	return(y)
}

real matrix truncated(real matrix x)
{
	y = abs(x):<=1
	return(y)
}

real matrix bartlett(real matrix x)
{
	y = (1-abs(x)):*(abs(x):<=1)
	return(y)
}

real matrix gaussian(real matrix x)
{
	y = normalden(x)
	return(y)
}


real matrix parzen(real matrix x)
{
	y = (1-6*x^2+6*abs(x)^3):*(abs(x):<=1/2)+2*((1-abs(x))^3):*(abs(x):>1/2):*(abs(x):<=1)
	return(y)
}


real matrix tukey_hanning(real matrix x)
{
	y = (abs(x):<=1):*(1+cos(pi*x))/2
	return(y)
}

real matrix quadratic_spectral(real matrix x)
{
	y = 25*(sin(6*pi*x/5):/(6*pi*x/5) - cos(6*pi*x/5)):/(12*x^2*pi^2) + (abs(x):<=1e-10)
	return(y)
}

real matrix dif(real matrix y)
{
	dy = y[(2::rows(y)),.] - y[(1::rows(y)-1),.]
	return(dy)
}

real matrix dif0(real matrix y)
{
	dy = y-(J(1,cols(y),0)\y[(1::rows(y)-1),.])
	return(dy)
}

real matrix lgm0(real matrix x, real scalar p)
{
	k = cols(x)
	sq = seqa(1,1,p)#J(k,1,1)
	xx = J(p,1,1)#x'
	lx1p = shiftr(xx,sq,0)'
	return(lx1p)
}

void output_on_screen(string matrix bs, string matrix test, real scalar lvl, real matrix ts, real matrix cv, real matrix pv)
{
	for (k=1;k<=rows(bs);k++)
	{
		if ((k == 1 && bs[k] == "bwb") || (k == 2 && bs[k] == "dwb") || (k == 3 && bs[k] == "awb"))
		{
			printf("\n")
			printf("%s",strupper(bs[k]))
			printf("\n")
			printf("{txt} Test {space 2} {c |} {space 1} Statistic \t {space 2} 5%% c.v. \t p-value \n")
			printf("{hline 9}{c +}{hline 60}\n")

			for (i=1;i<=rows(test);i++)
			{
				if ((i==1 && test[i] == "av") || (i==2 && test[i] == "med"))
				{
					if (i==1 && test[i] == "av") printf("{txt}Mean \t {c |} {res}")
					if (i==2 && test[i] == "med") printf("{txt}Median \t {c |} {res}")
					
					printf(" %9.3fc \t ", ts[i])
					printf(" %9.3fc ", cv[i,k])
					printf(" %9.3fc ", pv[i,k])
					printf("\n")
				}
			}
			printf("\n")
		}
	}
}


real matrix shiftr(real matrix x, real matrix s, real scalar f)
{
	y = x
	for (i=1;i<=rows(x);i++)
	{
		if (s[i]>=1)
		{
			for (j=1;j<=cols(x);j++)
			{
				if (j<=s[i]) y[i,j] = f

				else y[i,j] = x[i,j-s[i]]	
			}
		}
		
		else if (s[i] <=-1)
		{
			for (j=1;j<=cols(x);j++)
			{
				if (j<=-s[i]) y[i,j] = x[i,j-s[i]]
				else y[i,j] = f
			}
		}
	}
	return(y)
}

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
		for (k=1;k<=p1;k++) result[j,.] = result[j,.]+a[k,.]:*result[j-k,.]
	}
	return(result)
}
end

exit
