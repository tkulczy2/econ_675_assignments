clear all
set more off

global basedir "/home/theodor/Desktop/Courses/ECON_675/Assignments/Assignment 1"

import delimited "$basedir/LaLonde_1986.csv", clear

gen educ2 = educ * educ
gen blackearn74 = black * earn74
gen const = 1

* 2.4.i - standard inverse
mata:
	y = st_data(., "earn78")
	X = st_data(., "const treat black age educ educ2 earn74 blackearn74 u74 u75")
	n = rows(X)
	d = cols(X)
	XpX = quadcross(X, X)
	XXinv = invsym(XpX)
	bhat = XXinv * quadcross(X, y)
	e = y - X * bhat
	vcvhat = XXinv * X' * diag(e:^2) * X * XXinv
	se = (diagonal(vcvhat)*n/(n-d-1)):^(1/2)
	tstat = bhat :/ se
	pval = t(n-d, tstat)
end


* 2.4.ii - cholesky inverse
mata:
	y = st_data(., "earn78")
	X = st_data(., "const treat black age educ educ2 earn74 blackearn74 u74 u75")
	n = rows(X)
	d = cols(X)
	XpX = quadcross(X, X)
	XXinv = cholinv(XpX)
	bhat = XXinv * quadcross(X, y)
	e = y - X * bhat
	vcvhat = XXinv * X' * diag(e:^2) * X * XXinv
	se = (diagonal(vcvhat)*n/(n-d-1)):^(1/2)
	tstat = bhat :/ se
	pval = t(n-d, tstat)
end


* 2.5.b
reg earn78 const treat black age educ educ2 earn74 blackearn74 u74 u75, vce(robust)


mata:	
	y = st_data(., "earn78")
	t = st_data(., "treat")
	N1 = sum(t:==1); N0 = rows(t)-N1
	mu1 = mean(select(y,t:==1)); mu1
	mu0 = mean(select(y,t:==0)); mu0
	v1 = variance(select(y,t:==1))/N1;
	v0 = variance(select(y,t:==0))/N0;
	
	T = (mu1-mu0)/sqrt(v1+v0)
	pval = 2*normal(-abs(T))
	
	st_numscalar("T", T); st_numscalar("pval", pval)
end

display "t-test = " T
display "pval   = " pval



