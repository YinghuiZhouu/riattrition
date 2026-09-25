*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Replication for: Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	-------
*	Recruitment Experiment Figures 
* 	Figure 9: Share of Women Choosing and Showing Up for Job Outside the Home 
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	use "${datadir}/clean_data/05_recruitment_exp_clean.dta", clear 
	
* ------------------------------------------------------------------------------
* 	Recruitment Experiment - Figure 9: Share of Women Choosing and Showing Up for 
* 										Job Outside the Home 
* ------------------------------------------------------------------------------
	
	* 9(a): Job Sign-up
	
	bysort treated: egen outside_job_takeup_mean = mean(outside_job_takeup)
	
	* change to percentage
	replace outside_job_takeup_mean = outside_job_takeup_mean * 100
	
	gen outside_job_takeup_hi = .
	gen outside_job_takeup_lo = .
	
	forval i = 0/1 {
	sum outside_job_takeup if treated == `i'
	local mean`i': di %4.2f `r(mean)'*100
	ci proportions outside_job_takeup if treated == `i'
	replace outside_job_takeup_hi =`r(ub)'*100 if treated == `i'
	replace outside_job_takeup_lo =`r(lb)'*100 if treated == `i'
}
	
	prtest outside_job_takeup, by(treated) 
	local p2 = 2*(normprob(-abs(`r(z)')))
	local pvalue = `p2'/2
	local pvalue: di %4.3f `pvalue'
	display `pvalue'

	twoway (bar outside_job_takeup_mean treated, color(gs13) lcolor(black) lwidth(medium) barwidth(0.7)) ///
			(rcap outside_job_takeup_hi outside_job_takeup_lo treated, lcolor(black)), legend(off) ///
			graphregion(color(white)) ///
			title("Panel A. Outside-home job take-up", position(11)) ///
			ytitle("Percent sign up for outside-home job", margin(small)) yscale(range(0 40)) ylabel(0(10)40) ///
			xlabel(0 "Control" 1 "Treatment", labgap(small) nogrid) ///
			xtitle(" ", margin(small)) aspectratio(0.7) ///
			text(5 0 "`mean0'%", size(med)) text(5 1 "`mean1'%", size(med)) ///
			text(40 0.5 "{it:p}-value = `pvalue'", size(med)) ///
			ylabel(,glpattern(solid )glcolor(gs15) glwidth(thin) grid) ///
			name(job_takeup, replace)
	

	* 9(b) Job Take-up

	bysort treated: egen outside_job_takeupXshow_up_mean = mean(outside_job_takeupXshow_up)
	
	* change to percentage
	replace outside_job_takeupXshow_up_mean =outside_job_takeupXshow_up_mean * 100
	
	gen outside_job_takeupXshow_up_hi = .
	gen outside_job_takeupXshow_up_lo = .
	
	forval i = 0/1 {
	sum outside_job_takeupXshow_up if treated == `i'
	local mean`i': di %4.2f `r(mean)'*100
	ci proportions outside_job_takeupXshow_up if treated == `i'
	replace outside_job_takeupXshow_up_hi =`r(ub)'*100 if treated == `i'
	replace outside_job_takeupXshow_up_lo =`r(lb)'*100 if treated == `i'
}

	prtest outside_job_takeupXshow_up, by(treated) 
	local p2 = 2*(normprob(-abs(`r(z)')))
	local pvalue = `p2'/2
	local pvalue: di %4.3f `pvalue'

	twoway (bar outside_job_takeupXshow_up_mean treated, color(gs13) lcolor(black) lwidth(medium) barwidth(0.7)) ///
			(rcap outside_job_takeupXshow_up_hi outside_job_takeupXshow_up_lo treated, lcolor(black)), legend(off) ///
			graphregion(color(white)) ///
			title("Panel B. Outside-home job show up", position(11)) ///
			ytitle("Percent show up for outside-home job", margin(small)) yscale(range(0 35)) ylabel(0(10)35) ///
			xlabel(0 "Control" 1 "Treatment", labgap(small) nogrid) ///
			xtitle(" ", margin(small)) aspectratio(0.7) ///
			text(5 0 "`mean0'%", size(med)) text(5 1 "`mean1'%", size(med)) ///
			text(32 0.5 "{it:p}-value = `pvalue'", size(med)) ///
			ylabel(,glpattern(solid )glcolor(gs15) glwidth(thin) grid) ///
			name(job_showup, replace)
			
* ------------------------------------------------------------------------------
* Export Figures
* ------------------------------------------------------------------------------
	graph combine job_takeup job_showup, col(2) altshrink
	graph export "${outdir}/figure9.pdf", replace
