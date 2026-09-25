*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Replication for: Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	-------
*	2nd National Survey Figures 
* 	Figure 8a: Frequency of WWOH Discussion
*   Figure 8b: Perceived Support & Frequency of WWOH Discussion
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	use "${datadir}/clean_data/04_2nd_online_svy_clean.dta", clear 
	
* ------------------------------------------------------------------------------
*	2nd National Survey - Figure 8a: Frequency of WWOH Discussion 
* ------------------------------------------------------------------------------
	preserve 
	
	drop if discuss_freq == .  

	local N = _N 

	bysort discuss_freq: egen discuss = count(discuss_freq) 
	
	* proportion of each category 
	replace discuss = (discuss / `N') * 100

	forval i = 1/5 {
		sum discuss if discuss_freq == `i' 
		local mean`i': di %4.2f `r(mean)'
	}

	twoway (bar discuss discuss_freq, color(gs13) lcolor(black) lwidth(medium) barwidth(0.9)), ///
			graphregion(color(white)) ///
			title("Panel A. Frequency of discussion on WWOH", position(11)) ///
			ytitle("Percentage", margin(small)) yscale(range(0 40)) ylabel(0(10)40) ///
			xlabel(1 "Very often" 2 "Often" 3 "Sometimes" 4 "Rarely" 5 "Very Rarely", labgap(small) nogrid) ///
			text(5 1 "`mean1'%", size(med)) text(5 2 "`mean2'%", size(med)) text(5 3 "`mean3'%", size(med)) ///
			text(5 4 "`mean4'%", size(med)) text(5 5 "`mean5'%", size(med)) ///
			xtitle(" ", margin(small)) ///
			ylabel(,glpattern(solid )glcolor(gs15) glwidth(thin) grid) ///
			name(discuss_freq, replace)
	
	restore
	
* ------------------------------------------------------------------------------
* 	2nd National Survey - Figure 8b: Perceived Support & Frequency of WWOH Discussion 
* ------------------------------------------------------------------------------
	
	sum belief_fob, d
	local fob_mean = r(mean)

	* average perception by each category 
	bysort discuss_freq: egen sob_mean = mean(belief_sob)

	gen sob_hi = .
	gen sob_lo = . 

	* for the figure label
	gen discuss_label = discuss_freq - 1

	forval i = 0/4 {
		sum belief_sob if discuss_label == `i'
		local mean`i': di %4.2f `r(mean)'
		ci means belief_sob if discuss_label == `i'
		replace sob_hi =`r(ub)' if discuss_label == `i'
		replace sob_lo =`r(lb)' if discuss_label == `i'
	}

	twoway (bar sob_mean discuss_label, yline(`fob_mean', lpattern(dash) lcolor(black)) color(gs13) lcolor(black) lwidth(medium) barwidth(0.9)) ///
			(rcap sob_hi sob_lo discuss_label, lcolor(black)), legend(off) ///
			graphregion(color(white)) ///
			title("Panel B. Perceived support for WWOH", position(11)) ///
			ytitle("Perceived support for WWOH", margin(small)) yscale(range(0 100)) ylabel(0(20)100) ///
			xlabel(0 "Very often" 1 "Often" 2 "Sometimes" 3 "Rarely" 4 "Very rarely", labgap(small) nogrid) ///
			text(20 0 "`mean0'%", size(med)) text(20 1 "`mean1'%", size(med)) text(20 2 "`mean2'%", size(med)) ///
			text(20 3 "`mean3'%", size(med)) text(20 4 "`mean4'%", size(med)) ///
			xtitle(" ", margin(small)) ///
			ylabel(,glpattern(solid )glcolor(gs15) glwidth(thin) grid) ///
			name(perceived_support_discuss_freq, replace)

			
* ------------------------------------------------------------------------------
* Export Figures
* ------------------------------------------------------------------------------

	graph combine discuss_freq perceived_support_discuss_freq, col(2) altshrink
	graph export "${outdir}/figure8.pdf", replace
