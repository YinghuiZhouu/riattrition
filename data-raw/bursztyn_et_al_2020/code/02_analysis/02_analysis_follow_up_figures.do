*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Replication for: Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	-------
*	Follow-up Experiment Figures 
* 	Figure 4: Long-Term Labor Supply Outcomes
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	use "${datadir}/clean_data/02_follow_up_clean.dta", clear 
	
* ------------------------------------------------------------------------------
* 	Follow-up - Figure 4: Long-Term Labor Supply Outcomes 
* ------------------------------------------------------------------------------

	//drop all unmatched responses (!!)
	keep if matched == 1
	
	//(for looping purposes)
	rename driving_fl driving_out_fl

	bys condition2: egen applied_rate = mean(applied_out_fl)
	bys condition2: egen interviewed_rate = mean(interviewed_out_fl)
	bys condition2: egen driving_rate = mean(driving_out_fl)
	bys condition2: egen employed_now_rate = mean(employed_now_out_fl)
	replace applied_rate = applied_rate*100
	replace interviewed_rate = interviewed_rate*100
	replace driving_rate = driving_rate*100
	replace employed_now_rate = employed_now_rate*100
	
	//y coord. for p-values
	local applied_p 22
	local interviewed_p 9
	local driving_p 82
	local employed_now_p 13.5

	foreach var in applied interviewed employed_now driving {
	gen `var'_rate_lo = .
	gen `var'_rate_hi = .
	
	forval i=0/1{
		su `var'_rate if condition2==`i'
		local `var'_rate_`i': di %4.2f `r(mean)'
		ci proportions `var'_out_fl if condition2==`i'
		replace `var'_rate_hi =`r(ub)'*100 if condition2 == `i'
		replace `var'_rate_lo =`r(lb)'*100 if condition2 == `i'
	}

	prtest `var'_out_fl, by(condition2)
	local p2 = 2*(normprob(-abs(`r(z)')))
	local p1 = `p2'/2
	local p1_`var': di %4.3f `p1'
	}
	
	* 4(a)
	twoway (bar applied_rate condition2, color(gs13) lcolor(black) lwidth(medium) barwidth(0.7) base(0)) ///
		(rcap applied_rate_hi applied_rate_lo condition2, lcolor(black)), legend(off) ///
		graphregion(color(white)) ytitle("Percent", margin(small)) ///
		title("Panel A. Applied for job", position(11)) ///
		xlabel(0 "Control" 1 "Treatment", labgap(small) nogrid) ///
		xtitle(" ", margin(small)) aspectratio(0.7) ///
		text(1 0.2 "`applied_rate_0'%", size(median)) text(1 1.2 "`applied_rate_1'%", size(median)) ///
		text(`applied_p' 0.5 "{it:p}-value = `p1_applied'", size(median)) ///
		ylabel(0(5)25,glpattern(solid )glcolor(gs15) glwidth(thin) grid) ///
		name("applied", replace)
	
	* 4(b)
	twoway (bar interviewed_rate condition2, color(gs13) lcolor(black) lwidth(medium) barwidth(0.7) base(0)) ///
		(rcap interviewed_rate_hi interviewed_rate_lo condition2, lcolor(black)), legend(off) ///
		graphregion(color(white)) ytitle("Percent", margin(small)) ///
		title("Panel B. Interviewed for job", position(11)) ///
		xlabel(0 "Control" 1 "Treatment", labgap(small) nogrid) ///
		xtitle(" ", margin(small)) aspectratio(0.7) ///
		text(0.5 0.2 "`interviewed_rate_0'%", size(median)) text(0.5 1.2 "`interviewed_rate_1'%", size(median)) ///
		text(`interviewed_p' 0.5 "{it:p}-value = `p1_interviewed'", size(med)) ///
		ylabel(0(2)10,glpattern(solid )glcolor(gs15) glwidth(thin) grid) ///
		name("interviewed", replace)
	
	* 4(c)
	twoway (bar employed_now_rate condition2, color(gs13) lcolor(black) lwidth(medium) barwidth(0.7) base(0)) ///
		(rcap employed_now_rate_hi employed_now_rate_lo condition2, lcolor(black)), legend(off) ///
		graphregion(color(white)) ytitle("Percent", margin(small)) ///
		title("Panel C. Employed", position(11)) ///
		xlabel(0 "Control" 1 "Treatment", labgap(small) nogrid) ///
		xtitle(" ", margin(small)) aspectratio(0.7) ///
		text(1 0.2 "`employed_now_rate_0'%", size(median)) text(1 1.2 "`employed_now_rate_1'%", size(median)) ///
		text(`employed_now_p' 0.5 "{it:p}-value = `p1_employed_now'", size(med)) ///
		ylabel(0(2)14,glpattern(solid )glcolor(gs15) glwidth(thin) grid) ///
		name("employed_now", replace)
	
	* 4(d)
	twoway (bar driving_rate condition2, color(gs13) lcolor(black) lwidth(medium) barwidth(0.7) base(0)) ///
		(rcap driving_rate_hi driving_rate_lo condition2, lcolor(black)), legend(off) ///
		graphregion(color(white)) ytitle("Percent", margin(small)) ///
		title("Panel D. Driving lessons", position(11)) ///
		xlabel(0 "Control" 1 "Treatment", labgap(small) nogrid) ///
		xtitle(" ", margin(small)) aspectratio(0.7) ///
		text(5 0.2 "`driving_rate_0'%", size(median)) text(5 1.2 "`driving_rate_1'%", size(median)) ///
		text(`driving_p' 0.5 "{it:p}-value = `p1_driving'", size(med)) ///
		ylabel(0(10)90,glpattern(solid )glcolor(gs15) glwidth(thin) grid) ///
		name("driving", replace)
		
		
* ------------------------------------------------------------------------------
* Export Figures
* ------------------------------------------------------------------------------

    graph combine applied interviewed employed_now driving, col(2) altshrink
	graph export "${outdir}/figure4.pdf", replace

	
	
