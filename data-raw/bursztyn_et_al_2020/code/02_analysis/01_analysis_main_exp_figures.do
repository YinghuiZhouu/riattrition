*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Replication for: Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	-------
*	Main Experiment Figures 
*	Figure 2: Wedges in Perceptions of Others’ Beliefs (WWOH)
* 	Figure 3: Job-Matching Service Sign-Up 
* 	Figure 5: Job-Matching Service Sign-Up: Heterogeneity by Wedge 
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	use "${datadir}/clean_data/01_main_exp_clean.dta", clear 
	
* ------------------------------------------------------------------------------
* 	Main Experiment - Figure 2: Wedges in Perceptions of Others’ Beliefs 
* ------------------------------------------------------------------------------

	local q "outside"
	
	gen `q'_guess_share = `q'_guess / 30
	gen `q'_others_truth = `q'_objective / 30 if stag == 1
	su `q'_self
	local `q'_mean `r(mean)'

	gen `q'_wedge_share = (`q'_wedge / 30) * 100

	twoway (hist `q'_wedge_share, lwidth(med) lcolor(gs0) fcolor(gs13) width(10) barwidth(10.5) start(-100)), ///
		xtitle("Wedge (guess % - objective %)") ///
		ytitle("Density", height(5)) ///
		xlabel(-100(20)100, nogrid) ///
		graphregion(color(white)) ylabel(,glcolor(gs15) glpattern(solid) glwidth(thin) grid) ///
		xline(0, lpattern(_) lcolor(cranberry) lwidth(medthin)) ///
		aspect(0.6) ///
		name("`q'_wedges", replace)
	
	
		
* ------------------------------------------------------------------------------
* 	Main Experiment - Figure 3: Job-Matching Service Sign-Up 
* ------------------------------------------------------------------------------
	
	bys condition2: egen sign_up_rate = mean(signed_up_number)
	replace sign_up_rate = sign_up_rate*100

	gen rate_lo = .
	gen rate_hi = .

	forval i=0/1{
		su sign_up_rate if condition2==`i'
		local rate_`i': di %4.2f `r(mean)'
		ci proportions signed_up_number if condition2==`i'
		replace rate_hi =`r(ub)'*100 if condition2 == `i'
		replace rate_lo =`r(lb)'*100 if condition2 == `i'
	}

	prtest signed_up_number, by(condition2)
	local p2 = 2*(normprob(-abs(`r(z)')))
	local p1 = `p2'/2
	local p1: di %4.3f `p1'

	twoway (bar sign_up_rate condition2, color(gs13) lcolor(black) lwidth(medium) barwidth(0.7)) ///
		(rcap rate_hi rate_lo condition2, lcolor(black)), legend(off) ///
		graphregion(color(white)) ///
		ytitle("Percent sign up", margin(small)) yscale(range(0 40)) ylabel(0(10)40) ///
		xlabel(0 "Control" 1 "Treatment", labgap(small) nogrid) ///
		xtitle(" ", margin(small)) aspectratio(0.7) ///
		text(5 0 "`rate_0'%", size(small)) text(5 1 "`rate_1'%", size(small)) ///
		text(38 0.5 "{it:p}-value = `p1'", size(small)) ///
		ylabel(,glpattern(solid) glcolor(gs15) glwidth(thin) grid) ///
		name("signup", replace)	
		
* ------------------------------------------------------------------------------
* 	Main Experiment - Figure 5: Job-Matching Service Sign-Up: Heterogeneity by Wedge 
* ------------------------------------------------------------------------------
	
	local q "outside"

	gen rate_lo2_`q' = .
	gen rate_hi2_`q' = .

	gen condition_grade_`q' = condition2 if `q'_wedge <= 0
	replace condition_grade_`q' = 3+condition2 if `q'_wedge > 0

	bys condition_grade_`q': egen sign_up_rate_`q' = mean(signed_up_number)
	replace sign_up_rate_`q' = sign_up_rate_`q' * 100

	forval i=0/1{
		//sign up by treatment, by wedge
		su signed_up_number if condition2==`i' & `q'_wedge <=0
		local rate_pes_`i': di %4.2f `r(mean)'*100
		su signed_up_number if condition2==`i' & `q'_wedge > 0
		local rate_opt_`i': di %4.2f `r(mean)'*100
		//confidence intervals
		ci proportions signed_up_number if condition2==`i' & `q'_wedge <=0
		replace rate_hi2_`q' =`r(ub)'*100 if condition2 == `i' & `q'_wedge <=0
		replace rate_lo2_`q' =`r(lb)'*100 if condition2 == `i' & `q'_wedge <=0
		ci proportions signed_up_number if condition2==`i' & `q'_wedge >0
		replace rate_hi2_`q' =`r(ub)'*100 if condition2 == `i' & `q'_wedge >0
		replace rate_lo2_`q' =`r(lb)'*100 if condition2 == `i' & `q'_wedge >0
	}

	prtest signed_up_number if `q'_wedge<=0, by(condition2) 
	local p2 = 2*(normprob(-abs(`r(z)')))
	local p1 = `p2'/2
	local p1_pes: di %4.3f `p1'

	prtest signed_up_number if `q'_wedge>0, by(condition2) 
	local p2 = 2*(normprob(-abs(`r(z)')))
	local p1 = `p2'/2
	local p1_opt: di %4.3f `p1'

	twoway (bar sign_up_rate_`q' condition_grade_`q', color(gs13) lcolor(black) lwidth(medium) barwidth(0.9)) ///
		(rcap rate_hi2_`q' rate_lo2_`q' condition_grade_`q', lcolor(black)), legend(off) ///
		graphregion(color(white)) ///
		ytitle("Percent sign up", height(5)) yscale(range(0 60)) ylabel(0(10)50) ///
		xlabel(0 "Control" 1 "Treatment" 3 "Control" 4 "Treatment", labsize(small) labgap(vsmall) nogrid) ///
		xtitle(" ", margin(small)) aspectratio(0.7) xline( 2, lpattern(-) lcolor(gs8)) xscale(range(-1 5)) ///
		text(10 0 "`rate_pes_0'%", size(small)) text(10 1 "`rate_pes_1'%", size(small)) ///
		text(10 3 "`rate_opt_0'%", size(small)) text(10 4 "`rate_opt_1'%", size(small)) ///
		text(55 0.5 "Wedge {&le} 0", box margin(small) bcolor(gs14%60)) text(55 3.5 "Wedge > 0", box margin(small) bcolor(gs14%60)) ///
		text(50 0.5 "{it:p}-value = `p1_pes'", size(small)) text(50 3.5 "{it:p}-value = `p1_opt'", size(small)) ///
		ylabel(,glpattern(solid) glcolor(gs15) glwidth(thin) grid) ///
		name("signup_`q'_het", replace)
	
	
* ------------------------------------------------------------------------------
* Export Figures
* ------------------------------------------------------------------------------

	graph export "${outdir}/figure2.pdf",    name(outside_wedges) replace
	graph export "${outdir}/figure3.pdf",    name(signup) replace
	graph export "${outdir}/figure5.pdf",    name(signup_outside_het) replace

