*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Replication for: Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	-------
*	1st Online Survey Figures 
* 	Figure 6: Misperceptions about Others’ Beliefs
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	use "${datadir}/clean_data/03_1st_online_svy_clean.dta", clear 
		
* ------------------------------------------------------------------------------
* 	1st National Survey - Figure 6: Misperceptions about Others’ Beliefs 
* ------------------------------------------------------------------------------
	
	*	wedges
	su c_outside_self
	local c_outside_mean `r(mean)'
	
	su t_list
	local t_list_mean `r(mean)'
	su c_list
	local c_list_mean `r(mean)'
	local t_outside_mean = `t_list_mean' - `c_list_mean'
	
	preserve
	
	cumul c_outside_guess_frac, gen(c_outside_cdf)
	cumul t_outside_guess_frac, gen(t_outside_cdf)
	stack c_outside_cdf c_outside_guess_frac  t_outside_cdf t_outside_guess_frac , ///
		into(cdf resp_share) wide clear
	
	sort resp_share c_outside_cdf t_outside_cdf 
	
	twoway (line c_outside_cdf t_outside_cdf resp_share, lpattern("l" "-") lwidth(medthin medthin) lcolor(gs0 maroon) c(L) sort), ///
		xtitle("Share of others agreeing") ///
		ytitle("Cumulative Probability") ///
		legend(label(1 "Control (perceptions about others' {it:answers})") label(2 "Treatment (perceptions about others' {it:beliefs})") region(lcolor(black)) cols(1) position(11) ring(0)) ///
		graphregion(color(white)) ylabel(,glcolor(gs15) glwidth(vthin)) xlabel(,grid glcolor(gs15) glwidth(vthin)) ///
		xline(`c_outside_mean', lpattern(l) lcolor(gs0) lwidth(thin)) ///
		xline(`t_outside_mean', lpattern(-) lcolor(maroon) lwidth(thin)) ///
		ylabel(,glpattern(solid )glcolor(gs15) glwidth(thin) grid) ///
		text(0.1 0.83 "True proportion (control)",size(vsmall) placement(right) color(gs0)) ///
		text(0.1 0.79 "True proportion (treatment)",size(vsmall) placement(left) color(maroon)) ///
		aspect(0.6) ///
		name("misp_beliefs", replace)
		
	restore


* ------------------------------------------------------------------------------
* Export Figures
* ------------------------------------------------------------------------------

	graph export "${outdir}/figure6.pdf",    name(misp_beliefs) replace	
