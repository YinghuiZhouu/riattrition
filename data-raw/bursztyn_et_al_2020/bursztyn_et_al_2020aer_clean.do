*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	replication for: Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	prepared by: Z. Yu
*	6/19/2026
*	-------
*	Table 3: Effect of Treatment on Labor Supply Outcomes: Follow-Up
* treatment: condition2
* outcome: klk_index employed_now_out_fl applied_out_fl interviewed_out_fl driving_fl outside_others_fl2
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off
	
* ------------------------------------------------------------------------------
* change main directory
* ------------------------------------------------------------------------------

	global maindirectory = "/Users/peizan/github/riattrition/data-raw/bursztyn_et_al_2020" 
	cd "${maindirectory}"
	global datadir  = "${maindirectory}/data"

* ------------------------------------------------------------------------------
* partially replicate Panel A in Table 3
* ------------------------------------------------------------------------------

	use "${datadir}/clean_data/02_follow_up_clean.dta", clear 
	
	
	// generate Kling, Liebman, Katz (2007) index
	local index employed_now_out_fl ///
				applied_out_fl ///
				interviewed_out_fl interview_sched_out_fl ///
				driving_fl outside_others_fl2
	
	foreach var of varlist `index'{
		egen std_`var' = std(`var')
	}
	
	gen klk_index = std_employed_now_out_fl+std_applied_out_fl+std_interviewed_out_fl+ ///
		std_interview_sched_out_fl + std_driving_fl + std_outside_others_fl2
	replace klk_index = klk_index / 6
	
	*log using table_3_lee_bounds.log, replace
	
	// local outcome variables in Table 3
	local indvars klk_index employed_now_out_fl applied_out_fl interviewed_out_fl driving_fl outside_others_fl2
	
	foreach var in `indvars'{
		
		// regression estimates for ATE, outcomes are `indvars', treatment is condition2
		reg `var' condition2 if matched == 1, robust
		
		// estimated Lee bounds and inferential results
		leebounds `var' condition2, cie
	
	}
	
	*log close

* ------------------------------------------------------------------------------
* diagnose missing pattern: treatment groups has more missing
* ------------------------------------------------------------------------------

foreach var in klk_index employed_now_out_fl applied_out_fl interviewed_out_fl driving_fl outside_others_fl2{
	
	// gen missing indicator
	gen `var'_missing = 0
	replace `var'_missing = 1 if `var' == .
	
	// sum stats on missing indicator
	tabstat `var'_missing, by(condition2) statistics(n mean sd min max)
	ttest `var'_missing, by(condition2) unequal
	
	
}
