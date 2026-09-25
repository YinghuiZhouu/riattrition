*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Replication for: Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	-------
*	Arab Barometer Figure
* 	Figure 7: Support for WWOH Across Samples 
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
* ------------------------------------------------------------------------------
* 	Arab Barometer - Figure 7: Support for WWOH Across Samples 
* ------------------------------------------------------------------------------

***	arab barometer

	use "${datadir}/clean_data/06_arab_barometer_clean.dta", clear
	
	//sample mean: all ages
	su work_outside_ab [aw=wt]
	local ab_full_outside_mean = `r(mean)'*100
	ci means work_outside_ab [aw=wt]
		local ab_full_hi =`r(ub)'*100
		local ab_full_lo =`r(lb)'*100
	
	//ages 18-35
	su work_outside_ab if q1001<=35 & q1001>=18 [aw=wt]
	local ab_yng_outside_mean = `r(mean)'*100
	ci means work_outside_ab if q1001<=35 & q1001>=18 [aw=wt]
		local ab_yng_hi =`r(ub)'*100
		local ab_yng_lo =`r(lb)'*100
	
	//ages 36+
	su work_outside_ab if q1001>35 & !mi(q1001) [aw=wt]
	local ab_old_outside_mean = `r(mean)'*100
	ci means work_outside_ab if q1001>35 [aw=wt] 
		local ab_old_hi =`r(ub)'*100
		local ab_old_lo =`r(lb)'*100
	
***	1st Online Survey

	use "${datadir}/clean_data/03_1st_online_svy_clean.dta", clear 
	
*	wedges

	su c_outside_self
	local c_outside_mean `r(mean)'
	
	su t_list, d
	local t_list_mean `r(mean)'
	local t_list_mvar=`r(Var)'/`r(N)'
	su c_list, d
	local c_list_mean `r(mean)'
	local c_list_mvar=`r(Var)'/`r(N)'
	local t_outside_mean = `t_list_mean' - `c_list_mean'
	
	local t_outside_mean = `t_outside_mean'*100
	local c_outside_mean = `c_outside_mean'*100
	
	local natsvy_t_msd = sqrt(`t_list_mvar'+`c_list_mvar')*100
	local natsvy_t_hi = `t_outside_mean'+1.96*`natsvy_t_msd'
	local natsvy_t_lo = `t_outside_mean'-1.96*`natsvy_t_msd'
	
	ci proportions c_outside_self
		local natsvy_c_hi =`r(ub)'*100
		local natsvy_c_lo =`r(lb)'*100
		
*** 2nd Online Survey

	use "${datadir}/clean_data/04_2nd_online_svy_clean.dta", clear 
	
	sum belief_fob
	local tgm_mean = `r(mean)'
	ci means belief_fob
		local tgm_hi = `r(ub)'
		local tgm_lo = `r(lb)'
	
	
***	main experiment
	
	clear
	
	use "${datadir}/clean_data/01_main_exp_clean.dta", clear 
	
	keep outside_self 
	
	ci proportions outside_self
		local main_hi =`r(ub)'*100
		local main_lo =`r(lb)'*100
	collapse outside_self
	replace outside_self = outside_self * 100
	gen high_bound = `main_hi'
	gen low_bound = `main_lo'
	set obs `=7'
	su outside_self if _n==1
	local outside_mean `r(mean)'
	replace outside_self = `c_outside_mean' if _n==2
	replace high_bound = `natsvy_c_hi' if _n==2
	replace low_bound = `natsvy_c_lo' if _n==2
	replace outside_self = `t_outside_mean' if _n==3
	replace high_bound = `natsvy_t_hi' if _n==3
	replace low_bound = `natsvy_t_lo' if _n==3
	replace outside_self = `tgm_mean' if _n==4
	replace high_bound = `tgm_hi' if _n==4
	replace low_bound = `tgm_lo' if _n==4
	
	replace outside_self = `ab_full_outside_mean' if _n==5
	replace high_bound = `ab_full_hi' if _n==5
	replace low_bound = `ab_full_lo' if _n==5
	replace outside_self = `ab_yng_outside_mean' if _n==6
	replace high_bound = `ab_yng_hi' if _n==6
	replace low_bound = `ab_yng_lo' if _n==6
	replace outside_self = `ab_old_outside_mean' if _n==7
	replace high_bound = `ab_old_hi' if _n==7
	replace low_bound = `ab_old_lo' if _n==7
	
	gen exp = _n if _n <= 4
	replace exp = _n + 1 if _n > 4
	gen exp2 = _n if _n <= 4
	replace exp2 = _n - 4 if _n > 4
	
	gen exp_label = "Main Experiment" if exp == 1
	replace exp_label = "Online Survey (C)" if exp == 2
	replace exp_label = "Online Survey (T)" if exp == 3
	replace exp_label = "Follow-up Survey" if exp == 4
	replace exp_label = "Full Sample" if exp == 6
	replace exp_label = "Ages 18-35" if exp == 7
	replace exp_label = "Ages 36+" if exp == 8
	
	gen source_label = "Experiment" if exp <=4
	replace source_label = "Arab Barometer" if exp >=5
	
	gen source = 0 if exp <= 4
	replace source = 1 if exp >= 5
	
	format outside_self %5.2f
	
	twoway (bar outside_self exp if exp <=4, color(gs13) lcolor(black) lwidth(medium) barwidth(0.9)) ///
		(bar outside_self exp if exp >= 5, color(gs5) lcolor(black) lwidth(medium) barwidth(0.9)) ///
		(rcap high_bound low_bound exp, lcolor(black)), ///
		graphregion(color(white)) aspectratio(0.55) ///
		ylabel(0(10)100, glpattern(solid) glcolor(gs15) glwidth(thin) labsize(vsmall) format(%2.0f) grid) ///
		ytitle("Support for WWOH (percent)", size(small)) xtitle("") ///
		plotregion(margin(b = 1)) yscale(range(0 100)) ///
		xlabel(1 "Main Exp." 2 "Ntnl. Svy. (Direct)" 3 "Ntnl. Svy. (List)"  4 "Online Svy. 2" ///
		6 "All Males" 7 "Males 18-35" 8 "Males 36+", angle(45) labsize(vsmall) notick nogrid) ///
		legend(order(1 "Experimental Data" 2 "Arab Barometer") size(small) position(6) region(lcolor(black)) col(2)) ///
		text(93 1 "`:di %4.2f `outside_mean''", size(vsmall)) text(93 2 "`:di %4.2f `c_outside_mean''", size(vsmall)) text(93 3 "`:di %4.2f `t_outside_mean''", size(vsmall)) ///
		text(93 4 "`:di %4.2f `tgm_mean''", size(vsmall)) text(93 6 "`:di %4.2f `ab_full_outside_mean''", size(vsmall)) text(93 7 "`:di %4.2f `ab_yng_outside_mean''", size(vsmall)) ///
		text(93 8 "`:di %4.2f `ab_old_outside_mean''", size(vsmall)) ///
		name(flfp_compare,replace)

* ------------------------------------------------------------------------------
* Export Figure
* ------------------------------------------------------------------------------

	graph export "${outdir}/figure7.pdf",    name(flfp_compare) replace
	
	
