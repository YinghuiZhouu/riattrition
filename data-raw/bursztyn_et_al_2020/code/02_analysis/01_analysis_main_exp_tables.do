*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Replication for: Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	-------
*	Main Experiment Tables
* 	Table 1: Summary Statistics
*	Table 2: Job-Matching Service Sign-Up
* 	Table 4: Job-Matching Service Sign-Up: Heterogeneity by Wedge
* 	Table 5: Effect of Belief Update on Job-Matching Service Sign-Up
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	set seed 5052018
	
	//permutations and bootstraps
	local permutations 1000
	local bootstraps 1000
	
* ------------------------------------------------------------------------------
* 	Main Experiment - Table 1: Summary Statistics
* ------------------------------------------------------------------------------

	* retrospective follow up: wife working outside the home
	use "${datadir}/clean_data/02_follow_up_clean.dta", clear
	
	* variable of interest: employed_3mos_out_fl

	keep if matched == 1

	forval i=0/1{
			su employed_3mos_out_fl if condition2==`i'
				local employed_3mos_out_fl_`i': di %3.2f `r(mean)'*100
			}
			
	su employed_3mos_out_fl 
	local employed_3mos_out_fl_all: di %3.2f `r(mean)'*100
	
	* main experiment 
	
	use "${datadir}/clean_data/01_main_exp_clean.dta", clear 
	
	local sumvars "age college_deg children employed_now employed_wife num_know_per num_mfs_per"
		
	local age_label "Age"
	local children_label "Number of children"
	local college_deg_label "College degree (\%)"
	local employed_now_label "Employed (\%)"
	local employed_wife_label "Wife employed (\%)"
	local employed_3mos_out_fl_label "Wife working outside the home (\% retrospective follow-up)"
	local num_know_per_label "Other participants known (\%)"
	local num_mfs_per_label "Other participants with mutual friends (\%)"
	
	
	foreach var of varlist `sumvars'{
		//stats by treatment and control
		forval i=0/1{
			su `var' if condition2==`i'
			if inlist(`var', age, children){
				local `var'_`i': di %3.2f `r(mean)'
				local `var'_`i'_sd : di %3.2f `r(sd)'
			}
			else{
				local `var'_`i': di %3.2f `r(mean)'*100
				local `var'_`i'_sd : di %3.2f `r(sd)'*100
			}
		}

		//full sample
		su `var'
		if inlist(`var', age, children){
			local `var'_all: di %3.2f `r(mean)'
			local `var'_all_sd : di %3.2f `r(sd)'
		}
		else{
			local `var'_all: di %3.2f `r(mean)'*100
			local `var'_all_sd : di %3.2f `r(sd)'*100
		}
	}	
	
	//number of obs.
	forval i=0/1{
		count if condition2==`i'
		local count_`i' `r(N)'
	}

	local count_all = _N
	
	file open table using "${outdir}/table1.tex", write replace

	local fwt "file write table"
	local cols 3
	local header "\begin{tabular}{@{\extracolsep{0.1cm}}l*{`cols'}{c}} \toprule \midrule"
	local footer "\end{tabular}"
	
	`fwt' "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" _n
	`fwt' "`header'" _n
	`fwt' "& All & Control & Treatment \\" _n
	`fwt' "\midrule" _n

	`fwt' "Observations  & `count_all' & `count_0' & `count_1' \\" _n
	
	`fwt' "\midrule" _n

	foreach var of varlist age children{
		`fwt' "``var'_label' & ``var'_all' & ``var'_0' & ``var'_1' \\" _n
		`fwt' " & (``var'_all_sd') & (``var'_0_sd') & (``var'_1_sd') \\" _n
	}
	
	foreach var of varlist college_deg employed_now employed_wife{
		`fwt' "``var'_label' & ``var'_all' & ``var'_0' & ``var'_1' \\" _n
	}
	
	`fwt' "`employed_3mos_out_fl_label' & `employed_3mos_out_fl_all' & `employed_3mos_out_fl_0' & `employed_3mos_out_fl_1' \\" _n
	
	`fwt' "\midrule" _n
	
	foreach var of varlist num_know_per num_mfs_per{
		`fwt' "``var'_label' & ``var'_all' & ``var'_0' & ``var'_1' \\" _n
		`fwt' " & (``var'_all_sd') & (``var'_0_sd') & (``var'_1_sd') \\" _n
	}

	`fwt' "\bottomrule" _n
	`fwt' "`footer'" _n
	
	file close table
	
* ------------------------------------------------------------------------------
* 	Main Experiment - Table 2: Job-Matching Service Sign-Up 
* ------------------------------------------------------------------------------

	use "${datadir}/clean_data/01_main_exp_clean.dta", clear
	
	reg signed_up_number condition2 outside_wedge_pos interaction, robust
	boottest (interaction), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
	permute interaction conditionp=(abs(_b[interaction]/_se[interaction])), reps(`permutations'): reg signed_up_number condition2 outside_wedge_pos interaction, robust 
	
	reg signed_up_number condition2 outside_wedge_pos interaction i.session *_self *_others employed_wife employed_now i.education children num_know_per num_mfs_per age, robust
	boottest (interaction), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
	permute interaction conditionp=(abs(_b[interaction]/_se[interaction])), reps(`permutations'): reg signed_up_number condition2 outside_wedge_pos interaction i.session *_self *_others employed_wife employed_now i.education children num_know_per num_mfs_per age, robust 

	
	file open table using "${outdir}/table2.tex", write replace

	local fwt "file write table"
	
	local cols 4
	local header "\begin{tabular}{@{\extracolsep{0.1cm}}l*{`cols'}{c}} \toprule \midrule"
	local footer "\end{tabular}"
	`fwt' "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" _n
	`fwt' "`header'" _n
	
	file close table

	//Column 1
	eststo col1: reg signed_up_number condition2, robust
	estadd local inference ""
	local t = _b[condition2]/_se[condition2]
	local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
	estadd local p_robust `p_robust'

	qui boottest (condition2), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
	local p_wboot: di %4.3f `r(p)'
	estimates restore col1
	estadd local p_wboot `p_wboot'
	
	qui permute condition2 conditionp=(abs(_b[condition2]/_se[condition2])), reps(`permutations'): reg signed_up_number condition2, robust 
	mat p = r(p_upper)
	local p_permute: di %4.3f p[1,1]
	estimates restore col1
	estadd local p_permute `p_permute'
	
	estadd local fes ""
	estadd local baseline ""
	estadd local controls ""
	
	//Column 2
	eststo col2: reg signed_up_number condition2 i.session, robust
	estadd local inference ""
	local t = _b[condition2]/_se[condition2]
	local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
	estadd local p_robust `p_robust'
	
	qui boottest (condition2), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
	local p_wboot: di %4.3f `r(p)'
	estimates restore col2
	estadd local p_wboot `p_wboot'
	
	qui permute condition2 conditionp=(abs(_b[condition2]/_se[condition2])), reps(`permutations'): reg signed_up_number condition2 i.session, robust 
	mat p = r(p_upper)
	local p_permute: di %4.3f p[1,1]
	estimates restore col2
	estadd local p_permute `p_permute'
	
	estadd local fes "\checkmark"
	estadd local baseline ""
	estadd local controls ""
	
	//Column 3
	eststo col3: reg signed_up_number condition2 i.session *_self *_others, robust
	estadd local inference ""
	local t = _b[condition2]/_se[condition2]
	local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
	estadd local p_robust `p_robust'
	
	qui boottest (condition2), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
	local p_wboot: di %4.3f `r(p)'
	estimates restore col3
	estadd local p_wboot `p_wboot'
	
	qui permute condition2 conditionp=(abs(_b[condition2]/_se[condition2])), reps(`permutations'): reg signed_up_number condition2 i.session *_self *_others, robust 
	mat p = r(p_upper)
	local p_permute: di %4.3f p[1,1]
	estimates restore col3
	estadd local p_permute `p_permute'
	
	estadd local fes "\checkmark"
	estadd local baseline "\checkmark"
	estadd local controls ""
	
	//Column 4
	eststo col4: reg signed_up_number condition2 i.session *_self *_others employed_wife employed_now i.education children num_know_per num_mfs_per age, robust
	estadd local inference ""
	local t = _b[condition2]/_se[condition2]
	local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
	estadd local p_robust `p_robust'
	
	qui boottest (condition2), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
	local p_wboot: di %4.3f `r(p)'
	estimates restore col4
	estadd local p_wboot `p_wboot'
	
	qui permute condition2 conditionp=(abs(_b[condition2]/_se[condition2])), reps(`permutations'): reg signed_up_number condition2 i.session *_self *_others employed_wife employed_now i.education children num_know_per num_mfs_per age, robust 
	mat p = r(p_upper)
	local p_permute: di %4.3f p[1,1]
	estimates restore col4
	estadd local p_permute `p_permute'
	
	estadd local fes "\checkmark"
	estadd local baseline "\checkmark"
	estadd local controls "\checkmark"
	
	esttab using "${outdir}/table2.tex", append se noobs booktabs nostar /// 
	fragment drop(*session employed_* *education children age num_* *_self *_others) ///
	nobaselevels coeflabels(condition2 "Treatment (\$\beta\$)" _cons "Constant") ///
	scalars("inference \textit{Inference robustness} (\$\beta\$)" "p_robust \qquad \emph{p}-value: robust SE" ///
	"p_wboot \qquad \emph{p}-value: wild bootstrap" "p_permute \qquad \emph{p}-value: permutation test" ///
	"fes \midrule Session fixed effects" "baseline Baseline beliefs" "controls Controls" "N Observations" "r2 \$R^2\$") ///
	order(condition2 _cons) nomtitles
	
	eststo clear
	
	file open table using "${outdir}/table2.tex" , write append
	
	`fwt' "\bottomrule" _n
	`fwt' "`footer'" _n
	
	file close table
	
* ------------------------------------------------------------------------------
*   Main Experiment - Table 4: Job-Matching Service Sign-Up: Heterogeneity by Wedge
* ------------------------------------------------------------------------------
	
	file open table using "${outdir}/table4.tex", write replace

	local fwt "file write table"

	local cols 8
	local header "\begin{tabular}{@{\extracolsep{0.1cm}}l*{`cols'}{c}} \toprule \midrule"
	local footer "\end{tabular}"
	`fwt' "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" _n
	`fwt' "`header'" _n
	`fwt' "& \multicolumn{4}{c}{Wedge $\le$ 0} & \multicolumn{4}{c}{Wedge $>$ 0} \\ \cmidrule{2-5} \cmidrule{6-9}" _n

	file close table

	local j = 0
	
	forval i=0/1{
	
	//Column 1,5
	local j = `j'+1
	eststo col`j': reg signed_up_number condition2 if outside_wedge_pos == `i', robust
	estadd local inference ""
	local t = _b[condition2]/_se[condition2]
	local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
	estadd local p_robust `p_robust'

	qui boottest (condition2), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
	local p_wboot: di %4.3f `r(p)'
	estimates restore col`j'
	estadd local p_wboot `p_wboot'
	
	qui permute condition2 conditionp=(abs(_b[condition2]/_se[condition2])), reps(`permutations'): reg signed_up_number condition2 if outside_wedge_pos == `i', robust 
	mat p = r(p_upper)
	local p_permute: di %4.3f p[1,1]
	estimates restore col`j'
	estadd local p_permute `p_permute'
	
	estadd local fes ""
	estadd local baseline ""
	estadd local controls ""

	//Column 2,6
	local j = `j'+1
	eststo col`j': reg signed_up_number condition2 i.session if outside_wedge_pos == `i', robust
	estadd local inference ""
	local t = _b[condition2]/_se[condition2]
	local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
	estadd local p_robust `p_robust'
	
	qui boottest (condition2), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
	local p_wboot: di %4.3f `r(p)'
	estimates restore col`j'
	estadd local p_wboot `p_wboot'
	
	qui permute condition2 conditionp=(abs(_b[condition2]/_se[condition2])), reps(`permutations'): reg signed_up_number condition2 i.session if outside_wedge_pos == `i', robust 
	mat p = r(p_upper)
	local p_permute: di %4.3f p[1,1]
	estimates restore col`j'
	estadd local p_permute `p_permute'
	
	estadd local fes "\checkmark"
	estadd local baseline ""
	estadd local controls ""
	
	//Column 3,7
	local j = `j'+1
	eststo col`j': reg signed_up_number condition2 i.session *_self *_others if outside_wedge_pos == `i', robust
	estadd local inference ""
	local t = _b[condition2]/_se[condition2]
	local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
	estadd local p_robust `p_robust'
	
	qui boottest (condition2), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
	local p_wboot: di %4.3f `r(p)'
	estimates restore col`j'
	estadd local p_wboot `p_wboot'
	
	qui permute condition2 conditionp=(abs(_b[condition2]/_se[condition2])), reps(`permutations'): reg signed_up_number condition2 i.session *_self *_others if outside_wedge_pos == `i', robust 
	mat p = r(p_upper)
	local p_permute: di %4.3f p[1,1]
	estimates restore col`j'
	estadd local p_permute `p_permute'
	
	estadd local fes "\checkmark"
	estadd local baseline "\checkmark"
	estadd local controls ""
	
	//Column 4,8
	local j = `j'+1
	eststo col`j': reg signed_up_number condition2 i.session *_self *_others employed_wife employed_now i.education children num_know_per num_mfs_per age if outside_wedge_pos == `i', robust
	estadd local inference ""
	local t = _b[condition2]/_se[condition2]
	local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
	estadd local p_robust `p_robust'
	
	qui boottest (condition2), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
	local p_wboot: di %4.3f `r(p)'
	estimates restore col`j'
	estadd local p_wboot `p_wboot'
	
	qui permute condition2 conditionp=(abs(_b[condition2]/_se[condition2])), reps(`permutations'): reg signed_up_number condition2 i.session *_self *_others employed_wife employed_now i.education children num_know_per num_mfs_per age if outside_wedge_pos == `i', robust 
	mat p = r(p_upper)
	local p_permute: di %4.3f p[1,1]
	estimates restore col`j'
	estadd local p_permute `p_permute'
	
	estadd local fes "\checkmark"
	estadd local baseline "\checkmark"
	estadd local controls "\checkmark"
	}
	
	esttab using "${outdir}/table4.tex", append se noobs booktabs nostar ///
	fragment drop(*session employed_* *education children age num_* *_self *_others) ///
	nobaselevels coeflabels(condition2 "Treatment (\$\beta\$)" _cons "Constant") ///
	scalars("inference \textit{Inference robustness} (\$\beta\$)" "p_robust \qquad \emph{p}-value: robust SE" ///
	"p_wboot \qquad \emph{p}-value: wild bootstrap" "p_permute \qquad \emph{p}-value: permutation test" ///
	"fes \midrule Session fixed effects" "baseline Baseline beliefs" "controls Controls" "N Observations" "r2 \$R^2\$") ///
	order(condition2 _cons) nomtitles
	
	eststo clear
	
	file open table using "${outdir}/table4.tex" , write append
	
	`fwt' "\bottomrule" _n
	`fwt' "`footer'" _n
	
	file close table
	
* ------------------------------------------------------------------------------
* 	Main Experiment - Table 5: Effect of Belief Update on Job-Matching Service Sign-Up 
* ------------------------------------------------------------------------------	
	
	file open table using "${outdir}/table5.tex", write replace
	
	local fwt "file write table"
		
	local cols 3
	local header "\begin{tabular}{@{\extracolsep{0.1cm}}l*{`cols'}{c}} \toprule \midrule"
	local footer "\end{tabular}"
	`fwt' "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" _n
	`fwt' "`header'" _n
	
	file close table
	
	gen update = -outside_wedge if condition2 == 1
	replace update = 0 if condition2 == 0
	
	local spec1 "reg signed_up_number update outside_self outside_others outside_confidence, robust"
	local spec2 "reg signed_up_number update outside_self outside_others outside_confidence i.session, robust"
	local spec3 "reg signed_up_number update outside_self outside_others outside_confidence i.session employed_wife employed_now i.education children num_know_per num_mfs_per age, robust"

	forval i=1/3{
		eststo col`i': `spec`i''
		local t = _b[update]/_se[update]
		local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
		estadd local p_robust `p_robust'

		qui boottest (update), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
		local p_wboot: di %4.3f `r(p)'
		estimates restore col`i'
		estadd local p_wboot `p_wboot'
	
		qui permute update conditionp=(abs(_b[update]/_se[update])), reps(`permutations'): `spec`i'' 
		mat p = r(p_upper)
		local p_permute: di %4.3f p[1,1]
		estimates restore col`i'
		estadd local p_permute `p_permute'
	
		if `i' == 1{
		estadd local baseline "\checkmark"
		estadd local fes ""
		estadd local controls ""
		}
		if `i' == 2{
		estadd local baseline "\checkmark"
		estadd local fes "\checkmark"
		estadd local controls ""
		}
		if `i' == 3{
		estadd local baseline "\checkmark"
		estadd local fes "\checkmark"
		estadd local controls "\checkmark"
		}
	}
	
	esttab using "${outdir}/table5.tex", append se noobs booktabs nostar ///
	fragment drop(*session employed_* *education children age num_* *_self *_others *_confidence) ///
	nobaselevels coeflabels(update "Update (\$-\$Wedge\$\times\mathbf{1}_{\textit{Treatment}}\$;\$\beta\$)" _cons "Constant") ///
	scalars("inference \textit{Inference robustness} (\$\beta\$)" "p_robust \qquad \emph{p}-value: robust SE" ///
	"p_wboot \qquad \emph{p}-value: wild bootstrap" "p_permute \qquad \emph{p}-value: permutation test" ///
	"baseline \midrule Baseline beliefs and confidence" "fes Session fixed effects" "controls Controls" "N Observations" "r2 \$R^2\$") ///
	order(update _cons) nomtitles substitute (\_ _)
	
	eststo clear
	
	file open table using "${outdir}/table5.tex", write append
	
	`fwt' "\bottomrule" _n
	`fwt' "`footer'" _n
	
	file close table
