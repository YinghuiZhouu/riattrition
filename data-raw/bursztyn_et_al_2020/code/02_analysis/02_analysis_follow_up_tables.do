*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Replication for: Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	-------
*	Follow-up Experiment Tables
*	Table 3: Effect of Treatment on Labor Supply Outcomes
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	set seed 5052018
	
	//permutations and bootstraps
	local permutations 1000
	local bootstraps 1000
	
* ------------------------------------------------------------------------------
*   Follow-up - Table 3: Effect of Treatment on Labor Supply Outcomes
* ------------------------------------------------------------------------------
	
	use "${datadir}/clean_data/02_follow_up_clean.dta", clear 
	
	local index employed_now_out_fl ///
				applied_out_fl ///
				interviewed_out_fl interview_sched_out_fl ///
				driving_fl outside_others_fl2
	
	//Kling, Liebman, Katz (2007) Index
	foreach var of varlist `index'{
		egen std_`var' = std(`var')
	}
	gen klk_index = std_employed_now_out_fl+std_applied_out_fl+std_interviewed_out_fl+ ///
		std_interview_sched_out_fl + std_driving_fl + std_outside_others_fl2
	replace klk_index = klk_index / 6
		
	local employed_now_out_fl_label "Currently Employed (\%)"
	local applied_out_fl_label "Applied to a Job (\%)"
	local interviewed_out_fl_label "Interviewed for a Job (\%)"
	local interview_sched_out_fl_label "Interview Scheduled (\%)"
	local driving_fl_label "Driving lessons sign-up (\%)"
	local outside_others_2_label "Baseline (\%)"
	local outside_others_fl2_label "Follow-up (\%)"
	local klk_index_label "Kling-Liebman-Katz Index"
	
	
*	outcomes : regression

	local indvars klk_index employed_now_out_fl applied_out_fl interviewed_out_fl driving_fl outside_others_fl2
	
	file open table using "${outdir}/table3.tex", write replace

	local fwt "file write table"
	local cols 6
	local header "\begin{tabular}{@{\extracolsep{0.1cm}}l*{`cols'}{c}} \toprule \midrule"
	local footer "\end{tabular}"
	`fwt' "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" _n
	`fwt' "`header'" _n
	`fwt' "& K-L-K Index & Employed & Applied & Interviewed & Driving & Beliefs about \\" _n
	`fwt' "& & & & & Lessons & Neighbors \\ \midrule" _n
	`fwt' "\multicolumn{7}{@{}l}{\textit{Panel A. No controls}} \\" _n
	
	file close table

	//Panel A
	local i = 0
	
	//corrected p-values from mthexp (enter results manually since the provided command does not seem to store the results)
	preserve
	keep if matched == 1
	mhtexp `index' if matched==1, treatment(condition2)
	restore
	
	matrix mht = (0, 0.484, 0.003, 0.05, 0.039, 0.000)
	
	foreach var in `indvars'{
		eststo col`++i': reg `var' condition2 if matched==1, robust
		local t = _b[condition2]/_se[condition2]
		local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
		estadd local p_robust `p_robust'

		qui boottest (condition2), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
		local p_wboot: di %4.3f `r(p)'
		estimates restore col`i'
		estadd local p_wboot `p_wboot'
	
		qui permute condition2 conditionp=(abs(_b[condition2]/_se[condition2])), reps(`permutations'): reg `var' condition2 if matched==1, robust 
		mat p = r(p_upper)
		local p_permute: di %4.3f p[1,1]
		estimates restore col`i'
		estadd local p_permute `p_permute'
		
		local p_mhtexp: di %4.3f mht[1,`i']
		if `i' > 1 estadd local p_mhtexp `p_mhtexp'
		if `i' == 1 estadd local p_mhtexp "--"
		
		leebounds `var' condition2
		mat b = e(b)
		mat V = e(V)
		estimates restore col`i'
		local lower: di %4.3f b[1,1]
		local upper: di %4.3f b[1,2]
		local lower_se: di %5.4f sqrt(V[1,1])
		local upper_se: di %5.4f sqrt(V[2,2])

		estadd local lower `lower'
		estadd local upper `upper'
		estadd local lower_se (`lower_se')
		estadd local upper_se (`upper_se')
	}
	
	esttab using "${outdir}/table3.tex", append se noobs booktabs nostar /// 
	fragment sfmt(a3) ///
	nobaselevels coeflabels(condition2 "Treatment (\$\beta\$)" _cons "Constant") ///
	scalars("inference \textit{Inference robustness} (\$\beta\$)" "p_robust \qquad \emph{p}-value: robust SE" ///
	"p_wboot \qquad \emph{p}-value: wild bootstrap" "p_permute \qquad \emph{p}-value: permutation test" ///
	"p_mhtexp \qquad \emph{p}-value: L-S-X MHT correction" ///
	"lee \midrule \textit{Lee attrition bounds}" "lower \qquad Lower bound:" "lower_se \qquad" "upper \addlinespace \qquad Upper bound:" "upper_se \qquad" ///
	"N \midrule Observations" "r2 \$R^2\$") ///
	order(condition2 _cons) nonumbers nomtitles
	
	eststo clear
	file open table using "${outdir}/table3.tex" , write append
	`fwt' "\bottomrule" _n
	`fwt' "\noalign{\vskip 2mm} " _n
	`fwt' "\multicolumn{7}{@{}l}{\textit{Panel B. Session fixed effects, baseline beliefs and socioeconomic controls}} \\" _n
	file close table
	
	//Panel B
	local i = 0
	foreach var in `indvars'{
		eststo col`++i': reg `var' condition2 i.session *_self *_others employed_wife employed_now i.education ///
			children num_know_per num_mfs_per age if matched==1, robust
		local t = _b[condition2]/_se[condition2]
		local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
		estadd local p_robust `p_robust'

		qui boottest (condition2), boottype(wild) reps(`bootstraps') weighttype(webb) cluster(session) nograph
		local p_wboot: di %4.3f `r(p)'
		estimates restore col`i'
		estadd local p_wboot `p_wboot'
	
		qui permute condition2 conditionp=(abs(_b[condition2]/_se[condition2])), reps(`permutations'): reg `var' condition2 i.session *_self *_others employed_wife employed_now i.education children num_know_per num_mfs_per age if matched==1, robust
		mat p = r(p_upper)
		local p_permute: di %4.3f p[1,1]
		estimates restore col`i'
		estadd local p_permute `p_permute'
	}
	
	esttab using "${outdir}/table3.tex", append se noobs booktabs nostar ///
	fragment drop(*session employed_* *education children age num_* *_self *_others) ///
	nobaselevels coeflabels(condition2 "Treatment (\$\beta\$)" _cons "Constant") ///
	scalars("inference \textit{Inference robustness} (\$\beta\$)" "p_robust \qquad \emph{p}-value: robust SE" ///
	"p_wboot \qquad \emph{p}-value: wild bootstrap" "p_permute \qquad \emph{p}-value: permutation test" ///
	"N \midrule Observations" "r2 \$R^2\$") ///
	order(condition2 _cons) nonumbers nomtitles
	
	eststo clear
	file open table using "${outdir}/table3.tex" , write append
	
	`fwt' "\bottomrule" _n
	`fwt' "`footer'" _n
	file close table
	
	
	
