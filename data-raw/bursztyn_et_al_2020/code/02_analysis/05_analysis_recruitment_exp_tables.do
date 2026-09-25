*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Replication for: Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	-------
*	Recruitment Experiment Tables
* 	Table 6: Share of Women Choosing and Showing Up for Job Outside the Home 
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	set seed 960713
	local permutations 1000
	local bootstraps 1000
	
* ------------------------------------------------------------------------------
* 	Recruitment Experiment - Table 6: Share of Women Choosing and Showing Up for 
*									  Job Outside the Home 
* ------------------------------------------------------------------------------	
	
	use "${datadir}/clean_data/05_recruitment_exp_clean.dta", clear 
	
	file open table using "${outdir}/table6.tex", write replace
	
	local fwt "file write table"
	
	local cols 4
	local header "\begin{tabular}{@{\extracolsep{0.1cm}}l*{`cols'}{c}}\toprule"
	local footer "\end{tabular}"
	`fwt' "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" _n
	`fwt' "`header'" _n
	`fwt' "& Outside-home & Outside-home & Outside-home & Outside-home  \\" _n 
	`fwt' "& sign up & sign up & show up & show up \\" _n

	file close table
	
	local contlist age married num_child i.education employed_now 

	//Column 1
		eststo col1: reg outside_job_takeup treated, robust
		estadd local inference ""
		local t = _b[treated]/_se[treated]
		local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
		estadd local p_robust `p_robust'

		
		qui permute treated conditionp=(abs(_b[treated]/_se[treated])), reps(`permutations'): reg outside_job_takeup treated, robust 
		mat p = r(p_upper)
		local p_permute: di %4.3f p[1,1]
		estimates restore col1
		estadd local p_permute `p_permute'
		
		estadd local controls ""

	//Column 2
		eststo col2: reg outside_job_takeup treated `contlist', robust
		estadd local inference ""
		local t = _b[treated]/_se[treated]
		local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
		estadd local p_robust `p_robust'
		
		
		qui permute treated conditionp=(abs(_b[treated]/_se[treated])), reps(`permutations'): reg outside_job_takeup treated `contlist', robust 
		mat p = r(p_upper)
		local p_permute: di %4.3f p[1,1]
		estimates restore col2
		estadd local p_permute `p_permute'
		
		estadd local controls "\checkmark"

	//Column 3
		eststo col3: reg outside_job_takeupXshow_up treated, robust
		estadd local inference ""
		local t = _b[treated]/_se[treated]
		local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
		estadd local p_robust `p_robust'

		qui permute treated conditionp=(abs(_b[treated]/_se[treated])), reps(`permutations'): reg outside_job_takeupXshow_up treated, robust 
		mat p = r(p_upper)
		local p_permute: di %4.3f p[1,1]
		estimates restore col3
		estadd local p_permute `p_permute'
		
		estadd local controls ""

	//Column 4
		eststo col4: reg outside_job_takeupXshow_up treated `contlist', robust
		estadd local inference ""
		local t = _b[treated]/_se[treated]
		local p_robust: di %4.3f 2*ttail(e(df_r),abs(`t'))
		estadd local p_robust `p_robust'
		
		qui permute treated conditionp=(abs(_b[treated]/_se[treated])), reps(`permutations'): reg outside_job_takeupXshow_up treated `contlist', robust 
		mat p = r(p_upper)
		local p_permute: di %4.3f p[1,1]
		estimates restore col4
		estadd local p_permute `p_permute'
		
		estadd local controls "\checkmark"
		
		
		esttab using "${outdir}/table6.tex", append se noobs booktabs nostar ///
		fragment drop(age married num_child *education employed_now) ///
		nobaselevels coeflabels(treated "Treatment (\$\beta\$)" _cons "Constant") ///
		scalars("inference \textit{Inference robustness} (\$\beta\$)" "p_robust \qquad \emph{p}-value: robust SE" ///
		"p_permute \qquad \emph{p}-value: permutation test" ///
		"controls Controls" "N Observations" "r2 \$R^2\$") ///
		order(treated _cons) nomtitles 
		
		eststo clear
		
		file open table using "${outdir}/table6.tex" , write append

		`fwt' "\bottomrule" _n
		`fwt' "`footer'" _n
		
		file close table
	
