*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	Prepared by: R. Han & R. Wu
*	5/6/2020
*	-------
*	First National Survey Data Cleaning
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	use "${datadir}/raw_data/03_1st_online_svy_raw.dta" 
	
***	clean
	
*	drop respondents

	//36 did not finish
	drop if mi(c_outside_confidence_1) & mi(t_outside_confidence_1)

*	recode

	//treatment
	replace treatment = 0 if treatment == .
	
	//demographics
	
	// recode age (change arabic to numbers)
	replace age = "27" if age == "٢٧" 
	replace age = "28" if age == "٢٨"
	replace age = "29" if age == "٢٩"
	replace age = "30" if age == "٣٠"
	replace age = "31" if age == "٣١"
	replace age = "32" if age == "٣٢"
	replace age = "33" if age == "٣٣"
	replace age = "34" if age == "٣٤"
	replace age = "35" if age == "٣٥"
	replace age = "36" if age == "٣٦"
	replace age = "42" if age == "٤٢"
	
	replace age = "" if real(age) == .
	destring age, replace
	replace age=. if age<18 | age>36
	
	gen married = (marital == 2)
	
	gen college_deg = (education>=6)
	
	foreach var of varlist ever_employed-wife_employed_ft{
		replace `var' = 0 if `var' == 2
	}
	
	//correct wife employed outside of home 
	replace wife_employed_out = 0 if wife_employed == 0 
	
	//beliefs
	replace c_outside = 0 if c_outside == 2
	rename c_outside  c_outside_self
	rename c_outside_guess_1 c_outside_guess
	rename c_outside_confidence_1 c_outside_confidence
	rename t_outside_guess_1 t_outside_guess
	rename t_outside_confidence_1 t_outside_confidence
	
	//list
	replace t_list = t_list - 1
	replace c_list = c_list - 1
	
	//for consistency with main exp vars
	rename employed employed_now
	rename wife_employed employed_wife
	rename wife_employed_out employed_out_wife
	rename wife_employed_ft employed_ft_wife
	rename treatment condition
	
*	wedges

	su c_outside_self
	local c_outside_mean `r(mean)'
	
	su t_list
	local t_list_mean `r(mean)'
	su c_list
	local c_list_mean `r(mean)'
	local t_outside_mean = `t_list_mean' - `c_list_mean'
	
	gen t_outside_mean = `t_outside_mean'*100
	gen c_outside_mean = `c_outside_mean'*100
	
	gen c_outside_guess_frac = c_outside_guess / 100
	gen t_outside_guess_frac = t_outside_guess / 100
	gen c_outside_wedge = c_outside_guess - c_outside_mean
	gen t_outside_wedge = t_outside_guess - t_outside_mean
	
* 	keep relevant variables 
	keep age employed_now employed_wife employed_out_wife college_deg children condition ///
		 married c_list t_list c_outside_self c_outside_guess c_outside_confidence ///
		 t_outside_guess t_outside_confidence t_outside_mean c_outside_mean ///
		 c_outside_guess_frac t_outside_guess_frac c_outside_wedge t_outside_wedge
	
	label var age 					"Age"
	label var children 				"Number of Children in the hh"
	label var employed_now 			"Currently employed"
	label var employed_wife 		"Wife currently employed"
	label var employed_out_wife 	"Wife currently works outside the home"
	label var college_deg			"College degree or higher"
	label var married 				"Married"
	label var condition 			"Treatment"
	label var c_list 				"List of statements respondent agreed with in the control group"
	label var c_outside_self 		"Own belief of WWOH in the control group"
	label var c_outside_guess 		"Guess of others' beliefs of WWOH in the control group"
	label var c_outside_confidence  "Confidence about own guess in the control group"
	label var t_list 				"List of statements respondent agreed with in the treatment group"
	label var t_outside_guess 		"Guess of others' beliefs of WWOH in the treatment group"
	label var t_outside_confidence  "Confidence about own guess in the treatment group"
	label var t_outside_mean 		"True % of people support WWOH in the treatment group"
	label var c_outside_mean 		"True % of people support WWOH in the control group"
	label var c_outside_guess_frac  "Guess of fraction of people (out of 100) support WWOH in the control group"
	label var t_outside_guess_frac  "Guess of fraction of people (out of 100) support WWOH in the treatment group"
	label var c_outside_wedge 		"Difference between guess and truth in the control group"
	label var t_outside_wedge 		"Difference between guess and truth in the treatment group"
	
	save "${datadir}/clean_data/03_1st_online_svy_clean.dta", replace
