*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	Prepared by: R. Han & R. Wu
*	5/6/2020
*	-------
*	Follow-Up Experiment Data Cleaning
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	use "${datadir}/raw_data/01_main_exp_raw2.dta"
	
	*** merge

*	duplicates
	//people with same last 3 digits in same session in main exp
	//unique index that will allow the merge
	gen dup_index = _n + 1000
	replace session = dup_index if session == 1 	& digits == 136
	replace session = dup_index if session == 5 	& digits == 430
	replace session = dup_index if session == 6 	& digits == 651
	replace session = dup_index if session == 6 	& digits == 922
	replace session = dup_index if session == 7 	& digits == 646
	replace session = dup_index if session == 9 	& digits ==  16
	replace session = dup_index if session == 9 	& digits == 449
	replace session = dup_index if session == 12 	& digits == 282
	replace session = dup_index if session == 13 	& digits == 880
	replace session = dup_index if session == 15 	& digits ==  84
	replace session = dup_index if session == 17 	& digits == 341
	replace session = dup_index if session == 17 	& digits == 966
	
*	merge with follow-up

	//note: follow-up vars end with _fl

	merge 1:1 session digits using "${datadir}/raw_data/02_follow_up_raw1.dta"
	//(389 matches ; 111 no follow-up)
	
	gen matched = (_merge == 3)
	drop _merge 
	
	merge 1:1 session digits using "${datadir}/raw_data/02_follow_up_raw2.dta"
	
	replace matched = 1 if _merge == 3
	drop _merge
	
	replace work_outside_1_fl = work_outside_1 if mi(work_outside_1_fl)
	replace work_semiseg_1_fl = work_semiseg_1 if mi(work_semiseg_1_fl)
	drop work_outside_1 work_semiseg_1
	
	replace matched = 0 if session > 1000 & !mi(session)
	
	//note: 8 observations still missing 2-order beliefs in follow-up. They probably refused. Mark them as no follow-up.
	// total no-follow-up : 381
	replace matched = 0 if mi(work_semiseg_1_fl)
	
***	clean
	
**	recode

*	main exp
	
	//rename beliefs of labor demand
	rename q45_1_p2 labor_demand_guess

	//generating signup variable with wife's number
	gen signed_up_number=signed_up
	replace signed_up_number=0 if no_wife_number==1

	//moving no condition people to control and including them
	gen condition2=condition
	replace condition2=0 if condition==.

	//(from now on, I will use condition2 (results are robust to using condition))

	//replacing weird ages:
	replace age=. if age<18 | age>36
	
	gen college_deg = (education>=6)

	// PS: 3% of the sample reported to be single or divorced. We can drop them and things work just the same.

	//married
	gen married = (marital == 2)
	
	//absolute to percent
	gen num_know_per = num_know/30
	gen num_mfs_per = num_mutual_friends/30
	
	//wedge
	local statements outside semiseg mwage
	egen stag = tag(session)
	foreach q of local statements{
		gen `q'_guess = `q'_others + `q'_self
		bys session: egen `q'_objective = mean(`q'_self)
		replace `q'_objective = `q'_objective*30
		// wedge = guess - truth
		gen `q'_wedge = `q'_guess - `q'_objective
	}
	
*	follow-up

	foreach var of varlist employed_3mos_fl - driving_fl{
		replace `var' = 0 if `var' == 2
	}
	
	rename work_outside_1_fl outside_others_fl
	gen outside_others_fl2 = outside_others_fl / 30
	gen outside_others_2 = outside_others / 30
	gen min_wage_fl_per = min_wage_1_fl/30

* 	keep relevant variables 

	keep applied_out_fl driving_fl interviewed_out_fl employed_now_out_fl outside_objective ///
		 outside_guess outside_others_fl min_wage_fl_per children outside_others_fl2 ///
		 interview_sched_out_fl matched session *_self *_others employed_wife employed_now ///
		 education age college_deg num_know_per num_mfs_per signed_up_number outside_wedge ///
		 condition2 outside_confidence employed_3mos_out_fl
	
	label var age 						"Age"
	label var education 				"Highest education"
	label var employed_now 				"Currently employed"
	label var employed_wife				"Wife currently employed"
	label var employed_3mos_out_fl		"Wife employed outside the home three months before the follow-up"
	label var children 					"Number of children in the hh"
	label var haafez_self 				"Own belief of the unemployment insurance system"
	label var vjobs_self 				"Own belief of nationals' privileged access to job vacancies"
	label var mwage_self				"Own belief of the minimum wage"
	label var mwage_others 				"Guess of others' beliefs of the minimum wage"
	label var outside_self 				"Own belief of WWOH"
	label var outside_others 			"Guess of others' beliefs of WWOH"
	label var outside_confidence		"Confidence of guess of others' beliefs of WWOH" 
	label var semiseg_self 				"Own belief of women working in semi-segregated environment"
	label var semiseg_others			"Guess of others' belief of women working in semi-segregated environment"
	label var session 					"Session"
	label var signed_up_number 			"Signed up the job for wife"
	label var condition2 				"Treatment"
	label var college_deg 				"College degree or higher"
	label var num_know_per				"% participants known"
	label var num_mfs_per 				"% mututal friends"
	label var outside_guess 			"Guess of others' beliefs of WWOH (including self)"
	label var outside_objective 		"True # of people support WWOH"
	label var outside_wedge 			"Difference between guess and truth of WWOH"
	label var matched 					"Matched respondents"
	label var employed_now_out_fl		"Currently employed outside the home in the follow-up"
	label var applied_out_fl 			"Applied for job outside the home in the follow-up"
	label var interviewed_out_fl 		"Interviewed for job outside the home in the follow-up"
	label var interview_sched_out_fl 	"Interview scheduled for job outside the home in the follow-up"
	label var driving_fl 				"Driving lessons in the follow-up"
	label var outside_others_fl 		"Guess of others' beliefs of WWOH in the follow-up"
	label var outside_others_fl2		"Guess of % people support WWOH in the follow-up"
	label var min_wage_fl_per			"Guess of % people agree with the minimum wage statement in the follow-up"
	
	save "${datadir}/clean_data/02_follow_up_clean.dta", replace
