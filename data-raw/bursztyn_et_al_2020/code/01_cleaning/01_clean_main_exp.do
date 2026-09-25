*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Replication for: Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	-------
*	Main Experiment Data Cleaning
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	use "${datadir}/raw_data/01_main_exp_raw1.dta", clear 

*	recode

	//rename beliefs of labor demand
	rename q45_1_p2 labor_demand_guess
	replace labor_demand_guess=q45_4_p2 if mi(labor_demand_guess)
	replace labor_demand_guess=q46_1_p2 if mi(labor_demand_guess)
	replace labor_demand_guess=q43_1_p2 if mi(labor_demand_guess)

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
	
*	wedge
	
	local statements outside semiseg mwage
	
	egen stag = tag(session)
	
	foreach q of local statements{
		gen `q'_guess = `q'_others + `q'_self
		bys session: egen `q'_objective = mean(`q'_self)
		replace `q'_objective = `q'_objective*30
		// wedge = guess - truth
		gen `q'_wedge = `q'_guess - `q'_objective
	}
	
 	gen outside_wedge_pos = (outside_wedge > 0)
 	gen interaction = condition2*outside_wedge_pos
	
* 	drop irrelevant variables 
	drop condition marital employed_ever num_know num_mutual_friends q45_4_p2 ///
		q46_1_p2 q43_1_p2 count condition_txt glowork_choice second_order_total ///
		signed_up no_wife_number
		
	label var age 					"Age"
	label var education 			"Highest education"
	label var employed_now 			"Currently employed"
	label var employed_wife			"Wife currently employed"
	label var children 				"Number of children in the hh"
	label var haafez_self 			"Own belief of the unemployment insurance system"
	label var vjobs_self 			"Own belief of nationals' privileged access to job vacancies"
	label var mwage_self			"Own belief of the minimum wage"
	label var mwage_others 			"Guess of others' beliefs of the minimum wage"
	label var mwage_confidence 		"Confidence of guess of others' beliefs of the minimum wage"
	label var outside_self 			"Own belief of WWOH"
	label var outside_others 		"Guess of others' beliefs of WWOH"
	label var outside_confidence	"Confidence of guess of others' beliefs of WWOH" 
	label var semiseg_self 			"Own belief of women working in semi-segregated environment"
	label var semiseg_others		"Guess of others' belief of women working in semi-segregated environment"
	label var semiseg_confidence	"Confidence of guess of others' belief of women working in semi-segregated environment"
	label var labor_demand_guess 	"Guess the % firms with semi-segregated environment"
	label var session 				"Session"
	label var signed_up_number 		"Signed up the job for wife"
	label var condition2 			"Treatment"
	label var college_deg 			"College degree or higher"
	label var married 				"Married"
	label var num_know_per			"% participants known"
	label var num_mfs_per 			"% mututal friends"
	label var stag 					"Session tag"
	label var outside_guess 		"Guess of others' beliefs of WWOH (including self)"
	label var outside_objective 	"True # of people support WWOH"
	label var outside_wedge 		"Difference between guess and truth of WWOH"
	label var semiseg_guess 		"Guess of others' beliefs of semi-segregated environment (including self)"
	label var semiseg_objective 	"True # of people support semi-segregated environment"
	label var semiseg_wedge 		"Difference between guess and truth of semi-segregated environment"
	label var mwage_guess 		    "Guess of others' beliefs of minimum wage (including self)"
	label var mwage_objective 		"True # of people support minimum wage"
	label var mwage_wedge 			"Difference between guess and truth of minimum wage"
	label var outside_wedge_pos		"Positive wedge of WWOH"
	label var interaction 			"Interaction of treatment and positive wedge of WWOH"

	save "${datadir}/clean_data/01_main_exp_clean.dta", replace
