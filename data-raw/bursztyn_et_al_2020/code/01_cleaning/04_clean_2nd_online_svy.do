*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	Prepared by: R. Han & R. Wu
*	5/6/2020
*	-------
*	Second National Survey Data Cleaning
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	use "${datadir}/raw_data/04_2nd_online_svy_raw.dta" 
	
	compress 
	
* clean and recode variables 

	* q4 age 
	rename q4 age
	
	* q5 marital status
	gen married = q5 == 2 if !missing(q5)
	
	* q6 education 
	recode q6 (8=1) (2=2) (9=3) (3=4) (4=5) (10=6) (13=7) (14=8), gen(education)
	label define edu 1 "read and write" 2 "elementary" 3 "intermediate" 4 "secondary" ///
		5 "college" 6 "bachelor" 7 "master" 8 "doctor"
	label values education edu 
	
	gen higher_edu = education > 4 if !missing(education)
	
	* q7 ever been employed
	gen employed_before = q7 == 1 if !missing(q7)
	
	* q8 currently employed
	gen employed_now = q8 == 1 if !missing(q8)
	
	* q9 outside home or at home job 
	gen work_outside = q9 == 2 if !missing(q9)
	
	* q10 hours of work 
	gen work_morethan30hr = q10 == 1 if !missing(q10)
	
	* q11 wife currently employed
	gen wife_work = q11 == 1 if !missing(q11)
	
	* q12 wife outside home or at home job (recode: 0 = not working/working from home, 1 = working outside home)
	gen wife_work_outside = q12 == 2 
		
	label define work 0 "Wife Not Work Outside Home" 1 "Wife Work Outside Home"
	label values wife_work_outside work
	
	* q13 wife hours of work
	gen wife_work_morethan30hr = q13 == 1 if !missing(q13)
	
	* q14 number of children
	gen num_child = q14 - 1 if !missing(q14)
	
	* q15 city
	rename q15 city
	label define place 1 "Riyadh" 2 "Jeddah" 3 "Dammam" 4 "Other"
	label values city place 
	
	* q17 randomization 1 belief
	gen belief1_fob = q17 == 1 if !missing(q17)
		replace belief1_fob = 100 if belief1_fob == 1 
		
	* q18 randomization 1 perception
	gen belief1_sob = q18 
	
	* q19 randomization 2 belief
	gen belief2_fob = q19 == 1 if !missing(q19)
		replace belief2_fob = 100 if belief2_fob == 1 
		
	* q20 randomization 2 perception
	gen belief2_sob = q20 
	
	* q21 frequency of discussion 
	recode q21 (1=1) (2=2) (3=3) (6=4) (7=5), gen(discuss_freq)
	label define freq 1 "Very often" 2 "Often" 3 "Sometimes" 4 "Rarely" 5 "Very rarely"
		
	label values discuss_freq freq 
	
	gen talk_rarely = discuss_freq == 4 | discuss_freq == 5 if !missing(discuss_freq)

	* pool randomization 1 & 2 results 
	local list "fob sob"
	foreach x in `list' {
		gen belief_`x' = belief1_`x'
			replace belief_`x' = belief2_`x' if belief1_`x' == . 
	}
	
	gen randomization = 1 
		replace randomization = 2 if fl_65_do == "Belief2"
	
	* truth = average fob 
	sum belief_fob, d 
	gen truth = r(mean)

	* wedge = perception - truth
	gen wedge = belief_sob - truth
	
* keep only relevant variables
	
	keep age married education higher_edu employed_now wife_work wife_work_outside ///
		 work_outside num_child belief1_fob belief1_sob belief2_fob belief2_sob ///
		 discuss_freq talk_rarely belief_fob belief_sob randomization truth wedge
		 
	label var age 					"Age"
	label var married 				"Married"
	label var education 			"Highest education"
	label var higher_edu 			"College degree or higher"
	label var employed_now 			"Currently employed"
	label var work_outside 			"Currently work outside the home"
	label var wife_work 			"Wife currently employed"
	label var wife_work_outside 	"Wife works outside the home"
	label var num_child				"Number of Children"
	label var belief1_fob 			"Own belief - randomization arm 1"
	label var belief1_sob 			"Guess of others - randomization arm 1"
	label var belief2_fob 			"Own belief - randomization arm 2"
	label var belief2_sob 			"Guess of others - randomization arm 2"
	label var discuss_freq 			"Frequency of WWOH discussion"
	label var talk_rarely 			"Rarely/very rarely discuss WWOH"
	label var belief_fob 			"Own belief"
	label var belief_sob 			"Guess of others"
	label var randomization 		"Randomization arm"
	label var truth 				"True % of people support WWOH"
	label var wedge 				"Difference between guess and truth"
	
	save "${datadir}/clean_data/04_2nd_online_svy_clean.dta" , replace
	
	
	
