* ============================================================================ *
* data produced in this do file: chen_yang_2019aer_clean.dta
* treatment: treatment_vpnnl
* outcome: az_overall_a_w3_miss, az_overall_b_w3_miss, az_overall_c_w3_miss, 
* az_overall_d_w3_miss, az_overall_e_w3_miss
* ============================================================================ *

*** Do-file: replication of Chen and Yang (2018)
*** Part 1: based on student survey data
*** December 2018


clear all
set more off


*** Define root name for data folder
	local directory_path 	`"/Users/peizan/github/riattrition/data-raw/chen_and_yang_2019/raw_data"'
	cd "`directory_path'"

	
	
*** Program Anderson z-score generation

	cap program drop andersonz 
	program andersonz
	syntax varlist [if] [in], GENerate(name)

	* 	Convert to effect sizes
		local varliste ""
		foreach X of varlist `varlist' {
			qui sum `X' `if'
			gen `X'e = `X'/r(sd) `if'
			local csd = r(sd)
			qui sum `X'
			replace `X'e = `X'e - r(mean)/`csd'
			local varliste "`varliste' `X'e" 
			}
	
	*   Generate GLS weighting matrix
		matrix accum R = `varliste', nocons dev
		matrix R = R/r(N)
		matrix R = syminv(R)
		local counter1 = 1
		matrix J = J(colsof(R),1,1)
		while `counter1' <= colsof(R) {
			matrix T = R[`counter1',1..colsof(R)]
			matrix A = T*J
			local weight`counter1'f = A[1,1]
			local counter1 = `counter1' + 1
			}
	
	*   Generate total number of variables per obs
	*   As opposed to Anderson's code, do not replace missing with zeros
		cap drop sample
		gen sample = 0
		local counter1 = 1
		foreach X of varlist `varlist' {
			replace sample = sample + `weight`counter1'f' if `X'~=.
	*		replace `X'e = 0 if `X' ==.
			local counter1 = `counter1' + 1
			}
	
	*   Apply weights to outcomes
		quietly gen `generate' = 0 `if'
		
		local counter1 = 1
		foreach X of varlist `varlist' {
			replace `generate' = `X'e * (`weight`counter1'f') + `generate' `if'
			local counter1 = `counter1' + 1
			}
			
		replace `generate' = `generate' / sample `if'
		cap drop sample
		cap drop `varliste'
		
	* 	Re-standardize aggregate z-score
		egen z_`generate' = std(`generate')
		drop `generate'
		rename z_`generate' `generate'
	
	end	
	
	

*** Load merged raw data for student survey

	use panelsurvey_raw, clear
	
	

	
*** Define variables

	local var_identity_w1 					`"responseID_wave1"'
	local var_university 					`"university"'
	local var_treatment 					`"treatment_vpn treatment_newsletter vpn_current_paid_user treatment_vpnexpiration"'
	
	
* 	A: Media related beliefs, attitudes, and behaviors
	local var_info_ranking_reg_w1 			`"info_domestic_website_w1 info_foreign_website_w1 info_social_media_dom_w1 info_word_of_mouth_w1"'
	local var_info_ranking_reg_w2 			`"info_domestic_website_w2 info_foreign_website_w2 info_social_media_dom_w2 info_social_media_for_w2 info_word_of_mouth_w2"'
	local var_info_ranking_reg_w3 			`"info_domestic_website_w3 info_foreign_website_w3 info_social_media_dom_w3 info_social_media_for_w3 info_word_of_mouth_w3"'
	local var_info_ranking_fig_w23 			`"info_domestic_website info_foreign_website info_social_media_dom info_word_of_mouth info_social_media_for"'
	
	local var_info_freq_reg_w1 				`"info_freq_website_for_w1"'
	local var_info_freq_reg_w2 				`"info_freq_website_for_w2"'
	local var_info_freq_reg_w3 				`"info_freq_website_for_w3"'
	local var_info_freq_fig_w123 			`"info_freq_website_for"'
	
	local var_vpnusage_own 					`"vpn_usage vpn_neveruser vpn_payment_now vpn_stop_when vpn_stop_why vpn_stop_why_text"'
	
	local var_vpnusage_rm_all 				`"vpn_roommate_usage"'
	
	local var_vpn_purchase_w3 				`"vpn_purchase vpn_purchase_yes vpn_purchase_wmt vpn_purchase_premium"'
	local var_vpn_purchase_reg_w3 			`"vpn_purchase_wmt_record vpn_purchase_yes"'
		
	local var_media_valuation_reg_w1 		`"wtp_vpn_w1 added_value_foreign_media_w1"'
	local var_media_valuation_reg_w2 		`"wtp_vpn_w2 added_value_foreign_media_w2"'
	local var_media_valuation_reg_w3 		`"wtp_vpn_w3 added_value_foreign_media_w3"'
	local var_media_valuation_fig_w123 		`"wtp_vpn added_value_foreign_media az_belief_media_value"'
	
	local var_media_trust_reg_w1 			`"trust_media_dom_state_w1 trust_media_dom_private_w1 trust_media_foreign_w1"'
	local var_media_trust_reg_w2 			`"trust_media_dom_state_w2 trust_media_dom_private_w2 trust_media_foreign_w2"'
	local var_media_trust_reg_w3 			`"trust_media_dom_state_w3 trust_media_dom_private_w3 trust_media_foreign_w3"'
	local var_media_trust_fig_w123 			`"trust_media_dom_state trust_media_dom_private trust_media_foreign az_belief_media_trust"'
	
	local var_percmediabias_reg_di_cn_w2 	`"distneutral_cn_neg_cn_w2 distneutral_cn_pos_cn_w2 distneutral_us_neg_cn_w2 distneutral_us_pos_cn_w2"'
	local var_percmediabias_reg_di_us_w2 	`"distneutral_cn_neg_us_w2 distneutral_cn_pos_us_w2 distneutral_us_neg_us_w2 distneutral_us_pos_us_w2"'
	local var_percmediabias_reg_ce_cn_w2 	`"bias_cn_neg_cn_cens_w2 bias_cn_pos_cn_cens_w2 bias_us_neg_cn_cens_w2 bias_us_pos_cn_cens_w2"'
	local var_percmediabias_reg_ce_us_w2 	`"bias_cn_neg_us_cens_w2 bias_cn_pos_us_cens_w2 bias_us_neg_us_cens_w2 bias_us_pos_us_cens_w2"'
	local var_percmediabias_fig_w12 		`"distneutral_cn_neg_cn distneutral_cn_pos_cn distneutral_us_neg_cn distneutral_us_pos_cn distneutral_cn_neg_us distneutral_cn_pos_us distneutral_us_neg_us distneutral_us_pos_us bias_cn_neg_cn_cens bias_cn_pos_cn_cens bias_us_neg_cn_cens bias_us_pos_cn_cens bias_cn_neg_us_cens bias_cn_pos_us_cens bias_us_neg_us_cens bias_us_pos_us_cens"'
	
	local var_censor_justif_reg_w1 			`"censor_just_dom_economic_w1 censor_just_dom_political_w1 censor_just_dom_social_w1 censor_just_for_w1"'
	local var_censor_justif_reg_w2 			`"censor_just_dom_economic_w2 censor_just_dom_political_w2 censor_just_dom_social_w2 censor_just_for_w2 censor_just_porn_w2"'
	local var_censor_justif_reg_w3 			`"censor_just_dom_economic_w3 censor_just_dom_political_w3 censor_just_dom_social_w3 censor_just_for_w3 censor_just_porn_w3"'
	local var_censor_justif_reg_fp_w1 		`"censor_just_dom_economic_w1 censor_just_dom_political_w1 censor_just_dom_social_w1 censor_just_for_w1"'
	local var_censor_justif_reg_fp_w2 		`"censor_just_dom_economic_w2 censor_just_dom_political_w2 censor_just_dom_social_w2 censor_just_for_w2"'
	local var_censor_justif_reg_fp_w3 		`"censor_just_dom_economic_w3 censor_just_dom_political_w3 censor_just_dom_social_w3 censor_just_for_w3"'
	local var_censor_justif_fig_w123 		`"censor_just_dom_economic censor_just_dom_political censor_just_dom_social censor_just_for az_belief_media_justif_pf"'
	local var_censor_justif_fig_w23 		`"censor_just_porn"'
	
	local var_censor_level_reg_w1			`"bias_domestic_w1 bias_foreign_w1"'
	local var_censor_level_reg_w2			`"bias_domestic_w2 bias_foreign_w2"'
	local var_censor_level_reg_w3			`"bias_domestic_w3 bias_foreign_w3"'
	local var_censor_driver_reg_dom_w1 		`"bias_dom_govt_policy_t1_w1 bias_dom_firm_interest_t1_w1 bias_dom_media_pref_t1_w1 bias_dom_reader_demand_t1_w1"'
	local var_censor_driver_reg_dom_w2 		`"bias_dom_govt_policy_t1_w2 bias_dom_firm_interest_t1_w2 bias_dom_media_pref_t1_w2 bias_dom_reader_demand_t1_w2"'
	local var_censor_driver_reg_dom_w3 		`"bias_dom_govt_policy_t1_w3 bias_dom_firm_interest_t1_w3 bias_dom_media_pref_t1_w3 bias_dom_reader_demand_t1_w3"'
	local var_censor_driver_reg_for_w1		`"bias_for_govt_policy_t1_w1 bias_for_firm_interest_t1_w1 bias_for_media_pref_t1_w1 bias_for_reader_demand_t1_w1"'
	local var_censor_driver_reg_for_w2		`"bias_for_govt_policy_t1_w2 bias_for_firm_interest_t1_w2 bias_for_media_pref_t1_w2 bias_for_reader_demand_t1_w2"'
	local var_censor_driver_reg_for_w3		`"bias_for_govt_policy_t1_w3 bias_for_firm_interest_t1_w3 bias_for_media_pref_t1_w3 bias_for_reader_demand_t1_w3"'
	local var_censor_level_fig_w123 		`"bias_domestic bias_foreign"'
	local var_censor_driver_fig_w123 		`"bias_dom_govt_policy_t1 bias_dom_firm_interest_t1 bias_dom_media_pref_t1 bias_dom_reader_demand_t1 bias_for_govt_policy_t1 bias_for_firm_interest_t1 bias_for_media_pref_t1 bias_for_reader_demand_t1"'
	
	
* 	B: Knowledge
	local var_knowledge_news_reg_w1 		`"news_c_stock_rise news_c_exchange_depreciate news_c_train_brazil_peru news_c_army_reduce news_c_taiwan_election news_c_china_us_network_coop news_c_lijiacheng_china news_c_nanjing_memory_list"'
	local var_knowledge_news_reg_w2 		`"news_c_panamapapers news_c_stockcrash news_c_topinequality news_c_cpicensorship news_c_laborunrest news_c_tenyearshk news_c_economistcensor news_c_waterpollution news_c_applefbi news_c_taiwanelection news_c_yihehotel news_perccorrect_w2 news_perccor_cen_w2 news_perccor_unc_w2 news_perccor_qui"'
	local var_knowledge_news_reg_quiz 		`"news_c_topinequality news_c_cpicensorship news_c_laborunrest news_c_waterpollution"'
	local var_knowledge_news_reg_cen_w1 	`"news_c_stock_rise news_c_train_brazil_peru news_c_taiwan_election news_c_china_us_network_coop"'
	local var_knowledge_news_reg_cen_w2 	`"news_c_panamapapers news_c_tenyearshk news_c_stockcrash news_c_economistcensor"'
	local var_knowledge_news_reg_cen_w3 	`"news_c_coalprod news_c_trumpchina news_c_xiaojianhua news_c_xijiangcar news_c_chinanorway news_c_womenrights news_c_hkceelection"'
	local var_knowledge_news_reg_unc_w1 	`"news_c_exchange_depreciate news_c_army_reduce news_c_lijiacheng_china news_c_nanjing_memory_list"'
	local var_knowledge_news_reg_unc_w2 	`"news_c_applefbi news_c_taiwanelection news_c_yihehotel"'
	local var_knowledge_news_reg_unc_w3 	`"news_c_northkoreacoal news_c_birdflu news_c_ethiopiatrain news_c_foreignreserve"'
	local var_knowledge_news_fig_w123 		`"news_perccor_cen news_perccor_unc"'
	
	local var_knowledge_people_reg_tocens 	`"people_puzhiqiang_w2 people_renzhiqiang_w2 people_huangzhifeng_w2"'
	local var_knowledge_people_reg_censor 	`"people_lizehou_w2 people_chenguangcheng_w2 people_lixiaolin_w2"'
	local var_knowledge_people_reg_uncens 	`"people_maoyushi_w2 people_honghuang_w2 people_liuqiangdong_w2"'
	local var_knowledge_people_reg_fake 	`"people_jialequn_w2"'
	local var_knowledge_people_censor_w1 	`"people_lizehou_w1 people_huangzhifeng_w1 people_chenguangcheng_w1 people_lixiaolin_w1"'
	local var_knowledge_people_uncens_w1 	`"people_puzhiqiang_w1 people_renzhiqiang_w1 people_maoyushi_w1 people_honghuang_w1 people_liuqiangdong_w1"'
	local var_knowledge_people_fig_w12 		`"people_puzhiqiang people_lizehou people_huangzhifeng people_chenguangcheng people_lixiaolin people_renzhiqiang people_maoyushi people_honghuang people_liuqiangdong people_jialequn"'
	
	local var_knowledge_prot_reg_w2 		`"protest_2012_hk_curriculum_w2 protest_2014_umbrella_w2 protest_2016_mongkok_riot_w2 protest_2014_sun_flower_w2 protest_2014_europe_square_w2 protest_2010_arabic_spring_w2 protest_2014_crimea_vote_w2 protest_2010_catal_indep_w2 protest_2011_tmrw_parade_w2 protest_pcheard_total_w2 protest_pcheard_china_w2 protest_pcheard_foreign_w2"'
	local var_knowledge_prot_reg_w3 		`"protest_2012_hk_curriculum_w3 protest_2014_umbrella_w3 protest_2016_mongkok_riot_w3 protest_2014_sun_flower_w3 protest_2014_europe_square_w3 protest_2010_arabic_spring_w3 protest_2014_crimea_vote_w3 protest_2010_catal_indep_w3 protest_2011_tmrw_parade_w3 protest_2017_women_march_w3 protest_pcheard_total_w3 protest_pcheard_china_w3 protest_pcheard_foreign_w3"'
	local var_knowledge_prot_reg_chi_w1 	`"protest_2012_hk_curriculum_w1 protest_2014_umbrella_w1 protest_2014_sun_flower_w1"'
	local var_knowledge_prot_reg_chi_w2  	`"protest_2012_hk_curriculum_w2 protest_2014_umbrella_w2 protest_2016_mongkok_riot_w2 protest_2014_sun_flower_w2"'
	local var_knowledge_prot_reg_chi_w3  	`"protest_2012_hk_curriculum_w3 protest_2014_umbrella_w3 protest_2016_mongkok_riot_w3 protest_2014_sun_flower_w3"'
	local var_knowledge_prot_reg_for_w1 	`"protest_2014_europe_square_w1 protest_2010_arabic_spring_w1 protest_2014_crimea_vote_w1 protest_2010_catal_indep_w1"'
	local var_knowledge_prot_reg_for_w2 	`"protest_2014_europe_square_w2 protest_2010_arabic_spring_w2 protest_2014_crimea_vote_w2 protest_2010_catal_indep_w2"'
	local var_knowledge_prot_reg_for_w3 	`"protest_2014_europe_square_w3 protest_2010_arabic_spring_w3 protest_2014_crimea_vote_w3 protest_2010_catal_indep_w3 protest_2017_women_march_w3"'
	local var_knowledge_prot_reg_fak_w1 	`"protest_2011_tmrw_parade_w1"'
	local var_knowledge_prot_reg_fak_w2 	`"protest_2011_tmrw_parade_w2"'
	local var_knowledge_prot_reg_fak_w3 	`"protest_2011_tmrw_parade_w3"'
	local var_knowledge_prot_fig_w123 		`"protest_2014_europe_square protest_2014_sun_flower protest_2010_arabic_spring protest_2014_crimea_vote protest_2012_hk_curriculum protest_2010_catal_indep protest_2014_umbrella protest_2011_tmrw_parade protest_pcheard_china protest_pcheard_foreign"' 
	local var_knowledge_prot_fig_w23  		`"protest_2016_mongkok_riot"'
	local var_knowledge_prot_fig_w3  		`"protest_2017_women_march"'
	
	local var_knowledge_time_w2 			`"news_time_firstclick news_time_lastclick news_time_submit news_time_totalclick pppr_time_firstclick pppr_time_lastclick pppr_time_submit pppr_time_totalclick"'
	
	local var_knowledge_past_w1 			`"event_1994_beijing_jianguomen event_1994_kelamayi_fire event_2003_sunzhigang event_2008_guizhou_wengan event_2008_sanlu_milk_powder event_2010_hefei_wetland event_2011_guangdong_wukan event_2013_fudan_poison event_2014_caixin_xuxiao event_2015_sichuan_linshui_road"'
	
	local var_knowledge_meta_reg_w1 		`"familiar_china_issues_self_w1 familiar_china_others_w1"'
	local var_knowledge_meta_reg_w2 		`"familiar_china_issues_self_w2 familiar_china_others_w2"'
	local var_knowledge_meta_reg_w3 		`"familiar_china_issues_self_w3 familiar_china_others_w3"'
	local var_knowledge_meta_fig_w123 		`"familiar_china_issues_self familiar_china_others az_knowledge_meta"'

	
* 	C: Economic beliefs
	local var_econ_guess_reg_cn_perf_w1 	`"guess_gdp_growth_china_w1 guess_stock_index_sh_w1"'
	local var_econ_guess_reg_cn_perf_w2 	`"guess_gdp_growth_china_w2 guess_stock_index_sh_w2"'
	local var_econ_guess_reg_cn_perf_w3 	`"guess_gdp_growth_china_w3 guess_stock_index_sh_w3"'
	local var_econ_guess_reg_cn_conf_w1 	`"guess_gdp_growth_china_con_w1 guess_stock_index_sh_con_w1"'
	local var_econ_guess_reg_cn_conf_w2 	`"guess_gdp_growth_china_con_w2 guess_stock_index_sh_con_w2"'
	local var_econ_guess_reg_cn_conf_w3 	`"guess_gdp_growth_china_con_w3 guess_stock_index_sh_con_w3"'
	local var_econ_guess_reg_us_perf_w2 	`"guess_gdp_growth_us_w2 guess_stock_index_dj_w2"'
	local var_econ_guess_reg_us_perf_w3 	`"guess_gdp_growth_us_w3 guess_stock_index_dj_w3"'
	local var_econ_guess_reg_us_conf_w2 	`"guess_gdp_growth_us_con_w2 guess_stock_index_dj_con_w2"'
	local var_econ_guess_reg_us_conf_w3 	`"guess_gdp_growth_us_con_w3 guess_stock_index_dj_con_w3"'
	local var_econ_guess_fig_w123 			`"guess_gdp_growth_china guess_gdp_growth_china_con guess_stock_index_sh guess_stock_index_sh_con az_belief_econ_perf_cn az_belief_econ_conf_cn"'
	local var_econ_guess_fig_w23 			`"guess_gdp_growth_us guess_gdp_growth_us_con guess_stock_index_dj guess_stock_index_dj_con az_belief_econ_perf_us az_belief_econ_conf_us"'
	

* 	D: Political attitudes
	local var_demand_change_reg_w1 			`"inst_change_econ_w1 inst_change_poli_w1"'
	local var_demand_change_reg_w2 			`"inst_change_econ_w2 inst_change_poli_w2"'
	local var_demand_change_reg_w3 			`"inst_change_econ_w3 inst_change_poli_w3"'
	local var_demand_change_fig_w123 		`"inst_change_econ inst_change_poli az_belief_instchange"'
	
	local var_trust_inst_reg_govt_w1 		`"trust_central_govt_w1 trust_provincial_govt_w1 trust_local_govt_w1"'
	local var_trust_inst_reg_govt_w2 		`"trust_central_govt_w2 trust_provincial_govt_w2 trust_local_govt_w2"'
	local var_trust_inst_reg_govt_w3 		`"trust_central_govt_w3 trust_provincial_govt_w3 trust_local_govt_w3"'
	local var_trust_inst_reg_foreign_w1 	`"trust_japan_govt_w1 trust_us_govt_w1"'
	local var_trust_inst_reg_foreign_w2 	`"trust_japan_govt_w2 trust_us_govt_w2"'
	local var_trust_inst_reg_foreign_w3 	`"trust_japan_govt_w3 trust_us_govt_w3"'
	local var_trust_inst_reg_finance_w2 	`"trust_financial_domestic_w2 trust_financial_foreign_w2"'
	local var_trust_inst_reg_ngo_w2 		`"trust_ngo_w2"'
	local var_trust_inst_reg_copo_w2 		`"trust_court_w2 trust_police_w2"'
	local var_trust_inst_fig_w123  			`"trust_central_govt trust_provincial_govt trust_local_govt trust_japan_govt trust_us_govt az_belief_trust_govt az_belief_trust_foreign"'
	local var_trust_inst_fig_w12 			`"trust_court trust_police trust_ngo trust_financial_domestic trust_financial_foreign"'
	
	local var_eval_govt_reg_w1 				`"eval_govt_economic_w1 eval_govt_dom_politics_w1 eval_govt_for_relations_w1"'
	local var_eval_govt_reg_w2 				`"eval_govt_economic_w2 eval_govt_dom_politics_w2 eval_govt_for_relations_w2"'
	local var_eval_govt_reg_w3 				`"eval_govt_economic_w3 eval_govt_dom_politics_w3 eval_govt_for_relations_w3"'
	local var_eval_govt_fig_w123 			`"eval_govt_economic eval_govt_dom_politics eval_govt_for_relations az_belief_evalgovt"'
	
	local var_eval_criteria_reg_w2 			`"revalgovt_election_w2 revalgovt_economy_w2 revalgovt_equality_w2 revalgovt_ruleoflaw_w2 revalgovt_human_rights_w2 revalgovt_freedom_speech_w2 revalgovt_global_power_w2 revalgovt_fair_history_w2"'
	local var_eval_criteria_fig_w12 		`"revalgovt_election revalgovt_economy revalgovt_equality revalgovt_ruleoflaw revalgovt_human_rights revalgovt_freedom_speech revalgovt_global_power revalgovt_fair_history"'
	
	local var_severity_reg_w2 				`"severity_welfare_w2 severity_employment_w2 severity_pollution_w2 severity_inequality_w2 severity_corruption_w2 severity_dscrm_minority_w2"'
	local var_severity_fig_w12 				`"severity_welfare severity_employment severity_pollution severity_inequality severity_corruption severity_dscrm_minority"'
	
	local var_democracy_reg_w2 				`"china_interest_group_w2 china_rate_democracy_w2 china_rate_humanrights_w2 importance_live_in_demo_w2"'
	local var_democracy_reg_w3 				`"importance_live_in_demo_w3"'
	local var_democracy_reg_fp_w1 			`"importance_live_in_demo_w1"'
	local var_democracy_reg_fp_w2 			`"importance_live_in_demo_w2"'
	local var_democracy_reg_fp_w3 			`"importance_live_in_demo_w3"'
	local var_democracy_fig_w123 			`"importance_live_in_demo az_belief_democracy_fp"'
	local var_democracy_fig_w12 			`"china_interest_group china_rate_democracy china_rate_humanrights"'
	
	local var_contro_justi_reg_policy_w2 	`"justify_minority_policy_w2 justify_hukou_w2 justify_one_child_w2 justify_hongkong_policy_w2 justify_taiwan_policy_w2 justify_violence_stability_w2 justify_receive_refugee_w2 justify_reduce_pollution_w2 justify_gaokao_w2 justify_soe_privatize_w2"'
	local var_contro_justi_reg_liberal_w2 	`"justify_homo_marriage_w2 justify_legal_prostitute_w2 justify_abortion_w2 justify_exmaritalsex justify_transgene_w2 justify_soft_drug_w2"'
	local var_contro_justi_fig_w12 			`"justify_minority_policy justify_reduce_pollution justify_hukou justify_one_child justify_gaokao justify_hongkong_policy justify_taiwan_policy justify_transgene justify_receive_refugee justify_soe_privatize justify_homo_marriage justify_legal_prostitute justify_abortion justify_soft_drug justify_violence_stability"'
	
	local var_willing_fight_reg_w1 			`"willing_against_illi_govt_w1 willing_report_mis_w1 willing_protect_weak_w1"'
	local var_willing_fight_reg_w2 			`"willing_against_illi_govt_w2 willing_report_mis_w2 willing_protect_weak_w2"'
	local var_willing_fight_reg_w3 			`"willing_against_illi_govt_w3 willing_report_mis_w3 willing_protect_weak_w3"'
	local var_willing_fight_fig_123 		`"willing_against_illi_govt willing_report_mis willing_protect_weak az_belief_willing"'
	
	local var_interest_reg_w2 				`"interest_economic_w2 interest_politics_w2"'
	local var_interest_fig_w12 				`"interest_economic interest_politics"'
	
	local var_patriotism_reg_w2 			`"proud_being_chinese_w2"'
	local var_patriotism_fig_w12 			`"proud_being_chinese"'
	
	local var_fear_critgovt_reg_w2 			`"fear_critic_govt_self_w2"'
	local var_fear_critgovt_fig_w12 		`"fear_critic_govt_self"'
	
	
* 	E: Behaviors
	local var_socialinteract_reg_w1 		`"frequency_talk_politic_w1 frequency_persuade_friends_w1"'
	local var_socialinteract_reg_w2 		`"frequency_talk_politic_w2 frequency_persuade_friends_w2"'
	local var_socialinteract_reg_w3 		`"frequency_talk_politic_w3 frequency_persuade_friends_w3"'
	local var_socialinteract_fig_w123 		`"frequency_talk_politic frequency_persuade_friends az_var_socialinteract"'
	
	local var_socialknowledge_w1 			`"know_attitudes_relatives know_attitudes_schoolmates know_attitudes_outschool_friends"'
	
	local var_polparticipation_reg_w2 		`"participate_ngo_w2 participate_social_protest_w2 participate_plan_vote_w2 participate_complain_school_w2"'
	local var_polparticipation_reg_w3 		`"participate_social_protest_w3 participate_plan_vote_w3 participate_complain_school_w3"'
	local var_polparticipation_reg_pf_w1 	`"participate_social_protest_w1 participate_plan_vote_w1 participate_complain_school_w1"'
	local var_polparticipation_reg_pf_w2 	`"participate_social_protest_w2 participate_plan_vote_w2 participate_complain_school_w2"'
	local var_polparticipation_reg_pf_w3 	`"participate_social_protest_w3 participate_plan_vote_w3 participate_complain_school_w3"'
	local var_polparticipation_fig_w123 	`"participate_social_protest participate_plan_vote participate_complain_school az_var_polparticipation"'
	local var_polparticipation_fig_w12 		`"participate_ngo"'
	
	local var_polengagement_w1  			`"ccp_member ccp_join_year ccp_prep_member participate_stuunion participate_stuunion_ever participate_tuanwei participate_tuanwei_ever"'

	local var_planaftergrad_reg_w1 			`"plan_grad_gradschool_dom_w1 plan_grad_foreignmaster_w1 plan_grad_foreignphd_w1 plan_grad_military_w1 plan_grad_work_w1"'
	local var_planaftergrad_reg_w2 			`"plan_grad_gradschool_dom_w2 plan_grad_foreignmaster_w2 plan_grad_foreignphd_w2 plan_grad_military_w2 plan_grad_work_w2"'
	local var_planaftergrad_reg_w3 			`"plan_grad_gradschool_dom_w3 plan_grad_foreignmaster_w3 plan_grad_foreignphd_w3 plan_grad_military_w3 plan_grad_work_w3"'
	local var_planaftergrad_fig_w123 		`"plan_grad_gradschool_dom plan_grad_foreignmaster plan_grad_foreignphd plan_grad_military plan_grad_work"'
	
	local var_career_sector_reg_w1 			`"cp_t3_national_civil_w1 cp_t3_local_civil_w1 cp_t3_military_w1 cp_t3_chinese_private_w1 cp_t3_for_firm_w1 cp_t3_soe_w1 cp_t3_institutional_w1 cp_t3_entrepreneur_w1"'
	local var_career_sector_reg_w2 			`"cp_t3_national_civil_w2 cp_t3_local_civil_w2 cp_t3_military_w2 cp_t3_chinese_private_w2 cp_t3_for_firm_w2 cp_t3_soe_w2 cp_t3_institutional_w2 cp_t3_entrepreneur_w2"'
	local var_career_sector_reg_w3 			`"cp_t3_national_civil_w3 cp_t3_local_civil_w3 cp_t3_military_w3 cp_t3_chinese_private_w3 cp_t3_for_firm_w3 cp_t3_soe_w3 cp_t3_institutional_w3 cp_t3_entrepreneur_w3"'
	local var_career_sector_fig_w123 		`"cp_t1_national_civil cp_t1_local_civil cp_t1_military cp_t1_chinese_private cp_t1_for_firm cp_t1_soe cp_t1_institutional cp_t1_entrepreneur"'
	
	local var_career_loc_reg_w1 			`"cloc_beijing_w1 cloc_shanghai_w1 cloc_gzsz_w1 cloc_tjcq_w1 cloc_hkmc_w1 cloc_taiwan_w1 cloc_dom_w1 cloc_for_w1"'
	local var_career_loc_reg_w2 			`"cloc_beijing_w2 cloc_shanghai_w2 cloc_gzsz_w2 cloc_tjcq_w2 cloc_hkmc_w2 cloc_taiwan_w2 cloc_dom_w2 cloc_for_w2"'
	local var_career_loc_reg_w3 			`"cloc_beijing_w3 cloc_shanghai_w3 cloc_gzsz_w3 cloc_tjcq_w3 cloc_hkmc_w3 cloc_taiwan_w3 cloc_dom_w2 cloc_for_w3"'
	local var_career_loc_reg_w123 			`"cloc_beijing cloc_shanghai cloc_gzsz cloc_tjcq cloc_hkmc cloc_taiwan cloc_dom cloc_for"'
	
	local var_stock_invest_reg_w1  			`"stock_participation_w1"'
	local var_stock_invest_reg_w2  			`"stock_participation_w2"'
	local var_stock_invest_reg_w3  			`"stock_participation_w3"'
	local var_stock_invest_fig_w123 		`"stock_participation"'
	
	
* 	F: Demographics, background characteristics, and fundamental preferences
		
	// demographics
	local var_demog_reg 					`"gender birth_year ethnicity_han birthplace_coastal residence_coastal hukou_urban religion_religious ccp_member university_elite hs_track_science department_ssh siblings_total father_edu_hsabove father_ccp mother_edu_hsabove mother_ccp hh_income domestic_english_atleast4 foreign_english_yes travel_hktaiwan travel_foreign_yes"'
	local var_demog_reg_personal 			`"gender birth_year ethnicity_han birthplace_coastal residence_coastal hukou_urban religion_religious ccp_member"'
	local var_demog_reg_education 			`"university_elite hs_track_science department_ssh"'
	local var_demog_reg_english 			`"domestic_english_atleast4 foreign_english_yes"'
	local var_demog_reg_travel 				`"travel_hktaiwan travel_foreign_yes"'
	local var_demog_reg_household 			`"siblings_total father_edu_hsabove work_father_govt father_ccp mother_edu_hsabove work_mother_govt mother_ccp hh_income"'
	
	local var_demog_reg_imbalance 			`"residence_coastal hs_track_science father_ccp mother_ccp gift_amount risk_preference_ce"'
	
	// fundamental preferences
	local var_preference_risk 				`"willing_risk risk_preference_ce risk_lottery_choice"'
	local var_preference_time 				`"willing_future procrastinate"'
	local var_preference_altruism 			`"willing_goodcauses donate_amount"'
	local var_preference_reciprocity 		`"willing_returnfavor best_intentions gift_amount willing_punish_you willing_punish_others willing_revenge"'

		
* 	G: List experiment
	local var_listexp_all 					`"list_exp_direct list_trust_direct_count list_trust_direct_yes list_trust_veiled_count"'
	local var_listexp_reg 					`"list_count_trust"'
	
	
* 	Y: Overall effects
	local var_overall_w1 					`"info_freq_website_for_w1 az_belief_media_value_w1 az_belief_media_trust_w1 bias_domestic_w1 az_belief_media_justif_w1 bias_foreign_r_w1 news_perccor_cen_w1 news_perccor_unc_w1 protest_pcheard_china_w1 protest_pcheard_foreign_w1 az_knowledge_meta_w1 az_belief_econ_perf_cn_r_w1 az_belief_econ_conf_cn_w1 az_belief_instchange_w1 az_belief_trust_govt_r_w1 az_belief_trust_foreign_w1 az_belief_evalgovt_r_w1 importance_live_in_demo_r_w1 az_belief_willing_w1 az_var_socialinteract_w1 az_var_polparticipation_w1 plan_grad_gradschool_dom_w1 plan_grad_foreignmaster_w1 plan_grad_foreignphd_w1 plan_grad_military_w1 plan_grad_work_w1 cloc_beijing_w1 cloc_shanghai_w1 cloc_gzsz_w1 cloc_tjcq_w1 cloc_hkmc_w1 cloc_taiwan_w1 cloc_dom_w2 cloc_for_w1 stock_participation_w1"'
	local var_overall_w2 					`"info_freq_website_for_w2 az_belief_media_value_w2 az_belief_media_trust_w2 bias_domestic_w2 az_belief_media_justif_w2 bias_foreign_r_w2 news_perccor_cen_w2 news_perccor_unc_w2 protest_pcheard_china_w2 protest_pcheard_foreign_w2 az_knowledge_meta_w2 az_belief_econ_perf_cn_r_w2 az_belief_econ_conf_cn_w2 az_belief_instchange_w2 az_belief_trust_govt_r_w2 az_belief_trust_foreign_w2 az_belief_evalgovt_r_w2 importance_live_in_demo_r_w2 az_belief_willing_w2 az_var_socialinteract_w2 az_var_polparticipation_w2 plan_grad_gradschool_dom_w2 plan_grad_foreignmaster_w2 plan_grad_foreignphd_w2 plan_grad_military_w2 plan_grad_work_w2 cloc_beijing_w2 cloc_shanghai_w2 cloc_gzsz_w2 cloc_tjcq_w2 cloc_hkmc_w2 cloc_taiwan_w2 cloc_dom_w2 cloc_for_w2 stock_participation_w2"'
	local var_overall_w3 					`"info_freq_website_for_w3 az_belief_media_value_w3 az_belief_media_trust_w3 bias_domestic_w3 az_belief_media_justif_w3 bias_foreign_r_w3 news_perccor_cen_w3 news_perccor_unc_w3 protest_pcheard_china_w3 protest_pcheard_foreign_w3 az_knowledge_meta_w3 az_belief_econ_perf_cn_r_w3 az_belief_econ_conf_cn_w3 az_belief_instchange_w3 az_belief_trust_govt_r_w3 az_belief_trust_foreign_w3 az_belief_evalgovt_r_w3 importance_live_in_demo_r_w3 az_belief_willing_w3 az_var_socialinteract_w3 az_var_polparticipation_w3 plan_grad_gradschool_dom_w3 plan_grad_foreignmaster_w3 plan_grad_foreignphd_w3 plan_grad_military_w3 plan_grad_work_w3 cloc_beijing_w3 cloc_shanghai_w3 cloc_gzsz_w3 cloc_tjcq_w3 cloc_hkmc_w3 cloc_taiwan_w3 cloc_dom_w3 cloc_for_w3 stock_participation_w3"'

	
	
	
*** Prepare data: merge additional variables

* 	merge in vpn usage data: adoption dates
	merge 1:1 responseID_wave1 using vpn_date_adoption
	drop _merge
		

* 	merge in indicator of vpn active user
	merge 1:1 responseID_wave1 using vpn_browsing_active_user
	drop _merge
	replace active_user = 0 if active_user == .
	
	
* 	merge in vpn purchase data (from service provider server)
	merge 1:1 responseID_wave1 using vpn_purchase
	drop _merge
	
		


*** Prepare data: generate additional variables
	
* 	wave 1

	// foreach news, indicator correct or not
	gen news_c_stock_rise 				= (news_stock_rise == 1) if news_stock_rise != .
	gen news_c_exchange_depreciate 		= (news_exchange_depreciate == 1) if news_exchange_depreciate != .
	gen news_c_train_brazil_peru  		= (news_train_brazil_peru == 0) if news_train_brazil_peru != .
	gen news_c_army_reduce 				= (news_army_reduce == 1) if news_army_reduce != .
	gen news_c_taiwan_election 			= (news_taiwan_election == 0) if news_taiwan_election != .
	gen news_c_china_us_network_coop 	= (news_china_us_network_coop == 0) if news_china_us_network_coop != .
	gen news_c_lijiacheng_china 		= (news_lijiacheng_china == 0) if news_lijiacheng_china != .
	gen news_c_nanjing_memory_list 		= (news_nanjing_memory_list == 1) if news_nanjing_memory_list != .

	// generate count variable: # news answered correctly
	gen news_totalcorrect = news_c_stock_rise + news_c_exchange_depreciate + news_c_train_brazil_peru + news_c_army_reduce + news_c_taiwan_election + news_c_china_us_network_coop + news_c_lijiacheng_china + news_c_nanjing_memory_list
	
	// generate percentage correct
	gen news_perccorrect = news_totalcorrect / 8
	
	// correct count in percentage: censored news
	gen news_perccor_cen = (news_c_taiwan_election + news_c_lijiacheng_china) / 2
	
	// correct count in percentage: uncensored news
	gen news_perccor_unc = (news_c_stock_rise + news_c_exchange_depreciate + news_c_train_brazil_peru + news_c_army_reduce + news_c_china_us_network_coop + news_c_nanjing_memory_list) / 6

	// generate count variable: # people heard before
	gen people_totalheard = people_puzhiqiang_w1 + people_lizehou_w1 + people_huangzhifeng_w1 + people_chenguangcheng_w1 + people_lixiaolin_w1 + people_renzhiqiang_w1 + people_maoyushi_w1 + people_honghuang_w1 + people_liuqiangdong_w1
	gen people_pcheard = (people_puzhiqiang_w1 + people_lizehou_w1 + people_huangzhifeng_w1 + people_chenguangcheng_w1 + people_lixiaolin_w1 + people_renzhiqiang_w1 + people_maoyushi_w1 + people_honghuang_w1 + people_liuqiangdong_w1) / 9

	// generate count variable: # protests heard before
	gen protest_totalheard = protest_2014_europe_square_w1 + protest_2014_sun_flower_w1 + protest_2010_arabic_spring_w1 + protest_2014_crimea_vote_w1 + protest_2012_hk_curriculum_w1 + protest_2010_catal_indep_w1 + protest_2014_umbrella_w1
	gen protest_pcheard_total = (protest_2014_europe_square_w1 + protest_2014_sun_flower_w1 + protest_2010_arabic_spring_w1 + protest_2014_crimea_vote_w1 + protest_2012_hk_curriculum_w1 + protest_2010_catal_indep_w1 + protest_2014_umbrella_w1) / 9
	
	// generate percentage by category: full panel version
	gen protest_pcheard_china = (protest_2014_sun_flower_w1 + protest_2012_hk_curriculum_w1 + protest_2014_umbrella_w1) / 3
	gen protest_pcheard_foreign = (protest_2014_europe_square_w1 + protest_2010_arabic_spring_w1 + protest_2014_crimea_vote_w1 + protest_2010_catal_indep_w1) / 4

	// political participation: recode into ever participated
	gen participate_stuunion_ever = (participate_stuunion < 3) if participate_stuunion != .
	gen participate_tuanwei_ever = (participate_tuanwei < 3) if participate_tuanwei != .
	
	// generate distance from neutrality
	foreach var in us_pos_us us_pos_cn cn_pos_us cn_pos_cn us_neg_us us_neg_cn cn_neg_us cn_neg_cn {
		gen distneutral_`var' = abs(bias_`var'_media_w1 - 4) if bias_`var'_media_w1 != .
		}	
		
	// separate censorship and bias
	foreach var in us_pos_us us_pos_cn cn_pos_us cn_pos_cn us_neg_us us_neg_cn cn_neg_us cn_neg_cn {
		gen bias_`var'_cens = (bias_`var'_media_w1 == 1) if bias_`var'_media_w1 != .
		}
	
	// generate indicator for top categories: bias_for
	foreach var in bias_for_govt_policy bias_for_firm_interest bias_for_media_pref bias_for_reader_demand {
		gen `var'_t1 = (`var'_w1 == 1) if `var'_w1 != .
		}
		
	// generate indicator for top categories: bias_dom
	foreach var in bias_dom_govt_policy bias_dom_firm_interest bias_dom_media_pref bias_dom_reader_demand {
		gen `var'_t1 = (`var'_w1 == 1) if `var'_w1 != .
		}
		
	// generate relative weights for evaluating government performance
	gen evalgovt_total_w1 = evalgovt_election_w1 + evalgovt_economy_w1 + evalgovt_equality_w1 + evalgovt_ruleoflaw_w1 + evalgovt_human_rights_w1 + evalgovt_freedom_speech_w1 + evalgovt_global_power_w1 + evalgovt_fair_history_w1
	foreach var in evalgovt_election evalgovt_economy evalgovt_equality evalgovt_ruleoflaw evalgovt_human_rights evalgovt_freedom_speech evalgovt_global_power evalgovt_fair_history {
		gen r`var' = `var'_w1 / evalgovt_total_w1
		}

	// generate dummies for plan_graduation
	gen plan_grad_gradschool_dom 	= (plan_graduation_w1 == 1) if plan_graduation_w1 != .
	gen plan_grad_foreignmaster 	= (plan_graduation_w1 == 2) if plan_graduation_w1 != .
	gen plan_grad_foreignphd 		= (plan_graduation_w1 == 3) if plan_graduation_w1 != .
	gen plan_grad_military 			= (plan_graduation_w1 == 4) if plan_graduation_w1 != .
	gen plan_grad_work 				= (plan_graduation_w1 == 5) if plan_graduation_w1 != .
	gen plan_grad_dontknow 			= (plan_graduation_w1 == 6) if plan_graduation_w1 != .
	
	// generate top1 and top3 dummies for career preferences
	gen work_top3_missing = (work_top3_national_civil_w1 == . & work_top3_local_civil_w1 == . & work_top3_military_w1 == . & work_top3_chinese_private_w1 == . & work_top3_for_firm_w1 == . & work_top3_soe_w1 == . & work_top3_institutional_w1 == . & work_top3_entrepreneur_w1 == .)
	
	foreach var in national_civil local_civil military chinese_private for_firm soe institutional entrepreneur {
		gen cp_t1_`var' = (work_top3_`var'_w1 == 1) if work_top3_missing == 0
		gen cp_t3_`var' = (work_top3_`var'_w1 != .) if work_top3_missing == 0
		}
	drop work_top3_missing
	
	// generate dummies for career location preferences
	gen cloc_beijing 		= (place_top_w1 == 1) if place_top_w1 != .
	gen cloc_shanghai 		= (place_top_w1 == 2) if place_top_w1 != .
	gen cloc_gzsz 			= (place_top_w1 == 3 | place_top_w1 == 4) if place_top_w1 != .
	gen cloc_tjcq 			= (place_top_w1 == 5 | place_top_w1 == 6) if place_top_w1 != .
	gen cloc_hkmc 			= (place_top_w1 == 7 | place_top_w1 == 8) if place_top_w1 != .
	gen cloc_taiwan 		= (place_top_w1 == 9) if place_top_w1 != .
	gen cloc_dom 			= (place_top_w1 == 10) if place_top_w1 != .
	gen cloc_for  			= (place_top_w1 == 11) if place_top_w1 != .

	
	// add subscript: wave 1
	foreach var in news_totalcorrect news_perccorrect news_perccor_cen news_perccor_unc people_totalheard people_pcheard protest_totalheard protest_pcheard_total protest_pcheard_china protest_pcheard_foreign distneutral_us_pos_us distneutral_us_pos_cn distneutral_cn_pos_us distneutral_cn_pos_cn distneutral_us_neg_us distneutral_us_neg_cn distneutral_cn_neg_us distneutral_cn_neg_cn bias_us_pos_us_cens bias_us_pos_cn_cens bias_cn_pos_us_cens bias_cn_pos_cn_cens bias_us_neg_us_cens bias_us_neg_cn_cens bias_cn_neg_us_cens bias_cn_neg_cn_cens bias_for_govt_policy_t1 bias_for_firm_interest_t1 bias_for_media_pref_t1 bias_for_reader_demand_t1 bias_dom_govt_policy_t1 bias_dom_firm_interest_t1 bias_dom_media_pref_t1 bias_dom_reader_demand_t1 revalgovt_election revalgovt_economy revalgovt_equality revalgovt_ruleoflaw revalgovt_human_rights revalgovt_freedom_speech revalgovt_global_power revalgovt_fair_history plan_grad_gradschool_dom plan_grad_foreignmaster plan_grad_foreignphd plan_grad_military plan_grad_work plan_grad_dontknow cp_t1_national_civil cp_t3_national_civil cp_t1_local_civil cp_t3_local_civil cp_t1_military cp_t3_military cp_t1_chinese_private cp_t3_chinese_private cp_t1_for_firm cp_t3_for_firm cp_t1_soe cp_t3_soe cp_t1_institutional cp_t3_institutional cp_t1_entrepreneur cp_t3_entrepreneur cloc_beijing cloc_shanghai cloc_gzsz cloc_tjcq cloc_hkmc cloc_taiwan cloc_dom cloc_for {
		rename `var' 	`var'_w1
		}
	

* 	wave 2
	
	// foreach news, indicator correct or not
	gen news_c_topinequality 		= (news_topinequality == 0) if news_topinequality != .
	gen news_c_cpicensorship 	 	= (news_cpicensorship == 1) if news_cpicensorship != .
	gen news_c_laborunrest 			= (news_laborunrest == 1) if news_laborunrest != .
	gen news_c_waterpollution 		= (news_waterpollution == 0) if news_waterpollution != .
	gen news_c_stockcrash 			= (news_stockcrash == 1) if news_stockcrash != .
	gen news_c_taiwanelection 		= (news_taiwanelection == 0) if news_taiwanelection != .
	gen news_c_applefbi 			= (news_applefbi == 0) if news_applefbi != .
	gen news_c_tenyearshk 			= (news_tenyearshk == 1) if news_tenyearshk != .
	gen news_c_panamapapers 		= (news_panamapapers == 1) if news_panamapapers != .
	gen news_c_yihehotel 			= (news_yihehotel == 1) if news_yihehotel != .
	gen news_c_economistcensor 		= (news_economistcensor == 1) if news_economistcensor != .

	// generate count variable: # news answered correctly
	gen news_totalcorrect = news_c_topinequality + news_c_cpicensorship + news_c_laborunrest + news_c_waterpollution + news_c_stockcrash + news_c_taiwanelection + news_c_applefbi + news_c_tenyearshk + news_c_panamapapers + news_c_yihehotel + news_c_economistcensor

	// generate percentage correct
	gen news_perccorrect = news_totalcorrect / 11
	
	// correct count in percentage: quiz questions
	gen news_perccor_qui = (news_c_topinequality + news_c_cpicensorship + news_c_laborunrest + news_c_waterpollution) / 4
	
	// correct count in percentage: censored news
	gen news_perccor_cen = (news_c_stockcrash + news_c_panamapapers + news_c_tenyearshk + news_c_economistcensor) / 4
	
	// correct count in percentage: uncensored news
	gen news_perccor_unc = (news_c_taiwanelection + news_c_applefbi + news_c_yihehotel) / 3
		
	// generate count variable: # people heard before
	gen people_totalheard = people_puzhiqiang_w2 + people_lizehou_w2 + people_huangzhifeng_w2 + people_chenguangcheng_w2 + people_lixiaolin_w2 + people_renzhiqiang_w2 + people_maoyushi_w2 + people_honghuang_w2 + people_liuqiangdong_w2
	gen people_pcheard = (people_puzhiqiang_w2 + people_lizehou_w2 + people_huangzhifeng_w2 + people_chenguangcheng_w2 + people_lixiaolin_w2 + people_renzhiqiang_w2 + people_maoyushi_w2 + people_honghuang_w2 + people_liuqiangdong_w2) / 9
		
	// percentage of names heard by category
	gen people_perchd_censored = (people_lixiaolin_w2 + people_huangzhifeng_w2 + people_chenguangcheng_w2 + people_lizehou_w2) / 4
	gen people_perchd_uncensor = (people_maoyushi_w2 + people_liuqiangdong_w2 + people_honghuang_w2) / 3
	gen people_perchd_unctocen = (people_puzhiqiang_w2 + people_renzhiqiang_w2) / 2

	// generate count variable: # events heard before
	gen protest_totalheard = protest_2014_europe_square_w2 + protest_2014_sun_flower_w2 + protest_2010_arabic_spring_w2 + protest_2014_crimea_vote_w2 + protest_2012_hk_curriculum_w2 + protest_2010_catal_indep_w2 + protest_2014_umbrella_w2 + protest_2016_mongkok_riot_w2
	gen protest_pcheard_total = (protest_2014_europe_square_w2 + protest_2014_sun_flower_w2 + protest_2010_arabic_spring_w2 + protest_2014_crimea_vote_w2 + protest_2012_hk_curriculum_w2 + protest_2010_catal_indep_w2 + protest_2014_umbrella_w2 + protest_2016_mongkok_riot_w2) / 8
	
	// generate percentage by category: full panel version
	gen protest_pcheard_china = (protest_2014_sun_flower_w2 + protest_2012_hk_curriculum_w2 + protest_2014_umbrella_w2) / 3
	gen protest_pcheard_foreign = (protest_2014_europe_square_w2 + protest_2010_arabic_spring_w2 + protest_2014_crimea_vote_w2 + protest_2010_catal_indep_w2) / 4
		
	// generate distance from neutrality
	foreach var in us_pos_us us_pos_cn cn_pos_us cn_pos_cn us_neg_us us_neg_cn cn_neg_us cn_neg_cn {
		gen distneutral_`var' = abs(bias_`var'_media_w2 - 4) if bias_`var'_media_w2 != .
		}	
		
	// separate censorship and bias
	foreach var in us_pos_us us_pos_cn cn_pos_us cn_pos_cn us_neg_us us_neg_cn cn_neg_us cn_neg_cn {
		gen bias_`var'_cens = (bias_`var'_media_w2 == 1) if bias_`var'_media_w2 != .
		}
		
	// generate indicator for top categories: bias_for
	foreach var in bias_for_govt_policy bias_for_firm_interest bias_for_media_pref bias_for_reader_demand {
		gen `var'_t1 = (`var'_w2 == 1) if `var'_w2 != .
		}
		
	// generate indicator for top categories: bias_dom
	foreach var in bias_dom_govt_policy bias_dom_firm_interest bias_dom_media_pref bias_dom_reader_demand {
		gen `var'_t1 = (`var'_w2 == 1) if `var'_w2 != .
		}

	// generate relative weights for evaluating government performance
	gen evalgovt_total_w2 = evalgovt_election_w2 + evalgovt_economy_w2 + evalgovt_equality_w2 + evalgovt_ruleoflaw_w2 + evalgovt_human_rights_w2 + evalgovt_freedom_speech_w2 + evalgovt_global_power_w2 + evalgovt_fair_history_w2
	foreach var in evalgovt_election evalgovt_economy evalgovt_equality evalgovt_ruleoflaw evalgovt_human_rights evalgovt_freedom_speech evalgovt_global_power evalgovt_fair_history {
		gen r`var' = `var'_w2 / evalgovt_total_w2
		}

	// generate dummies for plan_graduation
	gen plan_grad_gradschool_dom 	= (plan_graduation_w2 == 1) if plan_graduation_w2 != .
	gen plan_grad_foreignmaster 	= (plan_graduation_w2 == 2) if plan_graduation_w2 != .
	gen plan_grad_foreignphd 		= (plan_graduation_w2 == 3) if plan_graduation_w2 != .
	gen plan_grad_military 			= (plan_graduation_w2 == 4) if plan_graduation_w2 != .
	gen plan_grad_work 				= (plan_graduation_w2 == 5) if plan_graduation_w2 != .
	gen plan_grad_dontknow 			= (plan_graduation_w2 == 6) if plan_graduation_w2 != .

	// generate top1 and top3 dummies for career preferences
	gen work_top3_missing = (work_top3_national_civil_w2 == . & work_top3_local_civil_w2 == . & work_top3_military_w2 == . & work_top3_chinese_private_w2 == . & work_top3_for_firm_w2 == . & work_top3_soe_w2 == . & work_top3_institutional_w2 == . & work_top3_entrepreneur_w2 == .)

	foreach var in national_civil local_civil military chinese_private for_firm soe institutional entrepreneur {
		gen cp_t1_`var' = (work_top3_`var'_w2 == 1) if work_top3_missing == 0
		gen cp_t3_`var' = (work_top3_`var'_w2 != .) if work_top3_missing == 0
		}
	drop work_top3_missing

	// generate dummies for career location preferences
	gen cloc_beijing 		= (place_top_w2 == 1) if place_top_w2 != .
	gen cloc_shanghai 		= (place_top_w2 == 2) if place_top_w2 != .
	gen cloc_gzsz 			= (place_top_w2 == 3 | place_top_w2 == 4) if place_top_w2 != .
	gen cloc_tjcq 			= (place_top_w2 == 5 | place_top_w2 == 6) if place_top_w2 != .
	gen cloc_hkmc 			= (place_top_w2 == 7 | place_top_w2 == 8) if place_top_w2 != .
	gen cloc_taiwan 		= (place_top_w2 == 9) if place_top_w2 != .
	gen cloc_dom 			= (place_top_w2 == 10) if place_top_w2 != .
	gen cloc_for  			= (place_top_w2 == 11) if place_top_w2 != .

	
	// add subscript: wave 2
	foreach var in news_totalcorrect news_perccorrect news_perccor_cen news_perccor_unc people_totalheard people_pcheard people_perchd_censored people_perchd_uncensor people_perchd_unctocen protest_totalheard protest_pcheard_total protest_pcheard_china protest_pcheard_foreign distneutral_us_pos_us distneutral_us_pos_cn distneutral_cn_pos_us distneutral_cn_pos_cn distneutral_us_neg_us distneutral_us_neg_cn distneutral_cn_neg_us distneutral_cn_neg_cn bias_us_pos_us_cens bias_us_pos_cn_cens bias_cn_pos_us_cens bias_cn_pos_cn_cens bias_us_neg_us_cens bias_us_neg_cn_cens bias_cn_neg_us_cens bias_cn_neg_cn_cens bias_for_govt_policy_t1 bias_for_firm_interest_t1 bias_for_media_pref_t1 bias_for_reader_demand_t1 bias_dom_govt_policy_t1 bias_dom_firm_interest_t1 bias_dom_media_pref_t1 bias_dom_reader_demand_t1 revalgovt_election revalgovt_economy revalgovt_equality revalgovt_ruleoflaw revalgovt_human_rights revalgovt_freedom_speech revalgovt_global_power revalgovt_fair_history plan_grad_gradschool_dom plan_grad_foreignmaster plan_grad_foreignphd plan_grad_military plan_grad_work plan_grad_dontknow cp_t1_national_civil cp_t3_national_civil cp_t1_local_civil cp_t3_local_civil cp_t1_military cp_t3_military cp_t1_chinese_private cp_t3_chinese_private cp_t1_for_firm cp_t3_for_firm cp_t1_soe cp_t3_soe cp_t1_institutional cp_t3_institutional cp_t1_entrepreneur cp_t3_entrepreneur cloc_beijing cloc_shanghai cloc_gzsz cloc_tjcq cloc_hkmc cloc_taiwan cloc_dom cloc_for {
		rename `var' 	`var'_w2
		}
	
	
* 	wave 3

	// foreach news, indicator correct or not
	gen news_c_coalprod  		= (news_coalprod == 0) if news_coalprod != .
	gen news_c_trumpchina 		= (news_trumpchina == 1) if news_trumpchina != .
	gen news_c_xiaojianhua 		= (news_xiaojianhua == 0) if news_xiaojianhua != .
	gen news_c_northkoreacoal  	= (news_northkoreacoal == 1) if news_northkoreacoal != .
	gen news_c_birdflu 			= (news_birdflu == 1) if news_birdflu != .
	gen news_c_ethiopiatrain 	= (news_ethiopiatrain == 0) if news_ethiopiatrain != .
	gen news_c_foreignreserve 	= (news_foreignreserve == 0) if news_foreignreserve != .
	gen news_c_xijiangcar 		= (news_xijiangcar == 0) if news_xijiangcar != .
	gen news_c_chinanorway 		= (news_chinanorway == 1) if news_chinanorway != .
	gen news_c_womenrights 		= (news_womenrights == 0) if news_womenrights != .
	gen news_c_hkceelection 	= (news_hkceelection == 0) if news_hkceelection != .
	
	// generate count variable: # news answered correctly
	gen news_totalcorrect = news_c_coalprod + news_c_trumpchina + news_c_xiaojianhua + news_c_northkoreacoal + news_c_birdflu + news_c_ethiopiatrain + news_c_foreignreserve + news_c_xijiangcar + news_c_chinanorway + news_c_womenrights + news_c_hkceelection

	// generate percentage correct
	gen news_perccorrect = news_totalcorrect / 11
	
	// correct count in percentage: censored news
	gen news_perccor_cen = (news_c_coalprod + news_c_trumpchina + news_c_xiaojianhua + news_c_xijiangcar + news_c_chinanorway + news_c_womenrights + news_c_hkceelection) / 7
	
	// correct count in percentage: uncensored news
	gen news_perccor_unc = (news_c_northkoreacoal + news_c_birdflu + news_c_ethiopiatrain + news_c_foreignreserve) / 4
	
	// generate count variable: # events heard before
	gen protest_totalheard = protest_2014_europe_square_w3 + protest_2014_sun_flower_w3 + protest_2010_arabic_spring_w3 + protest_2014_crimea_vote_w3 + protest_2012_hk_curriculum_w3 + protest_2010_catal_indep_w3 + protest_2014_umbrella_w3 + protest_2016_mongkok_riot_w3 + protest_2017_women_march_w3
	gen protest_pcheard_total = (protest_2014_europe_square_w3 + protest_2014_sun_flower_w3 + protest_2010_arabic_spring_w3 + protest_2014_crimea_vote_w3 + protest_2012_hk_curriculum_w3 + protest_2010_catal_indep_w3 + protest_2014_umbrella_w3 + protest_2016_mongkok_riot_w3 + protest_2017_women_march_w3) / 9
	
	// generate percentage by category: full panel version
	gen protest_pcheard_china = (protest_2014_sun_flower_w3 + protest_2012_hk_curriculum_w3 + protest_2014_umbrella_w3) / 3
	gen protest_pcheard_foreign = (protest_2014_europe_square_w3 + protest_2010_arabic_spring_w3 + protest_2014_crimea_vote_w3 + protest_2010_catal_indep_w3) / 4

	// generate indicator for top categories: bias_for
	foreach var in bias_for_govt_policy bias_for_firm_interest bias_for_media_pref bias_for_reader_demand {
		gen `var'_t1 = (`var'_w3 == 1) if `var'_w3 != .
		}
		
	// generate indicator for top categories: bias_dom
	foreach var in bias_dom_govt_policy bias_dom_firm_interest bias_dom_media_pref bias_dom_reader_demand {
		gen `var'_t1 = (`var'_w3 == 1) if `var'_w3 != .
		}

	// generate indicators for regression
	gen vpn_purchase_yes = (vpn_purchase < 6) if vpn_purchase != .
	
	gen vpn_purchase_wmt = ((vpn_purchase_yes == 1) & (vpn_purchase != 1)) if vpn_purchase != .
	
	gen vpn_purchase_premium = .
	replace vpn_purchase_premium = 0 	if (vpn_purchase == 6)
	replace vpn_purchase_premium = 1 	if (vpn_purchase == 1)
	replace vpn_purchase_premium = 2 	if (vpn_purchase == 4)
	replace vpn_purchase_premium = 3 	if (vpn_purchase == 3)
	replace vpn_purchase_premium = 4 	if (vpn_purchase == 5)	

	// generate dummies for plan_graduation
	gen plan_grad_gradschool_dom 	= (plan_graduation_w3 == 1) if plan_graduation_w3 != .
	gen plan_grad_foreignmaster 	= (plan_graduation_w3 == 2) if plan_graduation_w3 != .
	gen plan_grad_foreignphd 		= (plan_graduation_w3 == 3) if plan_graduation_w3 != .
	gen plan_grad_military 			= (plan_graduation_w3 == 4) if plan_graduation_w3 != .
	gen plan_grad_work 				= (plan_graduation_w3 == 5) if plan_graduation_w3 != .
	gen plan_grad_dontknow 			= (plan_graduation_w3 == 6) if plan_graduation_w3 != .

	// generate top1 and top3 dummies for career preferences
	gen work_top3_missing = (work_top3_national_civil_w3 == . & work_top3_local_civil_w3 == . & work_top3_military_w3 == . & work_top3_chinese_private_w3 == . & work_top3_for_firm_w3 == . & work_top3_soe_w3 == . & work_top3_institutional_w3 == . & work_top3_entrepreneur_w3 == .)

	foreach var in national_civil local_civil military chinese_private for_firm soe institutional entrepreneur {
		gen cp_t1_`var' = (work_top3_`var'_w3 == 1) if work_top3_missing == 0
		gen cp_t3_`var' = (work_top3_`var'_w3 != .) if work_top3_missing == 0
		}
	drop work_top3_missing

	// generate dummies for career location preferences
	gen cloc_beijing 		= (place_top_w3 == 1) if place_top_w3 != .
	gen cloc_shanghai 		= (place_top_w3 == 2) if place_top_w3 != .
	gen cloc_gzsz 			= (place_top_w3 == 3 | place_top_w3 == 4) if place_top_w3 != .
	gen cloc_tjcq 			= (place_top_w3 == 5 | place_top_w3 == 6) if place_top_w3 != .
	gen cloc_hkmc 			= (place_top_w3 == 7 | place_top_w3 == 8) if place_top_w3 != .
	gen cloc_taiwan 		= (place_top_w3 == 9) if place_top_w3 != .
	gen cloc_dom 			= (place_top_w3 == 10) if place_top_w3 != .
	gen cloc_for  			= (place_top_w3 == 11) if place_top_w3 != .

	
	// add subscript: wave 3
	foreach var in news_totalcorrect news_perccorrect news_perccor_cen news_perccor_unc protest_totalheard protest_pcheard_total protest_pcheard_china protest_pcheard_foreign bias_for_govt_policy_t1 bias_for_firm_interest_t1 bias_for_media_pref_t1 bias_for_reader_demand_t1 bias_dom_govt_policy_t1 bias_dom_firm_interest_t1 bias_dom_media_pref_t1 bias_dom_reader_demand_t1 plan_grad_gradschool_dom plan_grad_foreignmaster plan_grad_foreignphd plan_grad_military plan_grad_work plan_grad_dontknow cp_t1_national_civil cp_t3_national_civil cp_t1_local_civil cp_t3_local_civil cp_t1_military cp_t3_military cp_t1_chinese_private cp_t3_chinese_private cp_t1_for_firm cp_t3_for_firm cp_t1_soe cp_t3_soe cp_t1_institutional cp_t3_institutional cp_t1_entrepreneur cp_t3_entrepreneur cloc_beijing cloc_shanghai cloc_gzsz cloc_tjcq cloc_hkmc cloc_taiwan cloc_dom cloc_for {
		rename `var' 	`var'_w3
		}

	// generate percentage of total news quiz
	gen news_perccor_cen_all = (news_perccor_cen_w2 + news_perccor_cen_w3)/2
	

		
* 	generate certainty equivalent from Falk et al. 2015 staircase risk preference

	gen risk_preference_ce = .
	replace risk_preference_ce = 32 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 1 & risk_preference_29 == 1 & risk_preference_31 == 1)
	replace risk_preference_ce = 31 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 1 & risk_preference_29 == 1 & risk_preference_31 == 2)
	replace risk_preference_ce = 30 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 1 & risk_preference_29 == 2 & risk_preference_30 == 1)
	replace risk_preference_ce = 29 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 1 & risk_preference_29 == 2 & risk_preference_30 == 2)
	replace risk_preference_ce = 28 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 2 & risk_preference_26 == 1 & risk_preference_27 == 1)
	replace risk_preference_ce = 27 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 2 & risk_preference_26 == 1 & risk_preference_27 == 2)
	replace risk_preference_ce = 26 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 2 & risk_preference_26 == 2 & risk_preference_28 == 1)
	replace risk_preference_ce = 25 	if (risk_preference_1 == 1 & risk_preference_17 == 1 & risk_preference_25 == 2 & risk_preference_26 == 2 & risk_preference_28 == 2)
	replace risk_preference_ce = 24 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 1 & risk_preference_22 == 1 & risk_preference_23 == 1)
	replace risk_preference_ce = 23 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 1 & risk_preference_22 == 1 & risk_preference_23 == 2)
	replace risk_preference_ce = 22 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 1 & risk_preference_22 == 2 & risk_preference_24 == 1)
	replace risk_preference_ce = 21 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 1 & risk_preference_22 == 2 & risk_preference_24 == 2)
	replace risk_preference_ce = 20 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 2 & risk_preference_19 == 1 & risk_preference_20 == 1)
	replace risk_preference_ce = 19 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 2 & risk_preference_19 == 1 & risk_preference_20 == 2)
	replace risk_preference_ce = 18 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 2 & risk_preference_19 == 2 & risk_preference_21 == 1)
	replace risk_preference_ce = 17 	if (risk_preference_1 == 1 & risk_preference_17 == 2 & risk_preference_18 == 2 & risk_preference_19 == 2 & risk_preference_21 == 2)
	replace risk_preference_ce = 16 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 1 & risk_preference_14 == 1 & risk_preference_15 == 1)
	replace risk_preference_ce = 15 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 1 & risk_preference_14 == 1 & risk_preference_15 == 2)
	replace risk_preference_ce = 14 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 1 & risk_preference_14 == 2 & risk_preference_16 == 1)
	replace risk_preference_ce = 13 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 1 & risk_preference_14 == 2 & risk_preference_16 == 2)
	replace risk_preference_ce = 12 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 2 & risk_preference_11 == 1 & risk_preference_13 == 1)
	replace risk_preference_ce = 11 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 2 & risk_preference_11 == 1 & risk_preference_13 == 2)
	replace risk_preference_ce = 10 	if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 2 & risk_preference_11 == 2 & risk_preference_12 == 1)
	replace risk_preference_ce = 9 		if (risk_preference_1 == 2 & risk_preference_2  == 1 & risk_preference_10 == 2 & risk_preference_11 == 2 & risk_preference_12 == 2)
	replace risk_preference_ce = 8 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 1 & risk_preference_4  == 1 & risk_preference_5  == 1)
	replace risk_preference_ce = 7 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 1 & risk_preference_4  == 1 & risk_preference_5  == 2)
	replace risk_preference_ce = 6 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 1 & risk_preference_4  == 2 & risk_preference_6  == 1)
	replace risk_preference_ce = 5 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 1 & risk_preference_4  == 2 & risk_preference_6  == 2)
	replace risk_preference_ce = 4 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 2 & risk_preference_7  == 1 & risk_preference_8  == 1)
	replace risk_preference_ce = 3 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 2 & risk_preference_7  == 1 & risk_preference_8  == 2)
	replace risk_preference_ce = 2 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 2 & risk_preference_7  == 2 & risk_preference_9  == 1)
	replace risk_preference_ce = 1 		if (risk_preference_1 == 2 & risk_preference_2  == 2 & risk_preference_3  == 2 & risk_preference_7  == 2 & risk_preference_9  == 2)

	
* 	generate certainty equivalent from Falk et al. 2015 staircase risk preference

	gen time_preference_fe = .
	replace time_preference_fe = 32 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 1 & time_preference_22 == 1 & time_preference_23 == 1)
	replace time_preference_fe = 31 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 1 & time_preference_22 == 1 & time_preference_23 == 2)
	replace time_preference_fe = 30 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 1 & time_preference_22 == 2 & time_preference_24 == 1)
	replace time_preference_fe = 29 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 1 & time_preference_22 == 2 & time_preference_24 == 2)
	replace time_preference_fe = 28 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 2 & time_preference_19 == 1 & time_preference_20 == 1)
	replace time_preference_fe = 27 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 2 & time_preference_19 == 1 & time_preference_20 == 2)
	replace time_preference_fe = 26 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 2 & time_preference_19 == 2 & time_preference_21 == 1)
	replace time_preference_fe = 25 	if (time_preference_1 == 1 & time_preference_17 == 1 & time_preference_18 == 2 & time_preference_19 == 2 & time_preference_21 == 2)
	replace time_preference_fe = 24 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 1 & time_preference_29 == 1 & time_preference_31 == 1)
	replace time_preference_fe = 23 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 1 & time_preference_29 == 1 & time_preference_31 == 2)
	replace time_preference_fe = 22 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 1 & time_preference_29 == 2 & time_preference_30 == 1)
	replace time_preference_fe = 21 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 1 & time_preference_29 == 2 & time_preference_30 == 2)
	replace time_preference_fe = 20 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 2 & time_preference_26 == 1 & time_preference_28 == 1)
	replace time_preference_fe = 19 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 2 & time_preference_26 == 1 & time_preference_28 == 2)
	replace time_preference_fe = 18 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 2 & time_preference_26 == 2 & time_preference_27 == 1)
	replace time_preference_fe = 17 	if (time_preference_1 == 1 & time_preference_17 == 2 & time_preference_25 == 2 & time_preference_26 == 2 & time_preference_27 == 2)
	replace time_preference_fe = 16 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 1 & time_preference_14 == 1 & time_preference_16 == 1)
	replace time_preference_fe = 15 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 1 & time_preference_14 == 1 & time_preference_16 == 2)
	replace time_preference_fe = 14 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 1 & time_preference_14 == 2 & time_preference_15 == 1)
	replace time_preference_fe = 13 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 1 & time_preference_14 == 2 & time_preference_15 == 2)
	replace time_preference_fe = 12 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 2 & time_preference_11 == 1 & time_preference_13 == 1)
	replace time_preference_fe = 11 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 2 & time_preference_11 == 1 & time_preference_13 == 2)
	replace time_preference_fe = 10 	if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 2 & time_preference_11 == 2 & time_preference_12 == 1)
	replace time_preference_fe = 9 		if (time_preference_1 == 2 & time_preference_2  == 1 & time_preference_10 == 2 & time_preference_11 == 2 & time_preference_12 == 2)
	replace time_preference_fe = 8 		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 1 & time_preference_7  == 1 & time_preference_8  == 1)
	replace time_preference_fe = 7 		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 1 & time_preference_7  == 1 & time_preference_8  == 2)
	replace time_preference_fe = 6 		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 1 & time_preference_7  == 2 & time_preference_9  == 1)
	replace time_preference_fe = 5 		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 1 & time_preference_7  == 2 & time_preference_9  == 2)
	replace time_preference_fe = 4 		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 2 & time_preference_4  == 1 & time_preference_6  == 1)
	replace time_preference_fe = 3		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 2 & time_preference_4  == 1 & time_preference_6  == 2)
	replace time_preference_fe = 2		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 2 & time_preference_4  == 2 & time_preference_5  == 1)
	replace time_preference_fe = 1		if (time_preference_1 == 2 & time_preference_2  == 2 & time_preference_3  == 2 & time_preference_4  == 2 & time_preference_5  == 2)
		
	
	
* 	demographic and background characteristics

	gen ethnicity_han = (ethnicity == 1) if ethnicity != .
	gen hukou_urban = (hukou == 1) if hukou != .
	gen siblings_total = number_bro_younger + number_bro_older + number_sis_younger + number_sis_older
	gen father_edu_hsabove = (father_education > 3) if father_education != .
	gen mother_edu_hsabove = (mother_education > 3) if mother_education != .
	gen religion_religious = (religion > 1) if religion != .
	gen hs_track_science = (hs_track == 1) 	if hs_track != .

	// indicator of university
	gen university_elite = (university == 1) if university != .
	
	// recode birth_year
	recode birth_year (1=1990) (2=1991) (3=1992) (4=1993) (5=1994) (6=1995) (7=1996) (8=1997) (9=1998) (10=1999) (11=2000) (12=2001) (13=2002) (14=2003) (15=2004) (16=2005)

	// indicator of coastal provinces
	gen birthplace_coastal = (birthplace_province == "Beijing" | birthplace_province == "Fujian" | birthplace_province == "Guangdong" | birthplace_province == "Hainan" | birthplace_province == "Hebei" | birthplace_province == "Jiangsu" | birthplace_province == "Shandong" | birthplace_province == "Shanghai" | birthplace_province == "Tianjin" | birthplace_province == "Zhejiang" | birthplace_province == "Non-mainland")
	gen residence_coastal = (residence_province == "Beijing" | residence_province == "Fujian" | residence_province == "Guangdong" | residence_province == "Hainan" | residence_province == "Hebei" | residence_province == "Jiangsu" | residence_province == "Shandong" | residence_province == "Shanghai" | residence_province == "Tianjin" | residence_province == "Zhejiang" | residence_province == "Non-mainland")
	
	// indicator of english credentials
	gen domestic_english_atleast4 = (english_qual_domestic > 1) if english_qual_domestic != .
	gen foreign_english_yes = (english_qual_foreign < 3) if english_qual_foreign != .
	
	// indicator of parents' work sector
	gen work_father_govt = (father_work <= 4) if father_work != .
	gen work_mother_govt = (mother_work <= 4) if mother_work != .
	
	// indicator of travel experience
	gen travel_foreign_yes = (travel_foreign > 1) if travel_foreign != .
	
	
	
* 	treatment related indicators

	// indicators and dummies
	gen treatment_control 		= (treatment_newsletter == 0 & treatment_vpn == 0 & vpn_current_paid_user == 0)
	gen treatment_vpnonly 		= (treatment_newsletter == 0 & treatment_vpn == 1 & vpn_current_paid_user == 0)
	gen treatment_nlonly 		= (treatment_newsletter == 1 & treatment_vpn == 0 & vpn_current_paid_user == 0)
	gen treatment_vpnnl  		= (treatment_newsletter == 1 & treatment_vpn == 1 & vpn_current_paid_user == 0)
	gen treatment_user 			= (vpn_current_paid_user == 1)
	
	gen treatment_master = .
	replace treatment_master = 1 	if treatment_control == 1
	replace treatment_master = 2 	if treatment_vpnonly == 1
	replace treatment_master = 3 	if treatment_nlonly == 1
	replace treatment_master = 4 	if treatment_vpnnl == 1
	replace treatment_master = 5 	if treatment_user == 1
	label def treatment_master 1 "Group-C" 2 "Group-A" 3 "Group-CE" 4 "Group-AE" 5 "Existing users"
	label value treatment_master treatment_master
	
	gen treatment_main = .
	replace treatment_main = 1 	if (treatment_control == 1 | treatment_nlonly == 1)
	replace treatment_main = 2 	if treatment_vpnonly == 1
	replace treatment_main = 3 	if treatment_vpnnl == 1
	replace treatment_main = 4 	if treatment_user == 1
	label def treatment_main 1 "Control" 2 "Access" 3 "Access + Encour." 4 "Existing users"
	label value treatment_main treatment_main

	// treatment interactions	
	gen treatvpnXnl = treatment_vpn * treatment_newsletter

	
	
* 	vpn adoption status
	gen vpn_adopted = (vpn_date_adoption != .)	if vpn_date_adoption != .
	replace vpn_adopted = 0 					if vpn_adopted != 1 & (treatment_master == 2 | treatment_master == 4)

		
* 	roommates' vpn usage

	// generate social learning primitives
	gen vpn_roommate_existing = vpn_roommate_usage_w1
	gen vpn_roommate_new_w2 = vpn_roommate_usage_w2 - vpn_roommate_usage_w1
	gen vpn_roommate_new_w3 = vpn_roommate_usage_w3 - vpn_roommate_usage_w1
	
	// top code at 2
	replace vpn_roommate_existing = 2 		if vpn_roommate_existing > 2
	replace vpn_roommate_new_w2 = 0 		if vpn_roommate_new_w2 < 0
	replace vpn_roommate_new_w2 = 2 		if vpn_roommate_new_w2 > 2
	replace vpn_roommate_new_w3 = 0 		if vpn_roommate_new_w3 < 0
	replace vpn_roommate_new_w3 = 2 		if vpn_roommate_new_w3 > 2
	
	// finalize primitives
	gen soclearning_ownaccess = (treatment_master >= 4)		if treatment_master != .
	gen soclearning_rm_new_w2 = (vpn_roommate_new_w2) 		if vpn_roommate_new_w2 != .
	gen soclearning_rm_new_w3 = (vpn_roommate_new_w3) 		if vpn_roommate_new_w3 != .
	gen soclearning_rm_existing = (vpn_roommate_existing) 	if vpn_roommate_existing != .
	gen soclearning_ownXnew_w2 = soclearning_ownaccess * soclearning_rm_new_w2
	gen soclearning_ownXnew_w3 = soclearning_ownaccess * soclearning_rm_new_w3
		
		
		
* 	list experiment

	// list experiment: indicator for provision of veil
	gen list_exp_veiled = 1 - list_exp_direct

	// list experiment: calculated overall count for direct group
	foreach c in trust {
		gen list_countd_`c'_direct = list_`c'_direct_count + list_`c'_direct_yes
		gen list_count_`c' = .
		replace list_count_`c' = list_countd_`c'_direct 	if list_exp_direct == 1
		replace list_count_`c' = list_`c'_veiled_count 		if list_exp_direct == 0
		}
	

	
			
*** Prepare data: generate az-scores

* 	A: Beliefs and attitudes regarding media
	
	// A.2: Purchase of censorship circumvention tools
	local var_vpn_purchase_reg_w3 			`"vpn_purchase_wmt_record vpn_purchase_yes"'
	
	andersonz `var_vpn_purchase_reg_w3', gen(az_var_vpnpurchase)

	// A.3: Valuation of access to foreign media outlets
	local var_media_valuation_reg_w1 		`"wtp_vpn_w1 added_value_foreign_media_w1"'
	local var_media_valuation_reg_w2 		`"wtp_vpn_w2 added_value_foreign_media_w2"'
	local var_media_valuation_reg_w3 		`"wtp_vpn_w3 added_value_foreign_media_w3"'

	andersonz `var_media_valuation_reg_w1', gen(az_belief_media_value_w1)
	andersonz `var_media_valuation_reg_w2', gen(az_belief_media_value_w2)
	andersonz `var_media_valuation_reg_w3', gen(az_belief_media_value_w3)
	
	// A.4: Trust in media outlets
	local var_media_trust_reg_w1 			`"trust_media_dom_state_w1 trust_media_dom_private_w1 trust_media_foreign_w1"'
	local var_media_trust_reg_w2 			`"trust_media_dom_state_w2 trust_media_dom_private_w2 trust_media_foreign_w2"'
	local var_media_trust_reg_w3 			`"trust_media_dom_state_w3 trust_media_dom_private_w3 trust_media_foreign_w3"'

	andersonz `var_media_trust_reg_w1', gen(az_belief_media_trust_w1)
	andersonz `var_media_trust_reg_w2', gen(az_belief_media_trust_w2)
	andersonz `var_media_trust_reg_w3', gen(az_belief_media_trust_w3)
	
	// A.5: Calibration of news outlets' level of censorship and biases
	local var_percmediabias_reg_ce_cn_w2 	`"bias_cn_neg_cn_cens_w2 bias_cn_pos_cn_cens_w2 bias_us_neg_cn_cens_w2 bias_us_pos_cn_cens_w2"'
	local var_percmediabias_reg_ce_us_w2 	`"bias_cn_neg_us_cens_w2 bias_cn_pos_us_cens_w2 bias_us_neg_us_cens_w2 bias_us_pos_us_cens_w2"'
	local var_percmediabias_reg_di_cn_w2 	`"distneutral_cn_neg_cn_w2 distneutral_cn_pos_cn_w2 distneutral_us_neg_cn_w2 distneutral_us_pos_cn_w2"'
	local var_percmediabias_reg_di_us_w2 	`"distneutral_cn_neg_us_w2 distneutral_cn_pos_us_w2 distneutral_us_neg_us_w2 distneutral_us_pos_us_w2"'

	andersonz `var_percmediabias_reg_ce_cn_w2', gen(az_belief_media_cens_cn)
	andersonz `var_percmediabias_reg_ce_us_w2', gen(az_belief_media_cens_us)
	andersonz `var_percmediabias_reg_di_cn_w2', gen(az_belief_media_bias_cn)
	andersonz `var_percmediabias_reg_di_us_w2', gen(az_belief_media_bias_us)
	
	// A.6: Justification of media censorship
	local var_censor_justif_reg_w1 			`"censor_just_dom_economic_w1 censor_just_dom_political_w1 censor_just_dom_social_w1 censor_just_for_w1"'
	local var_censor_justif_reg_w2 			`"censor_just_dom_economic_w2 censor_just_dom_political_w2 censor_just_dom_social_w2 censor_just_for_w2 censor_just_porn_w2"'
	local var_censor_justif_reg_w3 			`"censor_just_dom_economic_w3 censor_just_dom_political_w3 censor_just_dom_social_w3 censor_just_for_w3 censor_just_porn_w3"'
	local var_censor_justif_reg_fp_w1 		`"censor_just_dom_economic_w1 censor_just_dom_political_w1 censor_just_dom_social_w1 censor_just_for_w1"'
	local var_censor_justif_reg_fp_w2 		`"censor_just_dom_economic_w2 censor_just_dom_political_w2 censor_just_dom_social_w2 censor_just_for_w2"'
	local var_censor_justif_reg_fp_w3 		`"censor_just_dom_economic_w3 censor_just_dom_political_w3 censor_just_dom_social_w3 censor_just_for_w3"'

	andersonz `var_censor_justif_reg_w1', gen(az_belief_media_justif_w1)
	andersonz `var_censor_justif_reg_w2', gen(az_belief_media_justif_w2)
	andersonz `var_censor_justif_reg_w3', gen(az_belief_media_justif_w3)
	andersonz `var_censor_justif_reg_fp_w1', gen(az_belief_media_justif_pf_w1)
	andersonz `var_censor_justif_reg_fp_w2', gen(az_belief_media_justif_pf_w2)
	andersonz `var_censor_justif_reg_fp_w3', gen(az_belief_media_justif_pf_w3)
	
	// A.7: Belief regarding drivers of media censorship
    local var_censor_driver_reg_dom_w2 		`"bias_dom_govt_policy_t1_w2 bias_dom_firm_interest_t1_w2 bias_dom_media_pref_t1_w2 bias_dom_reader_demand_t1_w2"'
	local var_censor_driver_reg_for_w2		`"bias_for_govt_policy_t1_w2 bias_for_firm_interest_t1_w2 bias_for_media_pref_t1_w2 bias_for_reader_demand_t1_w2"'

	andersonz `var_censor_driver_reg_dom_w2', gen(az_belief_media_dri_dom_w2)
	andersonz `var_censor_driver_reg_for_w2', gen(az_belief_media_dri_for_w2)
	
	
* 	B: Knowledge

	// B.1: Current news events covered in the demand treatment
	local var_knowledge_news_reg_quiz 		`"news_c_topinequality news_c_cpicensorship news_c_laborunrest news_c_waterpollution"'
	
	andersonz `var_knowledge_news_reg_quiz', gen(az_knowledge_news_quiz)
	
	// B.2: Current news events not covered in the demand treatment
	local var_knowledge_news_reg_cen_w2 	`"news_c_panamapapers news_c_tenyearshk news_c_stockcrash news_c_economistcensor"'
	local var_knowledge_news_reg_unc_w2 	`"news_c_applefbi news_c_taiwanelection news_c_yihehotel"'
	local var_knowledge_news_reg_cen_w1 	`"news_c_stock_rise news_c_train_brazil_peru news_c_taiwan_election news_c_china_us_network_coop"'
	local var_knowledge_news_reg_unc_w1 	`"news_c_exchange_depreciate news_c_army_reduce news_c_lijiacheng_china news_c_nanjing_memory_list"'

	andersonz `var_knowledge_news_reg_cen_w2', gen(az_knowledge_news_censored)
	andersonz `var_knowledge_news_reg_unc_w2', gen(az_knowledge_news_uncensor)
	andersonz `var_knowledge_news_reg_cen_w1', gen(az_knowledge_news_censored_w1)
	andersonz `var_knowledge_news_reg_unc_w1', gen(az_knowledge_news_uncensor_w1)
	
	// B.3: Awareness of notable figures
	local var_knowledge_people_reg_tocens 	`"people_puzhiqiang_w2 people_renzhiqiang_w2 people_huangzhifeng_w2"'
	local var_knowledge_people_reg_censor 	`"people_lizehou_w2 people_chenguangcheng_w2 people_lixiaolin_w2"'
	local var_knowledge_people_reg_uncens 	`"people_maoyushi_w2 people_honghuang_w2 people_liuqiangdong_w2"'
	local var_knowledge_people_censor_w1 	`"people_lizehou_w1 people_huangzhifeng_w1 people_chenguangcheng_w1 people_lixiaolin_w1"'
	local var_knowledge_people_uncens_w1 	`"people_puzhiqiang_w1 people_renzhiqiang_w1 people_maoyushi_w1 people_honghuang_w1 people_liuqiangdong_w1"'
	
	andersonz `var_knowledge_people_reg_tocens', gen(az_knowledge_people_tocens)
	andersonz `var_knowledge_people_reg_censor', gen(az_knowledge_people_censor)
	andersonz `var_knowledge_people_reg_uncens', gen(az_knowledge_people_uncens)
	andersonz `var_knowledge_people_censor_w1', gen(az_knowledge_people_censor_w1)
	andersonz `var_knowledge_people_uncens_w1', gen(az_knowledge_people_uncens_w1)

	// B.4: Awareness of protest events
	local var_knowledge_prot_reg_chi_w1 	`"protest_2012_hk_curriculum_w1 protest_2014_umbrella_w1 protest_2014_sun_flower_w1"'
	local var_knowledge_prot_reg_chi_w2  	`"protest_2012_hk_curriculum_w2 protest_2014_umbrella_w2 protest_2016_mongkok_riot_w2 protest_2014_sun_flower_w2"'
	local var_knowledge_prot_reg_for_w2 	`"protest_2014_europe_square_w2 protest_2010_arabic_spring_w2 protest_2014_crimea_vote_w2 protest_2010_catal_indep_w2"'
	
	andersonz `var_knowledge_prot_reg_chi_w2', gen(az_knowledge_protest_china)
	andersonz `var_knowledge_prot_reg_for_w2', gen(az_knowledge_protest_forei)
	andersonz `var_knowledge_prot_reg_chi_w1', gen(az_knowledge_protest_china_w1)
	
	// B.5: Meta-knowledge
	local var_knowledge_meta_reg_w1 		`"familiar_china_issues_self_w1 familiar_china_others_w1"'
	local var_knowledge_meta_reg_w2 		`"familiar_china_issues_self_w2 familiar_china_others_w2"'
	local var_knowledge_meta_reg_w3 		`"familiar_china_issues_self_w3 familiar_china_others_w3"'
	
	andersonz `var_knowledge_meta_reg_w1', gen(az_knowledge_meta_w1)
	andersonz `var_knowledge_meta_reg_w2', gen(az_knowledge_meta_w2)
	andersonz `var_knowledge_meta_reg_w3', gen(az_knowledge_meta_w3)
	

* 	C: Economic beliefs
	
	// C.1: Belief on economic performance in China
	local var_econ_guess_reg_cn_perf_w1 	`"guess_gdp_growth_china_w1 guess_stock_index_sh_w1"'
	local var_econ_guess_reg_cn_perf_w2 	`"guess_gdp_growth_china_w2 guess_stock_index_sh_w2"'
	local var_econ_guess_reg_cn_perf_w3 	`"guess_gdp_growth_china_w3 guess_stock_index_sh_w3"'
	
	andersonz `var_econ_guess_reg_cn_perf_w1', gen(az_belief_econ_perf_cn_w1)
	andersonz `var_econ_guess_reg_cn_perf_w2', gen(az_belief_econ_perf_cn_w2)
	andersonz `var_econ_guess_reg_cn_perf_w3', gen(az_belief_econ_perf_cn_w3)

	// C.2: Confidence on guesses regarding economic performance in China
	local var_econ_guess_reg_cn_conf_w1 	`"guess_gdp_growth_china_con_w1 guess_stock_index_sh_con_w1"'
	local var_econ_guess_reg_cn_conf_w2 	`"guess_gdp_growth_china_con_w2 guess_stock_index_sh_con_w2"'
	local var_econ_guess_reg_cn_conf_w3 	`"guess_gdp_growth_china_con_w3 guess_stock_index_sh_con_w3"'
	
	andersonz `var_econ_guess_reg_cn_conf_w1', gen(az_belief_econ_conf_cn_w1)
	andersonz `var_econ_guess_reg_cn_conf_w2', gen(az_belief_econ_conf_cn_w2)
	andersonz `var_econ_guess_reg_cn_conf_w3', gen(az_belief_econ_conf_cn_w3)

	// C.3: Belief on economic performance in the US
	local var_econ_guess_reg_us_perf_w2 	`"guess_gdp_growth_us_w2 guess_stock_index_dj_w2"'
	local var_econ_guess_reg_us_perf_w3 	`"guess_gdp_growth_us_w3 guess_stock_index_dj_w3"'
	
	andersonz `var_econ_guess_reg_us_perf_w2', gen(az_belief_econ_perf_us_w2)
	andersonz `var_econ_guess_reg_us_perf_w3', gen(az_belief_econ_perf_us_w3)
	
	// C.4: Confidence on guesses regarding economic performance in the US
	local var_econ_guess_reg_us_conf_w2 	`"guess_gdp_growth_us_con_w2 guess_stock_index_dj_con_w2"'
	local var_econ_guess_reg_us_conf_w3 	`"guess_gdp_growth_us_con_w3 guess_stock_index_dj_con_w3"'
	
	andersonz `var_econ_guess_reg_us_conf_w2', gen(az_belief_econ_conf_us_w2)
	andersonz `var_econ_guess_reg_us_conf_w3', gen(az_belief_econ_conf_us_w3)
	
	
* 	D: Political attitudes

	// D.1: Demand for institutional change
	local var_demand_change_reg_w1 			`"inst_change_econ_w1 inst_change_poli_w1"'
	local var_demand_change_reg_w2 			`"inst_change_econ_w2 inst_change_poli_w2"'
	local var_demand_change_reg_w3 			`"inst_change_econ_w3 inst_change_poli_w3"'
	
	andersonz `var_demand_change_reg_w1', gen(az_belief_instchange_w1)
	andersonz `var_demand_change_reg_w2', gen(az_belief_instchange_w2)
	andersonz `var_demand_change_reg_w3', gen(az_belief_instchange_w3)
				
	// D.2: Trust in institutions
	local var_trust_inst_reg_govt_w1 		`"trust_central_govt_w1 trust_provincial_govt_w1 trust_local_govt_w1"'
	local var_trust_inst_reg_govt_w2 		`"trust_central_govt_w2 trust_provincial_govt_w2 trust_local_govt_w2"'
	local var_trust_inst_reg_govt_w3 		`"trust_central_govt_w3 trust_provincial_govt_w3 trust_local_govt_w3"'
	local var_trust_inst_reg_foreign_w1 	`"trust_japan_govt_w1 trust_us_govt_w1"'
	local var_trust_inst_reg_foreign_w2 	`"trust_japan_govt_w2 trust_us_govt_w2"'
	local var_trust_inst_reg_foreign_w3 	`"trust_japan_govt_w3 trust_us_govt_w3"'
	local var_trust_inst_reg_copo_w2 		`"trust_court_w2 trust_police_w2"'
	
	andersonz `var_trust_inst_reg_govt_w1', gen(az_belief_trust_govt_w1)
	andersonz `var_trust_inst_reg_govt_w2', gen(az_belief_trust_govt_w2)
	andersonz `var_trust_inst_reg_govt_w3', gen(az_belief_trust_govt_w3)
	andersonz `var_trust_inst_reg_foreign_w1', gen(az_belief_trust_foreign_w1)
	andersonz `var_trust_inst_reg_foreign_w2', gen(az_belief_trust_foreign_w2)
	andersonz `var_trust_inst_reg_foreign_w3', gen(az_belief_trust_foreign_w3)
	andersonz `var_trust_inst_reg_copo_w2', gen(az_belief_trust_copo_w2)
	
	// D.3: Evaluation of government's performance
	local var_eval_govt_reg_w1 				`"eval_govt_economic_w1 eval_govt_dom_politics_w1 eval_govt_for_relations_w1"'
	local var_eval_govt_reg_w2 				`"eval_govt_economic_w2 eval_govt_dom_politics_w2 eval_govt_for_relations_w2"'
	local var_eval_govt_reg_w3 				`"eval_govt_economic_w3 eval_govt_dom_politics_w3 eval_govt_for_relations_w3"'
	
	andersonz `var_eval_govt_reg_w1', gen(az_belief_evalgovt_w1)
	andersonz `var_eval_govt_reg_w2', gen(az_belief_evalgovt_w2)
	andersonz `var_eval_govt_reg_w3', gen(az_belief_evalgovt_w3)

	// D.4: Performance evaluation criteria
	local var_eval_criteria_reg_w2 			`"revalgovt_election_w2 revalgovt_economy_w2 revalgovt_equality_w2 revalgovt_ruleoflaw_w2 revalgovt_human_rights_w2 revalgovt_freedom_speech_w2 revalgovt_global_power_w2 revalgovt_fair_history_w2"'

	andersonz `var_eval_criteria_reg_w2', gen(az_belief_evalcrit_w2)
	
	// D.5: Evaluation of severity of socioeconomic issues
	local var_severity_reg_w2 				`"severity_welfare_w2 severity_employment_w2 severity_pollution_w2 severity_inequality_w2 severity_corruption_w2 severity_dscrm_minority_w2"'

	andersonz `var_severity_reg_w2', gen(az_belief_severity_w2)

	// D.6: Evaluation of democracy and human rights protection in China
	local var_democracy_reg_fp_w1 			`"importance_live_in_demo_w1"'
	local var_democracy_reg_fp_w2 			`"importance_live_in_demo_w2"'
	local var_democracy_reg_fp_w3 			`"importance_live_in_demo_w3"'
	local var_democracy_reg_w2 				`"china_interest_group_w2 china_rate_democracy_w2 china_rate_humanrights_w2 importance_live_in_demo_w2"'
	
	andersonz `var_democracy_reg_fp_w1', gen(az_belief_democracy_fp_w1)
	andersonz `var_democracy_reg_fp_w2', gen(az_belief_democracy_fp_w2)
	andersonz `var_democracy_reg_fp_w3', gen(az_belief_democracy_fp_w3)
	andersonz `var_democracy_reg_w2', gen(az_belief_democracy_w2)
	
	// D.7: Justification of controversial policies and issues
	local var_contro_justi_reg_policy_w2 	`"justify_minority_policy_w2 justify_hukou_w2 justify_one_child_w2 justify_hongkong_policy_w2 justify_taiwan_policy_w2 justify_violence_stability_w2 justify_receive_refugee_w2 justify_reduce_pollution_w2 justify_gaokao_w2 justify_soe_privatize_w2"'
	local var_contro_justi_reg_liberal_w2 	`"justify_homo_marriage_w2 justify_legal_prostitute_w2 justify_abortion_w2 justify_exmaritalsex justify_transgene_w2 justify_soft_drug_w2"'

	andersonz `var_contro_justi_reg_policy_w2', gen(az_belief_justify_policy_w2)
	andersonz `var_contro_justi_reg_liberal_w2', gen(az_belief_justify_liberal_w2)
	
	// D.8: Willingness to act
	local var_willing_fight_reg_w1 			`"willing_against_illi_govt_w1 willing_report_mis_w1 willing_protect_weak_w1"'
	local var_willing_fight_reg_w2 			`"willing_against_illi_govt_w2 willing_report_mis_w2 willing_protect_weak_w2"'
	local var_willing_fight_reg_w3 			`"willing_against_illi_govt_w3 willing_report_mis_w3 willing_protect_weak_w3"'
	
	andersonz `var_willing_fight_reg_w1', gen(az_belief_willing_w1)
	andersonz `var_willing_fight_reg_w2', gen(az_belief_willing_w2)
	andersonz `var_willing_fight_reg_w3', gen(az_belief_willing_w3)

	// D.9: Interest in politics and economics
	local var_interest_reg_w2 				`"interest_economic_w2 interest_politics_w2"'

	andersonz `var_interest_reg_w2', gen(az_belief_interest_w2)
	
	
	
* 	E: behaviors

	// E.1: Social interactions
	local var_socialinteract_reg_w1 		`"frequency_talk_politic_w1 frequency_persuade_friends_w1"'
	local var_socialinteract_reg_w2 		`"frequency_talk_politic_w2 frequency_persuade_friends_w2"'
	local var_socialinteract_reg_w3 		`"frequency_talk_politic_w3 frequency_persuade_friends_w3"'
	
	andersonz `var_socialinteract_reg_w1', gen(az_var_socialinteract_w1)
	andersonz `var_socialinteract_reg_w2', gen(az_var_socialinteract_w2)
	andersonz `var_socialinteract_reg_w3', gen(az_var_socialinteract_w3)	
	
	// E.2: Participation
	local var_polparticipation_reg_pf_w1 	`"participate_social_protest_w1 participate_plan_vote_w1 participate_complain_school_w1"'
	local var_polparticipation_reg_pf_w2 	`"participate_social_protest_w2 participate_plan_vote_w2 participate_complain_school_w2"'
	local var_polparticipation_reg_pf_w3 	`"participate_social_protest_w3 participate_plan_vote_w3 participate_complain_school_w3"'
	
	andersonz `var_polparticipation_reg_pf_w1', gen(az_var_polparticipation_w1)
	andersonz `var_polparticipation_reg_pf_w2', gen(az_var_polparticipation_w2)
	andersonz `var_polparticipation_reg_pf_w3', gen(az_var_polparticipation_w3)

	
	
* 	F: Demographics, background characteristics, and fundamental preferences
	
	// F.1: Personal characteristics
	local var_demog_reg_personal 			`"gender birth_year ethnicity_han birthplace_coastal residence_coastal hukou_urban religion_religious ccp_member"'

	andersonz `var_demog_reg_personal', gen(az_demographics_personal)

	// F.2: Educational background
	local var_demog_reg_education 			`"university_elite hs_track_science department_ssh"'

	andersonz `var_demog_reg_education', gen(az_demographics_education)
	
	// F.3: English ability and oversea travel experiences
	local var_demog_reg_english 			`"domestic_english_atleast4 foreign_english_yes"'
	local var_demog_reg_travel 				`"travel_hktaiwan travel_foreign_yes"'
	
	andersonz `var_demog_reg_english', gen(az_demographics_english)
	andersonz `var_demog_reg_travel', gen(az_demographics_travel)
	
	// F.4: Household characteristics
	local var_demog_reg_household 			`"siblings_total father_edu_hsabove work_father_govt father_ccp mother_edu_hsabove work_mother_govt mother_ccp hh_income"'

	andersonz `var_demog_reg_household', gen(az_demographics_household)
	
	// F.5: Fundamental preferences
	local var_preference_risk 				`"willing_risk risk_preference_ce risk_lottery_choice"'
	local var_preference_time 				`"willing_future procrastinate"'
	local var_preference_altruism 			`"willing_goodcauses donate_amount"'
	local var_preference_reciprocity 		`"willing_returnfavor best_intentions gift_amount willing_punish_you willing_punish_others willing_revenge"'

	andersonz `var_preference_risk', gen(az_preference_risk)
	andersonz `var_preference_time', gen(az_preference_time)
	andersonz `var_preference_altruism', gen(az_preference_altruism)
	andersonz `var_preference_reciprocity', gen(az_preference_reciprocity)

	
	
* 	X: Overall effect index

	// flip variables to generate overall impact index
	forvalues i = 1/3 {
		gen bias_foreign_r_w`i' = 10 - bias_foreign_w`i'
		gen az_belief_econ_perf_cn_r_w`i' = 0 - az_belief_econ_perf_cn_w`i'
		gen az_belief_trust_govt_r_w`i' = 0 - az_belief_trust_govt_w`i'
		gen az_belief_evalgovt_r_w`i' = 0 - az_belief_evalgovt_w`i'
		gen importance_live_in_demo_r_w`i' = 10 - importance_live_in_demo_w`i'
		}

	local var_overall_w1 					`"info_freq_website_for_w1 az_belief_media_value_w1 az_belief_media_trust_w1 bias_domestic_w1 az_belief_media_justif_w1 bias_foreign_r_w1 news_perccor_cen_w1 news_perccor_unc_w1 protest_pcheard_china_w1 protest_pcheard_foreign_w1 az_knowledge_meta_w1 az_belief_econ_perf_cn_r_w1 az_belief_econ_conf_cn_w1 az_belief_instchange_w1 az_belief_trust_govt_r_w1 az_belief_trust_foreign_w1 az_belief_evalgovt_r_w1 importance_live_in_demo_r_w1 az_belief_willing_w1 az_var_socialinteract_w1 az_var_polparticipation_w1 plan_grad_gradschool_dom_w1 plan_grad_foreignmaster_w1 plan_grad_foreignphd_w1 plan_grad_military_w1 plan_grad_work_w1 cloc_beijing_w1 cloc_shanghai_w1 cloc_gzsz_w1 cloc_tjcq_w1 cloc_hkmc_w1 cloc_taiwan_w1 cloc_dom_w2 cloc_for_w1 stock_participation_w1"'
	local var_overall_w2 					`"info_freq_website_for_w2 az_belief_media_value_w2 az_belief_media_trust_w2 bias_domestic_w2 az_belief_media_justif_w2 bias_foreign_r_w2 news_perccor_cen_w2 news_perccor_unc_w2 protest_pcheard_china_w2 protest_pcheard_foreign_w2 az_knowledge_meta_w2 az_belief_econ_perf_cn_r_w2 az_belief_econ_conf_cn_w2 az_belief_instchange_w2 az_belief_trust_govt_r_w2 az_belief_trust_foreign_w2 az_belief_evalgovt_r_w2 importance_live_in_demo_r_w2 az_belief_willing_w2 az_var_socialinteract_w2 az_var_polparticipation_w2 plan_grad_gradschool_dom_w2 plan_grad_foreignmaster_w2 plan_grad_foreignphd_w2 plan_grad_military_w2 plan_grad_work_w2 cloc_beijing_w2 cloc_shanghai_w2 cloc_gzsz_w2 cloc_tjcq_w2 cloc_hkmc_w2 cloc_taiwan_w2 cloc_dom_w2 cloc_for_w2 stock_participation_w2"'
	local var_overall_w3 					`"info_freq_website_for_w3 az_belief_media_value_w3 az_belief_media_trust_w3 bias_domestic_w3 az_belief_media_justif_w3 bias_foreign_r_w3 news_perccor_cen_w3 news_perccor_unc_w3 protest_pcheard_china_w3 protest_pcheard_foreign_w3 az_knowledge_meta_w3 az_belief_econ_perf_cn_r_w3 az_belief_econ_conf_cn_w3 az_belief_instchange_w3 az_belief_trust_govt_r_w3 az_belief_trust_foreign_w3 az_belief_evalgovt_r_w3 importance_live_in_demo_r_w3 az_belief_willing_w3 az_var_socialinteract_w3 az_var_polparticipation_w3 plan_grad_gradschool_dom_w3 plan_grad_foreignmaster_w3 plan_grad_foreignphd_w3 plan_grad_military_w3 plan_grad_work_w3 cloc_beijing_w3 cloc_shanghai_w3 cloc_gzsz_w3 cloc_tjcq_w3 cloc_hkmc_w3 cloc_taiwan_w3 cloc_dom_w3 cloc_for_w3 stock_participation_w3"'

	andersonz `var_overall_w1', gen(az_overall_w1)
	andersonz `var_overall_w2', gen(az_overall_w2)
	andersonz `var_overall_w3', gen(az_overall_w3)
	
	// five main categories: wave 1
	andersonz info_freq_website_for_w1 az_belief_media_value_w1 az_belief_media_trust_w1 az_belief_media_justif_pf_w1 bias_domestic_w1 bias_foreign_r_w1, gen(az_overall_a_w1)
	andersonz news_perccor_cen_w1 news_perccor_unc_w1 protest_pcheard_china_w1 protest_pcheard_foreign_w1 az_knowledge_meta_w1, gen(az_overall_b_w1)
	andersonz az_belief_econ_perf_cn_r_w1 az_belief_econ_conf_cn_w1, gen(az_overall_c_w1)
	andersonz az_belief_instchange_w1 az_belief_trust_govt_r_w1 az_belief_trust_foreign_w1 az_belief_evalgovt_r_w1 importance_live_in_demo_r_w1 az_belief_willing_w1, gen(az_overall_d_w1)
	andersonz az_var_socialinteract_w1 az_var_polparticipation_w1 plan_grad_foreignmaster_w1 cloc_for_w1 stock_participation_w1, gen(az_overall_e_w1)
	
	// five main categories: wave 2
	andersonz info_freq_website_for_w2 az_belief_media_value_w2 az_belief_media_trust_w2 az_belief_media_justif_pf_w2 bias_domestic_w2 bias_foreign_r_w2, gen(az_overall_a_w2)
	andersonz news_perccor_cen_w2 news_perccor_unc_w2 protest_pcheard_china_w2 protest_pcheard_foreign_w2 az_knowledge_meta_w2, gen(az_overall_b_w2)
	andersonz az_belief_econ_perf_cn_r_w2 az_belief_econ_conf_cn_w2, gen(az_overall_c_w2)
	andersonz az_belief_instchange_w2 az_belief_trust_govt_r_w2 az_belief_trust_foreign_w2 az_belief_evalgovt_r_w2 importance_live_in_demo_r_w2 az_belief_willing_w2, gen(az_overall_d_w2)
	andersonz az_var_socialinteract_w2 az_var_polparticipation_w2 plan_grad_foreignmaster_w2 cloc_for_w2 stock_participation_w2, gen(az_overall_e_w2)

	// five main categories: wave 3
	andersonz info_freq_website_for_w3 az_belief_media_value_w3 az_belief_media_trust_w3 az_belief_media_justif_pf_w3 bias_domestic_w3 bias_foreign_r_w3, gen(az_overall_a_w3)
	andersonz news_perccor_cen_w3 news_perccor_unc_w3 protest_pcheard_china_w3 protest_pcheard_foreign_w3 az_knowledge_meta_w3, gen(az_overall_b_w3)
	andersonz az_belief_econ_perf_cn_r_w3 az_belief_econ_conf_cn_w3, gen(az_overall_c_w3)
	andersonz az_belief_instchange_w3 az_belief_trust_govt_r_w3 az_belief_trust_foreign_w3 az_belief_evalgovt_r_w3 importance_live_in_demo_r_w3 az_belief_willing_w3, gen(az_overall_d_w3)
	andersonz az_var_socialinteract_w3 az_var_polparticipation_w3 plan_grad_foreignmaster_w3 cloc_for_w3 stock_participation_w3, gen(az_overall_e_w3)
	
	// mega z-scores: wave 3
	andersonz az_overall_a_w3 az_overall_b_w3 az_overall_c_w3 az_overall_d_w3 az_overall_e_w3, gen(az_overall_all_w3)
	
		
	
	
*** Prepare data: generate heterogeneity cuts

	// demographics (category F)
	gen h_gender = (gender == 1) 							if gender != .
	gen h_birth_year = (birth_year < 1996) 					if birth_year != .
	gen h_ethnicity_han = (ethnicity_han == 1) 				if ethnicity_han != .
	gen h_birthplace_coastal = (birthplace_coastal == 1) 	if birthplace_coastal != .
	gen h_residence_coastal = (residence_coastal == 1) 		if residence_coastal != .
	gen h_hukou_urban = (hukou_urban == 1) 					if hukou_urban != .
	gen h_religion_religious = (religion_religious == 1) 	if religion_religious != .
	gen h_ccp_member = (ccp_member == 1) 					if ccp_member != .
	gen h_university_elite = (university_elite == 1) 		if university_elite != .
	gen h_hs_track_science = (hs_track_science == 1) 		if hs_track_science != .
	gen h_department_ssh = (department_ssh == 1) 			if department_ssh != .
	gen h_domestic_english_atleast4 = (domestic_english_atleast4 == 1) 	if domestic_english_atleast4 != .
	gen h_foreign_english_yes = (foreign_english_yes == 1) 	if foreign_english_yes != .
	gen h_travel_hktaiwan = (travel_hktaiwan == 1) 			if travel_hktaiwan != .
	gen h_travel_foreign_yes = (travel_foreign_yes == 1) 	if travel_foreign_yes != .
	gen h_siblings_total = (siblings_total > 0) 			if siblings_total != .
	gen h_father_edu_hsabove = (father_edu_hsabove == 1) 	if father_edu_hsabove != .
	gen h_work_father_govt = (work_father_govt == 1) 		if work_father_govt != .
	gen h_father_ccp = (father_ccp == 1) 					if father_ccp != .
	gen h_mother_edu_hsabove = (mother_edu_hsabove == 1) 	if mother_edu_hsabove != .
	gen h_work_mother_govt = (work_mother_govt == 1) 		if work_mother_govt != .
	gen h_mother_ccp = (mother_ccp == 1) 					if mother_ccp != .
	gen h_hh_income = (hh_income >= 75000) 					if hh_income != .
	gen h_az_preference_risk = (az_preference_risk > 0) 	if az_preference_risk != .
	gen h_az_preference_time = (az_preference_time > 0) 	if az_preference_risk != .
	gen h_az_preference_altruism = (az_preference_altruism > 0) 		if az_preference_altruism != .
	gen h_az_preference_reciprocity = (az_preference_reciprocity > 0) 	if az_preference_reciprocity != .
	
	// baseline outcomes
	gen h_az_overall_a_w1 = (az_overall_a_w1 > 0) 			if az_overall_a_w1 != .
	gen h_az_overall_b_w1 = (az_overall_b_w1 > 0) 			if az_overall_b_w1 != .
	gen h_az_overall_c_w1 = (az_overall_c_w1 > 0) 			if az_overall_c_w1 != .
	gen h_az_overall_d_w1 = (az_overall_d_w1 > 0) 			if az_overall_d_w1 != .
	gen h_az_overall_e_w1 = (az_overall_e_w1 > 0) 			if az_overall_e_w1 != .
	gen h_az_belief_media_value_w1 = (az_belief_media_value_w1 > 0) 		if az_belief_media_value_w1 != .
	gen h_az_belief_media_trust_w1 = (az_belief_media_trust_w1 > 0) 		if az_belief_media_trust_w1 != .
	gen h_az_knowledge_news_cens_w1 = (az_knowledge_news_censored_w1 > 0) 	if az_knowledge_news_censored_w1 != .
	gen h_az_knowledge_news_unce_w1 = (az_knowledge_news_uncensor_w1 > 0) 	if az_knowledge_news_uncensor_w1 != .
	gen h_az_knowledge_pp_censor_w1 = (az_knowledge_people_censor_w1 > 0) 	if az_knowledge_people_censor_w1 != .
	gen h_az_knowledge_pp_uncens_w1 = (az_knowledge_people_uncens_w1 > 0) 	if az_knowledge_people_uncens_w1 != .
	gen h_az_knowledge_pr_china_w1 = (az_knowledge_protest_china_w1 > 0) 	if az_knowledge_protest_china_w1 != .
	gen h_az_belief_trust_govt_w1 = (az_belief_trust_govt_w1 > 0) 			if az_belief_trust_govt_w1 != .
		
		
	
	
*** Prepare data: generate variable labels

	// A. Beleifs and attitudes regarding media
	local l_wtp_vpn 						`"WTP for uncensored Internet access ($/month)"'
	local l_added_value_foreign_media 		`"Value added of foreign media access"'
	local l_az_belief_media_value 			`"Valuation of access to foreign media outlets"'
	
	local l_trust_media_dom_state 			`"Distrust in domestic state-owned media"'
	local l_trust_media_dom_private 		`"Distrust in domestic privately-owned media"'
	local l_trust_media_foreign 			`"Trust in foreign media"'
	local l_az_belief_media_trust 			`"Trust in non-domestic media outlets"'
	
	local l_bias_domestic 					`"Degree of censorship on domestic news outlets"'
	local l_bias_foreign 					`"Degree of censorship on foreign news outlets"'
	local l_bias_dom_govt_policy_t1 		`"Domestic cens. driven by govt. policies"'
	local l_bias_dom_firm_interest_t1 		`"Domestic cens. driven by corp. interest"'
	local l_bias_dom_media_pref_t1 			`"Domestic cens. driven by media's ideology"'
	local l_bias_dom_reader_demand_t1		`"Domestic cens. driven by readers' demand"'
	local l_bias_for_govt_policy_t1 		`"Foreign cens. driven by govt. policies"'
	local l_bias_for_firm_interest_t1 		`"Foreign cens. driven by corp. interest"'
	local l_bias_for_media_pref_t1 			`"Foreign cens. driven by media's ideology"'
	local l_bias_for_reader_demand_t1 		`"Foreign cens. driven by readers' demand"'

	local l_distneutral_cn_neg_cn 			`"Bias: Chinese media on neg. news in China"'
	local l_distneutral_cn_pos_cn 			`"Bias: Chinese media on pos. news in China"'
	local l_distneutral_us_neg_cn 			`"Bias: Chinese media on neg. news in US"'
	local l_distneutral_us_pos_cn 			`"Bias: Chinese media on pos. news in US"'
	local l_distneutral_cn_neg_us 			`"Bias: US media on neg. news in China"'
	local l_distneutral_cn_pos_us 			`"Bias: US media on pos. news in China"'
	local l_distneutral_us_neg_us 			`"Bias: US media on neg. news in US"'
	local l_distneutral_us_pos_us 			`"Bias: US media on pos. news in US"'

	local l_bias_cn_neg_cn_cens 			`"Censorship: Chinese media on neg. news in China"'
	local l_bias_cn_pos_cn_cens 			`"Censorship: Chinese media on pos. news in China"'
	local l_bias_us_neg_cn_cens 			`"Censorship: Chinese media on neg. news in US"'
	local l_bias_us_pos_cn_cens 			`"Censorship: Chinese media on pos. news in US"'
	local l_bias_cn_neg_us_cens 			`"Censorship: US media on neg. news in China"'
	local l_bias_cn_pos_us_cens 			`"Censorship: US media on pos. news in China"'
	local l_bias_us_neg_us_cens 			`"Censorship: US media on neg. news in US"'
	local l_bias_us_pos_us_cens 			`"Censorship: US media on pos. news in US"'
	
	local l_censor_just_dom_economic 		`"Unjustified: censoring economic news"'
	local l_censor_just_dom_political 		`"Unjustified: censoring political news"'
	local l_censor_just_dom_social 			`"Unjustified: censoring social news"'
	local l_censor_just_for 				`"Unjustified: censoring foreign news"'
	local l_censor_just_porn 				`"Unjustified: censoring pornography"'
	local l_az_belief_media_justif_pf 		`"Censorship is unjustified"'
	
	
	// B. Knowledge
	local l_panamapapers 					`"Panama Papers"'
	local l_tenyearshk  					`"HK independence"'
	local l_stockcrash  					`"2016 stock mkt crash"'
	local l_economistcensor 				`"Censoring Economist"'
	local l_coalprod 						`"Steel prod. & pollution"'
	local l_trumpchina 						`"Trump trademark in China"'
	local l_xiaojianhua 					`"Jianhua Xiao kidnap"'
	local l_xijiangcar 						`"Tracking Xinjiang cars"'
	local l_chinanorway 					`"China-Norway relations"'
	local l_womenrights 					`"Women rights activisits"'
	local l_hkceelection 					`"HK CE election"'
	local l_news_perccor_cen 				`"% quizzes answered correctly: sensitive"'
	local l_news_perccor_unc 				`"% quizzes answered correctly: non-sensitive"'

	local l_people_puzhiqiang 				`"Aware of Zhiqiang Pu"'
	local l_people_lizehou 					`"Aware of Zehou Li"'
	local l_people_huangzhifeng 			`"Aware of Joshua Wong"'
	local l_people_chenguangcheng 			`"Aware of Guangcheng Cheng"'
	local l_people_lixiaolin 				`"Aware of Xiaolin Li"'
	local l_people_renzhiqiang 				`"Aware of Zhiqiang Ren"'
	local l_people_maoyushi 				`"Aware of Yushi Mao"'
	local l_people_honghuang 				`"Aware of Huang Hong"'
	local l_people_liuqiangdong 			`"Aware of Qiangdong Liu"'
	local l_people_jialequn 				`"Aware of Lequn Jia"'
	
	local l_protest_2014_europe_square 		`"Aware of 2014 Ukrainian Euromaidan Revolution"'
	local l_protest_2014_sun_flower 		`"Aware of 2014 Taiwan Sunflower Stud. Movement"'
	local l_protest_2010_arabic_spring 		`"Aware of 2010 Arab Spring"'
	local l_protest_2014_crimea_vote 		`"Aware of 2014 Crimean Status Referendum"'
	local l_protest_2012_hk_curriculum 		`"Aware of 2012 HK Anti-National Curr. Movement"'
	local l_protest_2010_catal_indep 		`"Aware of 2010 Catalonian Indep. Movement"'
	local l_protest_2014_umbrella 			`"Aware of 2014 HK Umbrella Revolution"'
	local l_protest_2011_tmrw_parade 		`"Aware of 2011 Tomorrow Revolution"'
	local l_protest_2016_mongkok_riot 		`"Aware of 2016 HK Mong Kok Revolution"'
	local l_protest_pcheard_china 			`"Awareness of protests in Greater China"'
	local l_protest_pcheard_foreign			`"Awareness of foreign protests"'
	
	local l_familiar_china_issues_self 		`"Informedness of issues in China"'
	local l_familiar_china_others 			`"Greater informedness than peers"'
	local l_az_knowledge_meta 				`"Self-assessment of knowledge level"'
	
	
	// C. Economic beliefs
	local l_guess_gdp_growth_china 			`"Guess on GDP growth rate in China"'
	local l_guess_stock_index_sh 			`"Guess on year-end SSCI"'
	local l_az_belief_econ_perf_cn 			`"Optimistic belief of Chinese economy"'
	local l_az_belief_econ_conf_cn 			`"Confidence of guesses on Chinese economy"'
	
	local l_guess_gdp_growth_china_con 		`"Confidence of China GDP guess"'
	local l_guess_stock_index_sh_con 		`"Confidence of SSCI guess"'

	local l_guess_gdp_growth_us 			`"Guess on GDP growth rate in US"'
	local l_guess_stock_index_dj			`"Guess on year-end DJI"'

	local l_guess_gdp_growth_us_con 		`"Confidence of US GDP guess"'
	local l_guess_stock_index_dj_con 		`"Confidence of DJI guess"'

	
	// D. Political attitudes
	local l_inst_change_econ 				`"Economic institution needs changes"'
	local l_inst_change_poli 				`"Political institution needs changes"'
	local l_az_belief_instchange 			`"Demand for institutional change"'
	
	local l_trust_central_govt 				`"Trust in central govt. of China"'
	local l_trust_provincial_govt 			`"Trust in provincial govt. of China"'
	local l_trust_local_govt 				`"Trust in local govt. of China"'
	local l_az_belief_trust_govt 			`"Trust in Chinese govt."'
	local l_trust_japan_govt 				`"Trust in central govt. of Japan"'
	local l_trust_us_govt 					`"Trust in federal govt. of US"'
	local l_az_belief_trust_foreign 		`"Trust in foreign govt."'
	local l_trust_court 					`"Trust in court"'
	local l_trust_police 					`"Trust in police"'
	local l_trust_ngo 						`"Trust in NGOs"'
	local l_trust_financial_domestic 		`"Trust in domestic financial inst."'
	local l_trust_financial_foreign 		`"Trust in foreign financial inst."'
	
	local l_eval_govt_economic 				`"Satisfaction of economic dev."'
	local l_eval_govt_dom_politics 			`"Satisfaction of domestic politics"'
	local l_eval_govt_for_relations 		`"Satisfaction of diplomatic affairs"'
	local l_az_belief_evalgovt 				`"Satisfaction of govt. performance"'
	
	local l_revalgovt_election 				`"Eval. importance: universal suffrage"'
	local l_revalgovt_economy 				`"Eval. importance: economic dev."'
	local l_revalgovt_equality 				`"Eval. importance: income and wealth equality"'
	local l_revalgovt_ruleoflaw 			`"Eval. importance: rule of law"'
	local l_revalgovt_human_rights 			`"Eval. importance: civil and human rights"'
	local l_revalgovt_freedom_speech 		`"Eval. importance: freedom of speech"'
	local l_revalgovt_global_power 			`"Eval. importance: intl. affairs"'
	local l_revalgovt_fair_history 			`"Eval. importance: handle history fairly"'

	local l_severity_welfare 				`"Severity: social security and welfare"'
	local l_severity_employment 			`"Severity: employments"'
	local l_severity_pollution 				`"Severity: environmental pollution"'
	local l_severity_inequality 			`"Severity: wealth inequality"'
	local l_severity_corruption 			`"Severity: govt. corruption"'
	local l_severity_dscrm_minority 		`"Severity: minority discrimination"'

	local l_importance_live_in_demo 		`"Living in democracy is not important"'
	local l_china_interest_group 			`"China cares interest for masses"'
	local l_china_rate_democracy 			`"Level of democracy in China"'
	local l_china_rate_humanrights 			`"Level of human rights protection"'
	local l_az_belief_democracy_fp 			`"Living in democracy is not important"'
	
	local l_justify_minority_policy 		`"Justified: minority policies"'
	local l_justify_reduce_pollution 		`"Justified: prod. cut to reduce pollution"'
	local l_justify_hukou 					`"Justified: migration restrictions"'
	local l_justify_one_child 				`"Justified: one-child policy"'
	local l_justify_gaokao 					`"Justified: college admission policies"'
	local l_justify_hongkong_policy 		`"Justified: policy towards HK"'
	local l_justify_taiwan_policy 			`"Justified: policy towards Taiwan"'
	local l_justify_transgene 				`"Justified: transgenetic food"'
	local l_justify_receive_refugee 		`"Justified: refusal of DPRK refugees"'
	local l_justify_soe_privatize 			`"Justified: privatization of SOEs"'
	local l_justify_homo_marriage 			`"Justified: legal. of homosexual marriages"'
	local l_justify_legal_prostitute 		`"Justified: legal. of prostitution"'
	local l_justify_abortion 				`"Justified: abortion"'
	local l_justify_soft_drug 				`"Justified: soft drugs usage"'
	local l_justify_violence_stability 		`"Justified: govt. use of violence"'

	local l_willing_against_illi_govt 		`"Willing to battle illegal govt. acts"'
	local l_willing_report_mis 				`"Willing to report govt. misconduct"'
	local l_willing_protect_weak 			`"Willing to stand up for the weak"'
	local l_az_belief_willing 				`"Willingness to act"'
	
	local l_interest_economic 				`"Interest in economics"'
	local l_interest_politics 				`"Interest in politics"'

	local l_proud_being_chinese 			`"Proud of being Chinese"'

	local l_fear_critic_govt_self 			`"Fear to criticize govt. in public"'
	
	
	// E. Behaviors
	local l_info_domestic_website 			`"Important info source: domestic websites"'
	local l_info_foreign_website 			`"Important info source: foreign websites"'
	local l_info_social_media_dom 			`"Important info source: domestic social media"'
	local l_info_word_of_mouth 				`"Important info source: word of mouth"'
	local l_info_social_media_for 			`"Important info source: foreign social media"'

	local l_info_freq_website_for 			`"Frequently visit foreign websites for info"'

	local l_participate_social_protest 		`"Participated in social protests"'
	local l_participate_plan_vote 			`"Plan to vote in PCR election"'
	local l_participate_complain_school 	`"Filed complaints to school"'
	local l_participate_ngo 				`"Participated in NGO activities"'
	local l_az_var_polparticipation 		`"Participation behaviors"'
	
	local l_frequency_talk_politic 			`"Frequency of discussing poli. with friends"'
	local l_frequency_persuade_friends 		`"Frequency of persuading others"'
	local l_az_var_socialinteract 			`"Social interaction in politics"'
	
	local l_plan_grad_gradschool_dom 		`"Plan: grad. school in China"'
	local l_plan_grad_foreignmaster 		`"Plan: master degree abroad"'
	local l_plan_grad_foreignphd 			`"Plan: PhD degree abroad"'
	local l_plan_grad_military 				`"Plan: military in China"'
	local l_plan_grad_work 					`"Plan: work right away"'

	local l_cp_t1_national_civil 			`"Sector pref.: national civil service"'
	local l_cp_t1_local_civil 				`"Sector pref.: local civil service"'
	local l_cp_t1_military 					`"Sector pref.: military"'
	local l_cp_t1_chinese_private 			`"Sector pref.: private firm in China"'
	local l_cp_t1_for_firm 					`"Sector pref.: private firm in China"'
	local l_cp_t1_soe 						`"Sector pref.: SOEs"'
	local l_cp_t1_institutional 			`"Sector pref.: inst. organizations"'
	local l_cp_t1_entrepreneur 				`"Sector pref.: entrepreneurship"'

	local l_cloc_beijing 					`"Location pref.: Beijing"'
	local l_cloc_shanghai 					`"Location pref.: Shanghai"'
	local l_cloc_gzsz 						`"Location pref.: Shenzhen/Guangzhou"'
	local l_cloc_tjcq 						`"Location pref.: tier 2 cities in central"'
	local l_cloc_hkmc 						`"Location pref.: HK and Macau"'
	local l_cloc_taiwan 					`"Location pref.: Taiwan"'
	local l_cloc_dom 						`"Location pref.: other cities in China"'
	local l_cloc_for 						`"Location pref.: foreign cities"'
	local l_stock_participation 			`"Currently invested in Chinese stock mkt."'

* ============================================================================ *
* reproduce Table 3 results partially
* ============================================================================ *

use "/Users/peizan/Library/CloudStorage/Dropbox/sens_lee_bound/data/chen_and_yang_2019aer_replication/chen_yang_2019aer_clean.dta", clear

* reproduce reduced form estimation in Panel A
foreach Y in a b c d e {
		
	// Panel A: reduced form
	reg az_overall_`Y'_w3 treatment_vpnonly treatment_nlonly treatment_vpnnl if treatment_user != 1, r
		
}

* ============================================================================ *
* keep only AE (Encouragement + Access) group
* ============================================================================ *

* drop units that are not used in regression sample
keep if treatment_user != 1

* keep only AE and control groups
* treatment_control: = 1 if unit is in control group; = 0 if unit is in one of three treatment groups
keep if treatment_vpnnl == 1 | treatment_control == 1

* double check: can reproduce reduced form estimation of AE effect in Panel A
foreach Y in a b c d e {
		
	// Panel A: reduced form
	reg az_overall_`Y'_w3 treatment_vpnnl, r
		
}

* ============================================================================ *
* use leebounds package to calculate estimated bounds and CI
* ============================================================================ *

foreach Y in a b c d e {
		
	// Panel A: reduced form
	leebounds az_overall_`Y'_w3 treatment_vpnnl, cie	
	
}
		
* ============================================================================ *
* check missingness pattern
* rmk: for all five outcome variables, there are more missing in the control
* ============================================================================ *

foreach Y in a b c d e {
	
  // remove old missingness indicator if it already exists
  capture drop az_overall_`Y'_w3_miss
  
  // generate missingness indicator
  gen az_overall_`Y'_w3_miss = 1 if az_overall_`Y'_w3 == .
  replace az_overall_`Y'_w3_miss = 0 if az_overall_`Y'_w3 != .
  
  // sum stats on missing indicator
  tabstat az_overall_`Y'_w3_miss, by(treatment_vpnnl) statistics(n mean sd min max)
  ttest az_overall_`Y'_w3_miss, by(treatment_vpnnl) unequal
  
  
}

* ============================================================================ *
* save data
* ============================================================================ *

save "/Users/peizan/github/riattrition/data-raw/chen_and_yang_2019/chen_yang_2019aer_clean.dta", replace
