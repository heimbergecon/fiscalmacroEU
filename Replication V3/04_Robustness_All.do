********************************************************************************
* Replication Code: 
* Macroeconomic and distributional effects of fiscal consolidation measures in EU countries
* Philipp Heimberger & Anna Matzner
* July 2026
* Figure 2, A1, A2, A3
********************************************************************************

sort country_id year
xtset country_id year, yearly

// IVREG without dependent ************************************************

* specs for local projection
local horizon = 5    // Impulse horizon
local estdiff = 2     // 0: level, 1: differences, 2: cumulative

local CI1 = 0.10      // Confidence level 1
local CI2 = 0.32	  // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock diff_STRUCBAL
local vars log_REALGDP STRUCBAL unemp log_PRICES gini_disp INFL
local varsnames " "Real GDP" "Structural balance" "Unemployment" "Prices"  "Income Inequality" "Inflation rate" " 
local labels " "%" "%-points" "%-points" "%" "index points" "%-points" "

local cntrls_log_REALGDP L(1/1).RYIELD L(1/1).OGAP L(1/1).REER
local cntrls_STRUCBAL  L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER
local cntrls_unemp  L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER
local cntrls_log_PRICES  L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER
local cntrls_gini_disp  L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER 
local cntrls_INFL L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER 

* options
global savefigs 1    
global verb qui      

cap gen t = _n 
cap gen h = t - 1 // h is the horizon for the irfs 

local ivar = 1
foreach var in `vars' { 
	
	* set controls
	local cntrls `cntrls_`var''
		
	* preallocate 
	qui gen biv`var' = .
	qui gen up90biv`var' = .
	qui gen lo90biv`var' = .
	qui gen up68biv`var' = .
	qui gen lo68biv`var' = .

	if `estdiff' > 0 {
		cap g d`var' = `var' - L.`var'
	}

	forvalues i = 0/`horizon' {

		local iname `i'
		cap gen d`iname'`var' = F`i'.`var' - L.`var' 
		
		if `estdiff' == 0 {
			`verb' ivreg2 F(`i').`var' (diff_STRUCBAL = TOTAL)   `cntrls' i.year i.country_id, dkraay(1)
		}
		else if `estdiff' == 1 {
			`verb' ivreg2 F(`i').d`var' (diff_STRUCBAL = TOTAL)  `cntrls' i.year i.country_id,  dkraay(1)
		}
		else if `estdiff' == 2 {
			`verb' ivreg2 d`iname'`var' (diff_STRUCBAL = TOTAL)   `cntrls'  i.year i.country_id,  dkraay(1)
		}
		
		cap gen biv`var'h`iname' = _b[`shock']
		cap gen seiv`var'h`iname' = _se[`shock']

		qui replace biv`var' = biv`var'h`iname' if h==`i'
		qui replace up90biv`var' = biv`var'h`iname' + `z1'*seiv`var'h`iname' if h==`i'
		qui replace lo90biv`var' = biv`var'h`iname' - `z1'*seiv`var'h`iname' if h==`i'
		qui replace up68biv`var' = biv`var'h`iname' + `z2'*seiv`var'h`iname' if h==`i'
		qui replace lo68biv`var' = biv`var'h`iname' - `z2'*seiv`var'h`iname' if h==`i'
		
	}

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	local labname : word `ivar' of `labels'
	tw (rarea up90biv`var' lo90biv`var' h, fcolor("$mblue%15") lcolor(mdred) lw(none) lpattern(solid)) ///
		(rarea up68biv`var' lo68biv`var' h, fcolor("$mblue%40") lcolor(mdred) lw(none) lpattern(solid)) ///
		(line biv`var' h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		title("`varname'", size(10) col(black) margin(b=2)) xtitle("Years", size(7)) ///
		ytitle("`labname'", margin(r=3) size(7))   ///
		xlabel(0(1)`horizon', labsize(7)) ///
		ylabel(,labsize(7)) /// 
		plotregion(color(white)) ///
		graphregion(color(white) )  legend(off) name("irfgs_`var'_ts", replace) ///
				saving("$FIGUREDIR/pirf_iv_wolagDep_`var'.gph", replace)

		graph export "$FIGUREDIR\pirf_iv_wolagDep_`var'.png", replace

			local ivar = `ivar'+1

		foreach stat in b se up90b lo90b up68b lo68b {
		capture drop `stat'
		capture drop `stat'*
		}
		}


// IVREG without output controls ***********************************************************************

* specs for local projection
local horizon = 5    // Impulse horizon
local estdiff = 2     // 0: level, 1: differences, 2: cumulative

local CI1 = 0.10      // Confidence level 1
local CI2 = 0.32	  // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock diff_STRUCBAL
local vars log_REALGDP STRUCBAL unemp log_PRICES gini_disp INFL
local varsnames " "Real GDP" "Structural balance" "Unemployment" "Prices"  "Income Inequality" "Inflation rate" " 
local labels " "%" "%-points" "%-points" "%" "index points" "%-points" "

local cntrls_log_REALGDP L(1/1).RYIELD L(1/1).REER
local cntrls_STRUCBAL L(1/1).RYIELD  L(1/1).REER
local cntrls_unemp  L(1/1).RYIELD L(1/1).REER L(1/1).unemp
local cntrls_log_PRICES L(1/1).RYIELD  L(1/1).REER L(1/1).diff_log_PRICES
local cntrls_gini_disp L(1/1).RYIELD L(1/1).REER L(1/1).gini_disp
local cntrls_INFL L(1/1).RYIELD  L(1/1).REER L(1/1).INFL

* options
global savefigs 1   
global verb qui     

cap gen t = _n 
cap gen h = t - 1 // h is the horizon for the irfs 

local ivar = 1
foreach var in `vars' { 
	
	* set controls
	local cntrls `cntrls_`var''
		
	* preallocate 
	qui gen biv`var' = .
	qui gen up90biv`var' = .
	qui gen lo90biv`var' = .
	qui gen up68biv`var' = .
	qui gen lo68biv`var' = .

	if `estdiff' > 0 {
		cap g d`var' = `var' - L.`var'
	}

	forvalues i = 0/`horizon' {

		local iname `i'
		cap gen d`iname'`var' = F`i'.`var' - L.`var' 
		
		if `estdiff' == 0 {
			`verb' ivreg2 F(`i').`var' (diff_STRUCBAL = TOTAL)   `cntrls' i.year i.country_id, dkraay(1)
		}
		else if `estdiff' == 1 {
			`verb' ivreg2 F(`i').d`var' (diff_STRUCBAL = TOTAL)  `cntrls' i.year i.country_id,  dkraay(1)
		}
		else if `estdiff' == 2 {
			`verb' ivreg2 d`iname'`var' (diff_STRUCBAL = TOTAL)   `cntrls'  i.year i.country_id,  dkraay(1)
		}
		
		cap gen biv`var'h`iname' = _b[`shock']
		cap gen seiv`var'h`iname' = _se[`shock']

		qui replace biv`var' = biv`var'h`iname' if h==`i'
		qui replace up90biv`var' = biv`var'h`iname' + `z1'*seiv`var'h`iname' if h==`i'
		qui replace lo90biv`var' = biv`var'h`iname' - `z1'*seiv`var'h`iname' if h==`i'
		qui replace up68biv`var' = biv`var'h`iname' + `z2'*seiv`var'h`iname' if h==`i'
		qui replace lo68biv`var' = biv`var'h`iname' - `z2'*seiv`var'h`iname' if h==`i'
		
	}

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	local labname : word `ivar' of `labels'
	tw (rarea up90biv`var' lo90biv`var' h, fcolor("$mblue%15") lcolor(mdred) lw(none) lpattern(solid)) ///
		(rarea up68biv`var' lo68biv`var' h, fcolor("$mblue%40") lcolor(mdred) lw(none) lpattern(solid)) ///
		(line biv`var' h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		title("`varname'", size(10) col(black) margin(b=2)) xtitle("Years", size(7)) ///
		ytitle("`labname'", margin(r=3) size(7))   ///
		xlabel(0(1)`horizon', labsize(7)) ///
		ylabel(,labsize(7)) /// 
		plotregion(color(white)) ///
		graphregion(color(white) )  legend(off) name("irfgs_`var'_ts", replace) ///
				saving("$FIGUREDIR/pirf_iv_wocontrols_`var'.gph", replace)

		graph export "$FIGUREDIR\pirf_iv_wocontrols_`var'.png",   replace

			local ivar = `ivar'+1

		foreach stat in b se up90b lo90b up68b lo68b {
		capture drop `stat'
		capture drop `stat'*
		}
		}


// IVREG Robustness - binary narrative variable *******************************

* specs for local projection
local horizon = 5    // Impulse horizon
local estdiff = 2     // 0: level, 1: differences, 2: cumulative

local CI1 = 0.10      // Confidence level 1
local CI2 = 0.32	  // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock diff_STRUCBAL
local vars log_REALGDP STRUCBAL unemp log_PRICES gini_disp INFL
local varsnames " "Real GDP" "Structural balance" "Unemployment" "Prices"  "Income Inequality" "Inflation rate"  " 
local labels " "%" "%-points" "%-points" "%" "index points" "%-points" "

local cntrls_log_REALGDP L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER
local cntrls_STRUCBAL L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER
local cntrls_unemp L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).unemp
local cntrls_log_PRICES L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).diff_log_PRICES
local cntrls_gini_disp L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).gini_disp
local cntrls_INFL L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).INFL

* options
global savefigs 1    
global verb qui      

cap gen t = _n 
cap gen h = t - 1 // h is the horizon for the irfs 

local ivar = 1
foreach var in `vars' { 
	
	* set controls
	local cntrls `cntrls_`var''
		
	* preallocate 
	qui gen biv`var' = .
	qui gen up90biv`var' = .
	qui gen lo90biv`var' = .
	qui gen up68biv`var' = .
	qui gen lo68biv`var' = .

	if `estdiff' > 0 {
		cap g d`var' = `var' - L.`var'
	}

	forvalues i = 0/`horizon' {

		local iname `i'
		cap gen d`iname'`var' = F`i'.`var' - L.`var' 
		
		if `estdiff' == 0 {
			`verb' ivreg2 F(`i').`var' (diff_STRUCBAL = TOTAL_binary)   `cntrls' i.year i.country_id, dkraay(1) partial(i.year i.country_id)
		}
		else if `estdiff' == 1 {
			`verb' ivreg2 F(`i').d`var' (diff_STRUCBAL = TOTAL_binary)  `cntrls' i.year i.country_id,  dkraay(1) partial(i.year i.country_id)
		}
		else if `estdiff' == 2 {
			`verb' ivreg2 d`iname'`var' (diff_STRUCBAL = TOTAL_binary)   `cntrls'  i.year i.country_id,  dkraay(1) partial(i.year i.country_id)
		}
		
		cap gen biv`var'h`iname' = _b[`shock']
		cap gen seiv`var'h`iname' = _se[`shock']

		qui replace biv`var' = biv`var'h`iname' if h==`i'
		qui replace up90biv`var' = biv`var'h`iname' + `z1'*seiv`var'h`iname' if h==`i'
		qui replace lo90biv`var' = biv`var'h`iname' - `z1'*seiv`var'h`iname' if h==`i'
		qui replace up68biv`var' = biv`var'h`iname' + `z2'*seiv`var'h`iname' if h==`i'
		qui replace lo68biv`var' = biv`var'h`iname' - `z2'*seiv`var'h`iname' if h==`i'
		
	}

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	local labname : word `ivar' of `labels'
	tw (rarea up90biv`var' lo90biv`var' h, fcolor("$mblue%15") lwidth(none) lpattern(solid)) ///
		(rarea up68biv`var' lo68biv`var' h, fcolor("$mblue%40") lwidth(none) lpattern(solid)) ///
		(line biv`var' h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		title("`varname'", size(10) col(black) margin(b=2)) xtitle("Years", size(7)) ///
		ytitle("`labname'", margin(r=3) size(7))   ///
		xlabel(0(1)`horizon', labsize(7)) ///
		ylabel(,labsize(7)) /// 
		plotregion(color(white)) ///
		graphregion(color(white) )  legend(off) ///
				saving("$FIGUREDIR/pirf_iv_binary_`var'.gph", replace)

		graph export "$FIGUREDIR\pirf_iv_binary_`var'.png",   replace

			local ivar = `ivar'+1

		foreach stat in b se up90b lo90b up68b lo68b {
		capture drop `stat'
		capture drop `stat'*
		}
		}


		
// IVREG without Ireland ********************************************************

drop if ccode == "IRE"

* specs for local projection
local horizon = 5    // Impulse horizon
local estdiff = 2     // 0: level, 1: differences, 2: cumulative

local CI1 = 0.10      // Confidence level 1
local CI2 = 0.32	  // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock diff_STRUCBAL
local vars log_REALGDP STRUCBAL unemp log_PRICES gini_disp INFL
local varsnames " "Real GDP" "Structural balance" "Unemployment" "Prices"  "Income Inequality" "Inflation rate"  " 
local labels " "%" "%-points" "%-points" "%" "index points" "%-points" "

local cntrls_log_REALGDP L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER
local cntrls_STRUCBAL L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER
local cntrls_unemp L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).unemp
local cntrls_log_PRICES L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).diff_log_PRICES
local cntrls_gini_disp L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).gini_disp

* options
global savefigs 1    
global verb qui     

cap gen t = _n 
cap gen h = t - 1 // h is the horizon for the irfs 


local ivar = 1
foreach var in `vars' { 
	
	* set controls
	local cntrls `cntrls_`var''
		
	* preallocate 
	qui gen biv`var' = .
	qui gen up90biv`var' = .
	qui gen lo90biv`var' = .
	qui gen up68biv`var' = .
	qui gen lo68biv`var' = .

	if `estdiff' > 0 {
		cap g d`var' = `var' - L.`var'
	}

	forvalues i = 0/`horizon' {

		local iname `i'
		cap gen d`iname'`var' = F`i'.`var' - L.`var' 
		
		if `estdiff' == 0 {
			`verb' ivreg2 F(`i').`var' (diff_STRUCBAL = TOTAL)   `cntrls' i.year i.country_id, dkraay(1) 
		}
		else if `estdiff' == 1 {
			`verb' ivreg2 F(`i').d`var' (diff_STRUCBAL = TOTAL)  `cntrls' i.year i.country_id,  dkraay(1) 
		}
		else if `estdiff' == 2 {
			`verb' ivreg2 d`iname'`var' (diff_STRUCBAL = TOTAL)   `cntrls'  i.year i.country_id,  dkraay(1) partial(i.year i.country_id)
			}
		
		cap gen biv`var'h`iname' = _b[`shock']
		cap gen seiv`var'h`iname' = _se[`shock']

		qui replace biv`var' = biv`var'h`iname' if h==`i'
		qui replace up90biv`var' = biv`var'h`iname' + `z1'*seiv`var'h`iname' if h==`i'
		qui replace lo90biv`var' = biv`var'h`iname' - `z1'*seiv`var'h`iname' if h==`i'
		qui replace up68biv`var' = biv`var'h`iname' + `z2'*seiv`var'h`iname' if h==`i'
		qui replace lo68biv`var' = biv`var'h`iname' - `z2'*seiv`var'h`iname' if h==`i'
		
	}

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	local labname : word `ivar' of `labels'
	tw (rarea up90biv`var' lo90biv`var' h, fcolor("$mblue%15") lwidth(none) lpattern(solid)) ///
		(rarea up68biv`var' lo68biv`var' h, fcolor("$mblue%40") lwidth(none) lpattern(solid)) ///
		(line biv`var' h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		title("`varname'", size(10) col(black) margin(b=2)) xtitle("Years", size(7)) ///
		ytitle("`labname'", margin(r=3) size(7))   ///
		xlabel(0(1)`horizon', labsize(7)) ///
		ylabel(,labsize(7)) /// 
		plotregion(color(white)) ///
		graphregion(color(white) )  legend(off) ///
				saving("$FIGUREDIR/pirf_iv_NOIRE_`var'.gph", replace)

		graph export "$FIGUREDIR\pirf_iv_NOIRE_`var'.png", replace

			local ivar = `ivar'+1

		foreach stat in b se up90b lo90b up68b lo68b {
		capture drop `stat'
		capture drop `stat'*
		}
		}
		

// JOIN PLOTS
local vars aipw_restricted iv_binary iv_wocontrols iv_NOIRE iv_wolagDep

foreach var in `vars' { 

graph use "$FIGUREDIR\pirf_`var'_log_REALGDP.gph", name(g1, replace)	
graph use "$FIGUREDIR\pirf_`var'_unemp.gph", name(g2, replace)	
graph use "$FIGUREDIR\pirf_`var'_INFL.gph", name(g3, replace)	
graph use "$FIGUREDIR\pirf_`var'_gini_disp.gph", name(g4, replace)	


grc1leg2 g1 g2 g3 g4, cols(4) ysize(3) xsize(12)  loff

graph export "$FIGUREDIRC\pirf_iv_robustness`var'.jpg", replace

}

