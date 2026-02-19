********************************************************************************
* Replication Code: 
* Macroeconomic and distributional effects of fiscal consolidation measures in EU countries
* Philipp Heimberger & Anna Matzner
* February 2026
* Figure 1
********************************************************************************

* packages
//ssc install ivreg2


sort country_id year
xtset country_id year, yearly


// OLS ************************************************************************

* specs for local projection
local horizon = 5    // Impulse horizon
local estdiff = 2     // 0: level, 1: differences, 2: cumulative

local CI1 = 0.10      // Confidence level 1
local CI2 = 0.32	  // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock diff_STRUCBAL
local vars log_REALGDP STRUCBAL unemp log_PRICES gini_disp INFL
local varsnames " "Real GDP" "Structural balance" "Unemployment rate" "Prices"  "Income inequality" "Inflation rate" " 
local labels " "%" "%-points" "%-points" "%" "index points" "%-points" "

local cntrls_log_REALGDP L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER
local cntrls_STRUCBAL L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER
local cntrls_unemp L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).unemp
local cntrls_log_PRICES L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).diff_log_PRICES
local cntrls_gini_disp L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).gini_disp

* options
global savefigs 1   
global verb qui     

gen t = _n 
gen h = t - 1 // h is the horizon for the irfs 

local ivar = 1
foreach var in `vars' { 
		
	* set controls
	local cntrls `cntrls_`var''
	
	* preallocate 
	qui gen b`var' = .
	qui gen up90b`var' = .
	qui gen lo90b`var' = .
	qui gen up68b`var' = .
	qui gen lo68b`var' = .

	if `estdiff' > 0 {
		g d`var' = `var' - L.`var'
	}

	forvalues i = 0/`horizon' {

		local iname `i'
		gen d`iname'`var' = F`i'.`var' - L.`var' 
		
		if `estdiff' == 0 {
			`verb' xtscc F(`i').`var' `shock' `cntrls'  i.year, fe 
		}
		else if `estdiff' == 1 { 
			`verb' xtscc F(`i').d`var' `shock'  `cntrls' i.year, fe 
		}
		else if `estdiff' == 2 {
			`verb' xtscc d`iname'`var' `shock' `cntrls' i.year, fe 
		}
		
		gen b`var'h`iname' = _b[`shock']
		gen se`var'h`iname' = _se[`shock']

		qui replace b`var' = b`var'h`iname' if h==`i'
		qui replace up90b`var' = b`var'h`iname' + `z1'*se`var'h`iname' if h==`i'
		qui replace lo90b`var' = b`var'h`iname' - `z1'*se`var'h`iname' if h==`i'
		qui replace up68b`var' = b`var'h`iname' + `z2'*se`var'h`iname' if h==`i'
		qui replace lo68b`var' = b`var'h`iname' - `z2'*se`var'h`iname' if h==`i'
		
	}

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	local labname : word `ivar' of `labels'
	tw (rarea up90b`var' lo90b`var' h, fcolor("$mdred%15") lcolor(mdred) lw(none) lpattern(solid)) ///
		(rarea up68b`var' lo68b`var' h, fcolor("$mdred%40") lcolor(mdred) lw(none) lpattern(solid)) ///
		(line b`var' h, lcolor("$mdred") lpattern(dash) lwidth(thick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		title("{bf:`varname'}", size(medsmall) col(black) margin(b=2)) xtitle("Years", size(medsmall)) ///
		ytitle("`labname'", margin(r=3) size(medsmall))   ///
		xlabel(0(1)`horizon') ///
		plotregion(color(white)) ///
		graphregion(color(white) )  legend(off) name("irfgs_`var'_ts", replace) ///
		saving("$FIGUREDIR/pirf_`var'.gph", replace)
		graph export "$FIGUREDIR\pirf_`var'.pdf", fontface($grfont) replace
		
	local ivar = `ivar'+1
	
		}


// IVREG ***********************************************************************

* specs for local projection
local horizon = 5    // Impulse horizon
local estdiff = 2     // 0: level, 1: differences, 2: cumulative

local CI1 = 0.10      // Confidence level 1
local CI2 = 0.32	  // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock diff_STRUCBAL
local vars log_REALGDP STRUCBAL unemp log_PRICES gini_disp INFL
local varsnames " "Real GDP" "Structural balance" "Unemployment rate" "Prices"  "Income inequality" "Inflation rate" " 
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
			`verb' ivreg2 F(`i').`var' (diff_STRUCBAL = TOTAL)   `cntrls' i.year i.country_id, dkraay(1
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
		title("{bf:Baseline IV results}", size(medsmall) col(black) margin(b=2)) xtitle("Years", size(medsmall)) ///
		ytitle("`labname'", margin(r=3) size(medsmall))   ///
		xlabel(0(1)`horizon') ///
		plotregion(color(white)) ///
		graphregion(color(white) )  legend(off) name("irfgs_`var'_ts", replace) ///
				saving("$FIGUREDIR/pirf_iv_`var'.gph", replace)

		graph export "$FIGUREDIR\pirf_iv_`var'.pdf", fontface($grfont) replace

			local ivar = `ivar'+1

		}

// JOINT PLOT ******************************************************************

* specs 
local horizon = 5    // Impulse horizon

local shock diff_STRUCBAL
local vars log_REALGDP STRUCBAL unemp log_PRICES gini_disp INFL
local varsnames " "Real GDP" "Structural balance" "Unemployment rate" "Prices"  "Income inequality" "Inflation rate" " 
local labels " "%" "%-points" "%-points" "%" "index points" "%-points" "

local ivar = 1

foreach var in `vars' { 

cap g zero = 0

	local varname : word `ivar' of `varsnames'
	local labname : word `ivar' of `labels'

	tw (rarea up90biv`var' lo90biv`var' h, fcolor("$mblue%15") lwidth(none) lpattern(solid)) ///
		(rarea up68biv`var' lo68biv`var' h, fcolor("$mblue%40") lwidth(none) lpattern(solid)) ///
		(line biv`var' h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
		(rarea up90b`var' lo90b`var' h, fcolor("$mdred%15") lwidth(none) lpattern(solid)) ///
		(rarea up68b`var' lo68b`var' h, fcolor("$mdred%40") lwidth(none) lpattern(solid)) ///
		(line b`var' h, lcolor("$mdred") lpattern(dash) lwidth(thick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		title("{bf:`varname'}", color(black) size(medsmall)) xtitle("Years", size(medsmall)) ///
		ytitle("`labname'", size(medsmall))   ///
		xlabel(0(1)`horizon') ///
		plotregion(color(white)) ///
		graphregion(color(white)) ///
		legend(order(3 6) label(3 "IV (narrative)") label(6 "OLS") size(small)) ///
		name("irfgs_`var'_ts", replace) ///
		saving("$FIGUREDIR/pirf_joint_`var'.gph", replace)
		
		graph export "$FIGUREDIR\pirf_joint_`var'.pdf", fontface($grfont) replace

local ivar = `ivar'+1

}


graph use "$FIGUREDIR\pirf_joint_log_REALGDP.gph", name(g1, replace)
graph use "$FIGUREDIR\pirf_joint_unemp.gph", name(g2, replace)
graph use "$FIGUREDIR\pirf_joint_log_PRICES.gph", name(g3, replace)
graph use "$FIGUREDIR\pirf_joint_gini_disp.gph", name(g4, replace)

grc1leg2 g1 g2 g3 g4, cols(2) ysize(6) xsize(7) 
graph export "$FIGUREDIR\iv_vs_ols.png", width(6000) height(4000) replace

// JOINT PLOT ******************************************************************


graph use "$FIGUREDIR\pirf_joint_log_REALGDP.gph", name(g1, replace)
graph use "$FIGUREDIR\pirf_joint_unemp.gph", name(g2, replace)
graph use "$FIGUREDIR\pirf_joint_INFL.gph", name(g3, replace)
graph use "$FIGUREDIR\pirf_joint_gini_disp.gph", name(g4, replace)

grc1leg2 g1 g2 g3 g4, cols(2) ysize(6) xsize(7) 
graph export "$FIGUREDIR\iv_vs_ols_CPI.pdf", replace

*******************************************************************************
foreach stat in b se up90b lo90b up68b lo68b {
		capture drop `stat'
		capture drop `stat'*
		}