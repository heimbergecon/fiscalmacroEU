********************************************************************************
* Replication Code: 
* Macroeconomic and distributional effects of fiscal consolidation measures in EU countries
* Philipp Heimberger & Anna Matzner
* February 2026
* Figure 2, A1, A2, A3: AIPW
* We follow Jorda, Taylor 2016
********************************************************************************

sort country_id year
xtset country_id year, yearly

cap	gen t = _n 
cap	gen h = t - 1 // h is the horizon for the ir

local aipwvars log_REALGDP unemp log_PRICES gini_disp INFL

local ivar = 1
foreach aipwvar in `aipwvars'{

	local varsnames " "Real GDP"  "Unemployment" "Prices"  "Income Inequality" "Inflation rate" " 
	local labels " "%"  "%-points" "%" "index points" "%-points" "

	local cntrls_log_REALGDP L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER
	local cntrls_STRUCBAL L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER
	local cntrls_unemp L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).unemp
	local cntrls_log_PRICES L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).diff_log_PRICES
	local cntrly_gini_disp L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).gini_disp
	local cntrls_INFL L(1/1).RGROWTH L(1/1).RYIELD L(1/1).OGAP L(1/1).REER L(1/1).INFL

	* set controls
	local cntrls `cntrls_`aipwvar''

pause on 
capture drop pihat pihat0
xi: probit TOTAL_binary `cntrls' debtgdp L.TOTAL_binary i.country_id i.year
outreg2 using "$TABLEDIR\probit_results_`aipwvar'.xls", excel replace ///
    se dec(3) label

* raw prscore, not truncated (pihat0)
predict pihat0

* truncate ipws at 10 (pihat)
gen pihat=pihat0
replace pihat = .9 if pihat>.9 & pihat~=.
replace pihat = .1 if pihat<.1 & pihat~=.



* sort again
sort country_id year
xtset country_id year

// AIPW 

capture drop a invwt
gen a = TOTAL_binary
gen invwt=a/pihat0 + (1-a)/(1-pihat0) if pihat~=. // invwt from Lunt et al.

* preallocate LPs
	qui gen b`aipwvar' = .
	qui gen up90b`aipwvar' = .
	qui gen lo90b`aipwvar' = .
	qui gen up68b`aipwvar' = .
	qui gen lo68b`aipwvar' = .

* LP Settings
	local horizon = 5    // Impulse horizon
	local CI1 = 0.10      // Confidence level 1
	local CI2 = 0.32	  // Confidence level 2
	local z1 = abs(invnormal(`CI1'/2))  
	local z2 = abs(invnormal(`CI2'/2))
	
forvalues i = 0/`horizon'{
		
		local iname `i'
		cap gen d`iname'`aipwvar' = F`i'.`aipwvar' - L.`aipwvar'

		reg d`iname'`aipwvar' TOTAL_binary `cntrls' i.country_id i.year [pweight=invwt], cluster(country_id) 
		gen samp=e(sample) // set sample 
		predict mu0 if samp==1 & TOTAL_binary==0 // actual
		predict mu1 if samp==1 & TOTAL_binary==1 // actual
		replace mu0 = mu1 - _b[TOTAL_binary] if samp==1 & TOTAL_binary==1 // ghost
		replace mu1 = mu0 + _b[TOTAL_binary] if samp==1 & TOTAL_binary==0 // ghost
		*from Lunt et al
		generate mdiff1=(-(a-pihat0)*mu1/pihat0)-((a-pihat0)*mu0/(1-pihat0))
		generate iptw=(2*a-1)*d`iname'`aipwvar'*invwt
		generate dr1=iptw+mdiff1
		qui gen ATE_IPWRA=1 // constant for convenience in next reg to get mean
		qui reg dr1 ATE_IPWRA, nocons cluster(country_id)
		eststo DR1_`i'
		
		gen b`aipwvar'h`iname' = _b[ATE_IPWRA]
		
		sum dr1
		local dr1m = r(mean)
		gen Isq  = (dr1-`dr1m')^2
		sum Isq
		estadd scalar RobustSE = sqrt(r(mean)/r(N))
		estadd scalar pvalue = normal(`dr1m'/sqrt(r(mean)/r(N)))	
		
		sum Isq
		gen se`aipwvar'h`iname' = sqrt(r(mean)/r(N))

		qui replace b`aipwvar' = b`aipwvar'h`iname' if h==`i'
		qui replace up90b`aipwvar' = b`aipwvar'h`iname' + `z1'*se`aipwvar'h`iname' if h==`i'
		qui replace lo90b`aipwvar' = b`aipwvar'h`iname' - `z1'*se`aipwvar'h`iname' if h==`i'
		qui replace up68b`aipwvar' = b`aipwvar'h`iname' + `z2'*se`aipwvar'h`iname' if h==`i'
		qui replace lo68b`aipwvar' = b`aipwvar'h`iname' - `z2'*se`aipwvar'h`iname' if h==`i'
	
		capture drop iptw Isq mdiff1 dr1 mu1 mu0 samp ATE_IPWRA
		capture scalar drop dr1m
	}
	
	// Plot LP
	cap g zero = 0
	local labname : word `ivar' of `labels'
	tw (rarea up90b`aipwvar' lo90b`aipwvar' h, fcolor("$mblue%15") lcolor(mdred) lw(none) lpattern(solid)) ///
		(rarea up68b`aipwvar' lo68b`aipwvar' h, fcolor("$mblue%40") lcolor(mdred) lw(none) lpattern(solid)) ///
		(line b`aipwvar' h, lcolor("$mblue") lpattern(dash) lwidth(thick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		title("{bf:AIPW}", size(medsmall) col(black) margin(b=2)) xtitle("Years", size(medsmall)) ///
		ytitle("`labname'", margin(r=3) size(medsmall))   ///
		xlabel(0(1)`horizon') ///
		plotregion(color(white)) ///
		graphregion(color(white) )  legend(off) ///
		saving("$FIGUREDIR/irf_aipw_restricted_`aipwvar'.gph", replace)
		graph export "$FIGUREDIR\irf_aipw_restricted_`aipwvar'.png", replace
		
		local ivar = `ivar'+1

foreach stat in b se up90b lo90b up68b lo68b {
		capture drop `stat'
		capture drop `stat'*
		}
	
	

}



*******************************************************************************
foreach stat in b se up90b lo90b up68b lo68b {
		capture drop `stat'
		capture drop `stat'*
		}