# clear all 
# cd ""
# *Get the data# set more off

rm(list=ls())

# use "nm.replication.dta", clear
library(foreign)
data <- read.dta("nm.replication.dta", convert.underscore=TRUE)
#*Generate Presidential 2-party Vote Share#	g dem_vs=kerry/(bush+kerry)

#Generate Presidential 2-party Vote Share
data$dem.vs <- data$kerry / (data$bush + data$kerry)
#*Generate Governor Vote#	g v_richardson=richard/(richard+dendahl)

#Generate Presidential 2-party Vote Share
data$v.richardson <- data$richard / (data$richard + data$dendahl)
#*Generate spending preferences Variable#	gen spending_pref = 1*q2full + 0.3*q2reduced + 0*q2none#	sum spending_pref, detail#		scalar spmed = r(p50)

#Generate spending preferences Variable
data$spending.pref <- (1* data$q2full) + (0.3* data$q2reduced) + (0* data$q2none)
spmed <- as.numeric(quantile(data$spending.pref, 0.5))
		#*Generate prior health bill vote#	gen prior_health = healthsolutions=="Y"

#Generate prior health bill vote
data$prior.health <- data$healthsolutions == "Y"


sort(names(data))
**************************************************************************************** Results in Paper ***********************************************************************************************#*Table 1
#tabstat rep q2full q1favor pres_2party self_2party runningforreelection unopposed prior_health, by(treatment)
#foreach var of varlist rep q2full q1favor pres_2party self_2party runningforreelection unopposed prior_health {
#        logit treatment `var'
#        }

####Table 1

table1 <- matrix(NA, 8, 3)
table1 <- as.data.frame(table1)
colnames(table1) <- c("Treatment", "Control", "p-value")
rownames(table1) <- c("Republican",
               "Constituent support for spending",
               "Constituent support for health care",
               "Bush vote-share 04",
               "Member vote-share 06",
               "Running for re-election",
               "Running unopposed",
               "Supported prior health care bill")

data.treat <- data[data$treatment==1,]
data.control <- data[data$treatment==0,]
vars <- cbind(data$rep, data$q2full, data$q1favor, data$pres.2party, data$self.2party, data$runningforreelection, data$unopposed, data$prior.health)
vars.treat <- cbind(data.treat$rep, data.treat$q2full, data.treat$q1favor, data.treat$pres.2party, data.treat$self.2party, data.treat$runningforreelection, data.treat$unopposed, data.treat$prior.health)
vars.control <- cbind(data.control$rep, data.control$q2full, data.control$q1favor, data.control$pres.2party, data.control$self.2party, data.control$runningforreelection, data.control$unopposed, data.control$prior.health)


#### First 2 columns
table.1.means <- function(var){
  percent(mean(var, na.rm=TRUE))
}

table1[,1] <- apply(vars.treat, 2, FUN=table.1.means)
table1[,2] <- apply(vars.control, 2, FUN=table.1.means)


#### P-values in the table
pvals <- function(var){
  value <- summary(glm(data$treatment ~ var, family="binomial"))$coef[2,4]
  format(value, digits=2)
}

table1[,3] <- apply(vars, 2, FUN=pvals)
table1


#### Table 2

#*Table 2: Results with Full Spending Preferences
#	sum q2full, detail
#		scalar q2med = r(p50)
#		scalar q2sd = r(sd)
#		scalar q225 = r(p25)
#	gen full_co = q2full<=scalar(q2med)
#	gen flt = treatment*full_co
#	** Column 1 - No Controls
#	probit sb24 treatment flt full_co, cluster(match_category)
#		test treatment=-flt
#		mfx 
#	** Column 2 - Controls for Partisanship & Vote Share
#	probit sb24 treatment flt full_co rep dem_vs, cluster(match_category)
#		test treatment=-flt
#		mfx

#### Table 2: Results with full spending preference

q2med <- as.numeric(quantile(data$q2full, 0.5))
q2sd <- as.numeric(sd(data$q2full))
q225 <- as.numeric(quantile(data$q2full, 0.25))

data$full.co <- 0
for (i in 1:length(data$full.co)){
  
  if (data$q2full[i] <= q2med){
    data$full.co[i] <- 1
  }
}

data$flt <- data$treatment * data$full.co
		
		
*Figure 2	reg q2full dem_vs	reg q2full v_richard	twoway (scatter q2full dem_vs) (lfit q2full dem_vs), yscale(range(.2 .6)) ylabel(.2(0.05).6) xscale(range(.1 .9)) xlabel(.1(0.2).9) ytitle(Percent of Constituents Prefering Spending Full Amount) xtitle("2004 Democratic Two-party" "Presidential Vote Share") legend(off) text(.57 .2 "R-squared=0.15", place(e)) saving(pres.gph, replace) plotregion(style(none)) plotregion(fcolor(white)) graphregion(fcolor(white))	twoway (scatter q2full v_richard) (lfit q2full v_richard), yscale(range(.2 .6)) ylabel(.2(0.05).6) xscale(range(.1 .9)) xlabel(.1(0.2).9) ytitle(Percent of Constituents Prefering Spending Full Amount) xtitle("2006 Democratic Two-party" "Gubernatorial Vote Share") legend(off) text(.57 .2 "R-squared=0.12", place(e)) saving(gov.gph, replace) plotregion(style(none)) plotregion(fcolor(white)) graphregion(fcolor(white))		gr combine pres.gph gov.gph,  rows(1) plotregion(style(none)) plotregion(fcolor(white)) graphregion(fcolor(white)) saving(figure2.q2full.gph, replace) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))

		*Figure 3	twoway (histogram q2full,  sort lcolor(gs8) fcolor(gs12)) (lowess sb24 q2full if treatment==1, bwidth(1) noweight yaxis(2) sort lcolor(black) lpattern(solid) lwidth(thick)) (lowess sb24 q2full if treatment==0, yaxis(2) sort lcolor(black) lpattern(dash) lwidth(thick)), ytitle(Density) xtitle(Percent Favoring Full Spending on Governor's Proposals)  ytitle(Probability of Voting in Favor of Bill, axis(2)) legend(off) text(2.15 .2 "Treatment", place(e)) text(7.7 .215 "Control", place(e)) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) ylabel(.4(.2)1, axis(2))*************************************************************************************** Supplementary Materials*******************************************************************************************Appendix B*Table 3 	ttest sb24, by(treatment)*Appendix C*Vary weight of the reduced category - Part or All of Tables 4, 6-8forvalues i = 0(1)7 {	gen cspref_`i' = 1*q2full + 0.`i'*q2reduced + 0*q2none	sum cspref_`i',detail	scalar t2 = r(p50)	gen byte co_`i' = cspref_`i'<=scalar(t2)	gen byte xtco_`i' = treatment*co_`i'	probit sb24 treatment co_`i' xtco_`i' rep dem_vs, cluster(match_category)		test treatment = -xtco_`i'		mfx	reg sb24 treatment co_`i' xtco_`i' rep dem_vs, cluster(match_category)	gen tXcsp_`i' = cspref_`i'*treatment	probit sb24 treatment cspref_`i' tXcsp_`i' rep dem_vs, cluster(match_category)		test treatment = -tXcsp_`i'		mfx	reg sb24 treatment cspref_`i' tXcsp_`i' rep dem_vs, cluster(match_category)	} * Results when varying the cut-off point: Results for Table 5**Varying cutoff for q2full	scalar q2co_low = scalar(q2med) - scalar(q2sd)	scalar q2co_high = scalar(q2med) + scalar(q2sd)		gen fcol = q2full<scalar(q2co_low)		gen fclt = treatment*fcol	probit sb24 treatment fclt fcol rep dem_vs, cluster(match_category)		test treatment=-fclt		mfx 	reg sb24 treatment fclt fcol rep dem_vs, cluster(match_category)		gen fcoh = q2full<scalar(q2co_high)		gen fcht = treatment*fcoh	probit sb24 treatment fcht fcoh rep dem_vs, cluster(match_category)	*Notes: Fails because of one empty cell among republicans		*1) Action is all on the low end**Varying cutoff for Constituency Spending Preferences*Move cut-off from median minus 1 sd	sum spending_pref, detail		scalar spsd = r(sd)		scalar col = scalar(spmed) - spsd		scalar coh = scalar(spmed) + spsd	gen byte cutoff_low = spending_pref <=scalar(col)	gen ctlow=treatment*cutoff_low	probit sb24 treatment cutoff_low ctlow rep dem_vs, cluster(match_category)		test treatment ct		mfx 	reg sb24 treatment cutoff_low ctlow rep dem_vs, cluster(match_category)*Move cut-off from median plus 1 sd	gen byte cutoff_high = spending_pref <=scalar(coh)		gen cthigh=treatment*cutoff_high			probit sb24 treatment cutoff_high cthigh rep dem_vs, cluster(match_category)		mfx *Using a Linear Interaction - "% Full Spending" Results: Table 7				gen Xfull = treatment*q2full	probit sb24 treatment q2full Xfull rep dem_vs, cluster(match_category)		test treatment=-Xfull		mfx 	reg sb24 treatment q2full Xfull rep dem_vs, cluster(match_category)