
use londondataset2002_2006, clear

**********************
*Descriptive analysis
**********************

*Raw data plot (Figure 1)
*(a) Create a list of dates to add as dashed guidelines in the plots
foreach year of numlist 2003/2006{
 local datejan1 = d(1/1/`year')
 local datelines `datelines' `datejan1'
 }
*(b) Plot the exposure and outcome data against time, and combine the plots
graph twoway scatter ozone date, msize(vsmall) mcol(gs5) name(ozone, replace) ytitle("Daily mean ozone level (ug/m3)") xtitle(Date) xline(`datelines', lp(dash) lc(gs10)) title(Ozone levels over time)
graph twoway scatter numdeaths date, msize(vsmall) mcol(black) name(deaths, replace) ytitle(Daily number of deaths) xtitle(Date) xline(`datelines', lp(dash) lc(gs10)) title(Daily deaths over time)
graph combine deaths ozone, cols(1) xcommon name(Fig1, replace)
graph drop deaths ozone

*Summary statistics
summ ozone temperature relative_humidity, d

*Correlation matrix for the explanatory variables
pwcorr ozone temperature relative_humidity

*Divide the ozone variable by 10 so that model estimates refer to a "per 10ug/m3 increase" (as per convention)
replace ozone = ozone/10
rename ozone ozone10
label var ozone10 "Ozone level in ug/m3 divided by 10"


******************************************************************************************************
*Modelling seasonality and long-term trend in the outcome data (leaving out the main exposure for now)
******************************************************************************************************

*Option 1 - Time-stratified model (simple indicator variables)
*(a) Generate a variable containing "elapsed month" (i.e. months since start of study period)
gen month=month(date)
gen year=year(date)
gen elapsedmonths=(year-1996)*12+month
drop month year

*(b) Fit elapsed month into the model with indicator variables ("i.")
glm numdeaths i.elapsedmonth, family(poisson) scale(x2) eform

*(c) Plot the raw outcome data and the fitted model over time
predict xb_monthstrata, xb
gen fitted_monthstrata =exp(xb_monthstrata)
graph twoway scatter numdeaths date, msize(tiny) mc(gs8) ylabel(100 200 300) ///
	|| line fitted_monthstrata date, lc(black) lw(thick) xtitle(Date) ytitle(N deaths) ///
	name(monthstrata, replace) legend(off) title("(a) Time stratified model (month strata)")
drop xb elapsedmonth fitted_monthstrata


*Option 2 - Periodic functions (Fourier terms) 
*(a) Generate sine and cosine functions of time with annual period, plus 3 harmonics
gen degrees=(date/365.25)*360
fourier degrees, n(4)
*(b) Fit in the model along with a linear term to capture long-term trend
* Note, coefficients of the sin/cos terms are estimated by maximum likelihood such that ...
*...the linear combination models the seasonal patterns in the outcome data as closely as possible
glm numdeaths cos* sin* date, family(poisson) scale(x2)
*(c) Plot the raw outcome data and the fitted model over time
predict xb_fourier, xb
gen fitted_fourier =exp(xb_fourier)
graph twoway scatter numdeaths date, msize(tiny) mc(gs8) ylabel(100 200 300) ///
	|| line fitted_fourier date, lc(black) lw(thick) xtitle(Date) ytitle(N deaths) ///
	name(fourier, replace) legend(off) title("(b) Sine/cosine functions (Fourier series)")
drop xb degrees cos* sin* fitted_fourier

*Option 3 - Flexible spline function
*(a) Create a basis for the spline with 34 equally spaced reference points
*Note that "frencurv" is not a core part of Stata, but can be downloaded (to find, type "net search frencurv")
frencurv, xvar(date) gen(spl) refpts(15341 15396 15451 15506 15562 15617 15672 15728 15783 15838 15894 15949 16004 16059 16115 16170 16225 16281 16336 16391 16447 16502 16557 16612 16668 16723 16778 16834 16889 16944 17000 17055 17110 17166) power(3)
*(b) Fit in the model - coefficients of the basis terms are estimated by maximum likelihood...
*...such that the linear combination models the seasonal patterns in the outcome data as closely as possible
glm numdeaths spl*, family(poisson) scale(x2) nocons
*(c) Plot the raw outcome data and the fitted model over time
predict xb_spline, xb
gen fitted_spline =exp(xb_spline)
graph twoway scatter numdeaths date, msize(tiny) mc(gs8) ylabel(100 200 300) ///
	|| line fitted_spline date, lc(black) lw(thick) xtitle(Date) ytitle(N deaths) ///
	name(spline, replace) legend(off) title("(c) Flexible cubic spline model")
drop xb_spline fitted_spline

*Combine the "fitted model" graphs from the 3 options
graph combine monthstrata fourier spline, cols(1) xcommon ysize(8) name(Fig2, replace)
graph drop monthstrata fourier spline

*Plot the raw residuals over tim in a model containing only a smooth function of time (from Option 3)
predict residuals_spline, r
graph twoway scatter residuals_spline date, msize(tiny) mc(black) yline(0, lc(gs7) lp(dash) lw(thick)) ///
	xtitle(Date) ytitle("Residual (Observed - Fitted*)") name(Fig3, replace)
drop residuals_spline

*******************************************************************************************************
*Model for ozone-mortality association, with and without adjustment for seasonality and long-term trend 
*******************************************************************************************************
*(a) without adjustment for long-term patterns
glm numdeaths ozone10, family(poisson) scale(x2) eform
*(b) with adjustment for long-term patterns (using flexible spline, option 3)
glm numdeaths ozone10 spl*, family(poisson) scale(x2) nocons eform

*Note: the scale(x2) adjusts the standard errors to account for overdispersion 

*********************************************************
*Model for ozone, further adjusting for daily temperature 
*********************************************************
*To allow for non-linearity in a simple way, we will fit temperature as a categorical variable in deciles
*However more sophisticated options are available - see Armstrong B., Epidemiology, 2006
xtile tempdecile=temperature, nq(10)
glm numdeaths ozone10 i.tempdecile spl* , family(poisson) scale(x2) nocons eform

******************************************************************************
*Exploring the lagged (delayed) effect of ozone (controlling for temperature)
******************************************************************************
*(a) Generate copies of the ozone and temperature variables, time-shifted by 0 to 7 days inclusive
sort date
for num 0/7: gen ozone10lagX = ozone10[_n-X]
for num 0/7: gen temperaturelagX = temperature[_n-X]
foreach lag of numlist 0/7{
	sort date
	xtile tempdecile_lag`lag'=temperaturelag`lag', nq(10)
	}
	drop temperaturelag*
	
*(b) Fit the individual lag models
foreach lag of numlist 0/7{
glm numdeaths ozone10lag`lag' i.tempdecile_lag`lag' spl*, family(poisson) scale(x2) nocons eform
}
*(c) Fit the unconstrained distributed lag model
glm numdeaths ozone10lag* i.tempdecile_lag* spl*, family(poisson) scale(x2) nocons eform

*(d) Fit the constrained ("lag stratified") distributed lag model
constraint define 1 ozone10lag1=ozone10lag2
constraint define 2 ozone10lag3=ozone10lag4
constraint define 3 ozone10lag3=ozone10lag5
constraint define 4 ozone10lag3=ozone10lag6
constraint define 5 ozone10lag3=ozone10lag7
glm numdeaths ozone10lag* i.tempdecile_lag* spl*, constraints(1/5) family(poisson) scale(x2) nocons eform

*Note: to plot Figure 4 displaying the estimates from the above models requires the model estimates 
* to be saved to a dataset for plotting; the coding for this is somewhat complex and is not relevant
* to the core aim of our tutorial paper, so we omit it here to avoid confusion

***************************
*Model checking (residuals)
***************************
*(a) Re-fit the unconstrained distributed lag model
glm numdeaths ozone10lag* i.tempdecile_lag* spl*, family(poisson) scale(x2) nocons eform

*(b) Generate the deviance residuals and plot vs time
predict dres, d
scatter dres date, name(FigA1, replace)

*(c) Partial autocorrelation function plot (PAF) of the deviance residuals - original model
tsset date
pac dres, title("(a) From original model") name(a2a, replace)

*(d) Include the 1-day lagged residuals in the model, and re-draw the PAF plot
gen dresL1=dres[_n-1]
noi glm numdeaths ozone10lag* i.tempdecile_lag*  spl* dresL1, family(poisson) scale(x2) nocons eform
predict dresV2, d
pac dresV2, title("(b) From model adjusted for residual autocorrelation") name(a2b, replace)

graph combine a2a a2b, cols(1) ysize(9) name(FigA2, replace)
graph drop a2a a2b

drop dres*


