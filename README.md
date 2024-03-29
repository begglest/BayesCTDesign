
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Overview

`BayesCTDesign` is a set of functions to help clinical trial researchers
calculate power and sample size for two-arm Bayesian randomized clinical
trials that do or do not incorporate historical control data. The power
and sample size calculations are always based on a two-sided test.
Historical control data is incorporated using a power prior. The primary
functions are:

-   `historic_sim()`
-   `simple_sim()`.

The `historic_sim()` function allows a clinical trialist to study trial
designs that incorporate historical control data. The `simple_sim()`
function allows a clinical trialist to study trial designs that are
simply two-arm Bayesian clinical trials. These two primary functions are
simulation functions and they output a list containing simulation
results. The first element of the returned list is a list that has at
least 2 arrays and depending on user requests may contain up to 5
arrays. The object created by `historic_sim()` and `simple_sim()` has
two methods available:

-   `print()`
-   `plot()`.

These methods allow the user to select a two-dimensional slice of an
array from `historic_sim()` or `simple_sim()` and print or plot the
slice. Both `historic_sim()` and `simple_sim()` allow the user to study
two-arm clinical trials with Gaussian, Poisson, Bernoulli, Lognormal,
Weibull, and Piecewise Exponential outcomes.

Note, since `historic_sim()` and `simple_sim()` are simulation programs,
a big simulation can take a long time to run. `historic_sim()` and
`simple_sim()` have a parameter called `quietly` that when it is set to
`FALSE` will make `historic_sim()` and `simple_sim()` print out a
progress note to indicate that the simulation is running and not hung-up
and to indicate where the simulation is within the set of user defined
scenarios. A nice feature, but the default for `quietly` is `TRUE`. If a
user sets `quietly` to `FALSE` and runs the code interactively in
RStudio or runs it in the R console, the printed information will always
be printed on the same line, and it can be a useful notice to the user
about how the simulation run is going. On the other hand, if a user sets
`quietly` to `FALSE` and runs the code in a Notebook or knitr document,
then each printed note is printed on a separate line within the output
section of the notebook or document. This result makes the output
section of the Notebook or knitr document unnecessarily long, so in
Notebooks and knitr documents it is best to set `quietly` to `TRUE` and
just expect the simulation will finish.

# Background

At some point during the design process, a clinical trial researcher who
is designing a two-arm Bayesian randomized clinical trial needs to make
decisions about power and sample size within the context of the
hypothesized treatment effects. Through simulation, the `simple_sim()`
function will estimate power and other user specified clinical trial
characteristics at desired sample sizes given user defined scenarios
about alpha, treatment effect, control group characteristics, and
outcome. Using simulation to study trial characteristics of simple two
arm studies may seem unnecessary, but `BayesCTDesign` was created for a
more complicated case where historical control data is incorporated into
the trial estimation process using a power prior. Once the primary
purpose of `BayesCTDesign` was accomplished, the more simpler case
represented by `simple_sim()` was a very easy addition to the package,
so it was included for completeness.

If the clinical trial researcher has access to historical control data,
then the researcher can use `BayesCTDesign` to design a two-arm Bayesian
randomized clinical trial that incorporates the historical data. In such
a case, the researcher needs to work through the potential consequences
of historical and randomized control differences on trial
characteristics, in addition to making a decision about power or sample
size. If a researcher designs a clinical trial that will incorporate
historical control data, the researcher needs the randomized controls to
be from the same population as the historical controls. What if this is
not the case when the designed trial is implemented? During the design
phase, the researcher needs to investigate the negative effects of
possible historic/randomized control differences on power, type one
error, and other trial characteristics. Using this information, the
researcher should design the trial to mitigate these negative effects.
Through simulation, the `historic_sim()` function will allow the
clinical trialist to investigate possible trial scenarios and make
decisions that will mitigate the negative effects of
historical/randomized control differences.

## Development Version Installation

You can install BayesCTDesign from github with:

``` r
# install.packages("devtools")
devtools::install_github("begglest/BayesCTDesign")
```

## Example Usage 1: A call to `historic_sim()`

If a user wants to use `historic_sim()`, then the user should have
access to historical control data with the required structure. See the
help page for `historic_sim()` for more details about the necessary
structure for historical control data. However, `BayesCTDesign` also has
several data generating functions that can be used for exploratory
purposes. For example, `genweibulldata()` can be used to simulate a
trial where the outcome is a Weibull time-to-event variable. See the
help pages for `BayesCTDesign` to read more about other data generating
functions. To generate a sample of historical control data with a
Weibull outcome, one can make a call to `genweibulldata()` and then keep
data from the control arm only.

For an example, consider a scenario where we have historical control
Weibull data from 60 subjects, which are right censored at 3.0. We want
to determine the sample size needed for a two-arm clinical trial that
will utilize information from these 60 control subjects to detect a
treatment hazard ratio of 0.7 with 80% power and a two-sided alpha of
0.05. Assume the outcome follows a Weibull distribution. Although the
targeted treatment hazard ratio is 0.7, assume the clinical trialist
believes the effect may range between 0.6 and 1.0. Initially, the
clinical trialist believes the required sample size per arm will be in
the range of 75 to 175 subjects per arm.

Now, since the clinical trialist might incorporate historical control
data into the trial design, the clinical trialist needs to assess the
risk of including these historical control data. The model will use
historical control and randomized control data as if they were samples
from the same control population. Yet, this historical control data
might be very different from the randomized control data. If the control
groups are different, including the historical control data may
significantly bias the results, because the model will produce a biased
estimate of the control group hazard. Assume the clinical trialist
believes the randomized and historical control data might differ in such
a way that the hazard ratio between the two controls groups will be in
the range of 0.8 to 1.2. Will such differences create unacceptable bias
in the final Bayesian estimate of the treatment hazard ratio? If so, how
might the clinical trialist mitigate against this possible bias?

The clinical trialist will need to determine if differences between
randomized and historical controls within this range will have a
significant effect on power and other trial characteristics. Now, one
way to control the influence of randomized/historical control
differences on the treatment effect estimate is to reduce the amount of
historical control information that is used in the estimation process.
`BayesCTDesign` uses a power prior to incorporate historical control
data information, and the amount of influence from historical control
data is defined via the `a0` parameter of the power prior.

Putting the scenario all together, the clinical trialist needs to study
power to detect treatment hazard ratios (experimental over control)
ranging from 0.6 to 1.0 while allowing sample size to range from 75 to
175 and including a sample of 60 historical controls. At the same time,
the clinical trialist also needs to study the effects of
randomized/control differences (control hazard ratio ranging from 0.8 to
1.2) when the power prior parameter, `a0`, ranges from 0(no historical
data included) to 1(all historical data included). Finally, assume the
outcome is a Weibull distributed time-to-event outcome, and the expected
treatment hazard ratio is 0.7. The test will be two-sided:

H0: Treatment hazard ratio is 1. Ha: Treatment hazard ratio is not equal
to 1.

With this set up information, the clinical trialist can run
`historic_sim()` and generate several arrays of simulated output. In
turn, these arrays of simulation output can be investigated to determine
what sample size the trial should use and make a decision on what value
of `a0` should be used.

In the following call to `historic_sim()`, we use the parameter
`subj_per_arm` to assign the samples sizes that will be investigated. We
use the parameter `effect_vals` to assign the hazard ratios that will be
investigated. The parameter `effect_vals` defines the treatment effect
among randomized subjects. We use the parameter `rand_control_diff` to
assign the possible randomized/control hazard ratios, and finally, we
use the parameter `a0_vals` to assign a set of historical control
information fractions to consider. The `trial_reps = 500` code states
that each combination of `subj_per_arm`, `effect_vals`,
`rand_control_diff`, and `a0_vals` will be simulated 500 times. The
parameter `hist_control_data` defines the dataset of historical
controls. The code, `censor_value = 3`, states that within each
replicated trial, the simulated trial outcome data will be right
censored at 3. The parameter `alpha` defines the acceptable Type 1 error
for making inference about the treatment hazard ratio with a two-sided
test. See the help page for `historic_sim()` for more details about how
inference is made about the treatment hazard ratio within each simulated
trial. The parameters `get_var`, `get_bias`, and `get_mse` are set to
`TRUE`, so in addition to getting an array of power and average effect
estimates, `historic_sim()` will also return arrays of variance, bias,
and mse. Finally, the `seedval` parameter is set for reproducibility.

On an i7 Windows machine the following call to `historic_sim()` took
about 1.5 hours to run, so it will take a while to complete; however,
within these 75 minutes the simulation made estimates under 300
scenarios and within each scenario the simulation ran 500 simulations.
If the same set of combinations were considered in a simulation where
MCMC techniques, such as a call to JAGS, were used within each trial
replication to make inference about the treatment hazard ratio, then the
length of time needed would probably take much longer. Although such a
longer simulation might give more exact results, for power and sample
size calculations the results from `BayesCTDesign` are sufficient for
design purposes.

``` r
#First, for illustrative purposes, generate the sample of historical data.
library(BayesCTDesign)
set.seed(2250)
SampleHistData <- genweibulldata(sample_size=60, scale1=2.82487,
                                 hazard_ratio=1.0, common_shape=3,
                                 censor_value=3)
histdata <- subset(SampleHistData, subset=(treatment==0))
histdata$id <- histdata$id+10000

#Next, run the simulation, using historic_sim().
weibull_test <- historic_sim(trial_reps = 500, outcome_type = "weibull",
                             subj_per_arm = c(75, 100, 125, 150, 175),
                             a0_vals = c(0, 0.25, 0.50, 0.75, 1),
                             effect_vals = c(0.6, 0.7, 0.8, 0.9, 1),
                             rand_control_diff = c(0.8, 1, 1.2),
                             hist_control_data = histdata, censor_value = 3, 
                             alpha = 0.05, get_var = TRUE, get_bias = TRUE, 
                             get_mse = TRUE, seedval=123)

#Tabulate the simulation results for power when no difference exists between
#historical and randomized controls.  For this first table, look at the case
#were experimental over control hazard ratio among randomized subjects is 0.6.
#The resulting table is a table of power with different sample sizes represented
#in the rows, and different a0 values represented in the columns.
test_table0_6 <- print(x=weibull_test, measure="power",
                       tab_type="WX|YZ", effect_val=0.6,
                       rand_control_diff_val=1.0)
print(test_table0_6)
#>         0  0.25   0.5  0.75     1
#> 75  0.670 0.728 0.756 0.796 0.780
#> 100 0.776 0.814 0.854 0.880 0.892
#> 125 0.876 0.906 0.940 0.918 0.940
#> 150 0.932 0.944 0.944 0.970 0.964
#> 175 0.960 0.972 0.980 0.980 0.988

#For this second table, look at the case were experimental over control hazard 
#ratio among randomized subjects is 0.7.  Look at the case where no differences
#exist between historical and randomized controls.
#The resulting table is a table of power with different sample sizes represented
#in the rows, and different a0 values represented in the columns.
test_table0_7 <- print(x=weibull_test, measure="power",
                       tab_type="WX|YZ", effect_val=0.7,
                       rand_control_diff_val=1.0)
print(test_table0_7)
#>         0  0.25   0.5  0.75     1
#> 75  0.434 0.394 0.458 0.488 0.520
#> 100 0.524 0.528 0.574 0.576 0.592
#> 125 0.602 0.654 0.642 0.662 0.730
#> 150 0.692 0.744 0.758 0.744 0.804
#> 175 0.788 0.742 0.816 0.822 0.834


#For the third table, look at the case were experimental over control hazard 
#ratio among randomized subjects is 1.0.  Look at the case where no differences
#exist between historical and randomized controls.
#The resulting table is a table of power with different sample sizes represented
#in the rows, and different a0 values represented in the columns.
test_table1_0 <- print(x=weibull_test, measure="power",
                       tab_type="WX|YZ", effect_val=1.0,
                       rand_control_diff_val=1.0)
print(test_table1_0)
#>         0  0.25   0.5  0.75     1
#> 75  0.044 0.048 0.044 0.030 0.028
#> 100 0.060 0.044 0.040 0.028 0.030
#> 125 0.052 0.052 0.032 0.036 0.028
#> 150 0.058 0.036 0.046 0.048 0.030
#> 175 0.060 0.036 0.028 0.044 0.042

#For the fourth table, look at the case were experimental over control hazard 
#ratio among randomized subjects is 0.7.  Look at the case where the difference
#between historical and randomized controls is represented by a hazard ratio
#of 1.2.  The resulting table is a table of power with different sample sizes 
#represented in the rows, and different a0 values represented in the columns.
test_table0_7_1_2 <- print(x=weibull_test, measure="power",
                           tab_type="WX|YZ", effect_val=0.7,
                           rand_control_diff_val=1.2)
print(test_table0_7_1_2)
#>         0  0.25   0.5  0.75     1
#> 75  0.430 0.402 0.362 0.364 0.344
#> 100 0.536 0.530 0.548 0.472 0.444
#> 125 0.700 0.638 0.630 0.584 0.574
#> 150 0.688 0.724 0.720 0.644 0.702
#> 175 0.798 0.782 0.804 0.750 0.728

#For the fifth table, look at the case were experimental over control hazard 
#ratio among randomized subjects is 0.7.  Look at the case where the difference
#between historical and randomized controls is represented by a hazard ratio
#of 0.8.  The resulting table is a table of power with different sample sizes 
#represented in the rows, and different a0 values represented in the columns.
test_table0_7_0_8 <- print(x=weibull_test, measure="power",
                           tab_type="WX|YZ", effect_val=0.7,
                           rand_control_diff_val=0.8)
print(test_table0_7_0_8)
#>         0  0.25   0.5  0.75     1
#> 75  0.380 0.470 0.530 0.606 0.702
#> 100 0.470 0.524 0.608 0.684 0.738
#> 125 0.556 0.636 0.686 0.776 0.842
#> 150 0.628 0.664 0.746 0.824 0.854
#> 175 0.694 0.758 0.818 0.830 0.862

#For the sixth table, look at the case were sample size is set to 150, which 
#appears to be equal to the required sample size for 0.80 power when all data
#from historical controls are included in trial, and effect is set to 0.7.  The 
#resulting table is a table of power with different historical/randomized control 
#differences in the rows, and different a0 values represented in the columns.
test_table150_0_7 <- print(x=weibull_test, measure="power",
                          tab_type="ZX|WY", subj_per_arm_val=150,
                          effect_val=0.7)
print(test_table150_0_7)
#>         0  0.25   0.5  0.75     1
#> 0.8 0.628 0.664 0.746 0.824 0.854
#> 1   0.692 0.744 0.758 0.744 0.804
#> 1.2 0.688 0.724 0.720 0.644 0.702

#For the seventh table, look at the power by control difference by a0 table again, 
#but this time set subj_per_arm_val equal to 175.
test_table175_0_7 <- print(x=weibull_test, measure="power",
                          tab_type="ZX|WY", subj_per_arm_val=175,
                          effect_val=0.7)
print(test_table175_0_7)
#>         0  0.25   0.5  0.75     1
#> 0.8 0.694 0.758 0.818 0.830 0.862
#> 1   0.788 0.742 0.816 0.822 0.834
#> 1.2 0.798 0.782 0.804 0.750 0.728

#Finally, for the eighth table, look at the power by control difference by a0 
#table again, but this time set subj_per_arm_val equal to 175 and the 
#effect value to 1.0.
test_table175_1_0 <- print(x=weibull_test, measure="power",
                          tab_type="ZX|WY", subj_per_arm_val=175,
                          effect_val=1.0)
print(test_table175_1_0)
#>         0  0.25   0.5  0.75     1
#> 0.8 0.076 0.048 0.040 0.064 0.074
#> 1   0.060 0.036 0.028 0.044 0.042
#> 1.2 0.052 0.036 0.046 0.058 0.056
```

From the output `test_table0_6` we see that if the true hazard ratio
between randomized experiment and control groups is 0.6, then a trial
that does not use any historical control data, a0 = 0 (first column),
will require a per arm sample size within the range of 100 to 125
subjects. Estimated power with 100 subjects is 0.776, but with 125
subjects the estimated power is 0.876. On the other hand, if the
historical and randomized controls are from the same population and all
information in the historical data is utilized in the estimation
process, the required per arm sample size is between 75 subjects and 100
subjects (power = 0.780 and 0.892 respectively).

Now, the clinical trialist primary focus is on the case where the true
hazard ratio is 0.7. From the output `test_table0_7` we see that if the
true hazard ratio between randomized experiment and control groups is
0.7, then a trial that does not use any historical control data, a0=0
(first column), will require a per arm sample size just greater than 175
subjects. Estimated power for 175 subjects when a0 is equal to 0.788. On
the other hand, if all the information in the historical data is
utilized in the estimation process, a0 = 1 (last column), the required
per arm sample size is about 150 subjects (power equal to 0.804).

Finally, `test_table1_0` shows us that when there is no difference
between historical and randomized controls, the Type 1 error is
maintained around 0.05 or below 0.05 for all combinations of sample size
and a0. Overall mean of all 25 scenarios represented in `test_table1_0`
is 0.04. Note, for the best estimates of Type 1 error, many more than
500 trial replicates are needed.

Based on these three tables, it seems that when all the information in
the historical control data is used, the trial needs to enroll about 150
subjects per arm for 80 percent power to detect a 0.7 hazard ratio
between experimental and controls groups. Yet, what are the effects of
randomized/historical control differences. A set of different tables can
be used to address this question. The table `test_table0_7_1_2` shows
that if the true randomized experiment to control group hazard ratio is
0.7, but the randomized and historical controls differ by a hazard ratio
of 1.2, then the power drops from 0.80 to 0.70 (fourth row, last
column). On the other hand, table `test_table0_7_0_8` shows that if the
true randomized experiment to control group hazard ratio is 0.7, but the
randomized and historical controls differ by a hazard ratio of 0.8, then
the power increases from 0.80 to 0.85 (fourth row, last column). In both
cases, the change in power is due to the resulting bias in experiment
over control hazard ratio estimation caused by the differences in
randomized and historical controls. This bias can go in either direction
depending on the true experimental over control group effect among
randomized in relation to the randomized and historical control group
differences.

To get a final view of the effect of historical/randomized control
differences, consider table `test_table150_0_7`. Notice that in this
table, we see that when a0 is equal to 0.25 (second column), the
historical/randomized control difference of 0.8 tended to have lower
power than the historical/randomized control difference of 1.2; however,
when a0 is equal to 0.5, then power is rather close regardless of the
historical/randomized control differences. If the randomized control
trial sample size could be increased so that this third column was
roughly 0.80, then a trial with a0 = 0.50 and that sample size per arm
would be rather robust to historical/randomized control differences
within the range of 0.8 to 1.2. If necessary, the clinical trialist
could rerun the simulations to look at larger sample sizes, but our
simulation did consider a larger sample size of 175 subjects per arm.
When a sample size of 175 per arm is considered, we find that at a0=0.5,
the estimated power is always 0.80 or greater. In conclusion, a clinical
trialist could conclude that at least 80% power will be obtained when
175 subjects are randomized into each of two arms (experimental and
control) and the true effect is equal to 0.7, and when data from 60
historical controls are included in the trial using a power prior a0
value of 0.5 and the historical controls do not differ from randomized
controls in terms of a hazard ratio which ranges from 0.8 to 1.2.

Earlier, we saw that a simple randomized trial of 175 subjects per arm
would not be sufficient to have 80% power to detect a true hazard ratio
of 0.7 with a two-sided test and an alpha of 0.05. With the addition of
60 historical controls using a power prior with a0 value of 0.5, we see
that the trial now has at least 80% power to detect the true hazard
ratio of 0.7 as long as the historical and randomized controls do not
differ more than what was explored in these simulations. Finally, one
look at `test_table175_1_0` shows that when the sample size is set to
175 and the true effect is 1.0, the power (which now is equal to Type 1
error) is always less than 0.05. It follows that now the clinical
trialist has a good idea of how to incorporate the historical controls
and gain some power without risking Type 1 error as long as the
differences between historical and randomized controls do not go beyond
what is expected.

In addition to tabulating the results as we have done, we can also plot
a slice of the results. For example, the following code will create a
plot of power as a function of sample size and a0 value, when the true
experiment to control hazard ratio is 0.7 and the randomized and
historical controls differ by a hazard ratio of 1.2.

``` r
#Create a plot of the power simulation results.
plot(x=weibull_test, measure="power", tab_type="WX|YZ",
     smooth=FALSE, plot_out=TRUE, effect_val=0.7,
     rand_control_diff_val=1.2)
```

![](README-example2-1.png)<!-- -->

This plot shows the uncertainty in the estimates and suggests that more
replicates are needed for an actual trial design study; however, the
plot does show how in general as sample size increases the design
increases in power and how at any sample size the inclusion of
historical control data that is different from the randomized controls
affects the power of the treatment effect test. As more replicates are
used in the simulation, this effect will be more obvious.

## Example Usage 2: A call to `simple_sim()`

Now, a clinical trialist may want to simply design a two-arm Bayesian
trial and not include historical data. Consider the case where a
clinical trialist wants to use a piecewise exponential model with 6 time
intervals (0 to 0.3 years, 0.3 to 0.9 years, 0.9 to 1.5 years, 1.5 to
2.1 years, 2.1 to 2.4 years, and 2.4 years and higher). Assume the
vector of constant hazards for this six piece pwe is
(0.19,0.35,0.56,0.47,0.38,0.34). Finally, assume the true hazard ratio
(experimental over control) within any time interval is 0.8, but the
clinical trialist wants to study hazard ratios ranging from 0.6 to 1.0.
We now are ready to use `simple_sim()` to determine the required sample
size if the clinical trialist wants the trial to have 80% power to
detect a hazard ratio of 0.8 using a two-sided test and alpha = 0.05.
The clinical trialist believes the required sample size will be in the
range of 400 to 550 subjects.

``` r
#Define control group pwe distribution parameters.
#Note: as with the call to historic_sim(), running this code takes a good while
#to run, on an i7 it took about 6 hours.  Also, notice that with results from simple_sim(), only the object from simple_sim() and the measure needs to be specified in the call to print_table().
set.seed(2250)
time.vec <- c(0.3,0.9,1.5,2.1,2.4)
lambdaHC.vec <- c(0.19,0.35,0.56,0.47,0.38,0.34)

#Run a pwe simulation, using historic_sim().
pwe_test <- simple_sim(trial_reps = 500, outcome_type = "pwe",
                        subj_per_arm = c(400,425,450,475,500,525,550),
                        effect_vals = c(0.6,0.7,0.8,0.9,1.0),
                        control_parms = lambdaHC.vec, time_vec = time.vec,
                        censor_value = 3, alpha = 0.05, get_var = TRUE,
                        get_bias = TRUE, get_mse = TRUE, seedval=123)

test_table <- print(x=pwe_test, measure="power")
#> [1] "Since simple_sim was used, tab_type was set to WY|XZ"
#> [1] "Values for tab_type, subj_per_arm_val, a0_val, effect_val, and rand_control_diff_val were ignored"
#> [1] "This works towards putting all results in a single table, effect by sample size"
print(test_table)
#>     0.6   0.7   0.8   0.9     1
#> 400   1 0.974 0.734 0.212 0.060
#> 425   1 0.978 0.754 0.230 0.052
#> 450   1 0.992 0.812 0.250 0.048
#> 475   1 0.990 0.796 0.260 0.054
#> 500   1 0.994 0.820 0.266 0.078
#> 525   1 0.988 0.828 0.278 0.068
#> 550   1 0.998 0.848 0.318 0.064
```

Based on these results, the trial will need to enroll between 475 and
500 subjects per arm to detect a hazard ratio of 0.8 with 80% power at
an alpha of 0.05.
