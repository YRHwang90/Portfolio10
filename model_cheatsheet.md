Model cheatsheet
================
Yoo Ri Hwang
2022-06-09

## Overview

This document presents eleven models that used in this project.

![](images/cheatsheet2-01.jpg)

![](images/cheatsheet3-01.jpg)

# Discordant kinship model

## sib0

= REG0

No interaction  
gender compositoin: mixed sex vs same sex

(y_diff \~ y_mean + con_mean + con_diff + gender_composition_two_eff,
data=dat3)

![y\_{diff} = \\beta\_{0} + \\beta\_{1}y\_{mean} + \\beta\_{3}Con\_{mean}+ \\beta\_{3}Con\_{diff} + \\beta\_{4}genderComposition\_{two eff}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bdiff%7D%20%3D%20%5Cbeta_%7B0%7D%20%2B%20%5Cbeta_%7B1%7Dy_%7Bmean%7D%20%2B%20%5Cbeta_%7B3%7DCon_%7Bmean%7D%2B%20%5Cbeta_%7B3%7DCon_%7Bdiff%7D%20%2B%20%5Cbeta_%7B4%7DgenderComposition_%7Btwo%20eff%7D "y_{diff} = \beta_{0} + \beta_{1}y_{mean} + \beta_{3}Con_{mean}+ \beta_{3}Con_{diff} + \beta_{4}genderComposition_{two eff}")

![health\_{diff} = \\beta\_{0} + \\beta\_{1}health\_{mean} + \\beta\_{3}edu\_{mean}+ \\beta\_{3}edu\_{diff} + \\beta\_{4}genderComposition\_{two level}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;health_%7Bdiff%7D%20%3D%20%5Cbeta_%7B0%7D%20%2B%20%5Cbeta_%7B1%7Dhealth_%7Bmean%7D%20%2B%20%5Cbeta_%7B3%7Dedu_%7Bmean%7D%2B%20%5Cbeta_%7B3%7Dedu_%7Bdiff%7D%20%2B%20%5Cbeta_%7B4%7DgenderComposition_%7Btwo%20level%7D "health_{diff} = \beta_{0} + \beta_{1}health_{mean} + \beta_{3}edu_{mean}+ \beta_{3}edu_{diff} + \beta_{4}genderComposition_{two level}")

The coefficient in plots : meaning

con: continuous : education  
con_mean: education mean in sibling pair  
con_diff: education difference score in sibling pair  
gender composition two eff : mixed vs same sex (eff=effect coding)

## sib1

=REG1

no intearction gender composition = three level (mm, ff, mixed) because
this have three levels, we have two variables.

lm(y_diff \~ y_mean + con_mean + con_diff + ev1 + ev2, data=dat3) // ev
= effect code variable

![y\_{diff} = \\beta\_{0} + \\beta\_{1}y\_{mean} + \\beta\_{3}Con\_{mean}+ \\beta\_{3}Con\_{diff} + \\beta\_{4}ev1\_{three eff}+\\beta\_{5}ev1\_{threeeff}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bdiff%7D%20%3D%20%5Cbeta_%7B0%7D%20%2B%20%5Cbeta_%7B1%7Dy_%7Bmean%7D%20%2B%20%5Cbeta_%7B3%7DCon_%7Bmean%7D%2B%20%5Cbeta_%7B3%7DCon_%7Bdiff%7D%20%2B%20%5Cbeta_%7B4%7Dev1_%7Bthree%20eff%7D%2B%5Cbeta_%7B5%7Dev1_%7Bthreeeff%7D "y_{diff} = \beta_{0} + \beta_{1}y_{mean} + \beta_{3}Con_{mean}+ \beta_{3}Con_{diff} + \beta_{4}ev1_{three eff}+\beta_{5}ev1_{threeeff}")

![health\_{diff} = \\beta\_{0} + \\beta\_{1}health\_{mean} + \\beta\_{3}edu\_{mean}+ \\beta\_{3}edu\_{diff} + \\beta\_{4}genderComposition1\_{three level} +\\beta\_{5}genderComposition2\_{three level} ](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;health_%7Bdiff%7D%20%3D%20%5Cbeta_%7B0%7D%20%2B%20%5Cbeta_%7B1%7Dhealth_%7Bmean%7D%20%2B%20%5Cbeta_%7B3%7Dedu_%7Bmean%7D%2B%20%5Cbeta_%7B3%7Dedu_%7Bdiff%7D%20%2B%20%5Cbeta_%7B4%7DgenderComposition1_%7Bthree%20level%7D%20%2B%5Cbeta_%7B5%7DgenderComposition2_%7Bthree%20level%7D%20 "health_{diff} = \beta_{0} + \beta_{1}health_{mean} + \beta_{3}edu_{mean}+ \beta_{3}edu_{diff} + \beta_{4}genderComposition1_{three level} +\beta_{5}genderComposition2_{three level} ")

The coefficient in plots : meaning

con: continuous : education  
con_mean: education mean in sibling pair  
con_diff: education difference score in sibling pair

gender composition variable (=ev1, ev2)

ev1 : effect code variable (mixed vs 11) 1/-1 is for male/femlae  
ev2 : effect code variable (mixed vs -1-1)

## sib2

=REG2 gender composition: two level (mixed vs same) interaction included

lm(y_diff \~ y_mean + con_mean + con_diff + gender_composition_two_eff +
gender_composition_two_eff:con_diff , data=data)

![y\_{diff} = \\beta\_{0} + \\beta\_{1}y\_{mean} + \\beta\_{3}Con\_{mean}+ \\beta\_{3}Con\_{diff} + \\beta\_{4}genderComposition\_{two eff} + \\beta\_{5}genderComposition\_{two eff} \* Con\_{diff}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bdiff%7D%20%3D%20%5Cbeta_%7B0%7D%20%2B%20%5Cbeta_%7B1%7Dy_%7Bmean%7D%20%2B%20%5Cbeta_%7B3%7DCon_%7Bmean%7D%2B%20%5Cbeta_%7B3%7DCon_%7Bdiff%7D%20%2B%20%5Cbeta_%7B4%7DgenderComposition_%7Btwo%20eff%7D%20%2B%20%5Cbeta_%7B5%7DgenderComposition_%7Btwo%20eff%7D%20%2A%20Con_%7Bdiff%7D "y_{diff} = \beta_{0} + \beta_{1}y_{mean} + \beta_{3}Con_{mean}+ \beta_{3}Con_{diff} + \beta_{4}genderComposition_{two eff} + \beta_{5}genderComposition_{two eff} * Con_{diff}")

![health\_{diff} = \\beta\_{0} + \\beta\_{1}health\_{mean} + \\beta\_{3}edu\_{mean}+ \\beta\_{3}edu\_{diff} + \\beta\_{4}genderComposition\_{two level} +\\beta\_{5}genderComposition\_{two eff} \* edu\_{diff}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;health_%7Bdiff%7D%20%3D%20%5Cbeta_%7B0%7D%20%2B%20%5Cbeta_%7B1%7Dhealth_%7Bmean%7D%20%2B%20%5Cbeta_%7B3%7Dedu_%7Bmean%7D%2B%20%5Cbeta_%7B3%7Dedu_%7Bdiff%7D%20%2B%20%5Cbeta_%7B4%7DgenderComposition_%7Btwo%20level%7D%20%2B%5Cbeta_%7B5%7DgenderComposition_%7Btwo%20eff%7D%20%2A%20edu_%7Bdiff%7D "health_{diff} = \beta_{0} + \beta_{1}health_{mean} + \beta_{3}edu_{mean}+ \beta_{3}edu_{diff} + \beta_{4}genderComposition_{two level} +\beta_{5}genderComposition_{two eff} * edu_{diff}")

The coefficient in plots : meaning

con: continuous : education  
con_mean: education mean in sibling pair  
con_diff: education difference score in sibling pair  
gender composition two eff : mixed vs same sex (eff=effect coding)  
con_diffxgender_composition_two_eff : interaction between con_diff and
gender_composition_two_level

## sib3

REG3

gender compisition: three level (mm ff mixed)  
interaction included because this have three levels, we have two
variables.

sib_model\<- lm(y_diff \~ y_mean + con_mean + con_diff + ev1 + ev2 +
ev1*con_diff + ev2*con_diff, data=data)

![y\_{diff} = \\beta\_{0} + \\beta\_{1}y\_{mean} + \\beta\_{3}Con\_{mean}+ \\beta\_{3}Con\_{diff} + \\beta\_{4}ev1\_{three eff}+\\beta\_{5}ev1\_{three eff} + \\beta\_{6}ev1\_{three eff}\*con\_{diff} +\\beta\_{7}ev2\_{three eff}\*con\_{diff} ](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bdiff%7D%20%3D%20%5Cbeta_%7B0%7D%20%2B%20%5Cbeta_%7B1%7Dy_%7Bmean%7D%20%2B%20%5Cbeta_%7B3%7DCon_%7Bmean%7D%2B%20%5Cbeta_%7B3%7DCon_%7Bdiff%7D%20%2B%20%5Cbeta_%7B4%7Dev1_%7Bthree%20eff%7D%2B%5Cbeta_%7B5%7Dev1_%7Bthree%20eff%7D%20%2B%20%5Cbeta_%7B6%7Dev1_%7Bthree%20eff%7D%2Acon_%7Bdiff%7D%20%2B%5Cbeta_%7B7%7Dev2_%7Bthree%20eff%7D%2Acon_%7Bdiff%7D%20 "y_{diff} = \beta_{0} + \beta_{1}y_{mean} + \beta_{3}Con_{mean}+ \beta_{3}Con_{diff} + \beta_{4}ev1_{three eff}+\beta_{5}ev1_{three eff} + \beta_{6}ev1_{three eff}*con_{diff} +\beta_{7}ev2_{three eff}*con_{diff} ")

![health\_{diff} = \\beta\_{0} + \\beta\_{1}health\_{mean} + \\beta\_{3}edu\_{mean}+ \\beta\_{3}edu\_{diff} + \\beta\_{4}genderComposition1\_{three level} +\\beta\_{5}genderComposition2\_{three level}+ \\beta\_{6}ev1\_{three eff}\*con\_{diff} +\\beta\_{7}ev2\_{three eff}\*con\_{diff} ](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;health_%7Bdiff%7D%20%3D%20%5Cbeta_%7B0%7D%20%2B%20%5Cbeta_%7B1%7Dhealth_%7Bmean%7D%20%2B%20%5Cbeta_%7B3%7Dedu_%7Bmean%7D%2B%20%5Cbeta_%7B3%7Dedu_%7Bdiff%7D%20%2B%20%5Cbeta_%7B4%7DgenderComposition1_%7Bthree%20level%7D%20%2B%5Cbeta_%7B5%7DgenderComposition2_%7Bthree%20level%7D%2B%20%5Cbeta_%7B6%7Dev1_%7Bthree%20eff%7D%2Acon_%7Bdiff%7D%20%2B%5Cbeta_%7B7%7Dev2_%7Bthree%20eff%7D%2Acon_%7Bdiff%7D%20 "health_{diff} = \beta_{0} + \beta_{1}health_{mean} + \beta_{3}edu_{mean}+ \beta_{3}edu_{diff} + \beta_{4}genderComposition1_{three level} +\beta_{5}genderComposition2_{three level}+ \beta_{6}ev1_{three eff}*con_{diff} +\beta_{7}ev2_{three eff}*con_{diff} ")

The coefficient in plots : meaning

con: continuous : education  
con_mean: education mean in sibling pair  
con_diff: education difference score in sibling pair

gender composition variable (=ev1, ev2)

ev1 : effect code variable (mixed vs 11) 1/-1 is for male/femlae ev2 :
effect code variable (mixed vs -1-1)

con_diffxev1 : interaction between ev1 and education difference
con_diffxev2 : intearction between ev2 and education difference

# MLM

## MLM0

null model

![y\_{ij} = \\beta\_{0j} + e\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bij%7D%20%3D%20%5Cbeta_%7B0j%7D%20%2B%20e_%7Bij%7D "y_{ij} = \beta_{0j} + e_{ij}")

![ \\beta\_{0j} = \\gamma\_{00} + u\_{0j}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%20%5Cbeta_%7B0j%7D%20%3D%20%5Cgamma_%7B00%7D%20%2B%20u_%7B0j%7D " \beta_{0j} = \gamma_{00} + u_{0j}")

## MLM1

MLM1: NO level 2 model without interaction

lmer(y\~1 + catwtn + con + (1\|pid), data=data)

y = health  
catwtn = cat within = gender within = individual gender  
con = education

## MLM2

MLM1: NO level 2 model with interaction

lmer(y\~1 + catwtn + con + catwtn\*con + (1\|pid), data=data)

extracted coefficient name version

y = health  
con = continuous = education  
catwtn = cat within = individual gender  
catwtnxcon = interaction between individual gender and education

## MLM3

level 2 variable: gender composition (mixed vs same)  
No interaction

lmer(y \~ 1 + con + catbtw + (1\|pid), data=data)

![y\_{ij} = \\beta\_{0j} + \\beta\_{1j}\*edu + e\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bij%7D%20%3D%20%5Cbeta_%7B0j%7D%20%2B%20%5Cbeta_%7B1j%7D%2Aedu%20%2B%20e_%7Bij%7D "y_{ij} = \beta_{0j} + \beta_{1j}*edu + e_{ij}")

![ \\beta\_{0j} = \\gamma\_{00} + \\gamma\_{01}\*genderComposition + u\_{0j}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%20%5Cbeta_%7B0j%7D%20%3D%20%5Cgamma_%7B00%7D%20%2B%20%5Cgamma_%7B01%7D%2AgenderComposition%20%2B%20u_%7B0j%7D " \beta_{0j} = \gamma_{00} + \gamma_{01}*genderComposition + u_{0j}")

![ \\beta\_{1j} = \\gamma\_{10}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%20%5Cbeta_%7B1j%7D%20%3D%20%5Cgamma_%7B10%7D " \beta_{1j} = \gamma_{10}")

catbtw : gender composition : categorical between

## MLM4

gender composition: two level : mixed vs same  
with interaction

lmer(y \~ 1 + con + catbtw + con:catbtw + (1\|pid), data=data)

con: education  
catbtw: gender composition  
conxcatbtw : interaction between education and gender composition

### MLM5

gender composition : three levels (mm ff mixed )  
without interaction

lmer(y \~ 1 + con + ev1 + ev2 + (1\|pid), data=data)

con: education

since gender composition has three levels,  
effect coded variable should be two : ev1 and ev2

## MLM6

gender composition: three levels (mm ff mixed )  
with interaction

lmer(y \~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1\|pid), data=data)

in extracted form,

conxev1 : interaction between continuous (edu) and ev1  
conxev2 : interaction between continuous (edu) and ev2

so, coef_fixed_conxev2_p\_MLM6

this variable means coeffcient, fixed effect, interaction between
continuous (edu) and ev2’s p-value in MLM6 model.

## 

SS = same sex sibling OS = opposite sex sibling mix = mixed sex sibling

# 
