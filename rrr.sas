OPTIONS nofmterr nocenter nonumber nodate;
libname cjdats 'C:\Users\Donna Coffman\Documents\Conferences\AHSR';

data cjdats.rrr;
run;


********************************************************************************************************;
*Continuous exposures;
********************************************************************************************************;

*estimate propensity scores and create weights;
*denominator model;
proc reg data=cjdats.rrr;
 model RRA_30 = COND CBLACK ABUSE6M suprob drgprb alcprb ownhouse nevermarried married 
                fulltime parttime LIVSP JAIL30D SATXIMP ARST30D MEDINS HS DGTXJLC ALCDYS FINSUP1 condom_intake;
 output out=cjdats.rrr student=rden;
run;
quit;

*numerator model;
proc reg data=cjdats.rrr;
 model RRA_30 = ;
 output out=cjdats.rrr student=rnum;
run;
quit;

data cjdats.rrr;
 set cjdats.rrr;
   pnum = exp(-.5*(rnum**2))/2.506;
   pden = exp(-.5*(rden**2))/2.506;
   wt_cont = pnum/pden;
run;

proc means data=cjdats.rrr;
 var wt_cont;
run;

*assess balance;
*before weighting;
proc corr data=cjdats.rrr;
 var RRA_30;
 with COND CBLACK ABUSE6M suprob drgprb alcprb ownhouse nevermarried married 
      fulltime parttime LIVSP JAIL30D SATXIMP ARST30D MEDINS HS DGTXJLC ALCDYS FINSUP1 condom_intake;
run;

*after weighting;
proc corr data=cjdats.rrr;
 var RRA_30;
 with COND CBLACK ABUSE6M suprob drgprb alcprb ownhouse nevermarried married 
      fulltime parttime LIVSP JAIL30D SATXIMP ARST30D MEDINS HS DGTXJLC ALCDYS FINSUP1 condom_intake;
 weight wt_cont;
run;

*outcome analysis;
proc genmod data=cjdats.rrr descending;
 model condomyn= RRA_30 /link=logit dist=binomial;
 weight wt_cont;
run;
quit;
