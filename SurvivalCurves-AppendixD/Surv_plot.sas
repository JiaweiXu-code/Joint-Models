
proc format;

   value afmt 
      1 = "(*ESC*){unicode '03B3'x} = (0.0,0.0,0.0,0.0), (*ESC*){unicode '03B1'x} = -0.2"
      2 = "(*ESC*){unicode '03B3'x} = (0.2,0.2,0.2,0.2), (*ESC*){unicode '03B1'x} = -0.2"
      3 = "(*ESC*){unicode '03B3'x} = (0.2,0.2,0.2,0.2), (*ESC*){unicode '03B1'x} = 0.0";

run;




data survp;
infile "\\Client\C$\D\paper1\revision\rev11\surv_curves\survp.csv" delimiter = ',';
input t s0 s1 alpha;
run;


%let root =  \\Client\C$\D\paper1\revision\rev11\surv_curves;


options Papersize=("6in","2.5in") nodate nonumber;
ODS PDF file = "&root.\surv_curves.pdf" dpi=400;
*ods escapechar='^';
ods html close;
ods graphics / height=2.5in width=6in noborder;

proc sgpanel data = survp;

 format alpha afmt.; 
 panelby alpha / layout=COLUMNLATTICE onepanel start=topleft novarname headerattrs=GraphUnicodeText(size=11);
 series x=t y = s1    / lineattrs=(color=black  thickness=1 pattern=1)   legendlabel = 'Treated Group';
 series x=t y = s0    / lineattrs=(color=black  thickness=1 pattern=2)   legendlabel = 'Control Group';
  rowaxis label='Survival Probability' min=0.7 max=1;
  colaxis label='Time (Years)';
  run;
  
ods layout end;
ods pdf close;



