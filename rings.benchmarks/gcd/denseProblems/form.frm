#-
Off statistics;

Symbols x1, x2, x3, x4, x5, x6, x7;

#$exp = 7;
#$a = (1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 + 15*x7)^$exp - 1;
#$b = (1 - 3*x1 - 5*x2 - 7*x3 + 9*x4 - 11*x5 - 13*x6 + 15*x7)^$exp + 1;
#$g = (1 + 3*x1 + 5*x2 + 7*x3 + 9*x4 + 11*x5 + 13*x6 - 15*x7)^$exp + 3;

#$ag = $a*$g;
#$bg = $b*$g;
#$agt = $ag + 1;

#$tstart = `timer_';
#$gcd = gcd_($ag, $bg);
#$tend = `timer_';

#$dt = ($tend - $tstart);
#write "Done in Z: %$" , $dt


#$tstart = `timer_';
#$gcd = gcd_(#$agt, $bg);
#$tend = `timer_';

#$dt = ($tend - $tstart);
#write "Done in Z (trivial): %$" , $dt

.end