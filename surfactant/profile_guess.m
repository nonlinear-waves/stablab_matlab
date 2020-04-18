function out = profile_guess(x,s,p)

% We are using a guess coming from the esp=0 case.

[frontnn, frontpp] = ep0integrated(p.hL,p.hR,p.D,p.be,s.I+.1,p.maxG);

tempp=frontpp(x);
tempn=frontnn(-x);
tempnd=(-frontnn(-x-0.01)+frontnn(-x))/0.01;
temppd=(frontpp(x+0.01)-frontpp(x))/0.01;
temppdd=(frontpp(x+2*0.01)-2*frontpp(x+0.01)+frontpp(x))/0.01;
tempndd=(frontnn(-x-2*0.01)-2*frontnn(-x-0.01)+frontnn(-x))/0.01;
                  
out=[tempp(1),temppd(1),temppdd(1),tempp(2),tempn(1),tempnd(1),tempndd(1),tempn(2)]; 