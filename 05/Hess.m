function h=Hess(s)
%通过syms 符号函数的 hessian函数直接得到
x=s(1);
m=s(2);
n=s(3);
q=s(4);
h=[                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        (2*((167*m)/1000 + 4019228480247545/144115188075855872)^2)/((167*n)/1000 + q + 4019228480247545/144115188075855872)^2 + (2*((357*m)/5000 + 2938773856812761/576460752303423488)^2)/((357*n)/5000 + q + 2938773856812761/576460752303423488)^2 + (2*(m + 1)^2)/(n + q + 1)^2 + (2*(2*m + 4)^2)/(2*n + q + 4)^2 + (2*(m/2 + 1/4)^2)/(n/2 + q + 1/4)^2 + (2*(4*m + 16)^2)/(4*n + q + 16)^2 + (2*(m/4 + 1/16)^2)/(n/4 + q + 1/16)^2 + (2*(m/8 + 1/64)^2)/(n/8 + q + 1/64)^2 + (2*(m/10 + 1/100)^2)/(n/10 + q + 1/100)^2 + (2*(m/16 + 1/256)^2)/(n/16 + q + 1/256)^2 + (2*((833*m)/10000 + 1999998874775351/288230376151711744)^2)/((833*n)/10000 + q + 1999998874775351/288230376151711744)^2,                                                                                                                                       (167*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000))/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)) + (357*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000))/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)) + (2*((x*(m + 1))/(n + q + 1) - 347/2000))/(n + q + 1) + ((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25)/(n/2 + q + 1/4) + ((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250)/(4*(n/8 + q + 1/64)) + ((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500)/(2*(n/4 + q + 1/16)) + ((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000)/(5*(n/10 + q + 1/100)) + ((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000)/(8*(n/16 + q + 1/256)) + (4*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4) + (8*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16) + (833*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)) + (4*x*(2*m + 4))/(2*n + q + 4)^2 + (x*(m/2 + 1/4))/(n/2 + q + 1/4)^2 + (8*x*(4*m + 16))/(4*n + q + 16)^2 + (x*(m/4 + 1/16))/(2*(n/4 + q + 1/16)^2) + (x*(m/8 + 1/64))/(4*(n/8 + q + 1/64)^2) + (x*(m/10 + 1/100))/(5*(n/10 + q + 1/100)^2) + (x*(m/16 + 1/256))/(8*(n/16 + q + 1/256)^2) + (833*x*((833*m)/10000 + 1999998874775351/288230376151711744))/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^2) + (167*x*((167*m)/1000 + 4019228480247545/144115188075855872))/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)^2) + (357*x*((357*m)/5000 + 2938773856812761/576460752303423488))/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)^2) + (2*x*(m + 1))/(n + q + 1)^2,                                                                                              - (833*((833*m)/10000 + 1999998874775351/288230376151711744)*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^2) - (167*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000)*((167*m)/1000 + 4019228480247545/144115188075855872))/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)^2) - (357*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000)*((357*m)/5000 + 2938773856812761/576460752303423488))/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)^2) - (4*x*(2*m + 4)^2)/(2*n + q + 4)^3 - (x*(m/2 + 1/4)^2)/(n/2 + q + 1/4)^3 - (8*x*(4*m + 16)^2)/(4*n + q + 16)^3 - (x*(m/4 + 1/16)^2)/(2*(n/4 + q + 1/16)^3) - (x*(m/8 + 1/64)^2)/(4*(n/8 + q + 1/64)^3) - (x*(m/10 + 1/100)^2)/(5*(n/10 + q + 1/100)^3) - (x*(m/16 + 1/256)^2)/(8*(n/16 + q + 1/256)^3) - (2*((x*(m + 1))/(n + q + 1) - 347/2000)*(m + 1))/(n + q + 1)^2 - (833*x*((833*m)/10000 + 1999998874775351/288230376151711744)^2)/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^3) - (167*x*((167*m)/1000 + 4019228480247545/144115188075855872)^2)/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)^3) - (357*x*((357*m)/5000 + 2938773856812761/576460752303423488)^2)/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)^3) - (2*x*(m + 1)^2)/(n + q + 1)^3 - ((m/2 + 1/4)*((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25))/(n/2 + q + 1/4)^2 - ((m/8 + 1/64)*((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250))/(4*(n/8 + q + 1/64)^2) - ((m/4 + 1/16)*((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500))/(2*(n/4 + q + 1/16)^2) - ((m/10 + 1/100)*((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000))/(5*(n/10 + q + 1/100)^2) - ((m/16 + 1/256)*((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000))/(8*(n/16 + q + 1/256)^2) - (4*(2*m + 4)*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4)^2 - (8*(4*m + 16)*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16)^2,                                                                                                            - (2*((833*m)/10000 + 1999998874775351/288230376151711744)*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/((833*n)/10000 + q + 1999998874775351/288230376151711744)^2 - (2*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000)*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872)^2 - (2*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000)*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488)^2 - (2*x*(2*m + 4)^2)/(2*n + q + 4)^3 - (2*x*(m/2 + 1/4)^2)/(n/2 + q + 1/4)^3 - (2*x*(4*m + 16)^2)/(4*n + q + 16)^3 - (2*x*(m/4 + 1/16)^2)/(n/4 + q + 1/16)^3 - (2*x*(m/8 + 1/64)^2)/(n/8 + q + 1/64)^3 - (2*x*(m/10 + 1/100)^2)/(n/10 + q + 1/100)^3 - (2*x*(m/16 + 1/256)^2)/(n/16 + q + 1/256)^3 - (2*((x*(m + 1))/(n + q + 1) - 347/2000)*(m + 1))/(n + q + 1)^2 - (2*x*((833*m)/10000 + 1999998874775351/288230376151711744)^2)/((833*n)/10000 + q + 1999998874775351/288230376151711744)^3 - (2*x*((167*m)/1000 + 4019228480247545/144115188075855872)^2)/((167*n)/1000 + q + 4019228480247545/144115188075855872)^3 - (2*x*((357*m)/5000 + 2938773856812761/576460752303423488)^2)/((357*n)/5000 + q + 2938773856812761/576460752303423488)^3 - (2*x*(m + 1)^2)/(n + q + 1)^3 - (2*(m/2 + 1/4)*((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25))/(n/2 + q + 1/4)^2 - (2*(m/8 + 1/64)*((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250))/(n/8 + q + 1/64)^2 - (2*(m/4 + 1/16)*((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500))/(n/4 + q + 1/16)^2 - (2*(m/10 + 1/100)*((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000))/(n/10 + q + 1/100)^2 - (2*(m/16 + 1/256)*((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000))/(n/16 + q + 1/256)^2 - (2*(2*m + 4)*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4)^2 - (2*(4*m + 16)*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16)^2;
(167*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000))/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)) + (357*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000))/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)) + (2*((x*(m + 1))/(n + q + 1) - 347/2000))/(n + q + 1) + ((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25)/(n/2 + q + 1/4) + ((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250)/(4*(n/8 + q + 1/64)) + ((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500)/(2*(n/4 + q + 1/16)) + ((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000)/(5*(n/10 + q + 1/100)) + ((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000)/(8*(n/16 + q + 1/256)) + (4*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4) + (8*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16) + (833*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)) + (4*x*(2*m + 4))/(2*n + q + 4)^2 + (x*(m/2 + 1/4))/(n/2 + q + 1/4)^2 + (8*x*(4*m + 16))/(4*n + q + 16)^2 + (x*(m/4 + 1/16))/(2*(n/4 + q + 1/16)^2) + (x*(m/8 + 1/64))/(4*(n/8 + q + 1/64)^2) + (x*(m/10 + 1/100))/(5*(n/10 + q + 1/100)^2) + (x*(m/16 + 1/256))/(8*(n/16 + q + 1/256)^2) + (833*x*((833*m)/10000 + 1999998874775351/288230376151711744))/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^2) + (167*x*((167*m)/1000 + 4019228480247545/144115188075855872))/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)^2) + (357*x*((357*m)/5000 + 2938773856812761/576460752303423488))/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)^2) + (2*x*(m + 1))/(n + q + 1)^2,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              (8*x^2)/(2*n + q + 4)^2 + x^2/(2*(n/2 + q + 1/4)^2) + (32*x^2)/(4*n + q + 16)^2 + x^2/(8*(n/4 + q + 1/16)^2) + x^2/(32*(n/8 + q + 1/64)^2) + x^2/(50*(n/10 + q + 1/100)^2) + x^2/(128*(n/16 + q + 1/256)^2) + (693889*x^2)/(50000000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^2) + (27889*x^2)/(500000*((167*n)/1000 + q + 4019228480247545/144115188075855872)^2) + (127449*x^2)/(12500000*((357*n)/5000 + q + 2938773856812761/576460752303423488)^2) + (2*x^2)/(n + q + 1)^2,                                                                                                                                                                                                                                                                                 - (127449*x*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000))/(12500000*((357*n)/5000 + q + 2938773856812761/576460752303423488)^2) - (2*x*((x*(m + 1))/(n + q + 1) - 347/2000))/(n + q + 1)^2 - (x*((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25))/(2*(n/2 + q + 1/4)^2) - (x*((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250))/(32*(n/8 + q + 1/64)^2) - (x*((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500))/(8*(n/4 + q + 1/16)^2) - (x*((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000))/(50*(n/10 + q + 1/100)^2) - (x*((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000))/(128*(n/16 + q + 1/256)^2) - (8*x*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4)^2 - (32*x*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16)^2 - (693889*x*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/(50000000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^2) - (8*x^2*(2*m + 4))/(2*n + q + 4)^3 - (x^2*(m/2 + 1/4))/(2*(n/2 + q + 1/4)^3) - (32*x^2*(4*m + 16))/(4*n + q + 16)^3 - (x^2*(m/4 + 1/16))/(8*(n/4 + q + 1/16)^3) - (x^2*(m/8 + 1/64))/(32*(n/8 + q + 1/64)^3) - (x^2*(m/10 + 1/100))/(50*(n/10 + q + 1/100)^3) - (x^2*(m/16 + 1/256))/(128*(n/16 + q + 1/256)^3) - (693889*x^2*((833*m)/10000 + 1999998874775351/288230376151711744))/(50000000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^3) - (27889*x^2*((167*m)/1000 + 4019228480247545/144115188075855872))/(500000*((167*n)/1000 + q + 4019228480247545/144115188075855872)^3) - (127449*x^2*((357*m)/5000 + 2938773856812761/576460752303423488))/(12500000*((357*n)/5000 + q + 2938773856812761/576460752303423488)^3) - (2*x^2*(m + 1))/(n + q + 1)^3 - (27889*x*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000))/(500000*((167*n)/1000 + q + 4019228480247545/144115188075855872)^2),                                                                                                                                                                                                                                                                                       - (357*x*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000))/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)^2) - (2*x*((x*(m + 1))/(n + q + 1) - 347/2000))/(n + q + 1)^2 - (x*((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25))/(n/2 + q + 1/4)^2 - (x*((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250))/(4*(n/8 + q + 1/64)^2) - (x*((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500))/(2*(n/4 + q + 1/16)^2) - (x*((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000))/(5*(n/10 + q + 1/100)^2) - (x*((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000))/(8*(n/16 + q + 1/256)^2) - (4*x*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4)^2 - (8*x*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16)^2 - (833*x*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^2) - (4*x^2*(2*m + 4))/(2*n + q + 4)^3 - (x^2*(m/2 + 1/4))/(n/2 + q + 1/4)^3 - (8*x^2*(4*m + 16))/(4*n + q + 16)^3 - (x^2*(m/4 + 1/16))/(2*(n/4 + q + 1/16)^3) - (x^2*(m/8 + 1/64))/(4*(n/8 + q + 1/64)^3) - (x^2*(m/10 + 1/100))/(5*(n/10 + q + 1/100)^3) - (x^2*(m/16 + 1/256))/(8*(n/16 + q + 1/256)^3) - (833*x^2*((833*m)/10000 + 1999998874775351/288230376151711744))/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^3) - (167*x^2*((167*m)/1000 + 4019228480247545/144115188075855872))/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)^3) - (357*x^2*((357*m)/5000 + 2938773856812761/576460752303423488))/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)^3) - (2*x^2*(m + 1))/(n + q + 1)^3 - (167*x*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000))/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)^2);
 - (833*((833*m)/10000 + 1999998874775351/288230376151711744)*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^2) - (167*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000)*((167*m)/1000 + 4019228480247545/144115188075855872))/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)^2) - (357*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000)*((357*m)/5000 + 2938773856812761/576460752303423488))/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)^2) - (4*x*(2*m + 4)^2)/(2*n + q + 4)^3 - (x*(m/2 + 1/4)^2)/(n/2 + q + 1/4)^3 - (8*x*(4*m + 16)^2)/(4*n + q + 16)^3 - (x*(m/4 + 1/16)^2)/(2*(n/4 + q + 1/16)^3) - (x*(m/8 + 1/64)^2)/(4*(n/8 + q + 1/64)^3) - (x*(m/10 + 1/100)^2)/(5*(n/10 + q + 1/100)^3) - (x*(m/16 + 1/256)^2)/(8*(n/16 + q + 1/256)^3) - (2*((x*(m + 1))/(n + q + 1) - 347/2000)*(m + 1))/(n + q + 1)^2 - (833*x*((833*m)/10000 + 1999998874775351/288230376151711744)^2)/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^3) - (167*x*((167*m)/1000 + 4019228480247545/144115188075855872)^2)/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)^3) - (357*x*((357*m)/5000 + 2938773856812761/576460752303423488)^2)/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)^3) - (2*x*(m + 1)^2)/(n + q + 1)^3 - ((m/2 + 1/4)*((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25))/(n/2 + q + 1/4)^2 - ((m/8 + 1/64)*((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250))/(4*(n/8 + q + 1/64)^2) - ((m/4 + 1/16)*((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500))/(2*(n/4 + q + 1/16)^2) - ((m/10 + 1/100)*((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000))/(5*(n/10 + q + 1/100)^2) - ((m/16 + 1/256)*((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000))/(8*(n/16 + q + 1/256)^2) - (4*(2*m + 4)*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4)^2 - (8*(4*m + 16)*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16)^2, - (127449*x*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000))/(12500000*((357*n)/5000 + q + 2938773856812761/576460752303423488)^2) - (2*x*((x*(m + 1))/(n + q + 1) - 347/2000))/(n + q + 1)^2 - (x*((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25))/(2*(n/2 + q + 1/4)^2) - (x*((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250))/(32*(n/8 + q + 1/64)^2) - (x*((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500))/(8*(n/4 + q + 1/16)^2) - (x*((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000))/(50*(n/10 + q + 1/100)^2) - (x*((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000))/(128*(n/16 + q + 1/256)^2) - (8*x*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4)^2 - (32*x*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16)^2 - (693889*x*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/(50000000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^2) - (8*x^2*(2*m + 4))/(2*n + q + 4)^3 - (x^2*(m/2 + 1/4))/(2*(n/2 + q + 1/4)^3) - (32*x^2*(4*m + 16))/(4*n + q + 16)^3 - (x^2*(m/4 + 1/16))/(8*(n/4 + q + 1/16)^3) - (x^2*(m/8 + 1/64))/(32*(n/8 + q + 1/64)^3) - (x^2*(m/10 + 1/100))/(50*(n/10 + q + 1/100)^3) - (x^2*(m/16 + 1/256))/(128*(n/16 + q + 1/256)^3) - (693889*x^2*((833*m)/10000 + 1999998874775351/288230376151711744))/(50000000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^3) - (27889*x^2*((167*m)/1000 + 4019228480247545/144115188075855872))/(500000*((167*n)/1000 + q + 4019228480247545/144115188075855872)^3) - (127449*x^2*((357*m)/5000 + 2938773856812761/576460752303423488))/(12500000*((357*n)/5000 + q + 2938773856812761/576460752303423488)^3) - (2*x^2*(m + 1))/(n + q + 1)^3 - (27889*x*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000))/(500000*((167*n)/1000 + q + 4019228480247545/144115188075855872)^2), (2*x^2*(m + 1)^2)/(n + q + 1)^4 + (8*x^2*(2*m + 4)^2)/(2*n + q + 4)^4 + (x^2*(m/2 + 1/4)^2)/(2*(n/2 + q + 1/4)^4) + (32*x^2*(4*m + 16)^2)/(4*n + q + 16)^4 + (x^2*(m/4 + 1/16)^2)/(8*(n/4 + q + 1/16)^4) + (x^2*(m/8 + 1/64)^2)/(32*(n/8 + q + 1/64)^4) + (x^2*(m/10 + 1/100)^2)/(50*(n/10 + q + 1/100)^4) + (x^2*(m/16 + 1/256)^2)/(128*(n/16 + q + 1/256)^4) + (693889*x^2*((833*m)/10000 + 1999998874775351/288230376151711744)^2)/(50000000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^4) + (27889*x^2*((167*m)/1000 + 4019228480247545/144115188075855872)^2)/(500000*((167*n)/1000 + q + 4019228480247545/144115188075855872)^4) + (127449*x^2*((357*m)/5000 + 2938773856812761/576460752303423488)^2)/(12500000*((357*n)/5000 + q + 2938773856812761/576460752303423488)^4) + (x*(m/2 + 1/4)*((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25))/(n/2 + q + 1/4)^3 + (x*(m/8 + 1/64)*((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250))/(16*(n/8 + q + 1/64)^3) + (x*(m/4 + 1/16)*((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500))/(4*(n/4 + q + 1/16)^3) + (x*(m/10 + 1/100)*((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000))/(25*(n/10 + q + 1/100)^3) + (x*(m/16 + 1/256)*((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000))/(64*(n/16 + q + 1/256)^3) + (16*x*(2*m + 4)*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4)^3 + (64*x*(4*m + 16)*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16)^3 + (693889*x*((833*m)/10000 + 1999998874775351/288230376151711744)*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/(25000000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^3) + (27889*x*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000)*((167*m)/1000 + 4019228480247545/144115188075855872))/(250000*((167*n)/1000 + q + 4019228480247545/144115188075855872)^3) + (127449*x*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000)*((357*m)/5000 + 2938773856812761/576460752303423488))/(6250000*((357*n)/5000 + q + 2938773856812761/576460752303423488)^3) + (4*x*((x*(m + 1))/(n + q + 1) - 347/2000)*(m + 1))/(n + q + 1)^3, (2*x^2*(m + 1)^2)/(n + q + 1)^4 + (4*x^2*(2*m + 4)^2)/(2*n + q + 4)^4 + (x^2*(m/2 + 1/4)^2)/(n/2 + q + 1/4)^4 + (8*x^2*(4*m + 16)^2)/(4*n + q + 16)^4 + (x^2*(m/4 + 1/16)^2)/(2*(n/4 + q + 1/16)^4) + (x^2*(m/8 + 1/64)^2)/(4*(n/8 + q + 1/64)^4) + (x^2*(m/10 + 1/100)^2)/(5*(n/10 + q + 1/100)^4) + (x^2*(m/16 + 1/256)^2)/(8*(n/16 + q + 1/256)^4) + (833*x^2*((833*m)/10000 + 1999998874775351/288230376151711744)^2)/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^4) + (167*x^2*((167*m)/1000 + 4019228480247545/144115188075855872)^2)/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)^4) + (357*x^2*((357*m)/5000 + 2938773856812761/576460752303423488)^2)/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)^4) + (2*x*(m/2 + 1/4)*((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25))/(n/2 + q + 1/4)^3 + (x*(m/8 + 1/64)*((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250))/(2*(n/8 + q + 1/64)^3) + (x*(m/4 + 1/16)*((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500))/(n/4 + q + 1/16)^3 + (2*x*(m/10 + 1/100)*((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000))/(5*(n/10 + q + 1/100)^3) + (x*(m/16 + 1/256)*((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000))/(4*(n/16 + q + 1/256)^3) + (8*x*(2*m + 4)*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4)^3 + (16*x*(4*m + 16)*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16)^3 + (833*x*((833*m)/10000 + 1999998874775351/288230376151711744)*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/(2500*((833*n)/10000 + q + 1999998874775351/288230376151711744)^3) + (167*x*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000)*((167*m)/1000 + 4019228480247545/144115188075855872))/(250*((167*n)/1000 + q + 4019228480247545/144115188075855872)^3) + (357*x*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000)*((357*m)/5000 + 2938773856812761/576460752303423488))/(1250*((357*n)/5000 + q + 2938773856812761/576460752303423488)^3) + (4*x*((x*(m + 1))/(n + q + 1) - 347/2000)*(m + 1))/(n + q + 1)^3;
 - (2*((833*m)/10000 + 1999998874775351/288230376151711744)*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/((833*n)/10000 + q + 1999998874775351/288230376151711744)^2 - (2*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000)*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872)^2 - (2*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000)*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488)^2 - (2*x*(2*m + 4)^2)/(2*n + q + 4)^3 - (2*x*(m/2 + 1/4)^2)/(n/2 + q + 1/4)^3 - (2*x*(4*m + 16)^2)/(4*n + q + 16)^3 - (2*x*(m/4 + 1/16)^2)/(n/4 + q + 1/16)^3 - (2*x*(m/8 + 1/64)^2)/(n/8 + q + 1/64)^3 - (2*x*(m/10 + 1/100)^2)/(n/10 + q + 1/100)^3 - (2*x*(m/16 + 1/256)^2)/(n/16 + q + 1/256)^3 - (2*((x*(m + 1))/(n + q + 1) - 347/2000)*(m + 1))/(n + q + 1)^2 - (2*x*((833*m)/10000 + 1999998874775351/288230376151711744)^2)/((833*n)/10000 + q + 1999998874775351/288230376151711744)^3 - (2*x*((167*m)/1000 + 4019228480247545/144115188075855872)^2)/((167*n)/1000 + q + 4019228480247545/144115188075855872)^3 - (2*x*((357*m)/5000 + 2938773856812761/576460752303423488)^2)/((357*n)/5000 + q + 2938773856812761/576460752303423488)^3 - (2*x*(m + 1)^2)/(n + q + 1)^3 - (2*(m/2 + 1/4)*((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25))/(n/2 + q + 1/4)^2 - (2*(m/8 + 1/64)*((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250))/(n/8 + q + 1/64)^2 - (2*(m/4 + 1/16)*((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500))/(n/4 + q + 1/16)^2 - (2*(m/10 + 1/100)*((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000))/(n/10 + q + 1/100)^2 - (2*(m/16 + 1/256)*((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000))/(n/16 + q + 1/256)^2 - (2*(2*m + 4)*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4)^2 - (2*(4*m + 16)*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16)^2,                                                         - (357*x*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000))/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)^2) - (2*x*((x*(m + 1))/(n + q + 1) - 347/2000))/(n + q + 1)^2 - (x*((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25))/(n/2 + q + 1/4)^2 - (x*((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250))/(4*(n/8 + q + 1/64)^2) - (x*((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500))/(2*(n/4 + q + 1/16)^2) - (x*((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000))/(5*(n/10 + q + 1/100)^2) - (x*((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000))/(8*(n/16 + q + 1/256)^2) - (4*x*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4)^2 - (8*x*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16)^2 - (833*x*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^2) - (4*x^2*(2*m + 4))/(2*n + q + 4)^3 - (x^2*(m/2 + 1/4))/(n/2 + q + 1/4)^3 - (8*x^2*(4*m + 16))/(4*n + q + 16)^3 - (x^2*(m/4 + 1/16))/(2*(n/4 + q + 1/16)^3) - (x^2*(m/8 + 1/64))/(4*(n/8 + q + 1/64)^3) - (x^2*(m/10 + 1/100))/(5*(n/10 + q + 1/100)^3) - (x^2*(m/16 + 1/256))/(8*(n/16 + q + 1/256)^3) - (833*x^2*((833*m)/10000 + 1999998874775351/288230376151711744))/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^3) - (167*x^2*((167*m)/1000 + 4019228480247545/144115188075855872))/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)^3) - (357*x^2*((357*m)/5000 + 2938773856812761/576460752303423488))/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)^3) - (2*x^2*(m + 1))/(n + q + 1)^3 - (167*x*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000))/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)^2),                                                   (2*x^2*(m + 1)^2)/(n + q + 1)^4 + (4*x^2*(2*m + 4)^2)/(2*n + q + 4)^4 + (x^2*(m/2 + 1/4)^2)/(n/2 + q + 1/4)^4 + (8*x^2*(4*m + 16)^2)/(4*n + q + 16)^4 + (x^2*(m/4 + 1/16)^2)/(2*(n/4 + q + 1/16)^4) + (x^2*(m/8 + 1/64)^2)/(4*(n/8 + q + 1/64)^4) + (x^2*(m/10 + 1/100)^2)/(5*(n/10 + q + 1/100)^4) + (x^2*(m/16 + 1/256)^2)/(8*(n/16 + q + 1/256)^4) + (833*x^2*((833*m)/10000 + 1999998874775351/288230376151711744)^2)/(5000*((833*n)/10000 + q + 1999998874775351/288230376151711744)^4) + (167*x^2*((167*m)/1000 + 4019228480247545/144115188075855872)^2)/(500*((167*n)/1000 + q + 4019228480247545/144115188075855872)^4) + (357*x^2*((357*m)/5000 + 2938773856812761/576460752303423488)^2)/(2500*((357*n)/5000 + q + 2938773856812761/576460752303423488)^4) + (2*x*(m/2 + 1/4)*((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25))/(n/2 + q + 1/4)^3 + (x*(m/8 + 1/64)*((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250))/(2*(n/8 + q + 1/64)^3) + (x*(m/4 + 1/16)*((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500))/(n/4 + q + 1/16)^3 + (2*x*(m/10 + 1/100)*((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000))/(5*(n/10 + q + 1/100)^3) + (x*(m/16 + 1/256)*((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000))/(4*(n/16 + q + 1/256)^3) + (8*x*(2*m + 4)*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4)^3 + (16*x*(4*m + 16)*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16)^3 + (833*x*((833*m)/10000 + 1999998874775351/288230376151711744)*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/(2500*((833*n)/10000 + q + 1999998874775351/288230376151711744)^3) + (167*x*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000)*((167*m)/1000 + 4019228480247545/144115188075855872))/(250*((167*n)/1000 + q + 4019228480247545/144115188075855872)^3) + (357*x*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000)*((357*m)/5000 + 2938773856812761/576460752303423488))/(1250*((357*n)/5000 + q + 2938773856812761/576460752303423488)^3) + (4*x*((x*(m + 1))/(n + q + 1) - 347/2000)*(m + 1))/(n + q + 1)^3,                                                                  (2*x^2*(m + 1)^2)/(n + q + 1)^4 + (2*x^2*(2*m + 4)^2)/(2*n + q + 4)^4 + (2*x^2*(m/2 + 1/4)^2)/(n/2 + q + 1/4)^4 + (2*x^2*(4*m + 16)^2)/(4*n + q + 16)^4 + (2*x^2*(m/4 + 1/16)^2)/(n/4 + q + 1/16)^4 + (2*x^2*(m/8 + 1/64)^2)/(n/8 + q + 1/64)^4 + (2*x^2*(m/10 + 1/100)^2)/(n/10 + q + 1/100)^4 + (2*x^2*(m/16 + 1/256)^2)/(n/16 + q + 1/256)^4 + (2*x^2*((833*m)/10000 + 1999998874775351/288230376151711744)^2)/((833*n)/10000 + q + 1999998874775351/288230376151711744)^4 + (2*x^2*((167*m)/1000 + 4019228480247545/144115188075855872)^2)/((167*n)/1000 + q + 4019228480247545/144115188075855872)^4 + (2*x^2*((357*m)/5000 + 2938773856812761/576460752303423488)^2)/((357*n)/5000 + q + 2938773856812761/576460752303423488)^4 + (4*x*(m/2 + 1/4)*((x*(m/2 + 1/4))/(n/2 + q + 1/4) - 4/25))/(n/2 + q + 1/4)^3 + (4*x*(m/8 + 1/64)*((x*(m/8 + 1/64))/(n/8 + q + 1/64) - 57/1250))/(n/8 + q + 1/64)^3 + (4*x*(m/4 + 1/16)*((x*(m/4 + 1/16))/(n/4 + q + 1/16) - 211/2500))/(n/4 + q + 1/16)^3 + (4*x*(m/10 + 1/100)*((x*(m/10 + 1/100))/(n/10 + q + 1/100) - 171/5000))/(n/10 + q + 1/100)^3 + (4*x*(m/16 + 1/256)*((x*(m/16 + 1/256))/(n/16 + q + 1/256) - 123/5000))/(n/16 + q + 1/256)^3 + (4*x*(2*m + 4)*((x*(2*m + 4))/(2*n + q + 4) - 1947/10000))/(2*n + q + 4)^3 + (4*x*(4*m + 16)*((x*(4*m + 16))/(4*n + q + 16) - 1957/10000))/(4*n + q + 16)^3 + (4*x*((833*m)/10000 + 1999998874775351/288230376151711744)*((x*((833*m)/10000 + 1999998874775351/288230376151711744))/((833*n)/10000 + q + 1999998874775351/288230376151711744) - 323/10000))/((833*n)/10000 + q + 1999998874775351/288230376151711744)^3 + (4*x*((x*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872) - 627/10000)*((167*m)/1000 + 4019228480247545/144115188075855872))/((167*n)/1000 + q + 4019228480247545/144115188075855872)^3 + (4*x*((x*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488) - 47/2000)*((357*m)/5000 + 2938773856812761/576460752303423488))/((357*n)/5000 + q + 2938773856812761/576460752303423488)^3 + (4*x*((x*(m + 1))/(n + q + 1) - 347/2000)*(m + 1))/(n + q + 1)^3];
 