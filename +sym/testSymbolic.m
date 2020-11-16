%% UNIT TESTS FOR THE pkg.sym LIBRARY

a = pkg.sym.Variable('a') ;
f = pkg.sym.Cos(a) + pkg.sym.Sin(a)

clf 
x = linspace(0,1)*2*pi ;
plot(x,cos(x)+sin(x)) ;
plot(x,f.evalAt(x)) ;

plot(x,-sin(x)+cos(x))
plot(x,f.diff.evalAt(x))





