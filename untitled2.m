syms r k y0 y t real;
f=( y0*k/(y0+(k-y0)*exp(-r*t)))
d=y-f;
Jsym=jacobian(d,[r k y0])