pro main

c = 3E8
H0 = 2.1916E-18
meters_to_kpc = double(3.086E19)
l_ratio = [0.001,0.005,0.01,0.03,0.06,0.1,0.2,0.5,0.6,0.7,1.0,1.5]

f_R = 0.84D
R_star = double((75.0)*meters_to_kpc)
beta = 0.23D
alpha = -1.67D
phi_star = double((0.000895)/((1000.0*meters_to_kpc)^3))

gamma_arg = alpha+(2*beta)+1
upper_gamma = GAMMA(gamma_arg)*(1 - IGAMMA(gamma_arg, l_ratio))
dNdX_nominal = (c*!PI/H0)*f_R*R_star*R_star*phi_star*upper_gamma

cgplot, l_ratio, dNdX_nominal,thick=2, xrange=[0.02,1.5],yrange=[0.0008,3],/YLOG, /XLOG, /Window



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



f_R = 0.82D
R_star = double((75.0-3.0)*meters_to_kpc)
beta = 0.24D
alpha = -1.65D
phi_star = double((0.000895-0.00006)/((1000.0*meters_to_kpc)^3))
gamma_arg = alpha+(2*beta)+1
upper_gamma = GAMMA(gamma_arg)*(1 - IGAMMA(gamma_arg, l_ratio))
dNdX_low = (c*!PI/H0)*f_R*R_star*R_star*phi_star*upper_gamma


f_R = 0.86D
R_star = double((75.0+20.0)*meters_to_kpc)
beta = 0.22D
alpha = -1.73D
phi_star = double((0.000895+0.0001)/((1000.0*meters_to_kpc)^3))
gamma_arg = alpha+(2*beta)+1
upper_gamma = GAMMA(gamma_arg)*(1 - IGAMMA(gamma_arg, l_ratio))
dNdX_high = (c*!PI/H0)*f_R*R_star*R_star*phi_star*upper_gamma



cgplot, l_ratio, dNdX_nominal,thick=2, xrange=[0.02,1.5],yrange=[0.0008,30],xstyle=1,ytitle='dN/dX',xtitle='LB/L*',/YLOG, /XLOG, /Window
cgplot, l_ratio, dNdX_low, thick=2,color=cgcolor('blue'),/Window, /Overplot
cgplot, l_ratio, dNdX_high, thick=2,color=cgcolor('blue'),/Window, /Overplot
cgplot, l_ratio, intarr(n_elements(l_ratio))+0.8-0.18, linestyle=2, /Window, /Overplot
cgplot, l_ratio, intarr(n_elements(l_ratio))+0.86+0.19, linestyle=2, /Window, /Overplot






































;;; at Lb.L_star = 1:

dNdX = (c*!PI/H0)*f_R*R_star*R_star*phi_star*exp(-gamma_arg)

print, dNdX


;;; different base units


c = 3E5
H0 = 67.74

f_R = 0.84D
R_star = double(75.0/1000.0)
beta = 0.23D
alpha = -1.76D
phi_star = double(0.79E-3)

l_ratio = [0.03,0.06,0.1,0.2,0.5,0.7,1.0]

gamma_arg = alpha+(2*beta)+1

upper_gamma = GAMMA(gamma_arg)*(1 - IGAMMA(gamma_arg, l_ratio))
dNdX = (c*!PI/H0)*f_R*R_star*R_star*phi_star*upper_gamma

cgplot, l_ratio, dNdX, xrange=[0.025,1.5],yrange=[0.0008,3],/YLOG, /XLOG, /Window






