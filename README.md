# read me

icassp_anon.pdf: anonymous version of our paper under review for ICASSP 2023

kx_noiseless_gaussian.m: code for strategic scalar quantization with \eta_E = (kx-y)^2, \eta_D = (x-y)^2, with source X ~ N(0,1) over a noiseless channel 

kx_noisy.m: code for strategic scalar quantization with \eta_E = (kx-y)^2, \eta_D = (x-y)^2, with source X ~ N(0,1) over a noisy channel with bit error rate pb

nonstr_xtheta_gaussian_noisy.m: code for non-strategic (X,\theta) quantization with \eta_E = (x+\theta-y)^2, \eta_D = (x-y)^2, with source X ~ N(0,1) over a noisy channel with bit error rate pb

plot_compare.m: code to plot figures comparing non-strategic (\eta_D=\eta_D=(x-y)^2) and strategic decoder distortion (\eta_E=(x+\theta-y)^2,\eta_D=(x-y)^2) for \rho=0.99

plot_figures_1pt5x_noiseless.m: code to plot figures for scalar strategic quantization distortions with \eta_E=(1.5x-y)^2,\eta_D=(x-y)^2 with source X ~ N(0,1) for a 
noiseless channel 

plot_figures_kx_noisy.m: code to plot figures for scalar strategic quantization distortions with \eta_E=(1.5x-y)^2,\eta_D=(x-y)^2 with source X ~ N(0,1) for a 
noisy channel with bit error rate pb

plot_figures_xcubed_noiseless.m: code to plot figures for scalar strategic quantization distortions with \eta_E=(x^3-y)^2,\eta_D=(x-y)^2 with source X ~ N(0,1) for a 
noiseless channel 

plot_figures_xcubed_noisy.m: code to plot figures for scalar strategic quantization distortions with \eta_E=(x^3-y)^2,\eta_D=(x-y)^2 with source X ~ N(0,1) for a 
noisy channel with bit error rate pb

plotfigures_xtheta_noisy_quadratic.m: code to plot figures for (X,\theta) strategic quantization distortions with \eta_E=(x+\theta-y)^2,\eta_D=(x-y)^2 with source 
X ~ N(0,1) for a noisy channel with bit error rate pb

xcubed_noiseless_gaussian.m: code for strategic scalar quantization with \eta_E = (x^3-y)^2, \eta_D = (x-y)^2, with source X ~ N(0,1) over a noiseless channel 

xcubed_noisy.m: code for strategic scalar quantization with \eta_E = (x^3-y)^2, \eta_D = (x-y)^2, with source X ~ N(0,1) over a noisy channel with bit error rate pb

xtheta_noisy_quadratic.m: code for (X,\theta) strategic quantization distortions with \eta_E=(x+\theta-y)^2,\eta_D=(x-y)^2 with source 
X ~ N(0,1) for a noisy channel with bit error rate pb 
