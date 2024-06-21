#%%
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

hep.styles.cms.CMS["lines.linewidth"] = 3

hep.styles.cms.CMS["legend.frameon"] = True
# hep.styles.cms.CMS["figure.autolayout"] = True


hep.style.use("CMS")

x = [
    1.44571e-09, 1.53714e-09, 1.62857e-09, 1.72e-09, 1.81143e-09, 1.90286e-09, 1.99429e-09, 2.08571e-09, 2.17714e-09, 2.26857e-09,
    2.36e-09, 2.45143e-09, 2.54286e-09, 2.63429e-09, 2.72571e-09, 2.81714e-09, 2.90857e-09, 3e-09, 3.09143e-09, 3.18286e-09,
    3.27429e-09, 3.36571e-09, 3.45714e-09, 3.54857e-09, 3.64e-09, 3.73143e-09, 3.82286e-09, 3.91429e-09, 4.00571e-09, 4.09714e-09,
    4.18857e-09, 4.28e-09, 4.37143e-09, 4.46286e-09, 4.55429e-09
]

NLL = [
    24.5725, 21.0262, 17.8138, 14.9176, 12.3217, 10.0115, 7.97303, 6.1936, 4.66124, 3.36479, 2.29382, 1.43856, 0.789864, 0.339138,
    0.0783152, 0., 0.0964723, 0.361575, 0.788763, 1.37282, 2.10623, 2.98496, 4.00363, 5.15747, 6.44186, 7.8525, 9.38522, 11.0361,
    12.8015, 14.6777, 16.6614, 18.7494, 20.9385, 23.2258, 25.6084
]

# Convert lists to numpy arrays
x = np.array(x)
NLL = np.array(NLL)

fig,ax=plt.subplots()

ax.plot(x*1e9,NLL)
ax.axhline(1, color='red', linestyle='--')
ax.grid()
ax.set_ylim(0,26)
ax.set_xlabel("$BR( B_{s} \\to \mu^+ \mu^-) \\times 10^{9}$")
ax.set_ylabel("-2$\Delta$ln(L)")

argmin=np.argmin(NLL)
best_fit=x[argmin]*1e9


sm=3.66
dsm=0.14
ax.axvline(sm, color='purple',label='SM prediction')
ax.axvspan(sm-dsm,sm+dsm, color='purple', alpha=0.3)

ax.axvline(best_fit, color='limegreen',label='Best fit')


upper_idx=np.argmin(np.abs(NLL[NLL>2.8]-1))
lower_idx=np.argmin(np.abs(NLL[NLL<2.8]-1))
upper_lim=x[NLL>2.8][upper_idx]*1e9-best_fit
lower_lim=x[NLL<2.8][lower_idx]*1e9-best_fit

SM_idx=np.argmin(np.abs(x-3.66e-09))
SM_x=x[SM_idx]*1e9
NLL_SM=NLL[SM_idx]

#find the number of signa using chi square cum function
from scipy.stats import chi2,norm
p_value=1-chi2.cdf(NLL_SM,df=1)
sigma = norm.ppf(1 - p_value/2)


arr_len=SM_x-best_fit
plt.arrow(best_fit, 15, arr_len, 0.0, color='red',
          head_length = 0.05, head_width = 0.4,
          length_includes_head = True)

plt.arrow(SM_x, 15, -arr_len, 0.0, color='red',
          head_length = 0.05, head_width = 0.4,
          length_includes_head = True)

ax.text(best_fit+arr_len/2-0.13, 15.5, f"{sigma:.2f} $\sigma$", color='red', fontsize=18)


ax.legend()

hep.cms.text("Private Work",ax=ax)
hep.cms.lumitext(f"${best_fit:.2f}^{{+{upper_lim:.2f}}}_{{{lower_lim:.2f}}}$",ax=ax)

fig.savefig("nll_plot.pdf")
fig.savefig("nll_plot.png")

#%%
