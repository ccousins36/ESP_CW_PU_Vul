# loads the right libraries in python
# these need to be installed in your system

# numpy allows for numerical calculations
import numpy as np
# matplotlib provides plotting capabilities
import matplotlib.pyplot as plt
# mesa_reader allows to get the data from the "./LOGS/" folder files easily
import mesa_reader as mr

# THIS IS TO GET NICE TEXT FOR THE CAPTIONS
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Load the "history.data" file and store it in "h" using Mesa_reader capabilities
h = mr.MesaData('LOGS/history.data')

# Retrieve the star age (in years) and save it in "time" array
time=h.data('star_age');

# You could opearate on this array to, say, transform it into thousands of years if you wanted to 
#time=time/1000;

# Retrieve the luminosity 
lumi=h.data('log_L');
temperature=h.data('log_center_T');
print(temperature)
# Labels for both x- and y-axis
plt.xlabel("Star age, "r"$\displaystyle \tau$ [yr]")
plt.ylabel("Luminosity, $\log L/L_{\odot}$");
#plt.xlim(12000,14000)
# Actual plotting
plt.plot(time, lumi)

# Save it on a pdf that can be reused later 
plt.savefig('luminosity.pdf')

# For on-screen matplotlib plot of data as well as a pdf saved copy, uncomment plt.show() below.
# plt.show()



# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 13:12:28 2021
@author: James
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
plt.rcParams['text.usetex'] = True
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams["font.size"] = 14


with open("PUvul.txt") as f:
    content = f.readlines()

magnitudes=[]
mag_unc=[]
JDs=[]
HJDs=[]
bands=[]
magnitude_type=[]


Vmagnitudes=[]
Vmag_unc=[]
VJDs=[]
VHJDs=[]
Vbands=[]
Vmagnitude_type=[]
test=[]

for i in range(1,len(content)):
    allsplit=content[i].split(",")
    magnitudes.append(allsplit[1])
    mag_unc.append(allsplit[2])
    JDs.append(allsplit[0])
    HJDs.append(allsplit[16])
    bands.append(allsplit[4])
    magnitude_type.append(allsplit[19])
    

    if allsplit[4]=="V" or allsplit[4]=="Vis.":#and allsplit[19]=="STD" and allsplit[1]!="" and  allsplit[0]!="":
        try:
            Vmagnitudes.append(float(allsplit[1]))
            VJDs.append(float(allsplit[0]))
        except:
            None
        try:
            
            Vmag_unc.append(float(allsplit[2]))
            
            Vbands.append(allsplit[4])
            Vmagnitude_type.append(allsplit[19])
            
            test.append([allsplit[1], allsplit[0]])
            
            
        except:
            None


### CONVERT MAGS TO LUMINOSITIES
d=3500 #distance in parsecs to PU Vul
Rv=3.1 # interstellar extinction coefficient, quoted to be 3.1+-0.2 but everyone uses 3.1, Fitzpatrick (late 1990s/2000s can't remember)
Ebv=0.4 # Hachisu & Kato 2019 for V1324 Sco, 1.32+-0.1
A=Rv*Ebv # Extinction coefficient

M_x_sol_V=4.8 #FOR THE V BAND, https://iopscience.iop.org/article/10.3847/1538-4365/aabfdf/pdf
L_V = np.log10((d/10)**2   /    (10**((np.asarray(Vmagnitudes)-M_x_sol_V-  A  )/2.5)))
print(len(L_V))
print(len(VJDs))

x_axis_scale=365.25 #change this to 1 for JD
x_axis_obs_offset=2230.2 #offset experimental data to line up with mesa plots

plt.scatter(x_axis_obs_offset+np.asarray(VJDs)/x_axis_scale, L_V, facecolors='orange', marker="x", edgecolors="orange", s=15, label="Vdata")
plt.xlabel("JD/years")
plt.ylabel(r"Log(Luminosity) (L$_\odot$)")


f = interpolate.interp1d(np.asarray(VJDs)/x_axis_scale + +x_axis_obs_offset, L_V,kind="linear")
all_arr_x=np.linspace(min(np.asarray(VJDs)/x_axis_scale +x_axis_obs_offset),max(np.asarray(VJDs)/x_axis_scale + x_axis_obs_offset),100000)
fit=f(all_arr_x)
plt.plot(all_arr_x,fit, linestyle="--", c="black")
plt.xlim(8920,8970) #adjust limit to fit burst within range
plt.savefig("luminosity_with_obs.pdf")



from scipy.interpolate import interp1d
from scipy import arange, array, exp

def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(list(map(pointwise, array(xs))))

    return ufunclike



mesa_interp_fit=interpolate.interp1d(np.asarray(time), np.asarray(lumi), kind="linear")
fit_mesa_extrapolated=extrap1d(mesa_interp_fit)

fit_mesa=fit_mesa_extrapolated(all_arr_x)

list_fit_for_chisq=[]
list_fit_mesa_for_chisq=[]
all_arr_x_filtered=[]
for i in range(0, len(all_arr_x)):
    if all_arr_x[i] > 8927 and all_arr_x[i] < 9999999999999999: # set lower limit to time of outburst in years
        list_fit_for_chisq.append(fit[i])
        list_fit_mesa_for_chisq.append(fit_mesa[i])
        all_arr_x_filtered.append(all_arr_x[i])
		
fit_for_chisq=np.asarray(list_fit_for_chisq)
fit_mesa_for_chisq=np.asarray(list_fit_mesa_for_chisq)
#equation for chi squared
chisq=np.sum((fit_mesa_for_chisq-fit_for_chisq)**2/(fit_for_chisq)) / len(list_fit_mesa_for_chisq)
print(chisq)
print(len(list_fit_mesa_for_chisq))



plt.plot(all_arr_x_filtered,fit_mesa_for_chisq)
plt.show()

