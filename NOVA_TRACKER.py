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
h=mr.MesaData('LOGS/history.data')

# Retrieve the star age (in years) and save it in "time" array
t=h.data('star_age');
# You could operate on this array to, say, transform it into thousands of years if you wanted to 
#t=t/1000;
# Or even days
#t=t*365.25;

# Retrieve the luminosity
L=h.data('log_L');

# Retrieve the nuclear energy generation rates of two high-burning zones: z1 and z2
# Where (R(0) <) R(z1) < R(z2)                                                                                                 ?
# z1:
#z1_epsilon_start=h.data('epsnuc_M_2');
#z1_epsilon_stop =h.data('epsnuc_M_3');
# z2:
#z2_epsilon_start=h.data('epsnuc_M_6');
#z2_epsilon_stop =h.data('epsnuc_M_7');

# Retrieve the ejecta mass fractions we care about individually
ejected_Li7 =h.data('ejected_li7');
ejected_Be7 =h.data('ejected_be7');
ejected_C13 =h.data('ejected_c13');
ejected_N15 =h.data('ejected_n15');
ejected_O17 =h.data('ejected_o17');
ejected_Al26=h.data('ejected_al26');
ejected_Mg26=h.data('ejected_mg26');

# Retrieve the ejecta mass fractions of other isotopes
# as to sum mass of ejecta
# and to compare C/N and O/N ratios
ejected_H1  =h.data('ejected_h1');
ejected_H2  =h.data('ejected_h2');
ejected_He3 =h.data('ejected_he3');
ejected_He4 =h.data('ejected_he4');
ejected_B8  =h.data('ejected_b8');
ejected_C12 =h.data('ejected_c12');
ejected_N13 =h.data('ejected_n13');
ejected_N14 =h.data('ejected_n14');
ejected_O14 =h.data('ejected_o14');
ejected_O15 =h.data('ejected_o15');
ejected_O16 =h.data('ejected_o16');
ejected_O18 =h.data('ejected_o18');
ejected_F17 =h.data('ejected_f17');
ejected_F18 =h.data('ejected_f18');
ejected_F19 =h.data('ejected_f19');
ejected_Ne18=h.data('ejected_ne18');
ejected_Ne19=h.data('ejected_ne19');
ejected_Ne20=h.data('ejected_ne20');
ejected_Ne22=h.data('ejected_ne22');
ejected_Na23=h.data('ejected_na23');
ejected_Mg22=h.data('ejected_mg22');
ejected_Mg24=h.data('ejected_mg24');
ejected_Mg25=h.data('ejected_mg25');
ejected_Al25=h.data('ejected_al25');

# -*- coding: utf-8 -*-
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
    

    if allsplit[4]=="V" or allsplit[4]=="Vis." and float(allsplit[0])>2.445E6:#and allsplit[19]=="STD" and allsplit[1]!="" and  allsplit[0]!="":
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
d=4700   # https://arxiv.org/pdf/1202.6171.pdf
Rv=3.1   # interstellar extinction coefficient, quoted to be 3.1+-0.2 but everyone uses 3.1, Fitzpatrick (late 1990s/2000s can't remember)
Ebv=0.3  # https://arxiv.org/pdf/1202.6171.pdf
A=Rv*Ebv # Extinction coefficient

M_x_sol_V=4.8 #FOR THE V BAND, https://iopscience.iop.org/article/10.3847/1538-4365/aabfdf/pdf
L_V = np.log10((d/10)**2   /    (10**((np.asarray(Vmagnitudes)-M_x_sol_V-  A  )/2.5)))
#print(len(L_V))
#print(len(VJDs))

x_axis_scale=365.25 #change this to 1 for JD

L_model_max=max(np.asarray(L))
L_exp_max =max(np.asarray(L_V))

for i in range(len(L)):
    if L[i] == L_model_max:
        t1 = t[i]

for i in range(len(L_V)):
    if L_V[i] == L_exp_max:
        t2 = VJDs[i]/365.25

#print(t1)
#print(t2)
#print(t1-t2)
x_axis_obs_offset=t1-t2

plt.scatter(x_axis_obs_offset+np.asarray(VJDs)/x_axis_scale, L_V, facecolors='orange', marker="x", edgecolors="orange", s=15, label="Vdata")
#plt.xlabel("JD [yr]")
plt.xlabel("Age of Model [yr]")
plt.ylabel(r"Log(L/L$_\odot$)")

f = interpolate.interp1d(np.asarray(VJDs)/x_axis_scale + +x_axis_obs_offset, L_V,kind="linear")
all_arr_x=np.linspace(min(np.asarray(VJDs)/x_axis_scale +x_axis_obs_offset),max(np.asarray(VJDs)/x_axis_scale + x_axis_obs_offset),100000)
fit=f(all_arr_x)
plt.plot(all_arr_x,fit, linestyle="--", c="black")

############################################################
flag = True
M_min = min(Vmagnitudes)
M_length = len(np.asarray(Vmagnitudes))
t_exp_1=-1
t_exp_2=-1
for i in range(M_length):
    if Vmagnitudes[i] == M_min:
        t_exp_1 = VJDs[i]
#        print("  ")
#        print(M_min)
#        print(np.asarray(L_V[i]))
#        print(t_exp_1)
#        print(i)
    if Vmagnitudes[i] > (M_min+2) and VJDs[i] > t_exp_1 and flag == True:
#        print("  ")
#        print(np.asarray(Vmagnitudes[i]))
#        print(np.asarray(L_V[i]))
#        print(np.asarray(VJDs[i]))
#        print(i)
#	print(t_exp_2)
        t_exp_2 = VJDs[i]-t_exp_1
	t_exp_2 = t_exp_2/365.25
#	print(t_exp_2)
        flag = False

print(" ")
print("-=============================-")
print("           EXP. DATA           ")
print("-=============================-")
print(" ")
print("log(L_peak) = "+ str(L_exp_max))
print("t_2         = "+ str(t_exp_2))
print(" ")

############################################################

eject_pos = 1E9
t_eject   = 1E9
L_eject   = -1
M_v_eject = -1
peak_pos  = -1
t_peak    = 1E9
t_2       = -1

prof_sens = 150 #this is the number of models that are required to write 1 profile. found in inlist_main_controls

flag   = True
flag_M = True

L_length = len(np.asarray(L))
L_peak = max(np.asarray(L))
ejected_O17_maximum = max(np.asarray(ejected_O17))

#define array for visual magnitudes as too calculate t_2
#initialise with same length as array of L
M_v=[None]*L_length

newflag = True

for i in range(L_length):

    M_v[i]  = 19.09 - 2.5*L[i] #this is using d = 4700 kpc E(B - V) = 0.3 from PU Vul data
    M_v_min = 19.09 - 2.5*L_peak

    #print(M_v[i])

    if ejected_O17[i] == ejected_O17_maximum and flag == True:
        eject_pos = i
    #if i==1123 and flag == True:
    if i > eject_pos and ejected_O17[i] == 0 and newflag == True:
        #eject_pos = i
        eject_pos = i-1
        t_eject   = t[eject_pos]
        L_eject   = L[eject_pos]
        M_v_eject = M_v[eject_pos]
        #
        ejected_Li7_max  = ejected_Li7[eject_pos]
        ejected_Be7_max  = ejected_Be7[eject_pos]
        ejected_C13_max  = ejected_C13[eject_pos]
        ejected_N15_max  = ejected_N15[eject_pos]
        ejected_O17_max  = ejected_O17[eject_pos]
        ejected_Al26_max = ejected_Al26[eject_pos]
        ejected_Mg26_max = ejected_Mg26[eject_pos]
        #
        ejected_C12_max  = ejected_C12[eject_pos]
        ejected_N13_max  = ejected_N13[eject_pos]
        ejected_N14_max  = ejected_N14[eject_pos]
        ejected_O14_max  = ejected_O14[eject_pos]
        ejected_O15_max  = ejected_O15[eject_pos]
        ejected_O16_max  = ejected_O16[eject_pos]
        ejected_O18_max  = ejected_O18[eject_pos]
	#
        ejected_H1_max   = ejected_H1[eject_pos]
	ejected_H2_max   = ejected_H2[eject_pos]
	ejected_He3_max  = ejected_He3[eject_pos]
        ejected_He4_max  = ejected_He4[eject_pos]
	ejected_B8_max   =  ejected_B8[eject_pos]
	ejected_C12_max  =  ejected_C12[eject_pos]
	ejected_N13_max  =  ejected_N13[eject_pos]
	ejected_N14_max  =  ejected_N14[eject_pos]
	ejected_O14_max  =  ejected_O14[eject_pos]
	ejected_O15_max  =  ejected_O15[eject_pos]
	ejected_O16_max  =  ejected_O16[eject_pos]
	ejected_O18_max  =  ejected_O18[eject_pos]
	ejected_F17_max  =  ejected_F17[eject_pos]
	ejected_F18_max  =  ejected_F18[eject_pos]
	ejected_F19_max  =  ejected_F19[eject_pos]
	ejected_Ne18_max = ejected_Ne18[eject_pos]
	ejected_Ne19_max = ejected_Ne19[eject_pos]
	ejected_Ne20_max = ejected_Ne20[eject_pos]
	ejected_Ne22_max = ejected_Ne22[eject_pos]
	ejected_Na23_max = ejected_Na23[eject_pos]
	ejected_Mg22_max = ejected_Mg22[eject_pos]
	ejected_Mg24_max = ejected_Mg24[eject_pos]
	ejected_Mg25_max = ejected_Mg25[eject_pos]
	ejected_Al25_max = ejected_Al25[eject_pos]
	#
        mass_tot=ejected_Li7_max+ejected_Be7_max+ejected_C13_max+ejected_N15_max+ejected_O17_max+\
	ejected_Al26_max+ejected_Mg26_max+ejected_C12_max+ejected_N13_max+ejected_N14_max+ejected_O14_max+\
	ejected_O15_max+ejected_O16_max+ejected_O18_max+ejected_H1_max+ejected_H2_max+ejected_He3_max+\
        ejected_He4_max+ejected_B8_max+ejected_C12_max+ejected_N13_max+ejected_N14_max+ejected_O14_max+\
	ejected_O15_max+ejected_O16_max+ejected_O18_max+ejected_F17_max+ejected_F18_max+ejected_F19_max+\
	ejected_Ne18_max+ejected_Ne19_max+ejected_Ne20_max+ejected_Ne22_max+ejected_Na23_max+ejected_Mg22_max+\
	ejected_Mg24_max+ejected_Mg25_max+ejected_Al25_max
	newflag = False

#    if i > eject_pos and ejected_O17[i] == 0 and newflag == True:
#	print(i-1)
#	newflag = False

    if L[i] == L_peak and t[i] > t_eject and flag == True:
        peak_pos = i
        t_peak   = t[peak_pos]
    if M_v[i] > (M_v_min+2) and t[i] > t_peak and flag_M == True:
        t_2 = t[i] - t_peak
        #print(t[i])
        flag_M = False

#calculating C/N and O/N ratios
C_tot=ejected_C12_max/12 + ejected_C13_max/13
N_tot=ejected_N13_max/13 + ejected_N14_max/14 + ejected_N15_max/15
O_tot=ejected_O14_max/14 + ejected_O15_max/15 + ejected_O16_max/16 + ejected_O17_max/17 + ejected_O18_max/18

C_N_ratio = C_tot/N_tot
O_N_ratio = O_tot/N_tot
#C_N_ratio = ejected_C13_max/ejected_N15_max
#O_N_ratio = ejected_O17_max/ejected_N15_max

print("-=============================-")
print("         MODEL OUTBURST        ")
print("-=============================-")
print(" ")
print("Model #         = "+ str(eject_pos+1))
eject_pos = eject_pos/prof_sens
print("Profile #       = "+ str(eject_pos+1))
print("log(L_outburst) = "+ str(L_eject))
print("M_v_outburst    = "+ str(M_v_eject))
print("t_outburst      = "+ str(t_eject))
print(" ")
print("Mass of ejecta from outburst for isotopes:")
print("H-1   = "+ str(format(ejected_H1_max,".8E")))
print("He-4  = "+ str(format(ejected_He4_max,".8E")))
print("Li-7  = "+ str(format(ejected_Li7_max,".8E")))
print("Be-7  = "+ str(format(ejected_Be7_max,".8E")))
print("C-13  = "+ str(format(ejected_C13_max,".8E")))
print("N-15  = "+ str(format(ejected_N15_max,".8E")))
print("O-17  = "+ str(format(ejected_O17_max,".8E")))
print("Al-26 = "+ str(format(ejected_Al26_max,".8E")))
print("Mg-26 = "+ str(format(ejected_Mg26_max,".8E")))
print(" ")
print("Ratios of elements within ejecta:")
print("C/N = "+ str(C_N_ratio))
print("O/N = "+ str(O_N_ratio))
print(" ")
print("Total mass of ejection = "+ str(format(mass_tot,".8E")))
print(" ")
print("-=============================-")
print("          MODEL PEAK           ")
print("-=============================-")
print(" ")
print("Model #     = "+ str(peak_pos+1))
peak_pos = peak_pos/prof_sens
print("Profile #   = "+ str(peak_pos+1))
print("log(L_peak) = "+ str(L_peak))
print("M_v_peak    = "+ str(M_v_min))
print("t_peak      = "+ str(t_peak))
print("t2          = "+ str(t_2))
print(" ")

############################################################

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

mesa_interp_fit=interpolate.interp1d(np.asarray(t), np.asarray(L), kind="linear")
fit_mesa_extrapolated=extrap1d(mesa_interp_fit)
fit_mesa=fit_mesa_extrapolated(all_arr_x)

list_fit_for_chisq=[]
list_fit_mesa_for_chisq=[]
all_arr_x_filtered=[]
for i in range(0, len(all_arr_x)):
	if all_arr_x[i] > 0 and all_arr_x[i] < 9999999999999999:
		list_fit_for_chisq.append(fit[i])
		list_fit_mesa_for_chisq.append(fit_mesa[i])
		all_arr_x_filtered.append(all_arr_x[i])

fit_for_chisq=np.asarray(list_fit_for_chisq)
fit_mesa_for_chisq=np.asarray(list_fit_mesa_for_chisq)
chisq=np.sum((fit_mesa_for_chisq-fit_for_chisq)**2/(fit_for_chisq)) / len(list_fit_mesa_for_chisq)

############################################################

delta_L   = abs(L_exp_max-L_peak)
delta_t_2 = abs(t_exp_2  -t_2   )

print("-=============================-")
print("           ANALYSIS            ")
print("-=============================-")
print(" ")
print("Chi-squared      = "+ str(chisq))
print("Change in L_peak = "+ str(delta_L))
print("Change in t_2    = "+ str(delta_t_2))
print("Let Phi equal the product, then")
print("Phi              = "+ str(chisq*delta_L*delta_t_2))
print(" ")

plt.plot(t, L)
plt.plot(all_arr_x,fit_mesa)
plt.xlim(t_peak-10,t_peak+50)
plt.ylim(1.0,5.0)
plt.savefig("FINAL.pdf")
#plt.plot(all_arr_x_filtered,fit_mesa_for_chisq)
#plt.show()
