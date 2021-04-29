# loads the right libraries in python
# these need to be installed in your system

# numpy allows for numerical calculations
import numpy as np
# matplotlib provides plotting capabilities
import matplotlib.pyplot as plt
# mesa_reader allows to get the data from the "./LOGS/" folder files easily
import mesa_reader as mr

# ==================================================================
# This program prints the novae's properites at maximum luminosities
# ==================================================================

# Load the "history.data" file and store it in "h" using Mesa_reader capabilities
h=mr.MesaData('LOGS/history.data')

# Retrieve the model number and save it in "counter" array
counter=h.data("model_number")

# Retrieve the star age (in years) and save it in "time" array
t=h.data('star_age');
# You could operate on this array to, say, transform it into thousands of years if you wanted to 
#t=t/1000;
# Or even days
t=t*365.25;

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
ejected_Li7=h.data('ejected_li7');
ejected_Be7=h.data('ejected_be7');
ejected_C13=h.data('ejected_c13');
ejected_N15=h.data('ejected_n15');
ejected_O17=h.data('ejected_o17');

# Retrieve the ejecta mass fractions of other CNO isotopes
# as to compare C/N and O/N ratios
ejected_C12=h.data('ejected_c12');
ejected_N13=h.data('ejected_n13');
ejected_N14=h.data('ejected_n14');
ejected_O14=h.data('ejected_o14');
ejected_O15=h.data('ejected_o15');
ejected_O16=h.data('ejected_o16');
ejected_O18=h.data('ejected_o18');
ejected_Al26=h.data('ejected_al26');
ejected_Mg26=h.data('ejected_mg26');

# ===================================
#Start of code

eject_pos = -1
#t_eject   = -1
t_eject   = 1E9
L_eject   = -1
M_v_eject = -1
#
peak_pos  = -1
#t_peak    = -1
t_peak    = 1E9
t_2       = -1
#
ejected_Li7_max=-1
ejected_Be7_max=-1
ejected_C13_max=-1
ejected_N15_max=-1
ejected_O17_max=-1
#
ejected_C12_max=-1
ejected_N13_max=-1
ejected_N14_max=-1
ejected_O14_max=-1
ejected_O15_max=-1
ejected_O16_max=-1
ejected_O18_max=-1
ejected_Al26_max=-1
ejected_Mg26_max=-1
#
prof_sens = 100 #this is the number of models that are required to write 1 profile. found in inlist_main_controls
#
flag   = True
flag_M = True
#
L_length = len(np.asarray(L))
L_peak = max(np.asarray(L))

ejected_O17_maximum = max(np.asarray(ejected_O17))

#define array for visual magnitudes as too calculate t_2
#initialise with same length as array of L
M_v=[None]*L_length

for i in range(L_length):

    M_v[i]  = 19.09 - 2.5*L[i] #this is using d = 4700 kpc E(B - V) = 0.3 from PU Vul data
    M_v_min = 19.09 - 2.5*L_peak

    #print(M_v[i])

    if ejected_O17[i] == ejected_O17_maximum and flag == True:
        eject_pos   = i
        t_eject     = t[eject_pos]
        L_eject     = L[eject_pos]
        M_v_eject   = M_v[eject_pos]
        #
        ejected_Li7_max = ejected_Li7[eject_pos]
        ejected_Be7_max = ejected_Be7[eject_pos]
        ejected_C13_max = ejected_C13[eject_pos]
        ejected_N15_max = ejected_N15[eject_pos]
        ejected_O17_max = ejected_O17[eject_pos]
        #
        ejected_C12_max = ejected_C12[eject_pos]
        ejected_N13_max = ejected_N13[eject_pos]
        ejected_N14_max = ejected_N14[eject_pos]
        ejected_O14_max = ejected_O14[eject_pos]
        ejected_O15_max = ejected_O15[eject_pos]
        ejected_O16_max = ejected_O16[eject_pos]
        ejected_O18_max = ejected_O18[eject_pos]
        ejected_Al26_max= ejected_Al26[eject_pos]
        ejected_Mg26_max= ejected_Mg26[eject_pos]
    if L[i] == L_peak and t[i] > t_eject and flag == True:
        peak_pos = i
        t_peak   = t[peak_pos]
    if M_v[i] > (M_v_min+2) and t[i] > t_peak and flag_M == True:
        t_2 = t[i] - t_peak
        #print(t[i])
        flag_M = False

#calculating C/N and O/N ratios
C_tot=ejected_C12_max + ejected_C13_max 
N_tot=ejected_N13_max + ejected_N14_max + ejected_N15_max
O_tot=ejected_O14_max + ejected_O15_max + ejected_O16_max + ejected_O17_max + ejected_O18_max 

C_N_ratio = C_tot/N_tot
O_N_ratio = O_tot/N_tot

print(" ")
print("-===================-")
print("    NOVA OUTBURST    ")
print("-===================-")
print(" ")
print("Model            = "+ str(eject_pos))
eject_pos = eject_pos/prof_sens
print("Profile Number   = "+ str(eject_pos))
print("Outburst log(L)  = "+ str(L_eject))
print("M_v  at Outburst = "+ str(M_v_eject))
print("Time at Outburst = "+ str(t_eject))
print(" ")
print("Mass fractions of ejecta from outburst:")
print("Li-7  = "+ str(format(ejected_Li7_max,".8E")))
print("Be-7  = "+ str(format(ejected_Be7_max,".8E")))
print("C-13  = "+ str(format(ejected_C13_max,".8E")))
print("N-15  = "+ str(format(ejected_N15_max,".8E")))
print("O-17  = "+ str(format(ejected_O17_max,".8E")))
print("Al-26 = "+ str(format(ejected_Al26_max,".8E")))
print("Mg-26 = "+ str(format(ejected_Mg26_max,".8E")))
#print(" ")
#print("Ratios of elements within ejecta:")
#print("C/N = "+ str(C_N_ratio))
#print("O/N = "+ str(O_N_ratio))
print(" ")
print("-===================-")
print("      t_2 CALC.      ")
print("-===================-")
print(" ")
print("Model          = "+ str(peak_pos))
peak_pos = peak_pos/prof_sens
print("Profile Number = "+ str(peak_pos))
print("Peak log(L)    = "+ str(L_peak))
print("M_v  at Peak   = "+ str(M_v_min))
print("Time at Peak   = "+ str(t_peak))
print("t2             = "+ str(t_2))
print(" ")
