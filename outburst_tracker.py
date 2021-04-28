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
# You could opearate on this array to, say, transform it into thousands of years if you wanted to 
#t=t/1000;
# Or even days
#t=t*365.25;

# Retrieve the luminosities 
L   =h.data('log_L');
L_H =h.data('log_LH');
L_He=h.data('log_LHe');

# Retrieve the temperature and density
T_eff=h.data('log_Teff');
rho_c=h.data('log_center_Rho');

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

peak_pos  = -1
t_max     = -1
t_2     = -1
L_H_max   = -1
L_He_max  = -1
T_eff_max = -1
rho_c_max = -1

ejected_Li7_max=-1
ejected_Be7_max=-1
ejected_C13_max=-1
ejected_N15_max=-1
ejected_O17_max=-1

ejected_C12_max=-1
ejected_N13_max=-1
ejected_N14_max=-1
ejected_O14_max=-1
ejected_O15_max=-1
ejected_O16_max=-1
ejected_O18_max=-1
ejected_Al26_max=-1
ejected_Mg26_max=-1

M_v_max = 10000
prof_sens = 100 #this is the number of models that are required to write 1 profile. found in inlist_main_controls

flag = True

L_length = len(np.asarray(L))

#ejected_C13_maximum = max(np.asarray(ejected_C13))
ejected_O17_maximum = max(np.asarray(ejected_O17))
#ejected_N15_maximum = max(np.asarray(ejected_N15))

#define array for visual magnitudes as too calculate t_2
#initialise with same length as array of L
M_v=[None]*L_length

for i in range(L_length):

    M_v[i] = 19.09 - 2.5*L[i] #this is using d = 4700 kpc E(B - V) = 0.3 from PU Vul data
    
    #print(M_v[i])

    #if ejected_C13[i] == ejected_C13_maximum and flag == True:
    if ejected_O17[i] == ejected_O17_maximum and flag == True:
    #if ejected_N15[i] == ejected_N15_maximum and flag == True:
        peak_pos  = i
        t_max     = t[peak_pos]
        L_max   = L[peak_pos]
        L_H_max   = L_H[peak_pos]
        L_He_max  = L_He[peak_pos]
        T_eff_max = T_eff[peak_pos]
        rho_c_max = rho_c[peak_pos]
        #
        ejected_Li7_max= ejected_Li7[peak_pos]
        ejected_Be7_max= ejected_Be7[peak_pos]
        ejected_C13_max= ejected_C13[peak_pos]
        ejected_N15_max= ejected_N15[peak_pos]
        ejected_O17_max= ejected_O17[peak_pos]
        #
        ejected_C12_max= ejected_C12[peak_pos]
        ejected_N13_max= ejected_N13[peak_pos]
        ejected_N14_max= ejected_N14[peak_pos]
        ejected_O14_max= ejected_O14[peak_pos]
        ejected_O15_max= ejected_O15[peak_pos]
        ejected_O16_max= ejected_O16[peak_pos]
        ejected_O18_max= ejected_O18[peak_pos]
        ejected_Al26_max= ejected_Al26[peak_pos]
        ejected_Mg26_max= ejected_Mg26[peak_pos]
        M_v_max = M_v[peak_pos]
    if M_v[i] > (M_v_max+2) and t[i] > t_max and flag == True:
        t_2 = t[i] - t_max
        flag = False

#calculating C/N and O/N ratios
C_tot=ejected_C12_max + ejected_C13_max 
#C_tot=ejected_C13_max 
N_tot=ejected_N13_max + ejected_N14_max + ejected_N15_max
O_tot=ejected_O14_max + ejected_O15_max + ejected_O16_max + ejected_O17_max + ejected_O18_max 
#O_tot=ejected_O17_max

C_N_ratio = C_tot/N_tot
O_N_ratio = O_tot/N_tot
#C_N_ratio = ejected_C13_max/ejected_N15_max
#O_N_ratio = ejected_O17_max/ejected_N15_max

print(" ")
print("-===================-")
print("    NOVA OUTBURST    ")
print("-===================-")
print(" ")
print("Model            = "+ str(peak_pos))
peak_pos = peak_pos/prof_sens
print("Profile Number   = "+ str(peak_pos))
print("Time of Outburst = "+ str(t_max))
print("t2               = "+ str(t_2))
print("log(L)           = "+ str(L_max))
print("log(L_H)         = "+ str(L_H_max))
print("log(L_He)        = "+ str(L_He_max))
print("log(T_eff)       = "+ str(T_eff_max))
print("log(rho_c)       = "+ str(rho_c_max))
print("M_v              = "+ str(M_v_max))
print(" ")
print("Mass fractions of ejecta from outburst:")
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
