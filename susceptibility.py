import numpy as n
import matplotlib.pyplot as plt
import matplotlib
import cmath
import time
import itertools
import pickle
from scipy import integrate, optimize
import sympy as sym
import scipy.integrate as integ

ao_ang = 0.529/1.26
lat_const = 3.0
boltz = 8.617333262145e-5

#%% DFT (FPLO) Orbital-Weighted Band Structure Plotter

fermi_band_number = 1
nb_bands = 14
nb_orbs = 14
orbs_array = [r"$1s$", r"$2s$", r"$3s$", r"$2p_y$", r"$2p_z$", r"$2p_x$", r"$3p_y$", r"$3p_z$", r"$3p_x$", r"$3d_{-2}$",r"$3d_{-1}$",r"$3d_{0}$",r"$3d_{1}$",r"$3d_{2}$",]

bandsw = n.loadtxt("fplo//band//eq_lat_var_x//data_files//3.00//+bweights_200", unpack=True)
weights = [[[]]]
temp_array = []
t_array=[]


# Parsing band weights output file to extract weights in arrays
for i in range(nb_bands):
    if i>0:
        weights += [[]]
    #Looping over kpoints + orbitals
    for j in range(nb_orbs + 2):
            k = 0
            #Looping over to separate the bands
            while (nb_bands*k + i)<len(bandsw[j]):
                temp_array += [bandsw[j][nb_bands*k + i]]
                k += 1
            weights[i] += [temp_array]
            temp_array=[]
#Removing first empty array in the first band array (not sure why it exists only there)
weights[0].pop(0)

#Setting up the color arrays with weight-dependent alpha transparency
weight_colors = np.zeros((len(weights[0][0]),4))
rgb_colors = rgb_colors = [[1,0,0], [0,1,0], [0,0,1], [1,0,1], [1,1,0], [0,1,1], [0.5,0,0], [0,0.5,0], [0,0,0.5],[0.5, 0.75, 0.25], [0.25,0.75,0.5],[0.75,0.5,0.25], [0.25,0.5,0.75],[0.5,0.25,0.75]]

#Plotting the bands themselves
for i in range(len(weights)):
    plt.plot(weights[i][0], weights[i][1], color = 'black')
    
    #Plotting the weights as color transparency scatterplots
    for j in range (2,len(weights[0])):
        weight_colors[:,0:3] = rgb_colors[j-2]
        weight_colors[:,3] = weights[i][j]
        if i == 0:
            plt.scatter(weights[i][0], weights[i][1], color = weight_colors)
        else:
            plt.scatter(weights[i][0], weights[i][1], color = weight_colors)
            
#Making x/y array for points for legend
x_points_legend = n.full(201,0)
y_points_legend = n.full(201,-30)
    
#Plotting random points for the legend
for j in range (2,len(weights[0])):
    weight_colors[:,0:3] = rgb_colors[j-2]
    weight_colors[:,3] = np.ones(len(weights[0][0]))
    plt.scatter(x_points_legend, y_points_legend, color = weight_colors, label = orbs_array[j-2])

#Adds vertical lines at high-symmetry k-points
i = 0
sympoints = []
while i  <= len(weights[0][0]):
    #plt.axvline(weights[0][0][i], c = 'black')
    sympoints += [weights[0][0][i]]
    i += 100


#Visual stuff
plt.legend(loc = 'upper right', fontsize = 20)
plt.xticks(np.sort(sympoints),kpath)
plt.tick_params(axis = 'x', which = 'both', bottom = False, top = False, labelsize = 20)
plt.tick_params(axis = 'y', labelsize = 20)
plt.ylabel(r'$E - E_F \; (eV)$',fontsize = 20)
plt.ylim(-1,2)
plt.xlim(min(weights[0][0]),max(weights[0][0]))
plt.grid(True)


#Saving the obtained spectra with the right resolution
figure = plt.gcf()
figure.set_size_inches(21,11.8)
plt.show()
plt.savefig("..//H//fplo//band//Li_eq_band_scale.png", bbox_inches='tight')
plt.close()



#%% Hamiltonian fit parameters from DFT (FPLO) calculation

fermi_band_number = 1
fermi_energy = 0.0
nb_bands = 14
nb_orbs = 14
bandsw = n.loadtxt("..//Li//fplo//band//eq_lat_var_x//data_files//3.00//+bweights")
weights = [[[]]]
temp_array = []

#Parsing band weights output file to extract weights in arrays
for i in range(nb_bands):
    if i>0:
        weights += [[]]
    #Looping over kpoints + orbitals
    for j in range(nb_orbs + 2):
            k = 0
            #Looping over to separate the bands
            if j == 1:
                while (nb_bands*k + i)<len(bandsw[j]):
                    temp_array += [bandsw[j][nb_bands*k + i] - fermi_energy]
                    k += 1
            else:
                while (nb_bands*k + i)<len(bandsw[j]):
                    temp_array += [bandsw[j][nb_bands*k + i]]
                    k += 1   
            weights[i] += [temp_array]
            temp_array=[]
#Removing first empty array in the first band array (not sure why it exists only there)
weights[0].pop(0)

#Zeroth order electron energy dispersion -2tcos(ka) fitting function
def dispertion(x,t,a):
    return 2*t*n.cos(x*2*n.pi/weights[fermi_band_number][0][-1]) + a

#Initializing k-points and q-points array
kpoints_disp_fit = weights[fermi_band_number][0]


######Calculating the energy parameters from fits of the band structure
#Contructing the s-band and p-band energies (removing the effect of crossing)
kpoints_disp_fit = weights[fermi_band_number][0]
s_band_energies = weights[fermi_band_number][1]
p_band_energies = []
i = 0
while i < len(kpoints_disp_fit) and weights[fermi_band_number + 1][1][i] != weights[fermi_band_number + 2][1][i]:
    p_band_energies += [weights[fermi_band_number+1][1][i]]
    i += 1
    
while i <= len(kpoints_disp_fit)//2:
    p_band_energies += [weights[fermi_band_number+3][1][i]]
    i += 1
  
if len(p_band_energies) < len(kpoints_disp_fit):
    for i in reversed(p_band_energies[:-1]):
        p_band_energies += [i]

#Fitting both s and p bands to extract tss and tpp
t_s,t_s_covariance = optimize.curve_fit(dispertion,kpoints_disp_fit,s_band_energies, p0 = [1,1])

tss = t_s[0]
delta = ( (p_band_energies[len(kpoints_disp_fit)//2] - s_band_energies[len(kpoints_disp_fit)//2]) - (p_band_energies[len(kpoints_disp_fit)//4] - s_band_energies[len(kpoints_disp_fit)//4]) )/2
tpp = ( (p_band_energies[len(kpoints_disp_fit)//2] - s_band_energies[len(kpoints_disp_fit)//2]) - delta )/2 - tss
tsp = n.sqrt((p_band_energies[len(kpoints_disp_fit)//4] - s_band_energies[len(kpoints_disp_fit)//4])**2 - delta**2)/4
print("delta = %f \ntss = %f \ntpp = %f \nts = %f" % (delta,tss,tpp,tsp))


#%% Chain 1D Li Tight Binding Model

#Hamiltonian parameter values
delta = 5.163
tss = 1.185 
tpp = 2.743 
tsp = 1.558

#Initializing the Hamiltonian and filling it for all k-points
chosen_dim = 2
k_array = n.linspace(-n.pi, n.pi, 501)
H = n.zeros((chosen_dim,chosen_dim),dtype=n.cdouble)

#eigenvalues and eigenvectors array initialization (eval indexing: [k][band] evec indexing: [k][orbital][band])
evals = n.zeros((len(k_array),chosen_dim), dtype=n.cdouble)
evecs = n.zeros((len(k_array),chosen_dim,chosen_dim), dtype=n.cdouble)

#Filling the Hamiltonian and diagonalizing for all k-points
for i in range(len(k_array)):
    k = k_array[i]

    #ss part
    H[0,0] = -2*tss*n.cos(k) 
    
    #pz-pz part
    H[1,1] = 2*tpp*n.cos(k) + delta
    
    #s-pz part
    H[0,1] = -2*tsp*sym.I*n.sin(k)
    H[1,0] = 2*tsp*sym.I*n.sin(k)
    
    #Diagonalizing the Hamiltonian and populating the evecs and evals arrays
    evals_temp,evecs_temp = n.linalg.eigh(H)
    evals_temp=evals_temp.real
    evals[i] = evals_temp
    evecs[i] = evecs_temp
            
    
### (defect_positions array + all_results)
with open("..//Li//chain_susc//chain_eigenstates_k_%d.pkl" % len(k_array), 'wb') as f:
    pickle.dump([evals] + [evecs], f)           
###          







#%% Fermi Function Definition And Plotting

fermi_energy = -2
temperature = 0.05

#Fermi Function (exponential)
def fermi_func(energy, fermi_energy, temp):
    
    #Fermi function at T>0
    if temp > 0.0:
        return 1/(1 + n.exp((energy - fermi_energy)/(temp)) ) 
    
    #Heaviside function at T=0
    else:
        if energy <= fermi_energy:
            return 1
        else:
            return 0
        

energy_try = n.linspace(-4,4,100)
fermi_plot_array = []
for i in energy_try:
    fermi_plot_array += [fermi_func(i,fermi_energy,temperature)]

plt.plot(energy_try,fermi_plot_array)





#%% Chain 1D Li X(q) 

#(eval indexing: [k][band] evec indexing: [k][orbital][band])

start_time = time.time()

#Chosen temperature and lifetime (imaginary part in the susceptibility
temperature = 0.0
temperature_comparison = 180*boltz
imag_part = 1E-6
ao_ang_modified_s_divider = 2.5
ao_ang_modified_p_divider = 1.6

#Overlap Integrals
omega_ss_1 = 0.5529
omega_pp_1 = -0.2757
omega_sp_1 = 0.3904

#Importing evals and evecs and defining the k and q arrays
k_density = 501
with open("..//Li//chain_susc//chain_eigenstates_k_%d.pkl" % k_density, 'rb') as f:
    diag_results = pickle.load(f)    
evals = n.array(diag_results[0])
evecs = n.array(diag_results[1])
k_points= n.linspace(-n.pi, n.pi, k_density)
q_points_indices = [i for i in range(0,2*k_density)]

 
#Importing the already-calculated integrals for each q
with open("..//Li//chain_susc//chain_integrals_k_%d_a0_s_%.2f_a0_p_%.2f.pkl" % (k_density, ao_ang_modified_s_divider, ao_ang_modified_p_divider), 'rb') as f:
    imported_integrals = pickle.load(f)[0]  

#Defining the fermi energy and the fermi band index
fermi_energy = evals[len(k_points)//4,0].real
fermi_band_index = 0

#Initializing necessary quantities to 0 so we can fill them up in the k-sum for each q point
chi_zero = n.zeros(len(q_points_indices))
chi_full = n.zeros(len(q_points_indices))
chi_ss = n.zeros(len(q_points_indices))
chi_pp = n.zeros(len(q_points_indices))
chi_sp = n.zeros(len(q_points_indices))
chi_pure_s = n.zeros(len(q_points_indices))
chi_temperature = n.zeros(len(q_points_indices))
q_sum_count = 0

###Determining all of the permutations of the orbital sum in n(q,k) in an array of lists (mu, nu, l, mu', nu', l')
nqk_sum_orbital_permutations = list(itertools.product(range(2), repeat=4))

#Calculating X(q) for each q
for q_index in q_points_indices:
    q_value = 2/k_density * n.pi/lat_const * q_index
     
    #Initializing the X(q) arrays that will need to be integrated over k for each q      
    chi_zero_array_to_integrate = n.zeros(len(k_points))
    chi_full_array_to_integrate = n.zeros(len(k_points))
    chi_ss_array_to_integrate = n.zeros(len(k_points))
    chi_pp_array_to_integrate = n.zeros(len(k_points))
    chi_sp_array_to_integrate = n.zeros(len(k_points))
    chi_pure_s_array_to_integrate = n.zeros(len(k_points))
    chi_temperature_array_to_integrate = n.zeros(len(k_points))
    
    #Summing over all k
    for k_index in range(len(k_points)):
        k_value = k_points[k_index]
        kplusq_index = (q_index + k_index)%k_density
        
        #Energy values from the chosen_band at k-point and at (k+q)-point
        energy_k = evals[k_index,fermi_band_index].real - fermi_energy
        energy_kplusq = evals[kplusq_index,fermi_band_index].real - fermi_energy
        
        #Occupation determined from the fermi function
        fermi_k = fermi_func(energy_k, 0, temperature)
        fermi_kplusq = fermi_func(energy_kplusq, 0, temperature)
        fermi_k_comparison = fermi_func(energy_k, 0, temperature_comparison)
        fermi_kplusq_comparison = fermi_func(energy_kplusq, 0, temperature_comparison)
        
        #Factor multiplying the Bloch state arising from the non-orthogonality of the basis functions
        ortho_factor = 1/n.sqrt( n.abs(evecs[k_index][0][fermi_band_index])**2 * (1 + 2*n.cos(k_value)*omega_ss_1) + n.abs(evecs[k_index][1][fermi_band_index])**2 * (1 + 2*n.cos(k_value)*omega_pp_1) + 2*n.sin(k_value)*omega_sp_1*( evecs[k_index][0][fermi_band_index]*n.conj(evecs[k_index][1][fermi_band_index]) - n.conj(evecs[k_index][0][fermi_band_index])*evecs[k_index][1][fermi_band_index])  )   
        
        #Calculating the n(k,q) factor by summing all orbital sum permutation contributions (permutation order: [mu, nu, mu', nu'])
        nqk = 0
        nqk_ss = 0
        nqk_pp = 0
        nqk_sp = 0
        for permutation in nqk_sum_orbital_permutations:
            
            coeff_part_ss = 0
            coeff_part_pp = 0
            coeff_part_sp = 0
            
            coeff_part = evecs[kplusq_index][permutation[0]][fermi_band_index] * n.conj(evecs[k_index][permutation[1]][fermi_band_index]) * n.conj(evecs[kplusq_index][permutation[2]][fermi_band_index]) * evecs[k_index][permutation[3]][fermi_band_index]
            
            #Defining the coefficient parts that depend only on ss, pp and sp
            if permutation == (0,0,0,0):
                coeff_part_ss = coeff_part
            elif permutation == (1,1,1,1):
                coeff_part_pp = coeff_part
            elif permutation == (0,0,1,1) or permutation == (1,1,0,0):
                coeff_part_sp = coeff_part


                
            #determining which integrals to consider for each term of the sum (Imported integral indexing: [ss onsite, pp onsite, sp onsite, ss nn, pp nn, sp nn])
            index_mu_nu = 20
            index_mu_nu_prime = 20
            
            if permutation[0] == 0 and permutation[1] == 0: 
                index_mu_nu = 0
            elif permutation[0] == 1 and permutation[1] == 1: 
                index_mu_nu  = 1
            elif permutation[0] == 0 and permutation[1] == 1:
                index_mu_nu  = 2
            elif permutation[0] == 1 and permutation[1] == 0:
                index_mu_nu  = 3
                
            if permutation[2] == 0 and permutation[3] == 0: 
                index_mu_nu_prime  = 0
            elif permutation[2] == 1 and permutation[3] == 1: 
                index_mu_nu_prime  = 1 
            elif permutation[2] == 0 and permutation[3] == 1:
                index_mu_nu_prime  = 2
            elif permutation[2] == 1 and permutation[3] == 0:
                index_mu_nu_prime  = 3

                
                
            integral_part_onsite = n.conj(imported_integrals[q_index][index_mu_nu]) * imported_integrals[q_index][index_mu_nu_prime]
            integral_part_nn = 0
            
     
            ###summing through the l and l' parts for the nn contributions
            #l=1, lprime=0
            integral_part_nn += cmath.exp(1j*((q_value * lat_const) + k_value)) * n.conj(imported_integrals[q_index][index_mu_nu+4]) * imported_integrals[q_index][index_mu_nu_prime]
            #l=-1, lprime=0
            integral_part_nn += cmath.exp(-1j*((q_value * lat_const) + k_value)) * n.conj(imported_integrals[q_index][index_mu_nu+4]) * imported_integrals[q_index][index_mu_nu_prime]
            #l=0, lprime=1
            integral_part_nn += cmath.exp(-1j*((q_value * lat_const) + k_value)) * n.conj(imported_integrals[q_index][index_mu_nu]) * imported_integrals[q_index][index_mu_nu_prime+4]
            #l=0, lprime=-1
            integral_part_nn += cmath.exp(1j*((q_value * lat_const) + k_value)) * n.conj(imported_integrals[q_index][index_mu_nu]) * imported_integrals[q_index][index_mu_nu_prime+4]
            #l=1, lprime=1
            integral_part_nn += n.conj(imported_integrals[q_index][index_mu_nu+3]) * imported_integrals[q_index][index_mu_nu_prime+3]
            #l=-1, lprime=-1
            integral_part_nn += n.conj(imported_integrals[q_index][index_mu_nu+3]) * imported_integrals[q_index][index_mu_nu_prime+3]
            
            
            ###summing through the l and l' parts for the nn2 contributions
            #l=2, lprime=0
            integral_part_nn2 = cmath.exp(1j*2*q_value*lat_const) * n.conj(imported_integrals[q_index][index_mu_nu+6]) * imported_integrals[q_index][index_mu_nu_prime]
            #l=-2, lprime=0
            integral_part_nn2 += cmath.exp(-1j*2*q_value*lat_const) * n.conj(imported_integrals[q_index][index_mu_nu+6]) * imported_integrals[q_index][index_mu_nu_prime]
            #l=0, lprime=2
            integral_part_nn2 += cmath.exp(-1j*2*q_value*lat_const) * n.conj(imported_integrals[q_index][index_mu_nu]) * imported_integrals[q_index][index_mu_nu_prime+6]
            #l=0, lprime=-2
            integral_part_nn2 += cmath.exp(1j*2*q_value*lat_const) * n.conj(imported_integrals[q_index][index_mu_nu]) * imported_integrals[q_index][index_mu_nu_prime+6]
            #l=2, lprime=1
            integral_part_nn2 += n.conj(imported_integrals[q_index][index_mu_nu+6]) * imported_integrals[q_index][index_mu_nu_prime+6]
            #l=-2, lprime=-1
            integral_part_nn2 += n.conj(imported_integrals[q_index][index_mu_nu+6]) * imported_integrals[q_index][index_mu_nu_prime+6]
            
            nqk +=   ortho_factor * coeff_part  * (integral_part_onsite + integral_part_nn) 
            nqk_ss +=   ortho_factor * coeff_part_ss * (integral_part_onsite + integral_part_nn)
            nqk_pp +=   ortho_factor * coeff_part_pp * (integral_part_onsite + integral_part_nn)
            nqk_sp +=   ortho_factor * coeff_part_sp * ( integral_part_onsite + integral_part_nn)

        
        chi_zero_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k - fermi_kplusq)/((energy_k - energy_kplusq)**2 + imag_part**2)
        chi_full_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k - fermi_kplusq)/((energy_k - energy_kplusq)**2 + imag_part**2) * nqk.real
        chi_ss_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k - fermi_kplusq)/((energy_k - energy_kplusq)**2 + imag_part**2) * nqk_ss.real
        chi_pp_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k - fermi_kplusq)/((energy_k - energy_kplusq)**2 + imag_part**2) * nqk_pp.real
        chi_sp_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k - fermi_kplusq)/((energy_k - energy_kplusq)**2 + imag_part**2) * nqk_sp.real
        chi_pure_s_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k - fermi_kplusq)/((energy_k - energy_kplusq)**2 + imag_part**2) * n.exp(-(q_value/1.2)**2)
        chi_temperature_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k_comparison - fermi_kplusq_comparison)/((energy_k - energy_kplusq)**2 + imag_part**2) * nqk.real

    chi_zero[q_sum_count] = integrate.trapezoid(chi_zero_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0)
    chi_full[q_sum_count] = integrate.trapezoid(chi_full_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0) 
    chi_ss[q_sum_count] = integrate.trapezoid(chi_ss_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0) 
    chi_pp[q_sum_count] = integrate.trapezoid(chi_pp_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0) 
    chi_sp[q_sum_count] = integrate.trapezoid(chi_sp_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0) 
    chi_pure_s[q_sum_count] = integrate.trapezoid(chi_pure_s_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0) 
    chi_temperature[q_sum_count] = integrate.trapezoid(chi_temperature_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0) 

    q_sum_count += 1


print("--- %s minutes ---" % ((time.time() - start_time)/60))

###Saving the susceptibility results
with open("..//Li//chain_susc//modified_chain_chi_nn2_%d.pkl" % len(k_points), 'wb') as f:
    pickle.dump([q_points_indices] + [chi_zero] + [chi_full] + [chi_ss] + [chi_pp] + [chi_sp] + [chi_pure_s] + [chi_temperature], f)           
###      
    


#%% Importing X(q) and plotting

### Importing evals and evecs
k_density =301
with open("..//Li//chain_susc//modified_chain_chi_nn2_%d.pkl" % k_density, 'rb') as f:
    chi_imported_results = pickle.load(f)    
q_path_plotting = n.linspace(0,100,len(chi_imported_results[0]))
q_path_plotting = [2/k_density * i for i in chi_imported_results[0]]
chi_zero = chi_imported_results[1]
chi_full = chi_imported_results[2]
chi_ss = chi_imported_results[3]
chi_pp = chi_imported_results[4]
chi_sp = chi_imported_results[5]
chi_pure_s = chi_imported_results[6]
chi_temperature = chi_imported_results[7]

#Plotting the normal chi fig
multiplier_res = 3
matplotlib.pyplot.figure(figsize=(1.5*multiplier_res, 1*multiplier_res), dpi=200)

plt.plot(q_path_plotting, -chi_zero, color = 'black',linewidth=3,  label=r"$\chi_0(q)$, n(q,k)=0")
plt.plot(q_path_plotting, -chi_pure_s, color = 'black',linewidth=3,  label=r"$\alpha_{2s}(k)=1, \alpha_{2p}(k)=0$")
plt.plot(q_path_plotting, -chi_full, color = 'red',linewidth=3,  label=r"Full $\chi(q)$ (T = 0 K)")
plt.plot(q_path_plotting, -chi_temperature, color = 'blue',linewidth=2.5, linestyle="dashed",  label=r"Full $\chi(q)$ (T = 180 K)")

plt.plot(q_path_plotting, -chi_ss, color = 'green',linewidth=5,  label=r"$X_{ss}(q)$")
plt.plot(q_path_plotting, -chi_pp, color = 'yellow',linewidth=5,  label=r"$X_{pp}(q)$")
plt.plot(q_path_plotting, -chi_sp, color = 'black',linewidth=5,  label=r"$X_{sp}(q)$")

plt.ylabel(r'$Re[\chi(q)] \; (eV^{-1})$',fontsize = 20)
plt.xlabel(r'$q \; (\frac{\pi}{a})$',fontsize = 20)
plt.tick_params(axis = 'x', labelsize = 16.4)
plt.tick_params(axis = 'y', labelsize = 16.4)
plt.xlim(min(q_path_plotting),max(q_path_plotting))
plt.xlim(0,4)
plt.legend(loc="upper right", fontsize=10)
plt.grid()
plt.show()
plt.close()





#%% Li Temperature multi-chain plotter

#(eval indexing: [k][band] evec indexing: [k][orbital][band])

temp_array_li = n.array([0.00, 0.005, 0.006, 0.006999999999999999, 0.008, 0.01, 0.0129259998932175, 0.01508033320875375, 0.01723466652429, 0.020000000000000004, 0.023266799807791498, 0.0267137331126495, 0.030000000000000006, 0.037, 0.04000000000000001, 0.05000000000000001, 0.06000000000000001, 0.07, 0.08, 0.09000000000000001, 0.1])

for temperature_comparison in temp_array_li:

    start_time = time.time()
    
    #Chosen temperature and lifetime (imaginary part in the susceptibility
    temperature = 0.0
    imag_part = 1E-6
    
    #Overlap Integrals
    omega_ss_1 = 0.5529
    omega_pp_1 = -0.2757
    omega_sp_1 = 0.3904
    
    
    #Importing evals and evecs and defining the k and q arrays
    k_density = 301
    with open("..//Li//chain_susc//chain_eigenstates_k_%d.pkl" % k_density, 'rb') as f:
        diag_results = pickle.load(f)    
    evals = n.array(diag_results[0])
    evecs = n.array(diag_results[1])
    k_points= n.linspace(-n.pi, n.pi, k_density)
    q_points_indices = [i for i in range(0,2*k_density)]
    
    #Importing the already-calculated integrals for each q
    with open("..//Li//chain_susc//chain_integrals_k_%d.pkl" % k_density, 'rb') as f:
        imported_integrals = pickle.load(f)[0]  
    
    
    #Defining the fermi energy and the fermi band index
    fermi_energy = evals[len(k_points)//4,0].real
    fermi_band_index = 0
    
    
    #Initializing necessary quantities to 0 so we can fill them up in the k-sum for each q point
    chi_zero = n.zeros(len(q_points_indices))
    chi_full = n.zeros(len(q_points_indices))
    chi_ss = n.zeros(len(q_points_indices))
    chi_pp = n.zeros(len(q_points_indices))
    chi_sp = n.zeros(len(q_points_indices))
    chi_pure_s = n.zeros(len(q_points_indices))
    chi_temperature = n.zeros(len(q_points_indices))
    q_sum_count = 0
    
    ###Determining all of the permutations of the orbital sum in n(q,k) in an array of lists (mu, nu, l, mu', nu', l')
    nqk_sum_orbital_permutations = list(itertools.product(range(2), repeat=4))
    
    
    #Calculating X(q) for each q
    for q_index in q_points_indices:
        q_value = 2/k_density * n.pi/lat_const * q_index
    
         
        #Initializing the X(q) arrays that will need to be integrated over k for each q      
        chi_zero_array_to_integrate = n.zeros(len(k_points))
        chi_full_array_to_integrate = n.zeros(len(k_points))
        chi_ss_array_to_integrate = n.zeros(len(k_points))
        chi_pp_array_to_integrate = n.zeros(len(k_points))
        chi_sp_array_to_integrate = n.zeros(len(k_points))
        chi_pure_s_array_to_integrate = n.zeros(len(k_points))
        chi_temperature_array_to_integrate = n.zeros(len(k_points))
        
        
        #Summing over all k
        for k_index in range(len(k_points)):
            k_value = k_points[k_index]
            kplusq_index = (q_index + k_index)%k_density
            
    
            #Energy values from the chosen_band at k-point and at (k+q)-point
            energy_k = evals[k_index,fermi_band_index].real - fermi_energy
            energy_kplusq = evals[kplusq_index,fermi_band_index].real - fermi_energy
            
            
            #Occupation determined from the fermi function
            fermi_k = fermi_func(energy_k, 0, temperature)
            fermi_kplusq = fermi_func(energy_kplusq, 0, temperature)
            fermi_k_comparison = fermi_func(energy_k, 0, temperature_comparison)
            fermi_kplusq_comparison = fermi_func(energy_kplusq, 0, temperature_comparison)
            
            #Factor multiplying the Bloch state arising from the non-orthogonality of the basis functions
            ortho_factor = 1/n.sqrt( n.abs(evecs[k_index][0][fermi_band_index])**2 * (1 + 2*n.cos(k_value)*omega_ss_1) + n.abs(evecs[k_index][1][fermi_band_index])**2 * (1 + 2*n.cos(k_value)*omega_pp_1) + 2*n.sin(k_value)*omega_sp_1*( evecs[k_index][0][fermi_band_index]*n.conj(evecs[k_index][1][fermi_band_index]) - n.conj(evecs[k_index][0][fermi_band_index])*evecs[k_index][1][fermi_band_index])  )   
            
            #Calculating the n(k,q) factor by summing all orbital sum permutation contributions (permutation order: [mu, nu, mu', nu'])
            nqk = 0
            nqk_ss = 0
            nqk_pp = 0
            nqk_sp = 0
            for permutation in nqk_sum_orbital_permutations:
                
                coeff_part_ss = 0
                coeff_part_pp = 0
                coeff_part_sp = 0
                
                coeff_part = evecs[kplusq_index][permutation[0]][fermi_band_index] * n.conj(evecs[k_index][permutation[1]][fermi_band_index]) * n.conj(evecs[kplusq_index][permutation[2]][fermi_band_index]) * evecs[k_index][permutation[3]][fermi_band_index]
                
                #Defining the coefficient parts that depend only on ss, pp and sp
                if permutation == (0,0,0,0):
                    coeff_part_ss = coeff_part
                elif permutation == (1,1,1,1):
                    coeff_part_pp = coeff_part
                elif permutation == (0,0,1,1) or permutation == (1,1,0,0):
                    coeff_part_sp = coeff_part
    
                    
                #determining which integrals to consider for each term of the sum (Imported integral indexing: [ss onsite, pp onsite, sp onsite, ss nn, pp nn, sp nn])
                index_mu_nu = 20
                index_mu_nu_prime = 20
                
                if permutation[0] == 0 and permutation[1] == 0: 
                    index_mu_nu = 0
                elif permutation[0] == 1 and permutation[1] == 1: 
                    index_mu_nu  = 1
                elif permutation[0] == 0 and permutation[1] == 1:
                    index_mu_nu  = 2
                elif permutation[0] == 1 and permutation[1] == 0:
                    index_mu_nu  = 3
                    
                if permutation[2] == 0 and permutation[3] == 0: 
                    index_mu_nu_prime  = 0
                elif permutation[2] == 1 and permutation[3] == 1: 
                    index_mu_nu_prime  = 1 
                elif permutation[2] == 0 and permutation[3] == 1:
                    index_mu_nu_prime  = 2
                elif permutation[2] == 1 and permutation[3] == 0:
                    index_mu_nu_prime  = 3
    
                    
                    
                integral_part_onsite = n.conj(imported_integrals[q_index][index_mu_nu]) * imported_integrals[q_index][index_mu_nu_prime]
                integral_part_nn = 0
                
         
                ###asumming through the l and l' parts for the nn contributions
                #l=1, lprime=0
                integral_part_nn += cmath.exp(1j*((q_value * lat_const) + k_value)) * n.conj(imported_integrals[q_index][index_mu_nu+4]) * imported_integrals[q_index][index_mu_nu_prime]
                #l=-1, lprime=0
                integral_part_nn += cmath.exp(-1j*((q_value * lat_const) + k_value)) * n.conj(imported_integrals[q_index][index_mu_nu+4]) * imported_integrals[q_index][index_mu_nu_prime]
                #l=0, lprime=1
                integral_part_nn += cmath.exp(-1j*((q_value * lat_const) + k_value)) * n.conj(imported_integrals[q_index][index_mu_nu]) * imported_integrals[q_index][index_mu_nu_prime+4]
                #l=0, lprime=-1
                integral_part_nn += cmath.exp(1j*((q_value * lat_const) + k_value)) * n.conj(imported_integrals[q_index][index_mu_nu]) * imported_integrals[q_index][index_mu_nu_prime+4]
                #l=1, lprime=1
                integral_part_nn += n.conj(imported_integrals[q_index][index_mu_nu+3]) * imported_integrals[q_index][index_mu_nu_prime+3]
                #l=-1, lprime=-1
                integral_part_nn += n.conj(imported_integrals[q_index][index_mu_nu+3]) * imported_integrals[q_index][index_mu_nu_prime+3]
                
                
                ###summing through the l and l' parts for the nn2 contributions
                #l=2, lprime=0
                integral_part_nn2 = cmath.exp(1j*2*q_value*lat_const) * n.conj(imported_integrals[q_index][index_mu_nu+6]) * imported_integrals[q_index][index_mu_nu_prime]
                #l=-2, lprime=0
                integral_part_nn2 += cmath.exp(-1j*2*q_value*lat_const) * n.conj(imported_integrals[q_index][index_mu_nu+6]) * imported_integrals[q_index][index_mu_nu_prime]
                #l=0, lprime=2
                integral_part_nn2 += cmath.exp(-1j*2*q_value*lat_const) * n.conj(imported_integrals[q_index][index_mu_nu]) * imported_integrals[q_index][index_mu_nu_prime+6]
                #l=0, lprime=-2
                integral_part_nn2 += cmath.exp(1j*2*q_value*lat_const) * n.conj(imported_integrals[q_index][index_mu_nu]) * imported_integrals[q_index][index_mu_nu_prime+6]
                #l=2, lprime=1
                integral_part_nn2 += n.conj(imported_integrals[q_index][index_mu_nu+6]) * imported_integrals[q_index][index_mu_nu_prime+6]
                #l=-2, lprime=-1
                integral_part_nn2 += n.conj(imported_integrals[q_index][index_mu_nu+6]) * imported_integrals[q_index][index_mu_nu_prime+6]
                
                
                
                #coeff_part = 1/4
                integral_part_nn = 0
                #integral_part_onsite = 1
                ortho_factor = 1
                #print(coeff_part)
                
                nqk +=   ortho_factor * coeff_part  * (integral_part_onsite + integral_part_nn) 
                nqk_ss +=   ortho_factor * coeff_part_ss * (integral_part_onsite + integral_part_nn)
                nqk_pp +=   ortho_factor * coeff_part_pp * (integral_part_onsite + integral_part_nn)
                nqk_sp +=   ortho_factor * coeff_part_sp * ( integral_part_onsite + integral_part_nn)
    
    
            #print(n.conj(imported_integrals[q_index][0]) * imported_integrals[q_index][0] - nqk_ss)
            
            chi_zero_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k - fermi_kplusq)/((energy_k - energy_kplusq)**2 + imag_part**2)
            chi_full_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k - fermi_kplusq)/((energy_k - energy_kplusq)**2 + imag_part**2) * nqk.real
            chi_ss_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k - fermi_kplusq)/((energy_k - energy_kplusq)**2 + imag_part**2) * nqk_ss.real
            chi_pp_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k - fermi_kplusq)/((energy_k - energy_kplusq)**2 + imag_part**2) * nqk_pp.real
            chi_sp_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k - fermi_kplusq)/((energy_k - energy_kplusq)**2 + imag_part**2) * nqk_sp.real
            chi_pure_s_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k - fermi_kplusq)/((energy_k - energy_kplusq)**2 + imag_part**2) * n.conj(imported_integrals[q_index][0]) * imported_integrals[q_index][0]
            chi_temperature_array_to_integrate[k_index] += (energy_k - energy_kplusq)*(fermi_k_comparison - fermi_kplusq_comparison)/((energy_k - energy_kplusq)**2 + imag_part**2) * nqk.real
    
        chi_zero[q_sum_count] = integrate.trapezoid(chi_zero_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0)
        chi_full[q_sum_count] = integrate.trapezoid(chi_full_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0) 
        chi_ss[q_sum_count] = integrate.trapezoid(chi_ss_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0) 
        chi_pp[q_sum_count] = integrate.trapezoid(chi_pp_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0) 
        chi_sp[q_sum_count] = integrate.trapezoid(chi_sp_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0) 
        chi_pure_s[q_sum_count] = integrate.trapezoid(chi_pure_s_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0) 
        chi_temperature[q_sum_count] = integrate.trapezoid(chi_temperature_array_to_integrate, dx=n.abs(k_points[1]-k_points[0])/lat_const, axis = 0) 
    
        q_sum_count += 1
    
    
    print("--- %s minutes ---" % ((time.time() - start_time)/60))
    
    print(temperature_comparison)

    ###Saving the susceptibility results
    with open("..//Li//chain_susc//erratum_results//temperature//chi_%d_temp_%.3f.pkl" % (len(k_points), temperature_comparison), 'wb') as f:
        pickle.dump([q_points_indices] + [chi_zero] + [chi_full] + [chi_ss] + [chi_pp] + [chi_sp] + [chi_pure_s] + [chi_temperature], f)           
    ###          
    
    
#%% Temperature Dependence Plotting

boltz = 8.617333262145e-5
temp_array_li = n.array([0.00, 0.005, 0.006, 0.006999999999999999, 0.008, 0.01, 0.0129259998932175, 0.01508033320875375, 0.01723466652429, 0.020000000000000004, 0.023266799807791498, 0.0267137331126495, 0.030000000000000006, 0.037, 0.04000000000000001, 0.05000000000000001, 0.06000000000000001, 0.07, 0.08, 0.09000000000000001, 0.1])


for temperature_comparison in temp_array_li:
#chosen_charge_transfer = charge_transfer_array[10]

    print(temperature_comparison)

    #Importing X(q) and plotting
    ### Importing evals and evecs
    k_density =301
    with open("..//Li//chain_susc//erratum_results//temperature//chi_%d_temp_%.3f.pkl" % (k_density, temperature_comparison), 'rb') as f:
        chi_imported_results = pickle.load(f)    
    q_path_plotting = n.linspace(0,100,len(chi_imported_results[0]))
    q_path_plotting = [2/k_density * i for i in chi_imported_results[0]]
    chi_zero = chi_imported_results[1]
    chi_full = chi_imported_results[2]
    chi_ss = chi_imported_results[3]
    chi_pp = chi_imported_results[4]
    chi_sp = chi_imported_results[6]
    chi_full_temperature = chi_imported_results[7]
    
    ###  
    
    #Plotting the normal chi fig
    multiplier_res = 3
    matplotlib.pyplot.figure(figsize=(1.5*multiplier_res, 1*multiplier_res), dpi=200)
    
    #plt.plot(q_path_plotting, -chi_zero, color = 'black',linewidth=3,  label=r"$\alpha_{2s}(k)=1, \alpha_{2p}(k)=0$")
    plt.plot(q_path_plotting, -chi_full, color = 'red',linewidth=3,  label=r"Full $\chi(q)$ (T = 0 K)")
    plt.plot(q_path_plotting, -chi_full_temperature, color = 'blue',linewidth=3, linestyle="dashed",  label=r"Full $\chi(q)$ (T = %.2f K)" % temperature_comparison)
    
    #plt.plot(q_path_plotting, -chi_ss, color = 'green',linewidth=5,  label=r"$X_{ss}(q)$")
    #plt.plot(q_path_plotting, -chi_pp, color = 'yellow',linewidth=5,  label=r"$X_{pp}(q)$")
    #plt.plot(q_path_plotting, -chi_sp, color = 'black',linewidth=5,  label=r"$X_{sp}(q)$")
    
    plt.ylabel(r'$Re[\chi(q)] \; (eV^{-1})$',fontsize = 20)
    plt.xlabel(r'$q \; (\frac{\pi}{a})$',fontsize = 20)
    plt.tick_params(axis = 'x', labelsize = 16.4)
    plt.tick_params(axis = 'y', labelsize = 16.4)
    plt.xlim(min(q_path_plotting),max(q_path_plotting))
    plt.legend(loc="upper right", fontsize=12)
    #plt.ylim(-6,14)
    plt.grid()
    plt.show()
    plt.close()


#%% Importing and Plotting Chain 1D Li Evals and Evecs

kpath = ( "X", r"$\Gamma$", "X") 
charge_transfer_array = [5.163*(0.1*i) for i in range(0,50)]
delta = charge_transfer_array[40]

### Importing evals and evecs
k_array = n.linspace(-n.pi, n.pi, 51)
#with open("..//Li//chain_susc//chain_eigenstates_k_%d.pkl" % len(k_array), 'rb') as f:
with open("..//Li//chain_susc//erratum_results//chain_eigenstates_k_%d_delta_%.2f.pkl" % (len(k_array), delta), 'rb') as f:
    diag_results = pickle.load(f)    
evals = diag_results[0]
evecs = diag_results[1]
###  

#Fermi enegy definition
fermi_energy = evals[len(k_array)//4,0]

'''
#Printing lengths for determining k-path array lengths
print(len(evals.T[i,len(k_array)//2,len(k_array)//2:len(k_array),len(k_array)//2].real))
print(len(evals.T[i,len(k_array)//2, len(k_array)-1, len(k_array)//2+1:len(k_array)].real))
print(len(evals.T[i,len(k_array)//2+1:len(k_array), len(k_array)-1, len(k_array)-1].real))
print(len([evals.T[i,len(k_array)-1 - bruh, len(k_array)-1 -bruh, len(k_array) -1 - bruh].real for bruh in range(1,len(k_array)//2)]))
print(len(k_array)//2)
'''

chosen_evecs = [0,1]
k_array_path = n.linspace(0,100,len(k_array))
sympoints = [k_array_path[0], k_array_path[len(k_array)//2], k_array_path[len(k_array)-1]]

#Plotting the energy bands
multiplier_res = 3
matplotlib.pyplot.figure(figsize=(1.87931317779*multiplier_res, 1*multiplier_res), dpi=200)
for i in chosen_evecs:
    plt.plot(k_array_path, evals.T[i] - fermi_energy, linewidth = 3, color='blue')

plt.ylabel(r"E (eV)", fontsize=17)
plt.xticks(n.sort(sympoints),kpath)
plt.tick_params(axis = 'x', which = 'both', bottom = False, top = False, labelsize = 17)
plt.tick_params(axis = 'y', labelsize = 16.4)
plt.xlim(min(k_array_path),max(k_array_path))
#plt.ylim(-6,14)
plt.grid()
plt.show()
plt.close()         



#%% Checking chosen band's orbital wavefunction 
   
#(eval indexing: [k][band] evec indexing: [k][orbital][band])

chosen_band = 0

legend_evec_array = [r"$2s$", r"$2p_z$"]

for i in range(2):
    plt.plot(k_array_path, evecs.T[chosen_band,i].imag, linewidth = 3, label = legend_evec_array[i])


plt.ylabel(r"E (eV)", fontsize=17)
plt.xticks(n.sort(sympoints),kpath)
plt.tick_params(axis = 'x', which = 'both', bottom = False, top = False, labelsize = 17)
plt.tick_params(axis = 'y', labelsize = 16.4)
plt.xlim(min(k_array_path),max(k_array_path))
#plt.ylim(-6,9)
plt.legend(loc="upper right")
plt.grid()
plt.show()
plt.close()             
            

print(evecs[1,7,24,1,1])
print(evecs.T[1,1,24,7,1])

brih = evals.T[1]



#%% Plotting real-space integrals

ao_ang = 0.529/1.26
ao_ang_modified_s = 0.529/2.5
ao_ang_modified_p = 0.529/1.6

def function_2s(r,theta,a_0,a):
    return (1/(4*a_0**(3/2)*n.sqrt(2*n.pi)))*n.exp(-(n.sqrt(r**2 + a**2 - 2*a*r*n.cos(theta)))/(2*a_0))*(2-n.sqrt(r**2 + a**2 - 2*a*r*n.cos(theta))/a_0)

def overlap_2s_2s(a_0,a):
    return 2*n.pi*integ.dblquad(lambda r, theta: function_2s(r,theta,a_0,0)*function_2s(r,theta,a_0,a)*r**2*n.sin(theta),0,n.pi,0,n.inf)[0]

def function_2p(r,theta,a_0,a):
    return (1/(4*a_0**(3/2)*n.sqrt(2*n.pi)))*n.exp(-n.sqrt(r**2 + a**2 - 2*a*r*n.cos(theta))/(2*a_0))*(n.sqrt(r**2 + a**2 - 2*a*r*n.cos(theta))/a_0)*n.cos(n.pi-n.arccos((a-r*n.cos(theta))/n.sqrt(r**2 + a**2 - 2*a*r*n.cos(theta))))

def overlap_2p_2p(a_0,a):
    return 2*n.pi*integ.dblquad(lambda r, theta: function_2p(r,theta,a_0,0)*function_2p(r,theta,a_0,a)*r**2*n.sin(theta),0,n.pi,0,n.inf)[0]

def overlap_2s_2p(a_0_s, a_0_p,a):
    return 2*n.pi*integ.dblquad(lambda r, theta: function_2s(r,theta,a_0_s,0)*function_2p(r,theta,a_0_p,a)*r**2*n.sin(theta),0,n.pi,0,n.inf)[0]

#Computing the norm to see if normalized
print("Norm ss= ", overlap_2s_2s(ao_ang_modified_s,0))
print("Norm pp= ", overlap_2p_2p(ao_ang_modified_p,0))
print("Overlap ss= ", overlap_2s_2s(ao_ang_modified_s,3))
print("Overlap pp= ", overlap_2p_2p(ao_ang_modified_p,3))
print("Overlap sp= ", overlap_2s_2p(ao_ang_modified_s,ao_ang_modified_p,3))

#Plotting s function
r_array = n.linspace(0, 3, 500)
plt.axhline(0, color='black')
plt.plot(r_array, n.vectorize(function_2s)(r_array, 0, ao_ang, 0), label="2s Original")
plt.plot(r_array, n.vectorize(function_2s)(r_array, 0, ao_ang_modified_s, 0), label="2s Compressed")
plt.legend(loc="upper right")
plt.grid()
plt.show()
plt.close()


#Plotting p function
r_array = n.linspace(0, 3, 500)
plt.axhline(0, color='black')
plt.plot(r_array, n.vectorize(function_2p)(r_array, 0, ao_ang, 0), label="2p Original")
plt.plot(r_array, n.vectorize(function_2p)(r_array, 0, ao_ang_modified_p, 0), label="2p Compressed")
plt.legend(loc="upper right")
plt.grid()
