import os, ctypes
from scipy import integrate, LowLevelCallable
import numpy as n
import time
import pickle
import sys

lib = ctypes.CDLL(os.path.abspath('integrals.so'))
start_time = time.time()

#Initial Parameters
k_density = 300
ao_ang = 0.529
lat_const = 1.0


#Importing all of the integral C functions and defining their input and outpup types
lib.real_1s1s.restype = ctypes.c_double
lib.real_1s1s.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)
lib.imag_1s1s.restype = ctypes.c_double
lib.imag_1s1s.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)

k_points= n.linspace(-n.pi, n.pi, k_density)
q_points_indices = [i for i in range(0,2*k_density)]

#Initializing the integral array
orb_integral_array = n.zeros((len(q_points_indices),3),dtype = n.cdouble)

#Calculating all of the orbital integrals for all q values, in order to export them for use in computing X(q)
for q_index in q_points_indices:
    qx_value = 0
    qy_value = 0
    qz_value = 2/k_density * n.pi/lat_const * q_index
     
    #Passing the q and a_0 parameters to the C functions 
    integrand_onsite_data_to_pass = ctypes.cast( (ctypes.c_double*7)(qx_value,qy_value,qz_value,ao_ang, 0,0,0) , ctypes.c_void_p)
    integrand_nn_data_to_pass = ctypes.cast( (ctypes.c_double*7)(qx_value,qy_value,qz_value,ao_ang, 0, 0, 1*lat_const) , ctypes.c_void_p)
    integrand_nn2_data_to_pass = ctypes.cast( (ctypes.c_double*7)(qx_value,qy_value,qz_value,ao_ang, 0, 0, 2*lat_const) , ctypes.c_void_p)
	
    #onsite integral functions
    onsite_1s1s_real = LowLevelCallable(lib.real_1s1s, integrand_onsite_data_to_pass)
    onsite_1s1s_imag = LowLevelCallable(lib.imag_1s1s, integrand_onsite_data_to_pass)

    #nn integral functions
    nn_1s1s_real = LowLevelCallable(lib.real_1s1s, integrand_nn_data_to_pass)
    nn_1s1s_imag = LowLevelCallable(lib.imag_1s1s, integrand_nn_data_to_pass)

    #nn2_integral functions
    nn2_1s1s_real = LowLevelCallable(lib.real_1s1s, integrand_nn2_data_to_pass)
    nn2_1s1s_imag = LowLevelCallable(lib.imag_1s1s, integrand_nn2_data_to_pass)

    #Computing all of the integrals, and populating the integral array
    int_1s1s_onsite = complex(integrate.nquad(onsite_1s1s_real, [[-n.inf,n.inf],[-n.inf,n.inf],[-n.inf,n.inf]])[0],  integrate.nquad(onsite_1s1s_imag, [[-n.inf,n.inf],[-n.inf,n.inf],[-n.inf,n.inf]])[0])
 
    int_1s1s_nn = complex(integrate.nquad(nn_1s1s_real, [[-n.inf,n.inf],[-n.inf,n.inf],[-n.inf,n.inf]])[0],  integrate.nquad(nn_1s1s_imag, [[-n.inf,n.inf],[-n.inf,n.inf],[-n.inf,n.inf]])[0])

    int_1s1s_nn2 = complex(integrate.nquad(nn2_1s1s_real, [[-n.inf,n.inf],[-n.inf,n.inf],[-n.inf,n.inf]])[0],  integrate.nquad(nn2_1s1s_imag, [[-n.inf,n.inf],[-n.inf,n.inf],[-n.inf,n.inf]])[0])


    orb_integral_array[q_index] = [int_1s1s_onsite, int_1s1s_nn, int_1s1s_nn2]
    print(q_index, int_1s1s_onsite, int_1s1s_nn, int_1s1s_nn2)

    
    

###Saving the integral results
with open("chain//chain_integrals_k_%d.pkl" % k_density, 'wb') as f:
    pickle.dump([orb_integral_array], f)           
###          





print("--- %s minutes ---" % ((time.time() - start_time)/60))
