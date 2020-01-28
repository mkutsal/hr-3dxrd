#########
### A script for comparing two peak list files (.flt's) 
### for 3DXRD data analysis and simulations
### M. Kutsal                                                        
### October-November 2019 - v0.7
#########


import numpy as np
from ImageD11.columnfile import columnfile as cl
from ImageD11 import transformer
import sys


################################
## flt inputs
################################ 

try:
    input_flt_ref = sys.argv[1]
    input_par_ref = sys.argv[2]
    input_flt_1   = sys.argv[3]
    #input_flt_2 = cl(sys.argv[4])

except:
    print("A script for comparing a peak list file (.flt or indexed, .flt.new) with respect to a supplied reference peak list file (.flt) for 3DXRD data analysis and simulations.")
    print("PolyXSim's output .flt file does not contain g-vectors. Thus reference should be supplied...")
    print("Supported formats: .flt, .flt.new, .hdf5(?) and .par for reference parameters file.")
    print("Usage flt_mining_for_two.py [flt_ref] [par_ref] [flt_1]")
    sys.exit()

# For testing...
#input_flt_ref  = "/users/kutsal/Desktop/MKutsal/HR_3DXRD_Simulations/compare_flts/sample_data1/simulation.flt"
#input_flt_1    = cl("/users/kutsal/Desktop/MKutsal/HR_3DXRD_Simulations/compare_flts/sample_data1/perfect.flt.new")
#input_flt_2    = cl("/users/kutsal/Desktop/MKutsal/HR_3DXRD_Simulations/compare_flts/sample_data1/harvested.flt.new")
#
#input_par_ref  = "/users/kutsal/Desktop/MKutsal/HR_3DXRD_Simulations/compare_flts/sample_data1/det_z_m3px_.par"

class mypeak:
    ''' peak entry having its detector-y (px), detector-z (px), omega (deg) positions, g-vector, hkl and intensity information'''
    def __init__(self,detector_y,detector_z,omega,num_of_px,avg_intensity,sum_intensity,gx,gy,gz,spot3d_id):
        self.detector_y = detector_y
        self.detector_z = detector_z
        self.omega = omega
        self.num_of_px = num_of_px
        self.avg_intensity = avg_intensity
        self.sum_intensity = sum_intensity
        self.gx = gx
        self.gy = gy
        self.gz = gz
        self.spot3d_id=spot3d_id
    def __str__(self):
        return str(self.spot3d_id)

def intorduce_gvectors(flt_file,par_file):
    
    obj = transformer.transformer()
    obj.loadfiltered( flt_file )
    obj.loadfileparameters( par_file )
    obj.compute_tth_eta()
    obj.addcellpeaks()
    obj.computegv()
    obj.write_colfile( str(flt_file.split(".")[0])+"_with_gvec.flt")
    
def get_attributes(flt_file, counter):
    peak_name = 'peak_no_' + str( flt_file.spot3d_id[counter] )    
    detector_y = flt_file.sc[counter]
    detector_z = flt_file.fc[counter]    
    omega = flt_file.omega[counter]
    num_of_px = flt_file.Number_of_pixels[counter]
    avg_intensity = flt_file.avg_intensity[counter]
    sum_intensity = flt_file.sum_intensity[counter]
    gx = flt_file.gx[counter]
    gy = flt_file.gy[counter]
    gz = flt_file.gz[counter]
    
    return detector_y,detector_z,omega,num_of_px,avg_intensity,sum_intensity,gx,gy,gz,peak_name

def calculate_dist_on_det(peak1, peak2):
    ''' A function calculating the difference of detector positions for two peaks (in pixel units) '''
    
    len_peak_1 = np.sqrt(peak1.detector_y * peak1.detector_y + peak1.detector_z * peak1.detector_z)
    len_peak_2 = np.sqrt(peak2.detector_y * peak2.detector_y + peak2.detector_z * peak2.detector_z)
    dist=np.subtract(len_peak_1,len_peak_2)
    return dist

    
def calculate_dist_on_gvector(peak1, peak2):
    ''' A function calculating the difference of g-vectors for two peaks '''
    len_peak_1 = np.sqrt(peak1.gx * peak1.gx + peak1.gy * peak1.gy + peak1.gz * peak1.gz)
    len_peak_2 = np.sqrt(peak2.gx * peak2.gx + peak2.gy * peak2.gy + peak2.gz * peak2.gz)
    dist=np.subtract(len_peak_1,len_peak_2)
    return dist

################################
## Peak matching
################################

##Dictionary holding the supplied peak objects
flt_ref = {}
flt_1 = {}
flt_2 = {}

intorduce_gvectors(input_flt_ref, input_par_ref)
intorduce_gvectors(input_flt_1, input_par_ref)
input_flt_ref = cl(str(input_flt_ref.split(".")[0])+"_with_gvec.flt")
input_flt_1 = cl(str(input_flt_1.split(".")[0])+"_with_gvec.flt")

for m in range(len(input_flt_ref.sc)):
    detector_y_temp,detector_z_temp,omega_temp,num_of_px_temp,avg_intensity_temp,sum_intensity_temp,gx_temp,gy_temp,gz_temp,peak_name_temp  =  get_attributes(input_flt_ref, m)   
    peak_temp  =  mypeak(detector_y_temp,detector_z_temp,omega_temp,num_of_px_temp,avg_intensity_temp,sum_intensity_temp,gx_temp,gy_temp,gz_temp,peak_name_temp)
    flt_ref[peak_name_temp]  =  peak_temp

for i in range(len(input_flt_1.sc)):
    detector_y_temp ,detector_z_temp ,omega_temp ,num_of_px_temp ,avg_intensity_temp ,sum_intensity_temp ,gx_temp ,gy_temp ,gz_temp , peak_name_temp  =  get_attributes(input_flt_1, i)   
    peak_temp  =  mypeak(detector_y_temp ,detector_z_temp ,omega_temp ,num_of_px_temp ,avg_intensity_temp ,sum_intensity_temp ,gx_temp ,gy_temp ,gz_temp , peak_name_temp)
    flt_1[peak_name_temp]  =  peak_temp
    
#for k in range(len(input_flt_2.sc)):
#    detector_y_temp,detector_z_temp,omega_temp,num_of_px_temp,avg_intensity_temp,sum_intensity_temp,gx_temp,gy_temp,gz_temp, h_temp, k_temp, l_temp, hr_temp,kr_temp,lr_temp,peak_name_temp  =  get_attributes(input_flt_2, k)   
#    peak_temp  =  mypeak(detector_y_temp,detector_z_temp,omega_temp,num_of_px_temp,avg_intensity_temp,sum_intensity_temp,gx_temp,gy_temp,gz_temp,h_temp, k_temp, l_temp, hr_temp,kr_temp,lr_temp,peak_name_temp)
#    flt_2[peak_name_temp]  =  peak_temp

##Dictionary holding the matched peaks
matched_peaks_ref_1 = {}
#matched_peaks_ref_2 = {}
#matched_peaks_1_2 = {}

for keys_ref in flt_ref:
    temp_list = []
    for keys_1 in flt_1:
        dist_diff = abs(calculate_dist_on_det ( flt_ref[keys_ref] , flt_1[keys_1] ) )
        #dist_diff = abs(calculate_dist_on_gvector ( flt_ref[keys_ref] , flt_1[keys_1] ) )
        temp_list.append(dist_diff)
        if min(temp_list) == dist_diff:
            matched_peaks_ref_1[ flt_ref[keys_ref] ] = flt_1[keys_1]
	    print("Hit!")



################################
## Analysis
################################

print("Comparison of two peak lists: ")
print(" ")

#####
## Number of matched peaks, listing unmatched peaks
#####

print("Number of peaks & number of matched peaks")
print
print("Number of peaks in the first flt: ",len(flt_1))
print
print("Unmatched peaks: ")
print

unmatched_flt_ref_1 = []


if len(flt_ref) == len(flt_1): # and len(flt_ref) == len(flt_2):
    print("All peaks are matched in both flt's!")
    dummy_equal = True
else:
    dummy_equal = False
    for keys in matched_peaks_ref_1:
        for keys_1 in flt_1:
            if matched_peaks_ref_1[keys].spot3d_id == flt_1[keys_1].spot3d_id:
                pass
            else:
                unmatched_flt_ref_1.append(keys_1)


    f5= "Found "+str(len(unmatched_flt_ref_1))+" peaks in the peak list "+str(input_flt_1.filename)+" and the reference"
    print(f5)


    
#####
## Peak position and intensity
#####
print
print("#######################")
print
print("Peak position & intensity")
print
    
det_y_ref_1_diff_temp=[]
det_z_ref_1_diff_temp=[]
omega_ref_1_diff_temp=[]
num_of_px_ref_1_diff_temp=[]
avg_intensity_ref_1_diff_temp=[]
sum_intensity_ref_1_diff_temp=[]

#det_y_ref_2_diff_temp=[]
#det_z_ref_2_diff_temp=[]
#omega_ref_2_diff_temp=[]
#num_of_px_ref_2_diff_temp=[]
#avg_intensity_ref_2_diff_temp=[]
#sum_intensity_ref_2_diff_temp=[]
#
#det_y_1_2_diff_temp=[]
#det_z_1_2_diff_temp=[]
#omega_1_2_diff_temp=[]
#num_of_px_1_2_diff_temp=[]
#avg_intensity_1_2_diff_temp=[]
#sum_intensity_1_2_diff_temp=[]

for keys in matched_peaks_ref_1:
    det_y_ref_1_diff_temp.append(float( (keys.detector_y - matched_peaks_ref_1[keys].detector_y) / keys.detector_y )*100)    
    det_z_ref_1_diff_temp.append(float( (keys.detector_z - matched_peaks_ref_1[keys].detector_z) / keys.detector_z )*100)
    omega_ref_1_diff_temp.append(float( (keys.omega - matched_peaks_ref_1[keys].omega) / keys.omega )*100)
    num_of_px_ref_1_diff_temp.append(float( (keys.num_of_px - matched_peaks_ref_1[keys].num_of_px) / keys.num_of_px )*100)
    avg_intensity_ref_1_diff_temp.append(float( (keys.avg_intensity - matched_peaks_ref_1[keys].avg_intensity) / keys.avg_intensity )*100)
    sum_intensity_ref_1_diff_temp.append(float( (keys.sum_intensity - matched_peaks_ref_1[keys].sum_intensity) / keys.sum_intensity )*100)
    
    
#for keys in matched_peaks_ref_2:
#    det_y_ref_2_diff_temp.append(float( (keys.detector_y - matched_peaks_ref_2[keys].detector_y) / keys.detector_y )*100)    
#    det_z_ref_2_diff_temp.append(float( (keys.detector_z - matched_peaks_ref_2[keys].detector_z) / keys.detector_z )*100)
#    omega_ref_2_diff_temp.append(float( (keys.omega - matched_peaks_ref_2[keys].omega) / keys.omega )*100)
#    num_of_px_ref_2_diff_temp.append(float( (keys.num_of_px - matched_peaks_ref_2[keys].num_of_px) / keys.num_of_px )*100)
#    avg_intensity_ref_2_diff_temp.append(float( (keys.avg_intensity - matched_peaks_ref_2[keys].avg_intensity) / keys.avg_intensity )*100)
#    sum_intensity_ref_2_diff_temp.append(float( (keys.sum_intensity - matched_peaks_ref_2[keys].sum_intensity) / keys.sum_intensity )*100)
#
#for keys in matched_peaks_1_2:
#    det_y_1_2_diff_temp.append(float( (keys.detector_y - matched_peaks_1_2[keys].detector_y) / keys.detector_y )*100)    
#    det_z_1_2_diff_temp.append(float( (keys.detector_z - matched_peaks_1_2[keys].detector_z) / keys.detector_z )*100)
#    omega_1_2_diff_temp.append(float( (keys.omega - matched_peaks_1_2[keys].omega) / keys.omega )*100)
#    num_of_px_1_2_diff_temp.append(float( (keys.num_of_px - matched_peaks_1_2[keys].num_of_px) / keys.num_of_px )*100)
#    avg_intensity_1_2_diff_temp.append(float( (keys.avg_intensity - matched_peaks_1_2[keys].avg_intensity) / keys.avg_intensity )*100)
#    sum_intensity_1_2_diff_temp.append(float( (keys.sum_intensity - matched_peaks_1_2[keys].sum_intensity) / keys.sum_intensity )*100)


det_y_ref_1_error_mean = np.mean(det_y_ref_1_diff_temp)
det_y_ref_1_error_stdev = np.std(det_y_ref_1_diff_temp)
#det_y_ref_2_error_mean = np.mean(det_y_ref_2_diff_temp)
#det_y_ref_2_error_stdev = np.std(det_y_ref_2_diff_temp)
#det_y_1_2_error_mean = np.mean(det_y_1_2_diff_temp)
#det_y_1_2_error_stdev = np.std(det_y_1_2_diff_temp)

det_z_ref_1_error_mean = np.mean(det_z_ref_1_diff_temp)
det_z_ref_1_error_stdev = np.std(det_z_ref_1_diff_temp)
#det_z_ref_2_error_mean = np.mean(det_z_ref_2_diff_temp)
#det_z_ref_2_error_stdev = np.std(det_z_ref_2_diff_temp)
#det_z_1_2_error_mean = np.mean(det_z_1_2_diff_temp)
#det_z_1_2_error_stdev = np.std(det_z_1_2_diff_temp)

omega_ref_1_error_mean = np.mean(omega_ref_1_diff_temp)
omega_ref_1_error_stdev = np.std(omega_ref_1_diff_temp)
#omega_ref_2_error_mean = np.mean(omega_ref_2_diff_temp)
#omega_ref_2_error_stdev = np.std(omega_ref_2_diff_temp)
#omega_1_2_error_mean = np.mean(omega_1_2_diff_temp)
#omega_1_2_error_stdev = np.std(omega_1_2_diff_temp)

num_of_px_ref_1_error_mean = np.mean(num_of_px_ref_1_diff_temp)
num_of_px_ref_1_error_stdev = np.std(num_of_px_ref_1_diff_temp)
#num_of_px_ref_2_error_mean = np.mean(num_of_px_ref_2_diff_temp)
#num_of_px_ref_2_error_stdev = np.std(num_of_px_ref_2_diff_temp)
#num_of_px_1_2_error_mean = np.mean(num_of_px_1_2_diff_temp)
#num_of_px_1_2_error_stdev = np.std(num_of_px_1_2_diff_temp)

avg_intensity_ref_1_error_mean = np.mean(avg_intensity_ref_1_diff_temp)
avg_intensity_ref_1_error_stdev = np.std(avg_intensity_ref_1_diff_temp)
#avg_intensity_ref_2_error_mean = np.mean(avg_intensity_ref_2_diff_temp)
#avg_intensity_ref_2_error_stdev = np.std(avg_intensity_ref_2_diff_temp)
#avg_intensity_1_2_error_mean = np.mean(avg_intensity_1_2_diff_temp)
#avg_intensity_1_2_error_stdev = np.std(avg_intensity_1_2_diff_temp)

sum_intensity_ref_1_error_mean = np.mean(sum_intensity_ref_1_diff_temp)
sum_intensity_ref_1_error_stdev = np.std(sum_intensity_ref_1_diff_temp)
#sum_intensity_ref_2_error_mean = np.mean(sum_intensity_ref_2_diff_temp)
#sum_intensity_ref_2_error_stdev = np.std(sum_intensity_ref_2_diff_temp)
#sum_intensity_1_2_error_mean = np.mean(sum_intensity_1_2_diff_temp)
#sum_intensity_1_2_error_stdev = np.std(sum_intensity_1_2_diff_temp)

pk1_1="Percent average error in detector-y position per peak for reference and peak list 1 is "+str(det_y_ref_1_error_mean)+"% with std. deviation of "+str(det_y_ref_1_error_stdev)
print(pk1_1)
#pk1_2="Percent average error in detector-y position per peak for reference and peak list 2 is "+str(det_y_ref_2_error_mean)+"% with std. deviation of "+str(det_y_ref_2_error_stdev)
#print(pk1_2)
#pk1_3="Percent average error in detector-y position per peak for peak list 1 and peak list 2 is "+str(det_y_1_2_error_mean)+"% with std. deviation of "+str(det_y_1_2_error_stdev)
#print(pk1_3)

pk2_1="Percent average error in detector-z position per peak for reference and peak list 1 is "+str(det_z_ref_1_error_mean)+"% with std. deviation of "+str(det_z_ref_1_error_stdev)
print(pk2_1)
#pk2_2="Percent average error in detector-z position per peak for reference and peak list 2 is "+str(det_z_ref_2_error_mean)+"% with std. deviation of "+str(det_z_ref_2_error_stdev)
#print(pk2_2)
#pk2_3="Percent average error in detector-z position per peak for peak list 1 and peak list 2 is "+str(det_z_1_2_error_mean)+"% with std. deviation of "+str(det_z_1_2_error_stdev)
#print(pk2_3)

pk3_1="Percent average error in omega position per peak for reference and peak list 1 is "+str(omega_ref_1_error_mean)+"% with std. deviation of "+str(omega_ref_1_error_stdev)
print(pk3_1)
#pk3_2="Percent average error in omega position per peak for reference and peak list 2 is "+str(omega_ref_2_error_mean)+"% with std. deviation of "+str(omega_ref_2_error_stdev)
#print(pk3_2)
#pk3_3="Percent average error in omega position per peak for peak list 1 and peak list 2 is "+str(omega_1_2_error_mean)+"% with std. deviation of "+str(omega_1_2_error_stdev)
#print(pk3_3)

pk4_1="Percent average error in number of pixels per peak for reference and peak list 1 is "+str(num_of_px_ref_1_error_mean)+"% with std. deviation of "+str(num_of_px_ref_1_error_stdev)
print(pk4_1)
#pk4_2="Percent average error in number of pixels per peak for reference and peak list 2 is "+str(num_of_px_ref_2_error_mean)+"% with std. deviation of "+str(num_of_px_ref_2_error_stdev)
#print(pk4_2)
#pk4_3="Percent average error in number of pixels per peak for peak list 1 and peak list 2 is "+str(num_of_px_1_2_error_mean)+"% with std. deviation of "+str(num_of_px_1_2_error_stdev)
#print(pk4_3)


pk5_1="Percent average error in average intensity per peak for reference and peak list 1 is "+str(avg_intensity_ref_1_error_mean)+"% with std. deviation of "+str(avg_intensity_ref_1_error_stdev)
print(pk5_1)
#pk5_2="Percent average error in average intensity per peak for reference and peak list 2 is "+str(avg_intensity_ref_2_error_mean)+"% with std. deviation of "+str(avg_intensity_ref_2_error_stdev)
#print(pk5_2)
#pk5_3="Percent average error in average intensity per peak for peak list1 and peak list 2 is "+str(avg_intensity_1_2_error_mean)+"% with std. deviation of "+str(avg_intensity_1_2_error_stdev)
#print(pk5_3)

pk6_1="Percent average error in sum intensity per peak for reference and peak list 1 is "+str(sum_intensity_ref_1_error_mean)+"% with std. deviation of "+str(sum_intensity_ref_1_error_stdev)
print(pk6_1)
#pk6_2="Percent average error in sum intensity per peak for reference and peak list 2 is "+str(sum_intensity_ref_2_error_mean)+"% with std. deviation of "+str(sum_intensity_ref_2_error_stdev)
#print(pk6_2)
#pk6_3="Percent average error in sum intensity per peak for peak list 1 and peak list 2 is "+str(sum_intensity_1_2_error_mean)+"% with std. deviation of "+str(sum_intensity_1_2_error_stdev)
#print(pk6_3)
#      
#####
## g-vector and hkl
#####       
print
print
print("#######################")
print
print("g-vector & hkl")
print


gx_diff_ref_1_temp=[]
gy_diff_ref_1_temp=[]
gz_diff_ref_1_temp=[]
h_diff_ref_1_temp=[]
k_diff_ref_1_temp=[]
l_diff_ref_1_temp=[]
#
#gx_diff_ref_2_temp=[]
#gy_diff_ref_2_temp=[]
#gz_diff_ref_2_temp=[]
#h_diff_ref_2_temp=[]
#k_diff_ref_2_temp=[]
#l_diff_ref_2_temp=[]
#
#gx_diff_1_2_temp=[]
#gy_diff_1_2_temp=[]
#gz_diff_1_2_temp=[]
#h_diff_1_2_temp=[]
#k_diff_1_2_temp=[]
#l_diff_1_2_temp=[]

for keys in matched_peaks_ref_1:
    gx_diff_ref_1_temp.append( float( (keys.gx - matched_peaks_ref_1[keys].gx) / keys.gx ) *100)
    gy_diff_ref_1_temp.append( float( (keys.gy - matched_peaks_ref_1[keys].gy) / keys.gy ) *100)
    gz_diff_ref_1_temp.append( float( (keys.gz - matched_peaks_ref_1[keys].gz) / keys.gz ) *100)
    
    
    

#for keys in matched_peaks_ref_2:
#    gx_diff_ref_2_temp.append( float( (keys.gx - matched_peaks_ref_2[keys].gx) / keys.gx ) *100)
#    gy_diff_ref_2_temp.append( float( (keys.gy - matched_peaks_ref_2[keys].gy) / keys.gy ) *100)
#    gz_diff_ref_2_temp.append( float( (keys.gz - matched_peaks_ref_2[keys].gz) / keys.gz ) *100)
#    
#    
#    h_err_1 = (keys.h - keys.hr)
#    h_err_1 = float("{:.6f}".format(h_err_1))
#    h_err_2 = (matched_peaks_ref_2[keys].h - matched_peaks_ref_2[keys].hr)
#    h_err_2 = float("{:.6f}".format(h_err_2))
#    k_err_1 = (keys.k - keys.kr)
#    k_err_1 = float("{:.6f}".format(k_err_1))
#    k_err_2 = (matched_peaks_ref_2[keys].k - matched_peaks_ref_2[keys].kr)
#    k_err_2 = float("{:.6f}".format(k_err_2))
#    l_err_1 = (keys.l - keys.lr)
#    l_err_1 = float("{:.6f}".format(h_err_1))
#    l_err_2 = (matched_peaks_ref_2[keys].l - matched_peaks_ref_2[keys].lr)
#    l_err_2 = float("{:.6f}".format(l_err_2))
#
#    h_diff_ref_2_temp.append( (h_err_1-h_err_2) ) 
#    k_diff_ref_2_temp.append( (k_err_1-k_err_2) )
#    l_diff_ref_2_temp.append( (l_err_1-l_err_2) )
#
#for keys in matched_peaks_1_2:
#    gx_diff_1_2_temp.append( float( (keys.gx - matched_peaks_1_2[keys].gx) / keys.gx ) *100)
#    gy_diff_1_2_temp.append( float( (keys.gy - matched_peaks_1_2[keys].gy) / keys.gy ) *100)
#    gz_diff_1_2_temp.append( float( (keys.gz - matched_peaks_1_2[keys].gz) / keys.gz ) *100)
#    
#    
#    h_err_1 = (keys.h - keys.hr)
#    h_err_1 = float("{:.6f}".format(h_err_1))
#    h_err_2 = (matched_peaks_1_2[keys].h - matched_peaks_1_2[keys].hr)
#    h_err_2 = float("{:.6f}".format(h_err_2))
#    k_err_1 = (keys.k - keys.kr)
#    k_err_1 = float("{:.6f}".format(k_err_1))
#    k_err_2 = (matched_peaks_1_2[keys].k - matched_peaks_1_2[keys].kr)
#    k_err_2 = float("{:.6f}".format(k_err_2))
#    l_err_1 = (keys.l - keys.lr)
#    l_err_1 = float("{:.6f}".format(h_err_1))
#    l_err_2 = (matched_peaks_1_2[keys].l - matched_peaks_1_2[keys].lr)
#    l_err_2 = float("{:.6f}".format(l_err_2))
#
#    h_diff_1_2_temp.append( (h_err_1-h_err_2) ) 
#    k_diff_1_2_temp.append( (k_err_1-k_err_2) )
#    l_diff_1_2_temp.append( (l_err_1-l_err_2) )
#    


gx_error_ref_1_mean = np.mean(gx_diff_ref_1_temp)
gx_error_ref_1_stdev = np.std(gx_diff_ref_1_temp)
#gx_error_ref_2_mean = np.mean(gx_diff_ref_2_temp)
#gx_error_ref_2_stdev = np.std(gx_diff_ref_2_temp)
#gx_error_1_2_mean = np.mean(gx_diff_1_2_temp)
#gx_error_1_2_stdev = np.std(gx_diff_1_2_temp)

gy_error_ref_1_mean = np.mean(gy_diff_ref_1_temp)
gy_error_ref_1_stdev = np.std(gy_diff_ref_1_temp)
#gy_error_ref_2_mean = np.mean(gy_diff_ref_2_temp)
#gy_error_ref_2_stdev = np.std(gy_diff_ref_2_temp)
#gy_error_1_2_mean = np.mean(gy_diff_1_2_temp)
#gy_error_1_2_stdev = np.std(gy_diff_1_2_temp)

gz_error_ref_1_mean = np.mean(gz_diff_ref_1_temp)
gz_error_ref_1_stdev = np.std(gz_diff_ref_1_temp)
#gz_error_ref_2_mean = np.mean(gz_diff_ref_2_temp)
#gz_error_ref_2_stdev = np.std(gz_diff_ref_2_temp)
#gz_error_1_2_mean = np.mean(gz_diff_1_2_temp)
#gz_error_1_2_stdev = np.std(gz_diff_1_2_temp)


print("###")
gv1_1="Percent average error in gx position per peak for peak list 1 w.r.t. reference is "+str(gx_error_ref_1_mean)+"% with std. deviation of "+str(gx_error_ref_1_stdev)
print(gv1_1)
#gv1_2="Percent average error in gx position per peak for peak list 2 w.r.t. reference is "+str(gx_error_ref_2_mean)+"% with std. deviation of "+str(gx_error_ref_2_stdev)
#print(gv1_2)
#gv1_3="Percent average error in gx position per peak for peak list 1 and peak list 2 is "+str(gx_error_1_2_mean)+"% with std. deviation of "+str(gx_error_1_2_stdev)
#print(gv1_3)
print("###")
gv2_1="Percent average error in gy position per peak for peak list 1 w.r.t. reference is "+str(gy_error_ref_1_mean)+"% with std. deviation of "+str(gy_error_ref_1_stdev)
print(gv2_1)
#gv2_2="Percent average error in gy position per peak for peak list 2 w.r.t. reference is "+str(gy_error_ref_2_mean)+"% with std. deviation of "+str(gy_error_ref_2_stdev)
#print(gv2_2)
#gv2_3="Percent average error in gy position per peak for peak list 1 and peak list 2 is "+str(gy_error_1_2_mean)+"% with std. deviation of "+str(gy_error_1_2_stdev)
#print(gv2_3)
print("###")
gv3_1="Percent average error in gz position per peak for peak list 1 w.r.t. reference is "+str(gz_error_ref_1_mean)+"% with std. deviation of "+str(gz_error_ref_1_stdev)
print(gv3_1)
#gv3_2="Percent average error in gz position per peak for peak list 2 w.r.t. reference is "+str(gz_error_ref_2_mean)+"% with std. deviation of "+str(gz_error_ref_2_stdev)
#print(gv3_2)
#gv3_3="Percent average error in gz position per peak for peak list 1 and peak list 2 is "+str(gz_error_1_2_mean)+"% with std. deviation of "+str(gz_error_1_2_stdev)
#print(gv3_3)
print("###")



print
print("#######################")
print




#####
## Output
##### 


#filename= "FLT_Comparison_"+str(input_flt_1.filename.split("/")[-1].split(".")[0])+"_and_"+str(str(input_flt_2.filename.split("/")[-1].split(".")[0]))+".txt"
filename= "FLT_Comparison_"+str(input_flt_ref.filename.split("/")[-1].split(".")[0])+"_and_"+str(str(input_flt_1.filename.split("/")[-1].split(".")[0]))+".txt"


f = open(filename,"w")
f.writelines("Comparison of two peak lists: "+'\n')
f.writelines(" "+'\n')
f.writelines("Input flt files: "+'\n')
f0="--- "+input_flt_ref.filename
f.writelines(f0+'\n')
f.writelines(" "+'\n')
f1="--- "+input_flt_1.filename
f.writelines(f1+'\n')
f.writelines(" "+'\n')
#f2="--- "+input_flt_2.filename
#f.writelines(f2+'\n')
f.writelines(" "+'\n')
f.writelines("#######################"+'\n')
f.writelines(" "+'\n')
f.writelines("Number of peaks & number of matched peaks"+'\n')
f.writelines(" "+'\n')

f3="Number of peaks in the first flt: "+str(len(flt_1))
#f4="Number of peaks in the second flt: "+str(len(flt_2))
f.writelines(f3+'\n')
#f.writelines(f4+'\n')
f.writelines(" "+'\n')
f.writelines("Unmatched peaks: "+'\n')
if dummy_equal==True:
    f.writelines("All peaks are matched in both flt's!"+'\n')
    f.writelines(" "+'\n')
elif dummy_equal==False:
    f.writelines(" "+'\n')
    f.writelines(f5+'\n')
    f.writelines(" "+'\n')
#    f.writelines(f6+'\n')
#    f.writelines(" "+'\n')
#    f.writelines(f7+'\n')
#    f.writelines(" "+'\n')
#    f.writelines(f8+'\n')
#    f.writelines(" "+'\n')
f.writelines("#######################"+'\n')
f.writelines(" "+'\n')
f.writelines("Peak position & intensity"+'\n')
f.writelines(" "+'\n')
f.writelines(pk1_1+'\n')
f.writelines(" "+'\n')
#f.writelines(pk1_2+'\n')
#f.writelines(" "+'\n')
#f.writelines(pk1_3+'\n')
#f.writelines(" "+'\n')
f.writelines("###"+'\n')
f.writelines(pk2_1+'\n')
f.writelines(" "+'\n')
#f.writelines(pk2_2+'\n')
#f.writelines(" "+'\n')
#f.writelines(pk2_3+'\n')
#f.writelines(" "+'\n')
f.writelines("###"+'\n')
f.writelines(pk3_1+'\n')
f.writelines(" "+'\n')
#f.writelines(pk3_2+'\n')
#f.writelines(" "+'\n')
#f.writelines(pk3_3+'\n')
#f.writelines(" "+'\n')
#f.writelines("###"+'\n')
f.writelines(pk4_1+'\n')
f.writelines(" "+'\n')
#f.writelines(pk4_2+'\n')
#f.writelines(" "+'\n')
#f.writelines(pk4_3+'\n')
#f.writelines(" "+'\n')
f.writelines("###"+'\n')
f.writelines(pk5_1+'\n')
f.writelines(" "+'\n')
#f.writelines(pk5_2+'\n')
#f.writelines(" "+'\n')
#f.writelines(pk5_3+'\n')
#f.writelines(" "+'\n')
#f.writelines("###"+'\n')
f.writelines(pk6_1+'\n')
f.writelines(" "+'\n')
#f.writelines(pk6_2+'\n')
#f.writelines(" "+'\n')
#f.writelines(pk6_3+'\n')
#f.writelines(" "+'\n')
f.writelines("#######################"+'\n')
f.writelines(" "+'\n')
f.writelines("g-vector & hkl"+'\n')
f.writelines(" "+'\n')
f.writelines(gv1_1+'\n')
f.writelines(" "+'\n')
#f.writelines(gv1_2+'\n')
#f.writelines(" "+'\n')
#f.writelines(gv1_3+'\n')
#f.writelines(" "+'\n')
f.writelines("###"+'\n')
f.writelines(gv2_1+'\n')
f.writelines(" "+'\n')
#f.writelines(gv2_2+'\n')
#f.writelines(" "+'\n')
#f.writelines(gv2_3+'\n')
#f.writelines(" "+'\n')
f.writelines("###"+'\n')
f.writelines(gv3_1+'\n')
f.writelines(" "+'\n')
#f.writelines(gv3_2+'\n')
#f.writelines(" "+'\n')
#f.writelines(gv3_3+'\n')
#f.writelines(" "+'\n')
f.writelines("#######################"+'\n')
f.close()

