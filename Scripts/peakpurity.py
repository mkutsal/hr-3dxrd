###############################################################################
### A script for checking the purity of peaks in peak lists
### with respect to a provided reference (e.g. PolyXSim output)
###
### The input requires a GFF file NOT a UBI or MAP file. You can use
### the polyxsim output gff file.
###
### M. Kutsal                                                        
### v0.1, April 2020
### DTU Physics & ESRF ID06-HXRM                                                
###############################################################################                                                        

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from ImageD11.columnfile import columnfile as cl
import ImageD11.columnfile as cl2
from ImageD11 import transformer
from sys import exit,argv
import time

start_time = time.time()

###############################################################################
## Inputs
###############################################################################
try:
    gff_ref = argv[1]
    flt_ref = argv[2]
    par_ref = argv[3]        
    ubi_1 = argv[4]
    flt_1 = argv[5]
    par_1 = argv[6]

except:
    print ("peakpurity.py [gff_ref] [flt_ref] [par_ref] [ubi_1] [flt_1] [par_1]")
    exit()
    
###############################################################################
## Few function and object definitions
###############################################################################

# Converting from ubi format to gff format using ubi_to_gff.py [code copy-pasted from 
# ImageD11 package]

def ubi_to_gff(ubi,par):
    from ImageD11 import grain as ig
    from ImageD11 import parameters as ip
    from string import split
    from xfab import tools
    from six.moves import range
    
    list_of_grains = ig.read_grain_file(ubi)
    p = ip.parameters()
    p.loadparameters(par)
    uc = [p.parameters['cell__a'],p.parameters['cell__b'],p.parameters['cell__c'],
          p.parameters['cell_alpha'],p.parameters['cell_beta'],p.parameters['cell_gamma']]
    #print(uc)
    
    grain_id = []
    x = []
    y = []
    z = []
    rodx = []
    rody = []
    rodz = []
    unitcell_a = []
    unitcell_b = []
    unitcell_c = []
    unitcell_alpha = []
    unitcell_beta = []
    unitcell_gamma = []
    U11 = []
    U12 = []
    U13 = []
    U21 = []
    U22 = []
    U23 = []
    U31 = []
    U32 = []
    U33 = []
    eps11 = []
    eps22 = []
    eps33 = []
    eps23 = []
    eps13 = []
    eps12 = []
    titles = ["grain_id","x","y","z",
              "rodx","rody","rodz",
              "U11","U12","U13","U21","U22","U23","U31","U32","U33",
              "unitcell_a","unitcell_b","unitcell_c","unitcell_alpha","unitcell_beta",
              "unitcell_gamma","eps11","eps22","eps33","eps23","eps13","eps12"]
    for i in range(len(list_of_grains)):
        grain_id.append(eval(split(list_of_grains[i].name,':')[0]))
        x.append(list_of_grains[i].translation[0]/1000.)
        y.append(list_of_grains[i].translation[1]/1000.)
        z.append(list_of_grains[i].translation[2]/1000.)
        ubi = list_of_grains[i].ubi
        (U,eps) = tools.ubi_to_u_and_eps(ubi,uc)
        uc_a, uc_b, uc_c, uc_alpha, uc_beta, uc_gamma = tools.ubi_to_cell(ubi)
        unitcell_a.append(uc_a)
        unitcell_b.append(uc_b)
        unitcell_c.append(uc_c)
        unitcell_alpha.append(uc_alpha)
        unitcell_beta.append(uc_beta)
        unitcell_gamma.append(uc_gamma)
        rod = tools.u_to_rod(U)
        rodx.append(rod[0])
        rody.append(rod[1])
        rodz.append(rod[2])
        U11.append(U[0,0])
        U12.append(U[0,1])
        U13.append(U[0,2])
        U21.append(U[1,0])
        U22.append(U[1,1])
        U23.append(U[1,2])
        U31.append(U[2,0])
        U32.append(U[2,1])
        U33.append(U[2,2])
        eps11.append(eps[0])
        eps12.append(eps[1])
        eps13.append(eps[2])
        eps22.append(eps[3])
        eps23.append(eps[4])
        eps33.append(eps[5])
    
    gff = cl2.newcolumnfile(titles)
    gff.ncols = len(titles)
    gff.nrows = len(grain_id)
    gff.bigarray = np.zeros((gff.ncols,gff.nrows))
    gff.set_attributes()
    gff.addcolumn(grain_id,"grain_id")
    gff.addcolumn(x,"x")
    gff.addcolumn(y,"y")
    gff.addcolumn(z,"z")
    gff.addcolumn(rodx,"rodx")
    gff.addcolumn(rody,"rody")
    gff.addcolumn(rodz,"rodz")
    gff.addcolumn(U11,"U11")
    gff.addcolumn(U12,"U12")
    gff.addcolumn(U13,"U13")
    gff.addcolumn(U21,"U21")
    gff.addcolumn(U22,"U22")
    gff.addcolumn(U23,"U23")
    gff.addcolumn(U31,"U31")
    gff.addcolumn(U32,"U32")
    gff.addcolumn(U33,"U33")
    gff.addcolumn(unitcell_a,"unitcell_a")
    gff.addcolumn(unitcell_b,"unitcell_b")
    gff.addcolumn(unitcell_c,"unitcell_c")
    gff.addcolumn(unitcell_alpha,"unitcell_alpha")
    gff.addcolumn(unitcell_beta,"unitcell_beta")
    gff.addcolumn(unitcell_gamma,"unitcell_gamma")    
    gff.addcolumn(eps11,"eps11")
    gff.addcolumn(eps22,"eps22")
    gff.addcolumn(eps33,"eps33")
    gff.addcolumn(eps23,"eps23")
    gff.addcolumn(eps13,"eps13")
    gff.addcolumn(eps12,"eps12")
    return gff, uc

class Umatrix(np.ndarray):
    ''' 3x3 orthogonal matrix corresponding to a rotation '''
    def __new__(cls, val=None):
        ''' 
        Create. If no argument given, create a completely random matrix, 
        otherwise use argument to initialize
        '''
        if val is None:
            self = np.zeros((3,3), dtype=np.float).view(cls)
           
        else:
            self = np.array(val, dtype=np.float).view(cls)
        return self

class mygrain:
    ''' grain object having its orientation, x,y,z positions, unit cell parameters,
        its (arbitrary) index and its assigned peaks'''

    def __init__(self,U,xyz,unitcell,grainname,peaks):
        self.U=U
        self.xyz=xyz
        self.unitcell=unitcell
        self.grainno=grainname
        self.peaks=peaks
    def __str__(self):
        return str(self.grainno)


def calc_dist(grain1,grain2):
	"""
	A function calculating the distance between two reconstructed grains
	"""
	dist=np.sqrt((grain1.xyz[0][0]-grain2.xyz[0][0])**2+(grain1.xyz[0][1]-grain2.xyz[0][1])**2
         +(grain1.xyz[0][2]-grain2.xyz[0][2])**2)
	return dist

class mypeak:
    ''' peak entry having its detector-y (px), detector-z (px), omega (deg) positions,
    g-vector, hkl and intensity information'''
    def __init__(self,detector_y,detector_z,omega,num_of_px,avg_intensity,sum_intensity,gx,gy,
                 gz,h,k,l,hr,kr,lr,diff,spot3d_id,labels,tth,eta,tth_per_grain,eta_per_grain):
        self.detector_y = detector_y
        self.detector_z = detector_z
        self.omega = omega
        self.num_of_px = num_of_px
        self.avg_intensity = avg_intensity
        self.sum_intensity = sum_intensity
        self.gx = gx
        self.gy = gy
        self.gz = gz
        self.diff = 0
        try:
            self.h = h
            self.k = k
            self.l = l
            self.hr = hr
            self.kr = kr
            self.lr = lr
        except:
            pass
        self.spot3d_id=spot3d_id
        self.labels=labels
        self.tth=tth
        self.eta=eta
        self.tth_per_grain=tth_per_grain
        self.eta_per_grain=eta_per_grain
    def __str__(self):
        return str(self.spot3d_id)

def introduce_gvectors(flt_file,par_file):
    
    obj = transformer.transformer()
    obj.loadfiltered( flt_file )
    obj.loadfileparameters( par_file )
    obj.colfile.sortby('omega')
    obj.compute_tth_eta()
    obj.addcellpeaks()
    obj.computegv()
    obj.write_colfile( str(flt_file.split(".")[0])+"_with_gvec.flt")

def get_attributes(flt_file, counter):
    #peak_name = 'peak_no_' + str( flt_file.spot3d_id[counter] )    
    detector_y = flt_file.sc[counter]
    detector_z = flt_file.fc[counter]    
    omega = flt_file.omega[counter]
    num_of_px = flt_file.Number_of_pixels[counter]
    avg_intensity = flt_file.avg_intensity[counter]
    sum_intensity = flt_file.sum_intensity[counter]
    gx = flt_file.gx[counter]
    gy = flt_file.gy[counter]
    gz = flt_file.gz[counter]
    if hasattr(flt_file,'h') == False:
        h = 1
        k = 1
        l = 1
    else:
        h = flt_file.h[counter]
        k = flt_file.k[counter]
        l = flt_file.l[counter] 
    if hasattr(flt_file,'hr') == False:
        hr = h
        kr = k
        lr = l
    else:        
        hr = flt_file.hr[counter]
        kr = flt_file.kr[counter]
        lr = flt_file.lr[counter]
    diff = 0
    spot3d_id = flt_file.spot3d_id[counter]
    labels = flt_file.labels[counter]
    tth=flt_file.tth[counter]
    eta=flt_file.eta[counter]
    if hasattr(flt_file,'tth_per_grain') == False:
        	tth_per_grain=flt_file.tth[counter]
		eta_per_grain=flt_file.eta[counter]
    else:        
        	tth_per_grain=flt_file.tth_per_grain[counter]
		eta_per_grain=flt_file.eta_per_grain[counter]
    
    return detector_y,detector_z,omega,num_of_px,avg_intensity,sum_intensity,gx
    return gy,gz,h,k,l,hr,kr,lr, diff,spot3d_id,labels,tth,eta,tth_per_grain
    return eta_per_grain

def getU_xyz_from_gff(gff_file,counter):

    grain_name = 'grain_no_' + str(	gff_file.grain_id[counter] ).split('.')[0]
    U_matrix = np.zeros((3,3))
    xyz_matrix = np.zeros((1,3))
       	
    U_matrix[0][0] = gff_file.U11[counter]
    U_matrix[0][1] = gff_file.U12[counter]
    U_matrix[0][2] = gff_file.U13[counter]
    U_matrix[1][0] = gff_file.U21[counter]
    U_matrix[1][1] = gff_file.U22[counter]
    U_matrix[1][2] = gff_file.U23[counter]
    U_matrix[2][0] = gff_file.U31[counter]
    U_matrix[2][1] = gff_file.U32[counter]
    U_matrix[2][2] = gff_file.U33[counter]
    
    xyz_matrix[0][0] = gff_file.x[counter]
    xyz_matrix[0][1] = gff_file.y[counter]
    xyz_matrix[0][2] = gff_file.z[counter]
    
    return U_matrix, xyz_matrix, grain_name

def calculate_dist_on_gvector(peak1, peak2):
    ''' A function calculating the difference of g-vectors for two peaks. 
        Idea based on CD's suggestions. '''

    dist = np.sqrt((peak1.gx-peak2.gx)**2 + 
                   (peak1.gy-peak2.gy)**2 + (peak1.gz-peak2.gz)**2 )
    return dist


def calculate_dist_on_detector(peak1, peak2):
    ''' A function calculating the difference of detector positions for 
        two peaks. MK's original idea... '''

    dist = np.sqrt((peak1.detector_y-peak2.detector_y)**2 + 
                   (peak1.detector_z-peak2.detector_z)**2  )
    return dist
 
   
def calculate_dist_on_gvector_HFP(peak1, peak2):
    ''' A function calculating the difference of g-vectors for two peaks. 
        Idea based on HFP's suggestions.'''

    g1 = np.array([peak1.gx, peak1.gy, peak1.gz])
    g2 = np.array([peak2.gx, peak2.gy, peak2.gz])
    
    gv_dot = float(g1[0]*g2[0]) + float(g1[1]*g2[1]) + 
             float(g1[2]*g2[2])
    
    norm_g1 = float(np.sqrt( float(g1[0]*g1[0]) + float(g1[1]*g1[1]) + 
                            float(g1[2]*g1[2]) ) )
    norm_g2 = float(np.sqrt( float(g2[0]*g2[0]) + float(g2[1]*g2[1]) + 
                            float(g2[2]*g2[2]) ) )
    
    
    dist =  float(np.degrees(float(np.arccos(float((gv_dot/(norm_g1*norm_g2)))))))
    
    return dist    

#####
## Open the input files
#####

gff_ref = cl(gff_ref)

if ubi_1.split(".")[-1] == "ubi" or ubi_1.split(".")[-1] == "map":
    gff_1, cellpars_1 = ubi_to_gff(ubi_1,par_1) 
elif ubi_1.split(".")[-1] == "gff":
    gff_1 = cl(ubi_1)

if flt_ref.split(".")[0].split("_")[-1] != "gvec":
	introduce_gvectors(flt_ref, par_ref)
	flt_ref = cl(str(flt_ref.split(".")[0])+"_with_gvec.flt")
else:
	flt_ref = cl(str(flt_ref))

if flt_1.split(".")[0].split("_")[-1] != "gvec":
	introduce_gvectors(flt_1, par_1)
	flt_1 = cl(str(flt_1.split(".")[0])+"_with_gvec.flt")
else:
	flt_1 = cl(str(flt_1))


#####
## Dictionaries holding the grains to be matched  and peak lists
#####

grains_ref = {}
grains_1 = {}
peaks_ref = {}
peaks_1 = {}

for i in range(len(gff_ref.U11)):
    u_temp, xyz_temp, name_temp = getU_xyz_from_gff(gff_ref,i)
    peaks_temp = []
    unitcell_temp = []
    grain_temp =  mygrain(u_temp, xyz_temp, unitcell_temp, name_temp, peaks_temp)
    grains_ref[name_temp]=grain_temp

for i in range(len(gff_1.U11)):
    u_temp, xyz_temp, name_temp = getU_xyz_from_gff(gff_1,i)
    peaks_temp = []
    unitcell_temp = []
    grain_temp = mygrain(u_temp, xyz_temp, unitcell_temp, name_temp, peaks_temp)
    grains_1[name_temp] = grain_temp

for m in range(len(flt_ref.sc)):
    detector_y_temp,detector_z_temp,omega_temp,num_of_px_temp,avg_intensity_temp,sum_intensity_temp,
    gx_temp,gy_temp,gz_temp,h_temp, k_temp, l_temp, hr_temp,kr_temp,lr_temp, diff_temp, name_temp, 
    labels_temp, tth_temp, eta_temp, tth_per_grain_temp, eta_per_grain_temp  =  get_attributes(flt_ref, m)  
    
    peak_temp  =  mypeak(detector_y_temp,detector_z_temp,omega_temp,num_of_px_temp,avg_intensity_temp,
                         sum_intensity_temp,gx_temp,gy_temp,gz_temp,h_temp, k_temp, l_temp, hr_temp,kr_temp,
                         lr_temp, diff_temp, name_temp, labels_temp, tth_temp, eta_temp, tth_per_grain_temp, 
                         eta_per_grain_temp )
    peaks_ref[name_temp]  =  peak_temp

for i in range(len(flt_1.sc)):
    detector_y_temp ,detector_z_temp ,omega_temp ,num_of_px_temp ,avg_intensity_temp ,sum_intensity_temp, 
    gx_temp ,gy_temp ,gz_temp , h_temp, k_temp, l_temp, hr_temp ,kr_temp ,lr_temp, diff_temp, name_temp, 
    labels_temp, tth_temp, eta_temp, tth_per_grain_temp, eta_per_grain_temp =  get_attributes(flt_1, i)   
    
    peak_temp  =  mypeak(detector_y_temp ,detector_z_temp ,omega_temp ,num_of_px_temp ,avg_intensity_temp,
                         sum_intensity_temp ,gx_temp ,gy_temp ,gz_temp , h_temp, k_temp, l_temp, hr_temp ,kr_temp ,
                         lr_temp, diff_temp, name_temp, labels_temp, tth_temp, eta_temp, tth_per_grain_temp, eta_per_grain_temp )
    peaks_1[name_temp]  =  peak_temp


for grain in grains_ref:
    for peak in peaks_ref:
        if peaks_ref[peak].labels == float(grains_ref[grain].grainno.split('_')[-1]):
            grains_ref[grain].peaks.append(peaks_ref[peak])

print("Reference grains are populated with their peaks")


for grain in grains_1:
    for peak in peaks_1:
        if peaks_1[peak].labels == float(grains_1[grain].grainno.split('_')[-1]):
            grains_1[grain].peaks.append(peaks_1[peak])

print("Analyzed grains are populated with their peaks")



###############################################################################
## Grain matching
###############################################################################
diff_dist_list = []
matched_grains_ref_1 = {}

from xfab.symmetry import Umis
for keys_ref in grains_ref:
        temp_listtt = []
        print()
        print()
        print("try to match grain ", keys_ref," of reference with grain list 1")
        for keys_1 in grains_1:
                    dist_diff = calc_dist(  grains_1[keys_1], grains_ref[keys_ref] )
                    test = Umis(grains_ref[keys_ref].U,grains_1[keys_1].U, 7)      #7 for cubic
                    or_diff = min(test[:,1])
                    diff_dist_list.append((or_diff,dist_diff))
                    ranking =  (dist_diff)/0.01  +  (or_diff)
                    temp_listtt.append(ranking)
                    if min(temp_listtt) == ranking  and or_diff < 0.2 and  dist_diff*1e3 < 0.5:
                        print("#####")
                        print("new minimum: ", ranking)
                        print("dist_diff: ",dist_diff*1e3)
                        print("or_diff: ",or_diff)
                        print("Ref grain: ",keys_ref,"  matched to grain list's ",keys_1)
                        print("#####")
                        matched_grains_ref_1[ grains_1[keys_1] ] = grains_ref[keys_ref]
print("#####")

completeness = []
purity_graindex = []
purity_grainspotter = []

dist_diff = 0

dummy_equal = True

all_purity_graindex_good=[]
all_purity_graindex_bad=[]
all_purity_grainspotter_good=[]
all_purity_grainspotter_bad=[]
all_completeness_good=[]
all_completeness_bad=[]


all_t_hkl_good=[]
all_t_hkl_bad=[]
all_t_tth_good=[]
all_t_tth_bad=[]
all_t_eta_good=[]
all_t_eta_bad=[]
all_t_omega_good=[]
all_t_omega_bad=[]

all_ref_peaks=[]

all_gx=[]
all_gy=[]
all_gz=[]


for keys in matched_grains_ref_1:
       
    matched_peaks = {}
    unmatched_peaks = []  
    
    
    for peak_1 in keys.peaks:
        temp_dist = []
        temp_list = []
        unmatched = True
        for peak_ref in matched_grains_ref_1[keys].peaks:
                        if peak_ref.omega < (peak_1.omega+0.2) and peak_ref.omega > (peak_1.omega-0.2):
                            dist_diff = abs(calculate_dist_on_gvector (peak_ref, peak_1))
#                            #dist_diff = abs(calculate_dist_on_detector (peak_ref, peak_1))
#                            #dist_diff = abs(calculate_dist_on_gvector_HFP(peak_ref, peak_1))
                            temp_list.append(dist_diff)
                            if min(temp_list) == dist_diff and float(dist_diff) < 0.1:
                                temp_dist.append(dist_diff)
                                matched_peaks[ peak_1 ] = peak_ref
                                unmatched = False
                    
        if unmatched == True:
            unmatched_peaks.append(peak_1)
            dummy_equal = False
            
        det_y_ref_1_diff_temp=[]
        det_z_ref_1_diff_temp=[]
        omega_ref_1_diff_temp=[]
        num_of_px_ref_1_diff_temp=[]
        avg_intensity_ref_1_diff_temp=[]
        sum_intensity_ref_1_diff_temp=[]
        gx_diff_ref_1_temp=[]
        gy_diff_ref_1_temp=[]
        gz_diff_ref_1_temp=[]
                           
        for peaks_per in matched_peaks:
            
            peaks_ref = matched_peaks[peaks_per]
            
            det_y_ref_1_diff_temp.append(float( (peaks_per.detector_y - peaks_ref.detector_y)  ))    
            det_z_ref_1_diff_temp.append(float( (peaks_per.detector_z - peaks_ref.detector_z) ))  
            omega_ref_1_diff_temp.append(float( (peaks_per.omega - peaks_ref.omega) ))  
            num_of_px_ref_1_diff_temp.append(float( (peaks_per.num_of_px - peaks_per.num_of_px) ))  
            avg_intensity_ref_1_diff_temp.append(float( (peaks_per.avg_intensity - peaks_ref.avg_intensity) ))  
            sum_intensity_ref_1_diff_temp.append(float( (peaks_per.sum_intensity - peaks_ref.sum_intensity) ))
            gx_diff_ref_1_temp.append(float( (peaks_per.gx - peaks_ref.gx) ))
            all_gx.append(float( (peaks_per.gx - peaks_ref.gx) ))
            gy_diff_ref_1_temp.append(float( (peaks_per.gy - peaks_ref.gy) ))
            gz_diff_ref_1_temp.append(float( (peaks_per.gz - peaks_ref.gz) ))

    print(np.mean(temp_dist))
    print(" ")
    print(str(len(keys.peaks)))
    print(str(len(matched_grains_ref_1[keys].peaks)))
    all_ref_peaks.append(len(matched_grains_ref_1[keys].peaks))
    print(str(len(matched_peaks)))
    print(str(len(unmatched_peaks)))
    print(" ***** ")
    
    p_gd = float( float(len(matched_peaks)) / float(len(keys.peaks)) )
    p_gs = float( float(len(matched_peaks)) / float(len(matched_grains_ref_1[keys].peaks)) )
    comp = float(float(len(keys.peaks))/len(matched_grains_ref_1[keys].peaks))
    
    purity_graindex.append( p_gd )
    purity_grainspotter.append( p_gs )    
    completeness.append(comp)

    print("purity_graindex: ",p_gd)
    print(" ")
    print("purity_grainspotter: ",p_gs)
    print(" ")
    print("completeness: ",comp)
    print(" ")
    
    if len(unmatched_peaks) < 100:
        all_purity_graindex_good.append(p_gd)
        all_purity_grainspotter_good.append(p_gs)
        all_completeness_good.append(comp)
        
    elif len(unmatched_peaks) > 100:
        all_purity_graindex_bad.append(p_gd)
        all_purity_grainspotter_bad.append(p_gs)
        all_completeness_bad.append(comp)

    def tol_hkl(peak_1, peak_ref):
        err = float(np.sqrt( float((peak_1.hr-peak_ref.h)**2) + float((peak_1.kr-peak_ref.k)**2) + float((peak_1.lr-peak_ref.l)**2) ) )
        return err

    def tol_tth(peak_1, peak_ref):
        err = float(peak_1.tth) - float(peak_ref.tth)
        return err
    
    def tol_eta(peak_1, peak_ref):
        err = float(peak_1.eta) - float(peak_ref.eta)
        return err
     
    def tol_omega(peak_1, peak_ref):
        err = float(peak_1.omega) - float(peak_ref.omega)
        return err   
    
    error_hkl=[]
    error_tth=[]
    error_eta=[]
    error_omega=[]
    
    for peak in matched_peaks:
        
        error_hkl.append( tol_hkl(peak, matched_peaks[peak]) )
        error_tth.append( tol_tth(peak, matched_peaks[peak]) )
        error_eta.append( tol_eta(peak, matched_peaks[peak]) )
        error_omega.append( tol_omega(peak, matched_peaks[peak]) )
        
        if len(unmatched_peaks) == 0 or len(unmatched_peaks) != 0:
            all_t_hkl_good.append(tol_hkl(peak, matched_peaks[peak]) )
            all_t_tth_good.append(tol_tth(peak, matched_peaks[peak]) )
            all_t_eta_good.append(tol_eta(peak, matched_peaks[peak]) )
            all_t_omega_good.append(tol_omega(peak, matched_peaks[peak]) )
       
    print("mean hkl_error: ",np.mean(error_hkl),"with stdev: ",np.std(error_hkl))
    print("mean tth_error: ",np.mean(error_tth),"with stdev: ",np.std(error_tth))
    print("mean eta_error: ",np.mean(error_eta),"with stdev: ",np.std(error_eta))
    print("mean omega_error: ",np.mean(error_omega),"with stdev: ",np.std(error_omega))
    print(" \n")
    print(" \n")
    print(" \n")
    print("****************************** \n")
    print(" \n")
    print(" \n")

print(" ")
print("Mean purity_graindex: ",np.mean(purity_graindex),
      "with stdev: ",np.std(purity_graindex))
print(" ")
print("Mean purity_grainspotter: ",np.mean(purity_grainspotter),"with stdev: ",np.std(purity_grainspotter))
print(" ")
print("Mean completeness: ",np.mean(completeness),"with stdev: ",np.std(completeness))
print(" \n")
print(" \n") 
print("Good grains:\n") 
print("Mean purity_graindex: ",np.mean(all_purity_graindex_good),"with stdev: ",np.std(all_purity_graindex_good))
print("Mean purity_grainspotter: ",np.mean(all_purity_grainspotter_good),"with stdev: ",np.std(all_purity_grainspotter_good))
print("Mean completeness: ",np.mean(all_completeness_good),"with stdev: ",np.std(all_completeness_good))
print("mean hkl_error: ",np.mean(all_t_hkl_good),"with stdev: ",np.std(all_t_hkl_good))
print("mean tth_error: ",np.mean(all_t_tth_good),"with stdev: ",np.std(all_t_tth_good))
print("mean eta_error: ",np.mean(all_t_eta_good),"with stdev: ",np.std(all_t_eta_good))
print("mean omega_error: ",np.mean(all_t_omega_good),"with stdev: ",np.std(all_t_omega_good))
print(" \n")
print(" \n")
print("Bad grains:\n") 
print("Mean purity_graindex: ",np.mean(all_purity_graindex_bad),"with stdev: ",np.std(all_purity_graindex_bad))
print("Mean purity_grainspotter: ",np.mean(all_purity_grainspotter_bad),"with stdev: ",np.std(all_purity_grainspotter_bad))
print("Mean completeness: ",np.mean(all_completeness_bad),"with stdev: ",np.std(all_completeness_bad))
print("mean hkl_error: ",np.mean(all_t_hkl_bad),"with stdev: ",np.std(all_t_hkl_bad))
print("mean tth_error: ",np.mean(all_t_tth_bad),"with stdev: ",np.std(all_t_tth_bad))
print("mean eta_error: ",np.mean(all_t_eta_bad),"with stdev: ",np.std(all_t_eta_bad))
print("mean omega_error: ",np.mean(all_t_omega_bad),"with stdev: ",np.std(all_t_omega_bad))
print(" \n")
print(" \n") 
print("Total execution time: %.3fs"%(time.time() - start_time)) 
print(" \n")
print(" \n") 
print(str(len(matched_grains_ref_1)))
print(np.mean(all_ref_peaks))


per_=[]
ref_=[]

for keys in matched_grains_ref_1:
    per_.append(int(keys.grainno.split("_")[-1]))
    ref_.append(int(matched_grains_ref_1[keys].grainno.split("_")[-1]))


print("Checking for possible duplicates")

def checkIfDuplicates_3(listOfElems):
    ''' Check if given list contains any duplicates '''    
    for elem in listOfElems:
        if listOfElems.count(elem) > 1:
            return True
    return False

result = checkIfDuplicates_3(per_)
 
if result:
    print('Yes, reference list contains duplicates')
else:
    print('No duplicates found in the reference list')  


result = checkIfDuplicates_3(ref_)
 
if result:
    print('Yes, analyzed list contains duplicates')
else:
    print('No duplicates found in the analyzed list')  

print("Duplicate check done!")


#########

filename = "matchinglist.flt"
f2 = open(filename, "w")
f2.write("\n")
f2.write("#  gs_grain ref_grain com_diff misori gs_x gs_y gs_z ref_x ref_y ref_z\n")
fmt = "%f "*4 + "%.9g "*6 + " \n"


for keys in matched_grains_ref_1.keys():
    
        com_diff = calc_dist( matched_grains_ref_1[keys], keys )*1e3
        or_diff = Umatrix()
        or_diff = or_diff.angle(keys.U,matched_grains_ref_1[keys].U)
        
        per_x = keys.xyz[0][0]*1e3
        per_y = keys.xyz[0][1]*1e3
        per_z = keys.xyz[0][2]*1e3
        
        ref_x = matched_grains_ref_1[keys].xyz[0][0]*1e3
        ref_y = matched_grains_ref_1[keys].xyz[0][1]*1e3
        ref_z = matched_grains_ref_1[keys].xyz[0][2]*1e3
        
        per_grain = float(keys.grainno.split("_")[-1])
        ref_grain = float(matched_grains_ref_1[keys].grainno.split("_")[-1])
               
        f2.write(fmt % ( per_grain, ref_grain, com_diff, or_diff, per_x, per_y, per_z, ref_x, ref_y , ref_z  ))
        
f2.close()


k=cl("matchinglist.flt")

k.sortby("gs_grain")
k.writefile("matchinglist.flt")

U_difference_ref_1 =[]
xyz_difference_ref_1 =[]

for keys in matched_grains_ref_1:
    test = Umis(keys.U,matched_grains_ref_1[keys].U, 7)      #7 for cubic
    misori = min(test[:,1])
    U_difference_ref_1.append(misori)

U_error_ref_1_mean = np.mean(U_difference_ref_1)
U_error_ref_1_stdev = np.std(U_difference_ref_1)

for keys in matched_grains_ref_1:
    xyz_difference_ref_1.append(calc_dist( keys, matched_grains_ref_1[keys] ))

xyz_error_ref_1_mean = np.mean(xyz_difference_ref_1)
xyz_error_ref_1_stdev = np.std(xyz_difference_ref_1)


U1="Average difference in U matrices per grain for reference and grain list 1 is "+str(U_error_ref_1_mean)+
    " (in degrees) with std. deviation of "+str(U_error_ref_1_stdev)
print(U1)
XYZ1 = "Average difference in CoM position of grains per grain for grain list 1 w.r.t reference is "+str(xyz_error_ref_1_mean*1000)+
    " (in micrometers) with std. deviation of "+str(xyz_error_ref_1_stdev*1000)
print(XYZ1)

print()
print()
print("Maximum U error (degrees): ", max(U_difference_ref_1))
print("Maximum CoM error (microns): ", max(xyz_difference_ref_1)*1000)


print("Unmatched grains: ")

temp_1 = set()
temp_2 = set()
temp_7 = set()
temp_8 = set()

if len(matched_grains_ref_1) == len(grains_1):
            print("All grains are matched in both grain lists!")
            dummy_equal = True
else:
    dummy_equal = False

    for keys in matched_grains_ref_1: 
        temp_1.add(keys.grainno)
        temp_2.add(matched_grains_ref_1[keys].grainno)

    for keys_ref in grains_ref:
        temp_7.add(grains_ref[keys_ref].grainno)        
    
    for keys_1 in grains_1:
        temp_8.add(grains_1[keys_1].grainno)
        
 
    unmatched_grains_1 = temp_7 - temp_1    #ref vs. grain list 1 -  ref
    unmatched_grains_2 = temp_8 - temp_2    #ref vs. grain list 1 - 1
        
    print("Reference vs. grain list 1: ")
    print("Found ",len(unmatched_grains_1)," unmatched grains in the grain list ", gff_ref.filename.split('/')[-1])
    print("Found ",len(unmatched_grains_2)," unmatched grains in the grain list ", gff_1.filename.split('/')[-1])
    print("#####")

    f8 = open("ref_not_found.gff", "w")
    f8.write("#  grain_id x y z U11 U12 U13 U21 U22 U23 U31 U32 U33\n")
    fmt = "%i " + "%9f "*12 + " \n"
    for i in unmatched_grains_1:
        grain_id = int(i.split("_")[-1])
#        i="grain_no_"+str(i)
        x=float(grains_ref[i].xyz[0][0])
        y=float(grains_ref[i].xyz[0][1])
        z=float(grains_ref[i].xyz[0][2])
        U11=float(grains_ref[i].U[0][0])
        U12=float(grains_ref[i].U[0][1])
        U13=float(grains_ref[i].U[0][2])
        U21=float(grains_ref[i].U[1][0])
        U22=float(grains_ref[i].U[1][1])
        U23=float(grains_ref[i].U[1][2])
        U31=float(grains_ref[i].U[2][0])
        U32=float(grains_ref[i].U[2][1])
        U33=float(grains_ref[i].U[2][2])
        f8.write(fmt % (grain_id, x, y, z, U11, U12, U13, U21, U22, U23, U31, U32, U33 ))
    f8.close() 
   
    print("##########################")