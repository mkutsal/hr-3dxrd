###############################################################################
### A script for comparing two 3DXRD outputs 
### with respect to a provided reference (e.g. PolyXSim output)
### M. Kutsal, C. Detlefs                                                        
### v0.5, October 2019
### DTU Physics & ESRF ID06-HXRM                                                
###############################################################################                                                        

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from ImageD11 import columnfile as cl
import sys

###############################################################################
## Inputs
###############################################################################
try:
    if str(sys.argv[1]).split(".")[1] == 'ubi':
        ubi_ref = sys.argv[1]
        par_ref = sys.argv[2]        
        ubi_1 = sys.argv[3]
        par_1 = sys.argv[4]
        ubi_2 = sys.argv[5]
        par_2 = sys.argv[6]
        file_format_ubi = True
        
    elif str(sys.argv[1]).split(".")[1] == 'map':
        ubi_ref = sys.argv[1]
        par_ref = sys.argv[2]        
        ubi_1 = sys.argv[3]
        par_1 = sys.argv[4]
        ubi_2 = sys.argv[5]
        par_2 = sys.argv[6]
        file_format_ubi = True
        
    elif str(sys.argv[1]).split(".")[1] == 'gff':
        input_gff_ref = sys.argv[1]        
        input_gff_1 = sys.argv[2]
        input_gff_2 = sys.argv[3]
        file_format_ubi = False
except:
    print(" Usage compare_ubis_3.py [ubi_ref] [par_ref] [ubi_1] [par_1] [ubi_2] [par_2] OR compare_ubis_2.py [gff_ref] [gff_1] [gff_2]")
    sys.exit()


###############################################################################
## Few function and object definitions
###############################################################################

#Converting from ubi format to gff format using ubi_to_gff.py [code copy-pasted from ImageD11 package]

def ubi_to_gff(ubi,par):
    from ImageD11 import grain as ig
    from ImageD11 import parameters as ip
    from string import split
    from xfab import tools
    from six.moves import range
    
    list_of_grains = ig.read_grain_file(ubi)
    p = ip.parameters()
    p.loadparameters(par)
    uc = [p.parameters['cell__a'],p.parameters['cell__b'],p.parameters['cell__c'],p.parameters['cell_alpha'],p.parameters['cell_beta'],p.parameters['cell_gamma']]
    #print(uc)
    
    grainno = []
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
    titles = ["grainno","x","y","z",
              "rodx","rody","rodz",
              "U11","U12","U13","U21","U22","U23","U31","U32","U33",
              "unitcell_a","unitcell_b","unitcell_c","unitcell_alpha","unitcell_beta","unitcell_gamma",
              "eps11","eps22","eps33","eps23","eps13","eps12"]
    for i in range(len(list_of_grains)):
        if hasattr(list_of_grains,'name') == False:
            grainno.append(i)
        elif hasattr(list_of_grains,'name') == True:
            grainno.append(eval(split(list_of_grains[i].name,':')[0]))
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
    
    gff = cl.newcolumnfile(titles)
    gff.ncols = len(titles)
    gff.nrows = len(grainno)
    gff.bigarray = np.zeros((gff.ncols,gff.nrows))
    gff.set_attributes()
    gff.addcolumn(grainno,"grainno")
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

# 24 cubic symmtry operations, point group 432 (no inversion symmetry)
CUBIC_SYMMETRIES = (
    np.array(((-1, 0, 0), (0, 0, 1), (0, 1, 0)), dtype=np.float),
    np.array(((0, 1, 0), (0, 0, -1), (-1, 0, 0)), dtype=np.float),
    np.array(((1, 0, 0), (0, -1, 0), (0, 0, -1)), dtype=np.float),
    np.array(((1, 0, 0), (0, 0, 1), (0, -1, 0)), dtype=np.float),
    np.array(((0, 1, 0), (1, 0, 0), (0, 0, -1)), dtype=np.float),
    np.array(((0, 0, -1), (0, 1, 0), (1, 0, 0)), dtype=np.float),
    np.array(((0, -1, 0), (-1, 0, 0), (0, 0, -1)), dtype=np.float),
    np.array(((-1, 0, 0), (0, 1, 0), (0, 0, -1)), dtype=np.float),
    np.array(((0, 1, 0), (-1, 0, 0), (0, 0, 1)), dtype=np.float),
    np.array(((0, 0, 1), (1, 0, 0), (0, 1, 0)), dtype=np.float),
    np.array(((-1, 0, 0), (0, 0, -1), (0, -1, 0)), dtype=np.float),
    np.array(((0, -1, 0), (1, 0, 0), (0, 0, 1)), dtype=np.float),
    np.array(((0, -1, 0), (0, 0, -1), (1, 0, 0)), dtype=np.float),
    np.array(((0, 0, 1), (-1, 0, 0), (0, -1, 0)), dtype=np.float),
    np.array(((0, 0, 1), (0, 1, 0), (-1, 0, 0)), dtype=np.float),
    np.array(((1, 0, 0), (0, 1, 0), (0, 0, 1)), dtype=np.float),
    np.array(((0, -1, 0), (0, 0, 1), (-1, 0, 0)), dtype=np.float),
    np.array(((1, 0, 0), (0, 0, -1), (0, 1, 0)), dtype=np.float),
    np.array(((0, 0, -1), (-1, 0, 0), (0, 1, 0)), dtype=np.float),
    np.array(((-1, 0, 0), (0, -1, 0), (0, 0, 1)), dtype=np.float),
    np.array(((0, 1, 0), (0, 0, 1), (1, 0, 0)), dtype=np.float),
    np.array(((0, 0, -1), (0, -1, 0), (-1, 0, 0)), dtype=np.float),
    np.array(((0, 0, -1), (1, 0, 0), (0, -1, 0)), dtype=np.float),
    np.array(((0, 0, 1), (0, -1, 0), (1, 0, 0)), dtype=np.float)
)

class Umatrix(np.ndarray):
    ''' 3x3 orthogonal matrix corresponding to a rotation '''
    def __new__(cls, val=None):
        ''' 
        Create. If no argument given, create a completely random matrix, 
        otherwise use argument to initializegff_1 = ubi_to_gff(ubi_1,par_1)
gff_2 = ubi_to_gff(ubi_2,par_2)
        '''
        if val is None:
            self = np.zeros((3,3), dtype=np.float).view(cls)
           
        else:
            self = np.array(val, dtype=np.float).view(cls)
        return self

    def angle(self, U1,U2): #This used to be angle(self,other)
        ''' calculate angle between two U-matrices, assuming cubic symmetry '''
        #x = self.dot(other.transpose());        
        x = U1.dot(U2.transpose());

        # find symmetry-equivalent with largest trace = smallest angle
        t_max = max([np.trace(x.dot(s)) for s in CUBIC_SYMMETRIES])

        if t_max > 3.0:
            # print("t_max = ", t_max)
            # rounding errors...
            t_max = 3.0
        # trace(matrix) = 1 + 2*cos(angle)
        angle = np.arccos(0.5*(t_max-1.0))
        # convert to degrees
        return angle*180.0/np.pi

class mygrain:
    ''' grain object having its orientation, x,y,z positions, unit cell parameters and its (arbitrary) index '''
    
    if file_format_ubi == True:
        def __init__(self,U,xyz,unitcell,grainname):
            self.U=U
            self.xyz=xyz
            self.unitcell=unitcell
            self.grainno=grainname
        def __str__(self):
            return str(self.grainno)
    
    if file_format_ubi == False:
        def __init__(self,U,xyz,grainname):
            self.U=U
            self.xyz=xyz
            self.grainno=grainname
        def __str__(self):
            return str(self.grainno)

def calc_dist(grain1,grain2):
	"""
	A function calculating the distance between two reconstructed grains
	"""
	dist=np.sqrt((grain1.xyz[0][0]-grain2.xyz[0][0])**2+(grain1.xyz[0][1]-grain2.xyz[0][1])**2+(grain1.xyz[0][2]-grain2.xyz[0][2])**2)
	return dist

def getU_xyz_cellpars_from_gff(gff_file,counter):

	grain_name = 'grain_no_' + str(	gff_file.grainno[counter] ).split('.')[0]
	U_matrix = np.zeros((3,3))
	xyz_matrix = np.zeros((1,3))
	unitcell_matrix = np.zeros((1,6))

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

	unitcell_matrix[0][0] = gff_file.unitcell_a[counter]
	unitcell_matrix[0][1] = gff_file.unitcell_b[counter]
	unitcell_matrix[0][2] = gff_file.unitcell_c[counter]
	unitcell_matrix[0][3] = gff_file.unitcell_alpha[counter]
	unitcell_matrix[0][4] = gff_file.unitcell_beta[counter]
	unitcell_matrix[0][5] = gff_file.unitcell_gamma[counter]
 
	return U_matrix, xyz_matrix, unitcell_matrix, grain_name
 
 
def getU_xyz_from_gff(gff_file,counter):

    grain_name = 'grain_no_' + str(	gff_file.grainno[counter] ).split('.')[0]
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

#####
## Dictionaries holding the grains to be matched
#####

grains_ref = {}
grains_1 = {}
grains_2 = {}

if file_format_ubi == True:
    gff_ref, cellpars_ref = ubi_to_gff(ubi_ref,par_ref)
    gff_1, cellpars_1 = ubi_to_gff(ubi_1,par_1)
    gff_2, cellpars_2 = ubi_to_gff(ubi_2,par_2)

    for i in range(len(gff_ref.U11)):
        u_temp, xyz_temp, unitcell_temp, name_temp=getU_xyz_cellpars_from_gff(gff_ref,i)
        grain_temp=mygrain(u_temp, xyz_temp, unitcell_temp, name_temp)
        grains_ref[name_temp]=grain_temp

    for i in range(len(gff_1.U11)):
        u_temp, xyz_temp, unitcell_temp, name_temp=getU_xyz_cellpars_from_gff(gff_1,i)
        grain_temp=mygrain(u_temp, xyz_temp, unitcell_temp, name_temp)
        grains_1[name_temp]=grain_temp

    for i in range(len(gff_2.U11)):
        u_temp, xyz_temp, unitcell_temp, name_temp=getU_xyz_cellpars_from_gff(gff_2,i)
        grain_temp=mygrain(u_temp, xyz_temp, unitcell_temp, name_temp)
        grains_2[name_temp]=grain_temp

    
elif file_format_ubi == False:
    gff_ref = cl.columnfile(input_gff_ref)
    gff_1 = cl.columnfile(input_gff_1)
    gff_2 = cl.columnfile(input_gff_2)
    
    for i in range(len(gff_ref.U11)):
        u_temp, xyz_temp, name_temp=getU_xyz_from_gff(gff_ref,i)
        grain_temp=mygrain(u_temp, xyz_temp, name_temp)
        grains_ref[name_temp]=grain_temp
    
    for i in range(len(gff_1.U11)):
        u_temp, xyz_temp, name_temp=getU_xyz_from_gff(gff_1,i)
        grain_temp=mygrain(u_temp, xyz_temp, name_temp)
        grains_1[name_temp]=grain_temp

    for i in range(len(gff_2.U11)):
        u_temp, xyz_temp, name_temp=getU_xyz_from_gff(gff_2,i)
        grain_temp=mygrain(u_temp, xyz_temp, name_temp)
        grains_2[name_temp]=grain_temp





###############################################################################
## Grain matching
###############################################################################


matched_grains_ref_1 = {}

for keys_ref in grains_ref:
    temp_listtt = []
    print()
    print()
    print("try to match grain ", keys_ref," of reference with grain list 1")
    for keys_1 in grains_1:
        dist_diff = calc_dist( grains_ref[keys_ref], grains_1[keys_1] )
        u_temp_1 = grains_ref[keys_ref].U
        u_temp_2 = grains_1[keys_1].U
        or_diff = Umatrix()
        or_diff = or_diff.angle(u_temp_1,u_temp_2)
        ranking =  (dist_diff/0.0003 ) +  (or_diff/0.05)
        temp_listtt.append(ranking)
        if min(temp_listtt) == ranking and or_diff < 1:
            print("#####")
            print("new minimum: ", ranking)
            print("dist_diff: ",dist_diff)
            print("or_diff: ",or_diff)
            print("#####")
            matched_grains_ref_1[ grains_ref[keys_ref] ] = grains_1[keys_1]

print("#####")

matched_grains_ref_2 = {}

for keys_ref in grains_ref:
    temp_listtt = []
    print()
    print()
    print("try to match grain ", keys_ref," of reference with grain list 2")
    for keys_2 in grains_2:
        dist_diff = calc_dist( grains_ref[keys_ref], grains_2[keys_2] )
        u_temp_1 = grains_ref[keys_ref].U
        u_temp_2 = grains_2[keys_2].U
        or_diff = Umatrix()
        or_diff = or_diff.angle(u_temp_1,u_temp_2)
        ranking =  (dist_diff/0.0003 ) +  (or_diff/0.05)
        temp_listtt.append(ranking)
        if min(temp_listtt) == ranking and or_diff < 1:
            print("#####")
            print("new minimum: ", ranking)
            print("dist_diff: ",dist_diff)
            print("or_diff: ",or_diff)
            print("#####")
            matched_grains_ref_2[ grains_ref[keys_ref] ] = grains_2[keys_2]

print("#####")

matched_grains_1_2 = {}

for keys_1 in grains_1:
    temp_listtt = []
    print()
    print()
    print("try to match grain ", keys_1," of grain list 1 with grain list 2")
    for keys_2 in grains_2:
        dist_diff = calc_dist( grains_1[keys_1], grains_2[keys_2] )
        u_temp_1 = grains_1[keys_1].U
        u_temp_2 = grains_2[keys_2].U
        or_diff = Umatrix()
        or_diff = or_diff.angle(u_temp_1,u_temp_2)
        ranking =  (dist_diff/0.0003 ) +  (or_diff/0.05)
        temp_listtt.append(ranking)
        if min(temp_listtt) == ranking and or_diff < 1:
            print("#####")
            print("new minimum: ", ranking)
            print("dist_diff: ",dist_diff)
            print("or_diff: ",or_diff)
            print("#####")
            matched_grains_1_2[ grains_1[keys_1] ] = grains_2[keys_2]



###############################################################################
## Analysis on U, B matrices and grain CoM positions
###############################################################################

#####
## Number of matched grains
#####
print(" "+'\n')
print("Number of grains & number of matched grains")
print("Number of grains in the reference grain list: ",len(gff_ref.U11))
print("Number of grains in the first grain list: ",len(gff_1.U11))
print("Number of grains in the second grain list: ",len(gff_2.U11))
print("#####")
print("Matched grains: ")
print("Number of matched grains between reference and grain list 1: ", len(matched_grains_ref_1))
print("Number of matched grains between reference and grain list 2: ", len(matched_grains_ref_2))
print("Number of matched grains between grain list 1 and grain list 2: ", len(matched_grains_1_2))
print("#####")
print("Unmatched grains: ")

temp_1 = set()
temp_2 = set()
temp_3 = set()
temp_4 = set()
temp_5 = set()
temp_6 = set()
temp_7 = set()
temp_8 = set()
temp_9 = set()

if len(matched_grains_ref_1) == len(grains_1):
    if len(matched_grains_ref_2) == len(grains_2):
        if len(matched_grains_1_2) == len(grains_1):
            print("All grains are matched in both grain lists!")
            dummy_equal = True
    else:
        dummy_equal = False
    
        for keys in matched_grains_ref_1: 
            temp_1.add(keys.grainno)
            temp_2.add(matched_grains_ref_1[keys].grainno)
    
        for keys in matched_grains_ref_2: 
            temp_3.add(keys.grainno)
            temp_4.add(matched_grains_ref_2[keys].grainno)
            
        for keys in matched_grains_1_2: 
            temp_5.add(keys.grainno)
            temp_6.add(matched_grains_1_2[keys].grainno)
        
        for keys_ref in grains_ref:
            temp_7.add(grains_ref[keys_ref].grainno)        
        
        for keys_1 in grains_1:
            temp_8.add(grains_1[keys_1].grainno)
            
        for keys_2 in grains_2:
            temp_9.add(grains_2[keys_2].grainno)
     
        unmatched_grains_1 = temp_7 - temp_1    #ref vs. grain list 1 -  ref
        unmatched_grains_2 = temp_8 - temp_2    #ref vs. grain list 1 - 1
        unmatched_grains_3 = temp_7 - temp_3    #ref vs. grain list 2 -  ref
        unmatched_grains_4 = temp_9 - temp_4    #ref vs. grain list 2 - 2
        unmatched_grains_5 = temp_8 - temp_5    #grain list 1 vs. grain list 2 - 1
        unmatched_grains_6 = temp_9 - temp_6    #grain list 1 vs. grain list 2 - 2
        
        print("Reference vs. grain list 1: ")
        print("Found ",len(unmatched_grains_1)," unmatched grains in the grain list ", ubi_ref.split('/')[-1])
        print("Found ",len(unmatched_grains_2)," unmatched grains in the grain list ", ubi_1.split('/')[-1])
        print("#####")
        print("Reference vs. grain list 2: ")
        print("Found ",len(unmatched_grains_3)," unmatched grains in the grain list ", ubi_ref.split('/')[-1])
        print("Found ",len(unmatched_grains_4)," unmatched grains in the grain list ", ubi_2.split('/')[-1])
        print("#####")
        print("Grain list 1 vs. grain list 2: ")
        print("Found ",len(unmatched_grains_5)," unmatched grains in the grain list ", ubi_1.split('/')[-1])
        print("Found ",len(unmatched_grains_6)," unmatched grains in the grain list ", ubi_2.split('/')[-1])
        print("#####")

print("#####")

#####
## U
#####

U_difference_ref_1 =[]
U_difference_ref_2 =[]
U_difference_1_2 =[]

for keys in matched_grains_ref_1:
    tmp_U1 = keys.U
    tmp_U2 = matched_grains_ref_1[keys].U
    or_temp = Umatrix()
    or_temp = or_temp.angle(tmp_U1,tmp_U2)
    U_difference_ref_1.append(or_temp)

for keys in matched_grains_ref_2:
    tmp_U1 = keys.U
    tmp_U2 = matched_grains_ref_2[keys].U
    or_temp = Umatrix()
    or_temp = or_temp.angle(tmp_U1,tmp_U2)
    U_difference_ref_2.append(or_temp)
    
for keys in matched_grains_1_2:
    tmp_U1 = keys.U
    tmp_U2 = matched_grains_1_2[keys].U
    or_temp = Umatrix()
    or_temp = or_temp.angle(tmp_U1,tmp_U2)
    U_difference_1_2.append(or_temp)

U_error_ref_1_mean = np.mean(U_difference_ref_1)
U_error_ref_1_stdev = np.std(U_difference_ref_1)
U_error_ref_2_mean = np.mean(U_difference_ref_2)
U_error_ref_2_stdev = np.std(U_difference_ref_2)
U_error_1_2_mean = np.mean(U_difference_1_2)
U_error_1_2_stdev = np.std(U_difference_1_2)

U1="Average difference in U matrices per grain for reference and grain list 1 is "+str(U_error_ref_1_mean)+"(in degrees) with std. deviation of "+str(U_error_ref_1_stdev)
print(U1)
U2="Average difference in U matrices per grain for reference and grain list 2 is "+str(U_error_ref_2_mean)+"(in degrees) with std. deviation of "+str(U_error_ref_2_stdev)
print(U2)
U3="Average difference in U matrices per grain for grain list 1 and grain list 2 is "+str(U_error_1_2_mean)+"(in degrees) with std. deviation of "+str(U_error_1_2_stdev)
print(U3)
    
#####
## B
#####
## Currently, this part only works for ubi&par input pairs.
## TODO: somehow add this part for the gff, as well...
#####

print("##########################")
      
uc_a_difference_ref_1 = []
uc_b_difference_ref_1 = []
uc_c_difference_ref_1 = []
uc_alpha_difference_ref_1 = []
uc_beta_difference_ref_1 = []
uc_gamma_difference_ref_1 = []

uc_a_difference_ref_2 = []
uc_b_difference_ref_2 = []
uc_c_difference_ref_2 = []
uc_alpha_difference_ref_2 = []
uc_beta_difference_ref_2 = []
uc_gamma_difference_ref_2 = []

uc_a_difference_1_2 = []
uc_b_difference_1_2 = []
uc_c_difference_1_2 = []
uc_alpha_difference_1_2 = []
uc_beta_difference_1_2 = []
uc_gamma_difference_1_2 = []
    
if file_format_ubi == False:
    pass
   
elif file_format_ubi == True:
    if cellpars_ref == cellpars_1:
        if cellpars_ref == cellpars_2:
            cell_start = "Both grain lists have the same starting unit cell parameters with respect to reference"
            print(cell_start)
    else:
        cell_start = "Both grain lists have different starting unit cell parameters with respect to reference"
        print(cell_start)
        tmp_a_diff_1 = cellpars_1[0] - cellpars_ref[0]
        tmp_a_diff_2 = cellpars_2[0] - cellpars_ref[0]
        tmp_b_diff_1 = cellpars_1[1] - cellpars_ref[1]
        tmp_b_diff_2 = cellpars_2[1] - cellpars_ref[1]
        tmp_c_diff_1 = cellpars_1[2] - cellpars_ref[2]
        tmp_c_diff_2 = cellpars_2[2] - cellpars_ref[2]
        tmp_alpha_diff_1 = cellpars_1[3] - cellpars_ref[3]
        tmp_alpha_diff_2 = cellpars_2[3] - cellpars_ref[3]
        tmp_beta_diff_1 = cellpars_1[4] - cellpars_ref[4]
        tmp_beta_diff_2 = cellpars_2[4] - cellpars_ref[4]
        tmp_gamma_diff_1 = cellpars_1[5] - cellpars_ref[5]
        tmp_gamma_diff_2 = cellpars_2[5] - cellpars_ref[5]
        
        B1_1="Difference in a (Angstroms): "+str(tmp_a_diff_1)
        B1_2="Difference in a (Angstroms): "+str(tmp_a_diff_2)
        B2_1="Difference in b (Angstroms): "+str(tmp_b_diff_1)
        B2_2="Difference in b (Angstroms): "+str(tmp_b_diff_2)
        B3_1="Difference in c (Angstroms): "+str(tmp_c_diff_1)
        B3_2="Difference in c (Angstroms): "+str(tmp_c_diff_2)
        B4_1="Difference in alpha (degrees): "+str(tmp_alpha_diff_1)
        B4_2="Difference in alpha (degrees): "+str(tmp_alpha_diff_2)
        B5_1="Difference in beta (degrees): "+str(tmp_beta_diff_1)
        B5_2="Difference in beta (degrees): "+str(tmp_beta_diff_2)
        B6_1="Difference in gamma (degrees): "+str(tmp_gamma_diff_1)
        B6_2="Difference in gamma (degrees): "+str(tmp_gamma_diff_2)
        
        print("Grain list 1 unit cell parameters with respect to reference")
        print(B1_1)
        print(B2_1)
        print(B3_1)
        print(B4_1)
        print(B5_1)
        print(B6_1)
        print("Grain list 2 unit cell parameters with respect to reference")
        print(B1_2)
        print(B2_2)
        print(B3_2)
        print(B4_2)
        print(B5_2)
        print(B6_2)
    
    for keys in matched_grains_ref_1:
        tmp_uc_1 = keys.unitcell
        tmp_uc_2 = matched_grains_ref_1[keys].unitcell
        
        uc_a_difference_ref_1.append( (tmp_uc_1[0][0] - tmp_uc_2[0][0] ) )
        uc_b_difference_ref_1.append( tmp_uc_1[0][1] - tmp_uc_2[0][1] )
        uc_c_difference_ref_1.append( tmp_uc_1[0][2] - tmp_uc_2[0][2] )
        uc_alpha_difference_ref_1.append( tmp_uc_1[0][3] - tmp_uc_2[0][3] )
        uc_beta_difference_ref_1.append( tmp_uc_1[0][4] - tmp_uc_2[0][4] )
        uc_gamma_difference_ref_1.append( tmp_uc_1[0][5] - tmp_uc_2[0][5] )
        
    for keys in matched_grains_ref_2:
        tmp_uc_1 = keys.unitcell
        tmp_uc_2 = matched_grains_ref_2[keys].unitcell
        
        uc_a_difference_ref_2.append( (tmp_uc_1[0][0] - tmp_uc_2[0][0] ) )
        uc_b_difference_ref_2.append( tmp_uc_1[0][1] - tmp_uc_2[0][1] )
        uc_c_difference_ref_2.append( tmp_uc_1[0][2] - tmp_uc_2[0][2] )
        uc_alpha_difference_ref_2.append( tmp_uc_1[0][3] - tmp_uc_2[0][3] )
        uc_beta_difference_ref_2.append( tmp_uc_1[0][4] - tmp_uc_2[0][4] )
        uc_gamma_difference_ref_2.append( tmp_uc_1[0][5] - tmp_uc_2[0][5] )
        
    for keys in matched_grains_1_2:
        tmp_uc_1 = keys.unitcell
        tmp_uc_2 = matched_grains_1_2[keys].unitcell
        
        uc_a_difference_1_2.append( (tmp_uc_1[0][0] - tmp_uc_2[0][0] ) )
        uc_b_difference_1_2.append( tmp_uc_1[0][1] - tmp_uc_2[0][1] )
        uc_c_difference_1_2.append( tmp_uc_1[0][2] - tmp_uc_2[0][2] )
        uc_alpha_difference_1_2.append( tmp_uc_1[0][3] - tmp_uc_2[0][3] )
        uc_beta_difference_1_2.append( tmp_uc_1[0][4] - tmp_uc_2[0][4] )
        uc_gamma_difference_1_2.append( tmp_uc_1[0][5] - tmp_uc_2[0][5] )

    a_diff_ref_1_mean = np.mean(uc_a_difference_ref_1)
    a_diff_ref_1_stdev = np.std(uc_a_difference_ref_1)
    B1_ref_1="Average difference in unit cell a-parameter per grain for grain list 1 w.r.t reference is "+str(a_diff_ref_1_mean)+" (in Angstroms) with std. deviation of "+str(float("{:.6f}".format(a_diff_ref_1_stdev)))
    a_diff_ref_2_mean = np.mean(uc_a_difference_ref_2)
    a_diff_ref_2_stdev = np.std(uc_a_difference_ref_2)
    B1_ref_2="Average difference in unit cell a-parameter per grain for grain list 2 w.r.t reference is "+str(a_diff_ref_2_mean)+" (in Angstroms) with std. deviation of "+str(float("{:.6f}".format(a_diff_ref_2_stdev)))
    a_diff_1_2_mean = np.mean(uc_a_difference_1_2)
    a_diff_1_2_stdev = np.std(uc_a_difference_1_2)
    B1_1_2="Average difference in unit cell a-parameter per grain for grain list 1 w.r.t grain list 2 is "+str(a_diff_1_2_mean)+" (in Angstroms) with std. deviation of "+str(float("{:.6f}".format(a_diff_1_2_stdev)))
    print(B1_ref_1)
    print(B1_ref_2)
    print(B1_1_2)
    
    b_diff_ref_1_mean = np.mean(uc_b_difference_ref_1)
    b_diff_ref_1_stdev = np.std(uc_b_difference_ref_1)
    B2_ref_1="Average difference in unit cell b-parameter per grain for grain list 1 w.r.t reference is "+str(b_diff_ref_1_mean)+" (in Angstroms) with std. deviation of "+str(float("{:.6f}".format(b_diff_ref_1_stdev)))
    b_diff_ref_2_mean = np.mean(uc_b_difference_ref_2)
    b_diff_ref_2_stdev = np.std(uc_b_difference_ref_2)
    B2_ref_2="Average difference in unit cell b-parameter per grain for grain list 2 w.r.t reference is "+str(b_diff_ref_2_mean)+" (in Angstroms) with std. deviation of "+str(float("{:.6f}".format(b_diff_ref_2_stdev)))
    b_diff_1_2_mean = np.mean(uc_b_difference_1_2)
    b_diff_1_2_stdev = np.std(uc_b_difference_1_2)
    B2_1_2="Average difference in unit cell b-parameter per grain for grain list 1 w.r.t grain list 2 is "+str(b_diff_1_2_mean)+" (in Angstroms) with std. deviation of "+str(float("{:.6f}".format(b_diff_1_2_stdev)))
    print(B2_ref_1)
    print(B2_ref_2)
    print(B2_1_2)
    
    c_diff_ref_1_mean = np.mean(uc_c_difference_ref_1)
    c_diff_ref_1_stdev = np.std(uc_c_difference_ref_1)
    B3_ref_1="Average difference in unit cell c-parameter per grain for grain list 1 w.r.t reference is "+str(c_diff_ref_1_mean)+" (in Angstroms) with std. deviation of "+str(float("{:.6f}".format(c_diff_ref_1_stdev)))
    c_diff_ref_2_mean = np.mean(uc_c_difference_ref_2)
    c_diff_ref_2_stdev = np.std(uc_c_difference_ref_2)
    B3_ref_2="Average difference in unit cell a-parameter per grain for grain list 2 w.r.t reference is "+str(c_diff_ref_2_mean)+" (in Angstroms) with std. deviation of "+str(float("{:.6f}".format(c_diff_ref_2_stdev)))
    c_diff_1_2_mean = np.mean(uc_c_difference_1_2)
    c_diff_1_2_stdev = np.std(uc_c_difference_1_2)
    B3_1_2="Average difference in unit cell c-parameter per grain for grain list 1 w.r.t grain list 2 is "+str(c_diff_1_2_mean)+" (in Angstroms) with std. deviation of "+str(float("{:.6f}".format(c_diff_1_2_stdev)))
    print(B3_ref_1)
    print(B3_ref_2)
    print(B3_1_2)
    
    alpha_diff_ref_1_mean = np.mean(uc_alpha_difference_ref_1)
    alpha_diff_ref_1_stdev = np.std(uc_alpha_difference_ref_1)
    B4_ref_1="Average difference in unit cell alpha-parameter per grain for grain list 1 w.r.t reference is "+str(alpha_diff_ref_1_mean)+" (in degrees) with std. deviation of "+str(float("{:.6f}".format(alpha_diff_ref_1_stdev)))
    alpha_diff_ref_2_mean = np.mean(uc_alpha_difference_ref_2)
    alpha_diff_ref_2_stdev = np.std(uc_alpha_difference_ref_2)
    B4_ref_2="Average difference in unit cell alpha-parameter per grain for grain list 2 w.r.t reference is "+str(alpha_diff_ref_2_mean)+" (in degrees) with std. deviation of "+str(float("{:.6f}".format(alpha_diff_ref_2_stdev)))
    alpha_diff_1_2_mean = np.mean(uc_alpha_difference_1_2)
    alpha_diff_1_2_stdev = np.std(uc_alpha_difference_1_2)
    B4_1_2="Average difference in unit cell alpha-parameter per grain for grain list 1 w.r.t grain list 2 is "+str(alpha_diff_1_2_mean)+" (in degrees) with std. deviation of "+str(float("{:.6f}".format(alpha_diff_1_2_stdev)))
    print(B4_ref_1)
    print(B4_ref_2)
    print(B4_1_2)
    
    beta_diff_ref_1_mean = np.mean(uc_beta_difference_ref_1)
    beta_diff_ref_1_stdev = np.std(uc_beta_difference_ref_1)
    B5_ref_1="Average difference in unit cell beta-parameter per grain for grain list 1 w.r.t reference is "+str(beta_diff_ref_1_mean)+" (in degrees) with std. deviation of "+str(float("{:.6f}".format(beta_diff_ref_1_stdev)))
    beta_diff_ref_2_mean = np.mean(uc_beta_difference_ref_2)
    beta_diff_ref_2_stdev = np.std(uc_beta_difference_ref_2)
    B5_ref_2="Average difference in unit cell beta-parameter per grain for grain list 2 w.r.t reference is "+str(beta_diff_ref_2_mean)+" (in degrees) with std. deviation of "+str(float("{:.6f}".format(beta_diff_ref_2_stdev)))
    beta_diff_1_2_mean = np.mean(uc_beta_difference_1_2)
    beta_diff_1_2_stdev = np.std(uc_beta_difference_1_2)
    B5_1_2="Average difference in unit cell beta-parameter per grain for grain list 1 w.r.t grain list 2 is "+str(beta_diff_1_2_mean)+" (in degrees) with std. deviation of "+str(float("{:.6f}".format(beta_diff_1_2_stdev)))
    print(B5_ref_1)
    print(B5_ref_2)
    print(B5_1_2)
    
    gamma_diff_ref_1_mean = np.mean(uc_gamma_difference_ref_1)
    gamma_diff_ref_1_stdev = np.std(uc_gamma_difference_ref_1)
    B6_ref_1="Average difference in unit cell gamma-parameter per grain for grain list 1 w.r.t reference is "+str(gamma_diff_ref_1_mean)+" (in degrees) with std. deviation of "+str(float("{:.6f}".format(gamma_diff_ref_1_stdev)))
    gamma_diff_ref_2_mean = np.mean(uc_a_difference_ref_2)
    gamma_diff_ref_2_stdev = np.std(uc_a_difference_ref_2)
    B6_ref_2="Average difference in unit cell gamma-parameter per grain for grain list 2 w.r.t reference is "+str(gamma_diff_ref_2_mean)+" (in degrees) with std. deviation of "+str(float("{:.6f}".format(gamma_diff_ref_2_stdev)))
    gamma_diff_1_2_mean = np.mean(uc_a_difference_1_2)
    gamma_diff_1_2_stdev = np.std(uc_a_difference_1_2)
    B6_1_2="Average difference in unit cell gamma-parameter per grain for grain list 1 w.r.t grain list 2 is "+str(gamma_diff_1_2_mean)+" (in degrees) with std. deviation of "+str(float("{:.6f}".format(gamma_diff_1_2_stdev)))
    print(B6_ref_1)
    print(B6_ref_2)
    print(B6_1_2)
 
print("##########################")
   
#####
## xyz - CoM position
#####    
    
xyz_difference_ref_1 =[]
xyz_difference_ref_2 =[]
xyz_difference_1_2 =[]

for keys in matched_grains_ref_1:
    xyz_difference_ref_1.append(calc_dist( keys, matched_grains_ref_1[keys] ))
    
for keys in matched_grains_ref_2:
    xyz_difference_ref_2.append(calc_dist( keys, matched_grains_ref_2[keys] ))
    
for keys in matched_grains_1_2:
    xyz_difference_1_2.append(calc_dist( keys, matched_grains_1_2[keys] ))
    
xyz_error_ref_1_mean = np.mean(xyz_difference_ref_1)
xyz_error_ref_1_stdev = np.std(xyz_difference_ref_1)

xyz_error_ref_2_mean = np.mean(xyz_difference_ref_2)
xyz_error_ref_2_stdev = np.std(xyz_difference_ref_2)

xyz_error_1_2_mean = np.mean(xyz_difference_1_2)
xyz_error_1_2_stdev = np.std(xyz_difference_1_2)

XYZ1 = "Average difference in CoM position of grains per grain for grain list 1 w.r.t reference is "+str(xyz_error_ref_1_mean*1000)+" (in micrometers) with std. deviation of "+str(xyz_error_ref_1_stdev*1000)
print(XYZ1)
XYZ2 = "Average difference in CoM position of grains per grain grain list 2 w.r.t reference is "+str(xyz_error_ref_2_mean*1000)+" (in micrometers) with std. deviation of "+str(xyz_error_ref_2_stdev*1000)
print(XYZ2)   
XYZ3 = "Average difference in CoM position of grains per grain for grain list 1 w.r.t grain list 2 is "+str(xyz_error_1_2_mean*1000)+" (in micrometers) with std. deviation of "+str(xyz_error_1_2_stdev*1000)
print(XYZ3)     
    
print("##########################")

#####
## Analysis on unmatched grains(s)
##### 

if dummy_equal == True:
    pass

elif dummy_equal == False:
    unmatched_grains_ref_1 = {}
    unmatched_grains_ref_2 = {}
    unmatched_grains_1_2 = {}
    
    for item_1 in unmatched_grains_1:
        temp_listttt = []
        for item_2 in unmatched_grains_2:
            dist_diff = calc_dist( grains_ref[item_1], grains_1[item_2] )
            u_temp_1 = grains_ref[item_1].U
            u_temp_2 = grains_1[item_2].U
            or_diff = Umatrix()
            or_diff = or_diff.angle(u_temp_1,u_temp_2)
            ranking =  (dist_diff/0.0003 ) +  (or_diff/0.05)
            temp_listttt.append(ranking)
            if min(temp_listttt) == ranking:
                unmatched_grains_ref_1[ grains_ref[item_1] ] = grains_1[item_2]
    
    for item_1 in unmatched_grains_3:
        temp_listttt = []
        for item_2 in unmatched_grains_4:
            dist_diff = calc_dist( grains_ref[item_1], grains_2[item_2] )
            u_temp_1 = grains_ref[item_1].U
            u_temp_2 = grains_2[item_2].U
            or_diff = Umatrix()
            or_diff = or_diff.angle(u_temp_1,u_temp_2)
            ranking =  (dist_diff/0.0003 ) +  (or_diff/0.05)
            temp_listttt.append(ranking)
            if min(temp_listttt) == ranking:
                unmatched_grains_ref_2[ grains_ref[item_1] ] = grains_2[item_2]
                
    for item_1 in unmatched_grains_5:
        temp_listttt = []
        for item_2 in unmatched_grains_6:
            dist_diff = calc_dist( grains_1[item_1], grains_2[item_2] )
            u_temp_1 = grains_1[item_1].U
            u_temp_2 = grains_2[item_2].U
            or_diff = Umatrix()
            or_diff = or_diff.angle(u_temp_1,u_temp_2)
            ranking =  (dist_diff/0.0003 ) +  (or_diff/0.05)
            temp_listttt.append(ranking)
            if min(temp_listttt) == ranking:
                unmatched_grains_1_2[ grains_1[item_1] ] = grains_2[item_2]            
                
    U_difference_unmatched_ref_1 =[]
    U_difference_unmatched_ref_2 =[]
    U_difference_unmatched_1_2 =[]
    xyz_difference_unmatched_ref_1 = []
    xyz_difference_unmatched_ref_2 = []
    xyz_difference_unmatched_1_2 = []
    
    for keys in unmatched_grains_ref_1:
        tmp_U1 = keys.U
        tmp_U2 = unmatched_grains_ref_1[keys].U
        or_temp = Umatrix()
        or_temp = or_temp.angle(tmp_U1,tmp_U2)
        U_difference_unmatched_ref_1.append(or_temp)
        xyz_difference_unmatched_ref_1.append(calc_dist( keys, unmatched_grains_ref_1[keys] ))
    
    for keys in unmatched_grains_ref_2:
        tmp_U1 = keys.U
        tmp_U2 = unmatched_grains_ref_2[keys].U
        or_temp = Umatrix()
        or_temp = or_temp.angle(tmp_U1,tmp_U2)
        U_difference_unmatched_ref_2.append(or_temp)
        xyz_difference_unmatched_ref_2.append(calc_dist( keys, unmatched_grains_ref_2[keys] ))
    
    for keys in unmatched_grains_1_2:
        tmp_U1 = keys.U
        tmp_U2 = unmatched_grains_1_2[keys].U
        or_temp = Umatrix()
        or_temp = or_temp.angle(tmp_U1,tmp_U2)
        U_difference_unmatched_1_2.append(or_temp)
        xyz_difference_unmatched_1_2.append(calc_dist( keys, unmatched_grains_1_2[keys] ))    
        
    
    U_unmatched_error_ref_1_mean = np.mean(U_difference_unmatched_ref_1)
    U_unmatched_error_ref_1_stdev = np.std(U_difference_unmatched_ref_1)
    U4="Average difference in U matrices between unmatched grain(s) for grain list 1 w.r.t reference is "+str(U_unmatched_error_ref_1_mean)+" (in degrees) with std. deviation of "+str(U_unmatched_error_ref_1_stdev)
    print(U4)
    
    U_unmatched_error_ref_2_mean = np.mean(U_difference_unmatched_ref_2)
    U_unmatched_error_ref_2_stdev = np.std(U_difference_unmatched_ref_2)
    U5="Average difference in U matrices between unmatched grain(s) for grain list 2 w.r.t reference is "+str(U_unmatched_error_ref_2_mean)+" (in degrees) with std. deviation of "+str(U_unmatched_error_ref_2_stdev)
    print(U5)
    
    U_unmatched_error_1_2_mean = np.mean(U_difference_unmatched_1_2)
    U_unmatched_error_1_2_stdev = np.std(U_difference_unmatched_1_2)
    U6="Average difference in U matrices between unmatched grain(s) for grain list 1 and grain list 2 is "+str(U_unmatched_error_1_2_mean)+" (in degrees) with std. deviation of "+str(U_unmatched_error_1_2_stdev)
    print(U6)
    
    xyz_unmatched_error_ref_1_mean = np.mean(xyz_difference_unmatched_ref_1)
    xyz_unmatched_error_ref_1_stdev = np.std(xyz_difference_unmatched_ref_1)
    XYZ4 = "Average difference in CoM position of unmatched grain(s) per grain for grain list 1 w.r.t reference is "+str(xyz_unmatched_error_ref_1_mean*1000)+" (in micrometers) with std. deviation of "+str(xyz_unmatched_error_ref_1_stdev)
    print(XYZ4)
    
    xyz_unmatched_error_ref_2_mean = np.mean(xyz_difference_unmatched_ref_2)
    xyz_unmatched_error_ref_2_stdev = np.std(xyz_difference_unmatched_ref_2)
    XYZ5 = "Average difference in CoM position of unmatched grain(s) per grain for grain list 2 w.r.t reference is "+str(xyz_unmatched_error_ref_2_mean*1000)+" (in micrometers) with std. deviation of "+str(xyz_unmatched_error_ref_2_stdev)
    print(XYZ5)
    
    xyz_unmatched_error_1_2_mean = np.mean(xyz_difference_unmatched_1_2)
    xyz_unmatched_error_1_2_stdev = np.std(xyz_difference_unmatched_1_2)
    XYZ6 = "Average difference in CoM position of unmatched grain(s) per grain  for grain list 1 w.r.t grain list 2 is "+str(xyz_unmatched_error_1_2_mean*1000)+" (in micrometers) with std. deviation of "+str(xyz_unmatched_error_1_2_stdev)
    print(XYZ6) 




print("##########################")
    
    
#####
## Output
##### 

if file_format_ubi == False:
    filename= "UBI_Comparison_"+str(input_gff_1.split('/')[-1].split('.')[0])+"_and_"+str(input_gff_2.split('/')[-1].split('.')[0])+".txt"
   
elif file_format_ubi == True:
    filename= "UBI_Comparison_"+str(ubi_1.split('/')[-1].split('.')[0])+"_and_"+str(ubi_2.split('/')[-1].split('.')[0])+".txt"

f = open(filename,"w")
f.writelines("Comparison of two grain lists: "+'\n')
f.writelines(" "+'\n')
f.writelines("Input files: "+'\n')

if file_format_ubi == True:
    f0="--- "+ubi_ref.split('/')[-1]
    f.writelines(f0+'\n')
    f.writelines(" "+'\n')    
    f1="--- "+ubi_1.split('/')[-1]
    f.writelines(f1+'\n')
    f.writelines(" "+'\n')
    f2="--- "+ubi_2.split('/')[-1]
    f.writelines(f2+'\n')
    f.writelines(" "+'\n')
elif file_format_ubi == False:
    f0="--- "+input_gff_ref.split('/')[-1]
    f.writelines(f0+'\n')
    f.writelines(" "+'\n')     
    f1="--- "+input_gff_1.split('/')[-1]
    f.writelines(f1+'\n')
    f.writelines(" "+'\n')
    f2="--- "+input_gff_2.split('/')[-1]
    f.writelines(f2+'\n')
    f.writelines(" "+'\n')

f.writelines("#######################"+'\n')
f.writelines(" "+'\n')
f.writelines("Number of grains & number of matched grains"+'\n')
f.writelines(" "+'\n')

f3="Number of grains in the first flt: "+str(len(gff_1.U11))
f4="Number of grains in the second flt: "+str(len(gff_1.U11))
f.writelines(f3+'\n')
f.writelines(f4+'\n')
f.writelines(" "+'\n')

f.writelines("Unmatched grains: "+'\n')
if dummy_equal==True:
    f.writelines("All grains are matched in both grain lists!"+'\n')
    f.writelines(" "+'\n')
elif dummy_equal==False:
    f.writelines(" "+'\n')
    
    f5 = "Reference vs. grain list 1: "
    f6 = "Found "+str(len(unmatched_grains_1))+" unmatched grains in the grain list "+ ubi_ref.split('/')[-1]
    f7 = "Found "+str(len(unmatched_grains_2))+" unmatched grains in the grain list "+ ubi_1.split('/')[-1]
    f8 = "#####"
    f9 = "Reference vs. grain list 2: "
    f10 = "Found "+str(len(unmatched_grains_3))+" unmatched grains in the grain list "+ ubi_ref.split('/')[-1]
    f11 = "Found "+str(len(unmatched_grains_4))+" unmatched grains in the grain list "+ ubi_2.split('/')[-1]
    f12 = "#####"
    f13 = "Grain list 1 vs. grain list 2: "
    f14 = "Found "+str(len(unmatched_grains_5))+" unmatched grains in the grain list "+ ubi_1.split('/')[-1]
    f15 = "Found "+str(len(unmatched_grains_6))+" unmatched grains in the grain list "+ ubi_2.split('/')[-1]

    f.writelines(f5+'\n')
    f.writelines(" "+'\n')
    f.writelines(f6+'\n')
    f.writelines(" "+'\n')
    f.writelines(f7+'\n')
    f.writelines(" "+'\n')
    f.writelines(f8+'\n')
    f.writelines(" "+'\n')
    f.writelines(f9+'\n')
    f.writelines(" "+'\n')
    f.writelines(f10+'\n')
    f.writelines(" "+'\n')
    f.writelines(f11+'\n')
    f.writelines(" "+'\n')
    f.writelines(f12+'\n')
    f.writelines(" "+'\n')
    f.writelines(f13+'\n')
    f.writelines(" "+'\n')
    f.writelines(f14+'\n')
    f.writelines(" "+'\n')
    f.writelines(f15+'\n')
    f.writelines(" "+'\n')

f.writelines("#######################"+'\n')
f.writelines(" "+'\n')
f.writelines("U-matrices"+'\n')
f.writelines(" "+'\n')
f.writelines(U1+'\n')
f.writelines(" "+'\n')
f.writelines(U2+'\n')
f.writelines(" "+'\n')
f.writelines(U3+'\n')
f.writelines(" "+'\n')

if file_format_ubi == True:
    f.writelines("Unit cell parameters"+'\n')
    f.writelines(" "+'\n')
    f.writelines(B1_ref_1+'\n')
    f.writelines(" "+'\n')
    f.writelines(B1_ref_2+'\n')
    f.writelines(" "+'\n')
    f.writelines(B1_1_2+'\n')
    f.writelines(" "+'\n')
    f.writelines(B2_ref_1+'\n')
    f.writelines(" "+'\n')
    f.writelines(B2_ref_2+'\n')
    f.writelines(" "+'\n')
    f.writelines(B2_1_2+'\n')
    f.writelines(" "+'\n')
    f.writelines(B3_ref_1+'\n')
    f.writelines(" "+'\n')
    f.writelines(B3_ref_2+'\n')
    f.writelines(" "+'\n')
    f.writelines(B3_1_2+'\n')
    f.writelines(" "+'\n')
    f.writelines(B4_ref_1+'\n')
    f.writelines(" "+'\n')
    f.writelines(B4_ref_2+'\n')
    f.writelines(" "+'\n')
    f.writelines(B4_1_2+'\n')
    f.writelines(" "+'\n')
    f.writelines(B5_ref_1+'\n')
    f.writelines(" "+'\n')
    f.writelines(B5_ref_2+'\n')
    f.writelines(" "+'\n')
    f.writelines(B5_1_2+'\n')
    f.writelines(" "+'\n')
    f.writelines(B6_ref_1+'\n')
    f.writelines(" "+'\n')
    f.writelines(B6_ref_2+'\n')
    f.writelines(" "+'\n')
    f.writelines(B6_1_2+'\n')
    f.writelines("#######################"+'\n')
    
elif file_format_ubi == False:
    pass    

f.writelines("#######################"+'\n')
f.writelines(" "+'\n')
f.writelines("Grain CoM positions"+'\n')
f.writelines(" "+'\n')
f.writelines(XYZ1+'\n')
f.writelines(" "+'\n')
f.writelines(XYZ2+'\n')
f.writelines(" "+'\n')
f.writelines(XYZ3+'\n')
f.writelines(" "+'\n')
f.writelines("#######################"+'\n')
f.writelines(" "+'\n')

if dummy_equal == True:
    pass    
elif dummy_equal == False:
    f.writelines("#######################"+'\n')
    f.writelines(" "+'\n')
    f.writelines("Unmatched grain(s) error on CoM positions and orientation"+'\n')
    f.writelines(" "+'\n')
    f.writelines(U4+'\n')
    f.writelines(" "+'\n')
    f.writelines(U5+'\n')
    f.writelines(" "+'\n')
    f.writelines(U6+'\n')
    f.writelines(" "+'\n')
    f.writelines(XYZ4+'\n')
    f.writelines(" "+'\n')
    f.writelines(XYZ5+'\n')
    f.writelines(" "+'\n')
    f.writelines(XYZ6+'\n')
    f.writelines(" "+'\n')
f.writelines("#######################"+'\n')
f.writelines(" "+'\n')
f.close()



   
    
