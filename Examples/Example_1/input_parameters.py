## Input file of source code vasp_input_hte_cif.py
## Author: Ilias Samathrakis

def calc_type():
        ## Type of Calculations you want to perform.

        SC = True     ## Self consistent
        DOS = False    ## Density of states
        WANN = False   ## Wannier functions construction
        BAND = False   ## Band structure
	
        SMART_SEARCH = False ## Searches and determines the progress of each directory. If set 'True', none of the previous tags is valid.

        CALCS = [SC,DOS,WANN,BAND]

        return CALCS, SMART_SEARCH

def running_parameters():
        cores = 24    ## number of cores (always integer)
        time = 24     ## in hours (always integer)
        memory = 0    ## in GB (if 0, it is set automatically)
        pr_id = ""    ## Project id (you can find it using csreport command)

        return cores, time, memory, pr_id

def directories():
        vasp_dir = ''   ## Path of VASP's executable
        wannier90_dir = '' ## Path of Wannier90's executable
        pseudopotentials_dir = ''  ## Path of pseudopotentials
        python_dir = '' ## Path of python#

        vasp_output = ''  ## Output file of VASP (only name, no path) eg.vasp.log
        wannier_output = '' ## Output file of wannier90 (only name, no path) eg. wannier90.wout

        return vasp_dir, wannier90_dir, pseudopotentials_dir, python_dir, vasp_output, wannier_output

def hte_tags():
        ## Input for INCAR and KPOINTS. Modify accordingly

        magnetism        = True   ## True if the system is magnetic, False otherwise.
        non_collinearity = True   ## True if non-collinear system, False otherwise.
        SOC              = True   ## True if you want to include SOC, False otherwise.
        constraint       = False   ## True if you want to constrain along a certain direction. Set MAGMOM in 'incar_tags' to specify the direction
        U_calc           = False   ## True if you want L(S)DA+U calculations, False otherwise

        ktype = "G"               ## M for Monkhorst Pack, G for Gamma

        kspace_density     = 50   ## KPOINTS
        kspace_density_wfs = 40   ## KPOINTS.wfs
        kspace_density_dos = 60   ## KPOINTS.dos

        ## NPAR is set automatically
        ## MAGMOM is parallel to x axis by default. Set it appropriately for different directions
        ## NBANDS is calculated automatically if it is not set (preferable for high throughput)
        ## ICHARG, LWANNIER90,  are set automatically

        incar_tags = {

        "prec" : "Accurate",
        "algo" : "NORMAL",
        "encut" : 500,
        "lorbit" : 10,
        "ediff" : 0.000001,
        "ediffg" : -0.001,
        "amix" : 0.001,
        "bmix" : 0.0001,
        "isym" : -1,
        "nelm" : 500,
        "istart" : 1,
        "ismear" : 1,
        "sigma" : 0.06,
        }

        return incar_tags, magnetism, non_collinearity, SOC, constraint, U_calc, ktype, kspace_density, kspace_density_wfs, kspace_density_dos

def pseudopotentials():
        pseudep = {
        "H" : "H", "He" : "He", "Li" : "Li", "Be" : "Be", "B" : "B", "C" : "C", "N" : "N", "O" : "O_s", "F" : "F", "Ne" : "Ne",
        "Na" : "Na", "Mg" : "Mg", "Al" : "Al", "Si" : "Si", "P" : "P", "S" : "S", "Cl" : "Cl", "Ar" : "Ar", "K" : "K_pv", "Ca" : "Ca_pv",
        "Sc" : "Sc", "Ti" : "Ti", "V" : "V", "Cr" : "Cr", "Mn" : "Mn", "Fe" : "Fe", "Co" : "Co", "Ni" : "Ni", "Cu" : "Cu", "Zn" : "Zn",
        "Ga" : "Ga_d", "Ge" : "Ge_d", "As" : "As_d", "Se" : "Se", "Br" : "Br", "Kr" : "Kr", "Rb" : "Rb_pv", "Sr" : "Sr_sv", "Y" : "Y_sv", "Zr" : "Zr_sv",
        "Nb" : "Nb_pv", "Mo" : "Mo", "Tc" : "Tc", "Ru" : "Ru", "Rh" : "Rh", "Pd" : "Pd", "Ag" : "Ag", "Cd" : "Cd", "In" : "In_d", "Sn" : "Sn_d",
        "Sb" : "Sb", "Te" : "Te", "I" : "I", "Xe" : "Xe", "Cs" : "Cs_sv", "Ba" : "Ba_sv", "La" : "La", "Ce" : "Ce_3", "Pr" : "Pr_3", "Nd" : "Nd_3",
        "Pm" : "Pm_3", "Sm" : "Sm_3", "Eu" : "Eu_3", "Gd" : "Gd_3", "Tb" : "Tb_3", "Dy" : "Dy_3", "Ho" : "Ho_3", "Er" : "Er_3", "Tm" : "Tm_3", "Yb" : "Yb_3",
        "Lu" : "Lu_3", "Hf" : "Hf", "Ta" : "Ta", "W" : "W", "Re" : "Re", "Os" : "Os", "Ir" : "Ir", "Pt" : "Pt", "Au" : "Au", "Hg" : "Hg",
        "Tl" : "Tl", "Pb" : "Pb", "Bi" : "Bi", "Po" : "Po", "At" : "At", "Rn" : "Rn", "Fr" : "Fr_sv", "Ra" : "Ra_sv", "Ac" : "Ac", "Th" : "Th",
        "Pa" : "Pa", "U" : "U", "Np" : "Np", "Pu" : "Pu", "Am" : "Am", "Cm" : "Cm", "Bk" : "NA", "Cf" : "Cf", "Es" : "NA", "Fm" : "NA",
        "Md" : "NA", "No" : "NA", "Lr" : "NA", "Rf" : "NA", "Db" : "NA", "Sg" : "NA", "Bh" : "NA", "Hs" : "NA", "Mt" : "NA", "Ds" : "NA",
        "Rg" : "NA", "Cn" : "NA", "Nh" : "NA", "Fl" : "NA", "Mc" : "NA", "Lv" : "NA", "Ts" : "NA", "Og" : "NA"
        }
        return pseudep

def ldau_values():
        u_val = {
        "La" : 9, "Ce" : 9, "Pr" : 9, "Nd" : 9, "Pm" : 9, "Sm" : 9, "Eu" : 9, "Gd" : 9, "Tb" : 9, "Dy" : 9, "Ho" : 9, "Er" : 9, "Tm" : 9, "Yb" : 9, "Lu" : 9, 
        "Ac" : 9, "Th" : 9, "Pa" : 9, "U"  : 9, "Np" : 9, "Pu" : 9, "Am" : 9, "Cm" : 9, "Bk" : 9, "Cf" : 9, "Es" : 9, "Fm" : 9, "Md" : 9, "No" : 9, "Lr" : 9
        }

        j_val = {
        "La" : 0.8, "Ce" : 0.8, "Pr" : 0.8, "Nd" : 0.8, "Pm" : 0.8, "Sm" : 0.8, "Eu" : 0.8, "Gd" : 0.8, "Tb" : 0.8, "Dy" : 0.8, "Ho" : 0.8, "Er" : 0.8, "Tm" : 0.8, "Yb" : 0.8, "Lu" : 0.8,
        "Ac" : 0.8, "Th" : 0.8, "Pa" : 0.8, "U"  : 0.8, "Np" : 0.8, "Pu" : 0.8, "Am" : 0.8, "Cm" : 0.8, "Bk" : 0.8, "Cf" : 0.8, "Es" : 0.8, "Fm" : 0.8, "Md" : 0.8, "No" : 0.8, "Lr" : 0.8
        }

        l_val = {
        "La" : 2, "Ce" : 2, "Pr" : 2, "Nd" : 2, "Pm" : 2, "Sm" : 2, "Eu" : 2, "Gd" : 2, "Tb" : 2, "Dy" : 2, "Ho" : 2, "Er" : 2, "Tm" : 2, "Yb" : 2, "Lu" : 2,
        "Ac" : 2, "Th" : 2, "Pa" : 2, "U"  : 2, "Np" : 2, "Pu" : 2, "Am" : 2, "Cm" : 2, "Bk" : 2, "Cf" : 2, "Es" : 2, "Fm" : 2, "Md" : 2, "No" : 2, "Lr" : 2
        }
        return u_val,j_val,l_val
	
def main():
        calc_type()
        hte_tags()
        pseudopotentials()
        directories()
        running_parameters()
        ldau_values()
