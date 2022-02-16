# Creates VASP input from a given POSCAR
# Author: Ilias Samathrakis

import subprocess
import os
import sys
import numpy as np
import re
import math
from fractions import Fraction
import input_parameters

def read_POSCAR():
	POSCAR_exists = False
	atoms, multiplicity, lat_con = [], [], []
	number_line, multiplier = 0, 0

	if os.path.exists('POSCAR')==True:
		POSCAR_exists = True

	if POSCAR_exists == True:
		with open("POSCAR","r") as rf:
			for i, line in enumerate(rf):
				line = line.split()
				if i == 5:
					for j in line:
						atoms.append(j)
				if i == 1:
					multiplier = float(line[0])
				if i == 2:
					line_info = []
					for j in line:
						line_info.append(multiplier*float(j))
					lat_con.append(line_info)
				if i == 3:
					line_info = []
					for j in line:
						line_info.append(multiplier*float(j))
					lat_con.append(line_info)
				if i == 4:
					line_info = []
					for j in line:
						line_info.append(multiplier*float(j))
					lat_con.append(line_info)
				if i == 6:
					if line[0].isdigit():
						number_line = 5
						for j in line:
							multiplicity.append(int(j))
			if not atoms:
				print("Information cannot be extracted from POSCAR, use cif instead")
				exit(1)
	else:
		print("POSCAR does not exist")
		exit(1)

	return atoms, multiplicity, lat_con

def estim_bands(atoms,multiplicity,in_tags,cores):
	nbands_index = search_index(in_tags,'NBANDS')

	if nbands_index != -1:
		minimum_bands = in_tags[nbands_index][1]
	else:
		number_of_atoms = 0
		for i in multiplicity:
			number_of_atoms = number_of_atoms + i
		minimum_bands = int(2 * 9 * number_of_atoms)

	num_bands = minimum_bands
	while True:
		if (num_bands/cores)%2 == 0:
			break
		num_bands = num_bands + 1

	return num_bands

def read_INCAR():
        data = []
        lorbit_index = 0
        with open("INCAR","r") as rf:
                for i, line in enumerate(rf):
                        line = line.split()
                        data.append(line)
        for i in range(len(data)):
                if data[i][0] == 'LORBIT':
                        lorbit_index = i
                        break
        return data, lorbit_index

def get_additional_tags(incar_tags,const_tag,number_bands,atoms,multiplicity,OXIDES,FELEM,u_calc,u_val,j_val,l_val,submitted,reason,magnetism,non_collinearity,soc_tag,cores):
	number_of_atoms = 0
	for i in multiplicity:
		number_of_atoms = number_of_atoms + i

	if number_of_atoms <= 50:
		submitted = True
	if number_of_atoms > 50:
		submitted = False
		reason.append("Too many atoms")

	incar_tags.append(['npar'.upper(),int(math.sqrt(cores))])
	incar_tags.append(['nbands'.upper(), number_bands])

	if magnetism:
		mag_ind = search_index(incar_tags,"MAGMOM")
		if mag_ind == -1:
			if non_collinearity:
				mag_init = '3 0 0 '
				incar_tags.append(['lnoncollinear'.upper(), ".TRUE."])
			else:
				mag_init = '3 '

			incar_tags.append(['magmom'.upper(), mag_init * number_of_atoms])
		else:
			if const_tag:
				if non_collinearity:
					incar_tags.append(['lnoncollinear'.upper(), ".TRUE."])
				rw = get_rws()
				incar_tags.append(['i_constrained_m'.upper(),1])
				incar_tags.append(['lambda'.upper(),10])
				incar_tags.append(['rwigs'.upper(),rw])
				incar_tags.append(['m_constr'.upper(),incar_tags[mag_ind][1]])	
		if soc_tag:
			incar_tags.append(['lsorbit'.upper(),".TRUE."])

	if OXIDES:
		submitted = False
		reason.append("3doxides")
	if u_calc and FELEM:
		u_values = "" 
		j_values = ""
		l_values = ""

		f_elements = ["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Ac","Th","Pa","U","Np","Pu","Am","Cm","Cf"]

		incar_tags.append(['ldau'.upper(),".TRUE."])
		incar_tags.append(['ldautype'.upper(),1])
		incar_tags.append(['ldauprint'.upper(),2])
		
		for i in atoms:
			if i in f_elements:
				u_values = u_values + str(u_val[i])
				j_values = j_values + str(j_val[i])
				l_values = l_values + str(l_val[i])
			else:
				u_values = u_values + ' 0 '
				j_values = j_values + ' 0 '
				l_values = l_values + ' -1 '

		incar_tags.append(["ldauu".upper(),u_values])
		incar_tags.append(["ldauj".upper(),j_values])
		incar_tags.append(["ldaul".upper(),l_values])

	return submitted, incar_tags

def get_INCARdos_tags(tags,lorbit_index):
	new_tags = []
	if lorbit_index != -1:
		tags[lorbit_index][1] = 10
	for i in range(len(tags)):
		new_tags.append([tags[i][0],tags[i][1]])
	new_tags.append(["ICHARG".upper(),"11"])
	return new_tags

def get_INCARwfs_tags(tags,lorbit_index):
	new_tags = []
	if lorbit_index != -1:
		tags[lorbit_index][1] = 10
	for i in range(len(tags)):
		new_tags.append([tags[i][0],tags[i][1]])
	new_tags.append(["LWANNIER90".upper(),".TRUE.".upper()])
	new_tags.append(["ICHARG".upper(),"11"])
	npar_index = search_index(new_tags,'NPAR')
	del new_tags[npar_index]
	return new_tags

def get_INCARband_tags(tags,lorbit_index):
	new_tags = []
	if lorbit_index != -1:
		tags[lorbit_index][1] = 10
	for i in range(len(tags)):
		new_tags.append([tags[i][0],tags[i][1]])
	new_tags.append(["ICHARG".upper(),"11"])
	return new_tags

def get_file(data,name):
	with open(name,'w') as wf:
		for i in range(len(data)):
			wf.write("{} = {}\n".format(data[i][0],data[i][1]))
		
	
def get_POTCAR(pseudo,atoms,submitted,reason,potcar_dir):
	full_dir = []
	elements_dir = []
	submitted = True

	for i in atoms:
		elements_dir.append(pseudo[i])
	if 'NA' in elements_dir:
		submitted = False
		reason.append("pseudopotential does not exist")
		return submitted
	for i in elements_dir:
		full_dir.append(potcar_dir + i + '/POTCAR')

	my_files, my_dirs = [], []
	my_files, my_dirs = search_files_and_dirs()

	if 'POTCAR' in my_files:
		os.system("rm POTCAR")

	with open('potgen.sh','w') as wf:
		for i in range(len(full_dir)):
			wf.write("{} {} {} {}\n".format('cat',full_dir[i],'>>','POTCAR'))
	
	os.system('sh potgen.sh')
	os.system('rm potgen.sh')
	return submitted

def get_rws():
	h = []
	with open('POTCAR','r') as rf:
		for i, line in enumerate(rf):
			line = line.split()
			if line and line[0] == 'RWIGS':
				for j in range(len(line)):
					line[j] = line[j].replace(";","")
				h.append(line[5])
	rws = " ".join(h)
	return rws

def get_KPOINTS(ktyp,lat_con,number,name):
	lat_magn, kpts = [], []
	for i in range(len(lat_con)):
		summation = 0
		for j in range(len(lat_con[0])):
			summation = summation + lat_con[i][j]*lat_con[i][j]
		summation = math.sqrt(summation)
		lat_magn.append(summation)
	for i in range(len(lat_magn)):
		kpts.append(int(round(number/lat_magn[i])))
	with open(name,"w")as wf:
		wf.write("Self-generated KPOINTS with latt*nkp product:\n")
		wf.write("0\n")
		if ktyp == "M":
			wf.write("Monkhorst-Pack\n")
		if ktyp == "G":
			wf.write("Gamma\n")
		wf.write("{} {} {}\n".format(kpts[0],kpts[1],kpts[2]))

def get_KPOINTS_dos(ktyp,lat_con,number,name):
	lat_magn, kpts = [], []
	for i in range(len(lat_con)):
		summation = 0
		for j in range(len(lat_con[0])):
			summation = summation + lat_con[i][j]*lat_con[i][j]
		summation = math.sqrt(summation)
		lat_magn.append(summation)
	for i in range(len(lat_magn)):
		kpts.append(int(round(number/lat_magn[i])))
	with open(name,"w")as wf:
		wf.write("Self-generated KPOINTS with latt*nkp product:\n")
		wf.write("0\n")
		if ktyp == "M":
			wf.write("Monkhorst-Pack\n")
		if ktyp == "G":
			wf.write("Gamma\n")
		wf.write("{} {} {}\n".format(kpts[0],kpts[1],kpts[2]))

def get_KPOINTS_wfs(ktyp,lat_con,number,name):
	lat_magn, kpts = [], []
	for i in range(len(lat_con)):
		summation = 0
		for j in range(len(lat_con[0])):
			summation = summation + lat_con[i][j]*lat_con[i][j]
		summation = math.sqrt(summation)
		lat_magn.append(summation)
	for i in range(len(lat_magn)):
		kpts.append(int(round(number/lat_magn[i])))
	with open(name,"w")as wf:
		wf.write("Self-generated KPOINTS with latt*nkp product:\n")
		wf.write("0\n")
		if ktyp == "M":
			wf.write("Monkhorst-Pack\n")
		if ktyp == "G":
			wf.write("Gamma\n")
		wf.write("{} {} {}\n".format(kpts[0],kpts[1],kpts[2]))

def get_restrictions(atoms):
	OXIDES = False
	FELEM = False

	oxides = ["O"]
	f_elements = ["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Ac","Th","Pa","U","Np","Pu","Am","Cm","Cf"]

	for i in atoms:
		if i in oxides:
			OXIDES = True

	for i in atoms:
		if i in f_elements:
			FELEM = True

	return OXIDES,FELEM

def get_run(CALC,atoms,multiplicity,names_dos,names_wfs,names_band,vasp_dir,wann_dir,cores,time,memory,increase_memory,pr_id,python_dir,vasp_output,wannier_output):
	multi_str = []
	name = ""
	GB = 1024
	for i in range(len(multiplicity)):
		multi_str.append(str(multiplicity[i]))
	for i in range(len(atoms)):
		name = name + atoms[i] + multi_str[i]
	total_atoms = 0
	for i in range(len(multiplicity)):
		total_atoms = total_atoms + multiplicity[i]

	memory = memory*GB

	if memory == 0:
		if total_atoms > 0 and total_atoms <= 10:
			memory = 4*GB
		elif total_atoms > 10 and total_atoms <= 20:
			memory = 6*GB
		elif total_atoms > 20 and total_atoms <= 30:
			memory = 8*GB
		elif total_atoms > 30 and total_atoms <= 40:
			memory = 10*GB
		elif total_atoms > 40 and total_atoms <= 50:
			memory = 12*GB
		else:
			memory = 16*GB

	if increase_memory:
		memory = memory + 2048

	time = int(time)
	if time > 168:
		time = 168
	if time == 0:
		time = 24

	with open("run.sh","w") as wf:
		wf.write("#!/bin/bash  -l\n")
		wf.write("#SBATCH -J {}\n".format(name))
		wf.write("#SBATCH -A {}\n".format(pr_id))
		wf.write("#SBATCH -n {}\n".format(cores))
		wf.write("#SBATCH --time={}:00:00\n".format(time))
		wf.write("#SBATCH --export=ALL\n")
		wf.write("#SBATCH --mem-per-cpu={}\n".format(memory))
		wf.write("#SBATCH -C avx512\n")
                wf.write("#SBATCH -p test24\n")
		wf.write("\n")
		wf.write("export SMPD_OPTION_NO_DYNAMIC_HOSTS=1\n")
		wf.write("export OMP_NUM_THREADS=1\n")
		wf.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/shared/apps/.intel/2020/mpich/3.4.1/lib/\n")
		wf.write("\n")
		wf.write("module purge\n")
		wf.write("module load intel\n")
		wf.write("module load intelmpi\n")
		wf.write("mudule load mpich\n")
		wf.write("\n")

		if CALC[0]:
			wf.write("srun -K {} >& {}\n".format(vasp_dir,vasp_output))
			wf.write("rm WAVECAR CHG CONTCAR DOSCAR EIGENVAL IBZKPT OSZICAR PCDAT slurm* XDATCAR core*\n")
			wf.write("\n")
		if CALC[1]:
			wf.write("mkdir dos\n")
			wf.write("mv {} dos/INCAR\n".format(names_dos[0]))
			wf.write("mv {} dos/KPOINTS\n".format(names_dos[1]))
			wf.write("cp CHGCAR POSCAR POTCAR dos/\n")
			wf.write("\n")
			wf.write("cd dos/\n")
			wf.write("srun -K {} >& {}\n".format(vasp_dir,vasp_output))
			wf.write("rm WAVECAR PROCAR CHG CONTCAR EIGENVAL IBZKPT OSZICAR PCDAT slurm* XDATCAR core*\n")
			wf.write("cd ../\n")
			wf.write("\n")
		if CALC[2]:
			wf.write("mkdir wann\n")
			wf.write("cp CHGCAR POSCAR POTCAR wann/\n")
			wf.write("mv {} wann/INCAR\n".format(names_wfs[0]))
			wf.write("mv {} wann/KPOINTS\n".format(names_wfs[1]))
			wf.write("mv autoconstruction.py wann/\n")
                        wf.write("\n")
			wf.write("cd wann/\n")
			wf.write("{} autoconstruction.py\n".format(python_dir))
			wf.write("srun -K {} >& {}\n".format(vasp_dir,vasp_output))
			if CALC[3]:
				wf.write("srun -K {} wannier90 >& {}\n".format(wann_dir,wannier_output))
			wf.write("rm WAVECAR CHG EIGENVAL OUTCAR CONTCAR DOSCAR IBZKPT OSZICAR PCDAT slurm* vasprun.xml XDATCAR wannier90_node* core*\n")
			wf.write("cd ../\n")
			wf.write("\n")
		if CALC[2] == False and CALC[3] == True:
			wf.write("cd wann\n")
			wf.write("srun -K {} wannier90 >& {}\n".format(wann_dir,wannier_output))
			wf.write("rm WAVECAR CHG EIGENVAL OUTCAR CONTCAR DOSCAR IBZKPT OSZICAR PCDAT slurm* vasprun.xml XDATCAR wannier90_node* core*\n")
			wf.write("cd ../\n")
			wf.write("\n")
		if CALC[4]:
			wf.write("mkdir band\n")
			wf.write("cp CHGCAR POSCAR POTCAR band/\n")
			wf.write("mv {} band/INCAR\n".format(names_band[0]))
			wf.write("mv highsymk.py band/\n")
			wf.write("mv plotband.py band/\n")
			wf.write("\n")
			wf.write("cd band\n")
			wf.write("{} highsymk.py\n".format(python_dir))
			wf.write("srun -K {} >& {}\n".format(vasp_dir,vasp_output))
			wf.write("rm WAVECAR CHG CHGCAR OUTCAR CONTCAR DOSCAR IBZKPT OSZICAR PCDAT slurm* vasprun.xml XDATCAR core*\n")
			wf.write("{} plotband.py\n".format(python_dir))
			wf.write("cd ../\n")
			wf.write("\n")
			wf.write("{} compareband.py\n".format(python_dir))
			wf.write("\n")

def search_files_and_dirs():
	my_files, my_dirs = [],[]
	for x in os.listdir("."):
		if os.path.isfile(x):
			my_files.append(x)
		elif os.path.isdir(x):
			my_dirs.append(x)
	return my_files, my_dirs

def check_vasp_output_file(my_files,filename,increase_memory):
	if filename in my_files:
		with open(filename,'r') as rf:
			for i, line in enumerate(rf):
				if increase_memory == False and "out-of-memory" in line:
					increase_memory = True
				line = line.split()
		if line and line[0] == "writing" and line[1] == "wavefunctions":
			restart = False
		else:
			restart = True
	else:
		restart = True

	return restart, increase_memory

def progress(CALC,vasp_output, wannier_output):
	my_dirs, my_files = [], []
	increase_memory = False

	my_files, my_dirs = search_files_and_dirs()
	CALC[0], increase_memory = check_vasp_output_file(my_files,vasp_output,increase_memory)
	cwd = os.getcwd()

	if "dos" in my_dirs:
		os.chdir("dos")
		my_files_dos, my_dirs_dos = search_files_and_dirs()
		CALC[1], increase_memory = check_vasp_output_file(my_files_dos,vasp_output,increase_memory)
		os.chdir(cwd)
	else:
		CALC[1] = True

	if "wann" in my_dirs:
		os.chdir("wann")
		my_files_wann, my_dirs_wann = search_files_and_dirs()
		CALC[2], increase_memory = check_vasp_output_file(my_files_wann,vasp_output,increase_memory)

		if wannier_output in my_files_wann:
			line = ''
			with open(wannier_output,'r') as rf:
				for i, line in enumerate(rf):
					line = line.split()
			if line and line[0] == "All" and line[1] == "done:" and line[2] == "wannier90" and line[3] == "exiting":
				CALC[3] = False
			elif line and line[0] == "Cycle:":
				CALC[3] = True
				number = int(line[1])
				w90_data, cycles, tag_exists = read_wann90_win('wannier90.win')
				modify_wann90_win(w90_data,tag_exists,cycles,number)
			else:
				CALC[3] = True
		else:
			CALC[3] = True
		os.chdir(cwd)
	else:
		CALC[2] = True
		CALC[3] = True

	if "band" in my_dirs:
		os.chdir("band")
		my_files_band, my_dirs_band = search_files_and_dirs()
		CALC[4], increase_memory = check_vasp_output_file(my_files_band,vasp_output,increase_memory)
		os.chdir(cwd)
	else:
		CALC[4] = True

	if CALC[0]: ### If scf is incomplete, everything has to be re-submitted
		if 'dos' in my_dirs:
			os.chdir("dos")
			os.system("cp INCAR ../INCAR.dos")
			os.system("cp KPOINTS ../KPOINTS.dos")
			os.chdir(cwd)
			os.system("rm -r dos")
		if 'wann' in my_dirs:
			os.chdir("wann")
			os.system("cp auto* ../")
			os.system("cp INCAR ../INCAR.wfs")
			os.system("cp KPOINTS ../KPOINTS.wfs")
			os.chdir(cwd)
			os.system("rm -r wann")
		if 'band' in my_dirs:
			os.chdir("band")
			os.system("cp plotband.py ../")
			os.system("cp high* ../")
			os.system("cp INCAR ../INCAR.band")
			os.chdir(cwd)
			os.system("rm -r band")
		files_to_del = ['CHG', 'CHGCAR', 'CONTCAR', 'DOSCAR', 'EIGENVAL', 'IBZKPT', 'OSZICAR', 'OUTCAR', 'PCDAT', 'slurm', 'PROCAR', 'SYM', 'vasprun.xml', 'WAVECAR', 'XDATCAR',vasp_output]
		string_del = ""
		for i in files_to_del:
			if i in my_files:
				string_del = i + " "
		if string_del != "":
			string_del = "rm " +  string_del
			os.system(string_del)
		CALC[1] = True
		CALC[2] = True
		CALC[3] = True
		CALC[4] = True
	if CALC[1] == True and CALC[0] == False: ### If dos is wrong and scf is correct, then wann has to be re-submitted
		if 'dos' in my_dirs:
			os.chdir('dos')
			os.system("cp INCAR ../INCAR.dos")
			os.chdir(cwd)
			os.system("rm -r dos")
		if 'wann' in my_dirs:
			os.chdir("wann")
			os.system("cp auto* ../")
			os.system("cp INCAR ../INCAR.wfs")
			os.system("cp KPOINTS ../KPOINTS.wfs")
			os.chdir(cwd)
			os.system("rm -r wann")
		CALC[2] = True
		CALC[3] = True
	if CALC[2] == True and CALC[1] == False and CALC[0] == False: ## If the first step of wann is wrong and dos and scf are correct, then the whole wann (1 & 2) have to be re-submitted
		if 'wann'in my_dirs:
			os.chdir("wann")
			os.system("cp auto* ../")
			os.system("cp INCAR ../INCAR.wfs")
			os.system("cp KPOINTS ../KPOINTS.wfs")
			os.chdir(cwd)
			os.system("rm -r wann")
	if CALC[4] == True and CALC[0] == False: ## If band is wrong and scf is correct, only band has to be re-submitted
		if 'band' in my_dirs:
			os.chdir('band')
			os.system("cp plotband.py ../")
			os.system("cp high* ../")
			os.system("cp INCAR ../INCAR.band")
			os.chdir(cwd)
			os.system("rm -r band")
	return CALC, increase_memory

def search_index(data,search_for):
	index = -1
	for i in range(len(data)):
		if data[i][0] == search_for:
			index = i
			break
	return index

def modify_wann90_win(data,tag_exists,cycles,number):
	with open('wannier90.win','w') as wf:
		for i in range(len(data)):
			if data[i] and data[i][0] == 'num_iter' and tag_exists == False:
				wf.write("restart = wannierise\n")
			if data[i] and data[i][0] == 'num_iter':
				data[i][2] = cycles-number+1
			for j in range(len(data[i])):
				wf.write("{} ".format(data[i][j]))
			wf.write("\n")

def read_wann90_win(file_name):
	data_w90 = []
	tag_exists = False
	with open(file_name,'r') as rf:
		for i, line in enumerate(rf):
			line = line.replace("="," = ")
			line = line.split()
			data_w90.append(line)
			if line and line[0] == 'num_iter':
				number_cycles = int(line[2])
			if line and line[0] == 'restart' and line[2] == 'wannierise':
				tag_exists = True
	for i in range(len(data_w90)):
		if data_w90[i] and data_w90[i][0] == "dis_win_min":
			data_w90[i] = ["".join(data_w90[i])]
		if data_w90[i] and data_w90[i][0] == "dis_win_max":
			data_w90[i] = ["".join(data_w90[i])]
		if data_w90[i] and data_w90[i][0] == "dis_froz_min":
			data_w90[i] = ["".join(data_w90[i])]
		if data_w90[i] and data_w90[i][0] == "dis_froz_max":
			data_w90[i] = ["".join(data_w90[i])]
	return data_w90, number_cycles, tag_exists

def modify_POSCAR():
	data = []

	with open('POSCAR',"r") as rf:
		for i, line in enumerate(rf):
			line = line.split()
			data.append(line)

	if data[5][0].isdigit():
        	with open("POSCAR","w") as wf:
                	for i in range(len(data)):
                        	for j in range(len(data[i])):
                                	wf.write("{} ".format(data[i][j]))
                        	wf.write("\n")
                        	if i == 4:
                                	for j in range(len(data[0])):
                                        	wf.write("{} ".format(data[0][j]))
                                	wf.write("\n")

def update_tags(incar_tags, magnetism, non_collinearity, soc_tag):

	if magnetism == False:
		non_collinearity = False
		incar_tags["ispin"] = 1

	if magnetism == True:
		incar_tags["saxis"] = "0 0 1"
		incar_tags["amix_mag"] = "0.0001"
		incar_tags["bmix_mag"] = "0.00001"
		incar_tags["lorbmom"] = ".TRUE."
		incar_tags["ispin"] = 2

	return incar_tags, non_collinearity, soc_tag

def main():
	###   Define files and tags

	names_scf = ['INCAR','POTCAR','KPOINTS']
	names_dos = ['INCAR.dos','KPOINTS.dos']
	names_wfs = ['INCAR.wfs','KPOINTS.wfs']
	names_band = ['INCAR.band']
	submitted = False ## Initial tag for submission (DO NOT MODIFY)
	increase_memory = False ## Initial tag for memory (DO NOT MODIFY)
	in_tags, reason = [], []

	###   Read input file

	CALC, SMART_SEARCH = input_parameters.calc_type()
	pseudo = input_parameters.pseudopotentials()
	incar_tags, magnetism, non_collinearity, soc_tag, const_tag, u_calc, ktyp, ksp, ksp_wfs, ksp_dos = input_parameters.hte_tags()
	vasp_dir, wann_dir, ps_dir, python_dir, vasp_output, wannier_output = input_parameters.directories()
	cores,time,memory,pr_id = input_parameters.running_parameters()
	u_val, j_val, l_val = input_parameters.ldau_values()
	incar_tags, non_collinearity, soc_tag = update_tags(incar_tags,magnetism, non_collinearity, soc_tag)

	for key, value in incar_tags.items():
		temp = [key.upper(),value]
		in_tags.append(temp)

	if CALC[2] == True:
		wann1 = True
		CALC.insert(3,wann1)
	elif CALC[2] == False:
		wann1 = False
		CALC.insert(3,wann1)

	if SMART_SEARCH:
		CALC,increase_memory = progress(CALC,vasp_output,wannier_output)

	modify_POSCAR()
	atoms, multiplicity, lat_con = read_POSCAR()
	number_bands = estim_bands(atoms,multiplicity,in_tags,cores)

	OXIDES,FELEM = get_restrictions(atoms)

	submitted = get_POTCAR(pseudo,atoms,submitted,reason,ps_dir)

	get_KPOINTS(ktyp,lat_con,ksp,names_scf[2])
	submitted, all_incar_tags = get_additional_tags(in_tags,const_tag,number_bands,atoms,multiplicity,OXIDES,FELEM,u_calc,u_val,j_val,l_val,submitted,reason,magnetism,non_collinearity,soc_tag,cores)

	lorbit_index = search_index(all_incar_tags,'lorbit'.upper())

	incar_dos_tags = get_INCARdos_tags(all_incar_tags,lorbit_index)
	incar_wfs_tags = get_INCARwfs_tags(all_incar_tags,lorbit_index)
	incar_band_tags = get_INCARband_tags(all_incar_tags,lorbit_index)

	if CALC[0]:
		get_file(all_incar_tags,names_scf[0])
	if CALC[1]:
		get_file(incar_dos_tags,names_dos[0])
		get_KPOINTS_dos(ktyp,lat_con,ksp_dos,names_dos[1])
	if CALC[2]:
		get_file(incar_wfs_tags,names_wfs[0])
		get_KPOINTS_wfs(ktyp,lat_con,ksp_wfs,names_wfs[1])
	if CALC[4]:
		get_file(incar_band_tags,names_band[0])
	
	get_run(CALC,atoms,multiplicity,names_dos,names_wfs,names_band,vasp_dir,wann_dir,cores,time,memory,increase_memory,pr_id,python_dir,vasp_output,wannier_output)

	with open('info.dat','w') as wf:
		wf.write("Submitted: {} ".format(submitted))
		for i in range(len(reason)):
			wf.write("{} ".format(reason[i]))
	os.system("rm -r __py*")

main()
