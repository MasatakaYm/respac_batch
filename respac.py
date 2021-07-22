#!/usr/bin/env python

import os
import re
import sys


# -------------------- Commands & Paths ---------------------
pqr_exe     = "pdb2pqr30 "
apbs_exe    = "apbs "
dxmath_exe  = "dxmath "
surface_exe = "surface "
pdcp_exe    = "pdcp "
env_ldlib_path = "export LD_LIBRARY_PATH=/usr/local/lib/:/usr/lib/:$LD_LIBRARY_PATH; "


# -------------------- Defaults of Basic Variables ----------
ionic_strength  = 0.15
apbs_box_margin = 20.0
apbs_grid_size  = 0.45
apbs_radius_A   = 3.0
apbs_radius_B   = 12.0


class Respac:
    
    def __init__(self, pro_name, pdb_dir = './pdb_protein', out_dir = ".", template_dir = "./lib/template"):

        # Set required commands
        self.pqr_exe     = pqr_exe
        self.apbs_exe    = apbs_exe
        self.dxmath_exe  = dxmath_exe
        self.surface_exe = surface_exe
        self.pdcp_exe    = pdcp_exe
        self.env_ldlib_path = env_ldlib_path

        # Set filenames
        self.pdb_dir         = pdb_dir
        self.out_dir         = out_dir
        self.pdb_name        = pdb_dir + '/'              + pro_name + '.pdb'
        self.pdb_tmp_name    = out_dir + '/run/pdb/'      + pro_name + '.pdb'
        self.pqr_name        = out_dir + '/run/pqr/'      + pro_name + '.pqr'
        self.apbs_name       = out_dir + '/run/apbs_in/'  + pro_name + ".in"
        self.apbs_vol_A_name = out_dir + '/run/apbs_in/'  + pro_name + "_vol_A.in"
        self.apbs_vol_B_name = out_dir + '/run/apbs_in/'  + pro_name + "_vol_B.in"
        self.apbs_out_name   = out_dir + '/run/apbs_out/' + pro_name + '_apbs_potential.dx'
        self.volm_out_name   = out_dir + '/run/apbs_out/' + pro_name + '_delta_volm.dx'
        self.apbs_io_mc      = out_dir + '/run/apbs_out/' + 'io.mc'
        self.surf_name       = out_dir + '/run/surf_in/'  + pro_name + '.surf'
        self.pdc_name        = out_dir + '/run/pdc_in/'   + pro_name + '.pdcin'
        self.charge_name     = out_dir + '/results/'      + pro_name + '.charge'

        # Template files
        self.apbs_in_template     = template_dir + '/apbs_in_template'
        self.apbs_vol_in_template = template_dir + '/apbs_vol_in_template'
        self.pdc_in_template      = template_dir + '/pdc_in_template'
        self.dxmath_template      = template_dir + '/dxmath.inp'

        # Conditions
        self.ionic_strength  = ionic_strength
        self.apbs_box_margin = apbs_box_margin
        self.apbs_grid_size  = apbs_grid_size
        self.apbs_radius_A   = apbs_radius_A
        self.apbs_radius_B   = apbs_radius_B

        self.verbose = False


    # --------------------------------------------------------------------------------
    # Utilities

    def init(self):
        # Make output directory
        os.system('mkdir -p ' + self.out_dir + '/run/pqr')
        os.system('mkdir -p ' + self.out_dir + '/run/apbs_in')
        os.system('mkdir -p ' + self.out_dir + '/run/apbs_out')
        os.system('mkdir -p ' + self.out_dir + '/run/surf_in')
        os.system('mkdir -p ' + self.out_dir + '/run/pdc_in')
        os.system('mkdir -p ' + self.out_dir + '/run/pdb')
        os.system('mkdir -p ' + self.out_dir + '/results')
    

    def show_basic_settings(self):
        print("============================================================")
        print(" Basic settings")
        print("============================================================")
        print(" Input PDB directory = {}".format(self.pdb_dir))
        print(" Output directory = {0}/run & {0}/results".format(self.out_dir))
        
        print(" Ionic strength = {}".format(self.ionic_strength))        
        print(" Box margin     = {}".format(self.apbs_box_margin))
        print(" APBS grid size = {}".format(self.apbs_grid_size))
        print(" APBS radius A  = {}".format(self.apbs_radius_A))
        print(" APBS radius B  = {}".format(self.apbs_radius_B))


    def is_available(self, filename):
        if not os.path.exists(filename):
            print(" !!! ERROR: {} is not found. ".format(filename))
            return False
        else:
            print(" Find {} ".format(filename))
            return True

        
    def processing_pdb(self):
        print("")
        print("============================================================")
        print(" Processing PDB file: {}".format(self.pdb_name))
        print("============================================================")

        # Check file
        self.is_available(self.pdb_name)

        aa_id = -1
        ch_id = '*'
        cg_id = 0
        fout_pdb_tmp = open(self.pdb_tmp_name, 'w')
        
        with open(self.pdb_name, 'r') as pdb_in:
            for line in pdb_in:
                if line.startswith('ATOM  ') or line.startswith('HETATM'):
                    res_name = line[17:20]
                    chain_id = line[21]
                    pdb_resid = int(line[22:26])
                    
                    if res_name == 'ZN ' or res_name == ' ZN':
                        res_name = 'ZN2'
                    
                    if pdb_resid > aa_id or chain_id != ch_id:
                        aa_id  = pdb_resid
                        ch_id  = chain_id
                        cg_id += 1
                        
                    newline = 'ATOM  ' + line[6:17] + res_name + line[20:22] + str(cg_id).rjust(4) + line.rstrip()[26:] + '\n'
                    fout_pdb_tmp.write(newline)
                    
            fout_pdb_tmp.write('END')
        fout_pdb_tmp.close()
        print(" Output PDB: {}".format(self.pdb_tmp_name))

        
    def measure_boxsize(self):
        
        print("")
        print("============================================================")
        print(" Getting box size from {}".format(self.pdb_name))
        print("============================================================")

        # Check file
        self.is_available(self.pdb_name)

        aa_id = -1
        ch_id = '*'
        cg_id = 0
        
        with open(self.pdb_name, 'r') as pdb_in:
            coord_x = []
            coord_y = []
            coord_z = []
            
            for line in pdb_in:
                if line.startswith('ATOM  ') or line.startswith('HETATM'):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coord_x.append(x)
                    coord_y.append(y)
                    coord_z.append(z)
                    
            x_min, x_max = min(coord_x), max(coord_x)
            y_min, y_max = min(coord_y), max(coord_y)
            z_min, z_max = min(coord_z), max(coord_z)
            length_x = round(x_max - x_min) + 2.0
            length_y = round(y_max - y_min) + 2.0
            length_z = round(z_max - z_min) + 2.0

        print(" Getting box size: {} * {} * {}".format(length_x, length_y, length_z))

        return length_x, length_y, length_z


    def generate_apbs_inputs(self):

        length_x, length_y, length_z = self.measure_boxsize()

        print("")
        print("============================================================")
        print(" Generating APBS inputs for: {}".format(self.pqr_name))
        print("============================================================")
              
        apbs_box_xlen = length_x + self.apbs_box_margin
        apbs_box_ylen = length_y + self.apbs_box_margin
        apbs_box_zlen = length_z + self.apbs_box_margin
        apbs_grid_n_x = int(apbs_box_xlen / self.apbs_grid_size)
        apbs_grid_n_y = int(apbs_box_ylen / self.apbs_grid_size)
        apbs_grid_n_z = int(apbs_box_zlen / self.apbs_grid_size)

        fout_apbs_in = open(self.apbs_name, 'w')
        with open(self.apbs_in_template) as fin_apbs_temp:
            for line in fin_apbs_temp:
                new_line = line
                new_line = re.sub('PQRFILE',        str(self.pqr_name),       new_line)
                new_line = re.sub('DIMX',           str(apbs_grid_n_x),       new_line)
                new_line = re.sub('DIMY',           str(apbs_grid_n_y),       new_line)
                new_line = re.sub('DIMZ',           str(apbs_grid_n_z),       new_line)
                new_line = re.sub('BOXLX',          str(apbs_box_xlen),       new_line)
                new_line = re.sub('BOXLY',          str(apbs_box_ylen),       new_line)
                new_line = re.sub('BOXLZ',          str(apbs_box_zlen),       new_line)
                new_line = re.sub('IONIC_STRENGTH', str(self.ionic_strength), new_line)
                fout_apbs_in.write(new_line)
        fout_apbs_in.close()
        print(" APBS input file {} created".format(self.apbs_name))

        fout_apbs_vol_A_in = open(self.apbs_vol_A_name, 'w')
        with open(self.apbs_vol_in_template) as fin_apbs_vol_temp:
            for line in fin_apbs_vol_temp:
                new_line = line
                new_line = re.sub('PQRFILE',        str(self.pqr_name),       new_line)
                new_line = re.sub('DIMX',           str(apbs_grid_n_x),       new_line)
                new_line = re.sub('DIMY',           str(apbs_grid_n_y),       new_line)
                new_line = re.sub('DIMZ',           str(apbs_grid_n_z),       new_line)
                new_line = re.sub('BOXLX',          str(apbs_box_xlen),       new_line)
                new_line = re.sub('BOXLY',          str(apbs_box_ylen),       new_line)
                new_line = re.sub('BOXLZ',          str(apbs_box_zlen),       new_line)
                new_line = re.sub('RADIUS',         str(self.apbs_radius_A),  new_line)
                new_line = re.sub('OUTPUT',         'vol_A',                  new_line)
                new_line = re.sub('IONIC_STRENGTH', str(self.ionic_strength), new_line)
                fout_apbs_vol_A_in.write(new_line)
        fout_apbs_vol_A_in.close()
        print(" APBS input file {} created".format(self.apbs_vol_A_name))

        fout_apbs_vol_B_in = open(self.apbs_vol_B_name, 'w')
        with open(self.apbs_vol_in_template) as fin_apbs_vol_temp:
            for line in fin_apbs_vol_temp:
                new_line = line
                new_line = re.sub('PQRFILE',        str(self.pqr_name),       new_line)
                new_line = re.sub('DIMX',           str(apbs_grid_n_x),       new_line)
                new_line = re.sub('DIMY',           str(apbs_grid_n_y),       new_line)
                new_line = re.sub('DIMZ',           str(apbs_grid_n_z),       new_line)
                new_line = re.sub('BOXLX',          str(apbs_box_xlen),       new_line)
                new_line = re.sub('BOXLY',          str(apbs_box_ylen),       new_line)
                new_line = re.sub('BOXLZ',          str(apbs_box_zlen),       new_line)
                new_line = re.sub('RADIUS',         str(self.apbs_radius_B),  new_line)
                new_line = re.sub('OUTPUT',         'vol_B',                  new_line)
                new_line = re.sub('IONIC_STRENGTH', str(self.ionic_strength), new_line)
                fout_apbs_vol_B_in.write(new_line)
        fout_apbs_vol_B_in.close()
        print(" APBS input file {} created".format(self.apbs_vol_B_name))
    

    def generate_pdc_input(self):
        print("")
        print("============================================================")
        print(" Generating PDC input file: {}".format(self.pdc_name))
        print("============================================================")

        # Check files
        self.is_available(self.apbs_io_mc)

        # Reading Debye Length
        debye_length = -1.0
        
        with open(self.apbs_io_mc, 'r') as fin_iomc:
            for line in fin_iomc:
                if "Debye length =" in line:
                    try:
                        words = line.split()
                        debye_length = float(words[4])
                        assert debye_length > 0
                    except:
                        print("!!! ERROR: Invalid debye_length (={})".format(debye_length))
                    break
        print(" Detected Debye Length = ", debye_length)

        fout_pdc_in = open(self.pdc_name, 'w')
        with open(self.pdc_in_template) as fin_pdc_temp:
            for line in fin_pdc_temp:
                new_line = line
                new_line = re.sub('DEBYE', str(debye_length), new_line)
                fout_pdc_in.write(new_line)
        fout_pdc_in.close()
        

    #--------------------------------------------------------------------------------
    # Exec commands
        
    def run_pdb2pqr(self):

        self.processing_pdb()
        
        print("")
        print("============================================================")
        print(" Generating pqr file {} from {}".format(self.pqr_name, self.pdb_tmp_name))
        print("============================================================")

        # Check file
        self.is_available(self.pdb_tmp_name)

        pqr_log = " > " + self.out_dir + "/run/PDB2PQR.log 2>&1"
        pqr_command_args = "--ff=CHARMM --whitespace " + self.pdb_tmp_name + " " + self.pqr_name + pqr_log
        pqr_command = self.pqr_exe + " " + pqr_command_args
        
        try:
            if self.verbose:
                print(pqr_command)
            os.system(pqr_command)
        except:
            print(" !!! ERROR: pdb2pqr failed!")
            return
        print(" Done... ")
        

    def run_apbs(self):

        print("")
        print("============================================================")
        print(" Calculating electrostatic potential... ")
        print("============================================================")

        # Check files
        self.is_available(self.apbs_name)
        self.is_available(self.apbs_vol_A_name)
        self.is_available(self.apbs_vol_B_name)


        print(" Step 1 of 3: APBS potentials...")
        apbs_log = " > " + self.out_dir + "/run/APBS1.log 2>&1"
        apbs_command = self.env_ldlib_path + self.apbs_exe + " " + self.apbs_name + apbs_log
        try:
            if self.verbose:
                print(apbs_command)
            os.system(apbs_command)
        except:
            print(" !!! ERROR: APBS calculation 1 failed!")
            return
        print(" Done... \n")

        print(" Step 2 of 3: APBS volume A...")
        apbs_log = " > " + self.out_dir + "/run/APBS2.log 2>&1"
        apbs_command = self.env_ldlib_path + self.apbs_exe + " " + self.apbs_vol_A_name + apbs_log
        try:
            if self.verbose:
                print(apbs_command)
            os.system(apbs_command)
        except:
            print(" !!! ERROR: APBS calculation 2 failed!")
            return
        print(" Done... \n")

        print(" Step 3 of 3: APBS volume B...")
        apbs_log = " > " + self.out_dir + "/run/APBS3.log 2>&1"
        apbs_command = self.env_ldlib_path + self.apbs_exe + " " + self.apbs_vol_B_name + apbs_log
        try:
            if self.verbose:
                print(apbs_command)
            os.system(apbs_command)
        except:
            print(" !!! ERROR: APBS calculation 3 failed!")
            return
        print(" Done... \n")

        print(" DXMATH calculating...")
        dxmath_log     = " > " + self.out_dir + "/run/DXMATH.log 2>&1"
        dxmath_command = self.env_ldlib_path + self.dxmath_exe + self.dxmath_template + dxmath_log
        try:
            if self.verbose:
                print(dxmath_command)
            os.system(dxmath_command)
        except:
            print(" !!! ERROR: dxmath calculation failed!")
            return
        print(" Done... ")

        # Move output files
        mv_apbs_out_command = "mv apbs_potential.dx " + self.apbs_out_name
        mv_volm_out_command = "mv delta_vol.dx "      + self.volm_out_name
        mv_io_mc_command    = "mv io.mc "             + self.apbs_io_mc
        try:
            os.system(mv_apbs_out_command)
            os.system(mv_volm_out_command)
            os.system(mv_io_mc_command)
            os.system("rm -f vol_*.dx")
        except:
            print(" Something wrong with APBS claculations...")
            return

        
    def run_surface(self):
        print("")
        print("============================================================")
        print(" Determining surface residues... ")
        print("============================================================")

        # Check files
        self.is_available(self.pqr_name)

        surface_log = " > " + self.out_dir + "/run/SURFACE.log 2>&1"
        surface_args1 = " --pqr " + self.pqr_name + " --ofname " + self.surf_name
        surface_args2 = " --dbox 6.0 --r_probe 4.0 "
        surface_command = self.env_ldlib_path + self.surface_exe + surface_args1 + surface_args2 + surface_log
        try:
            if self.verbose:
                print(surface_command)
            os.system(surface_command)
        except:
            print(" !!! ERROR: Program surface failed!")
            print(" Done...")


    def run_pdc(self):

        print("")
        print("============================================================")
        print(" Computing Potential Derived Charges (PDC)...")
        print("============================================================")

        # Check files
        self.is_available(self.pdc_name)
        self.is_available(self.pqr_name)
        self.is_available(self.apbs_out_name)
        self.is_available(self.volm_out_name)
        self.is_available(self.surf_name)

        pdcp_log   = " > " + self.out_dir + "/run/RESPAC.log 2>&1"
        pdcp_args1 = ' --ifname '   + self.pdc_name      + ' --pqr ' + self.pqr_name
        pdcp_args2 = ' --pot '      + self.apbs_out_name + ' --vol ' + self.volm_out_name
        pdcp_args3 = ' --site All ' + ' --residue ' + self.surf_name
        pdcp_args4 = ' --ofname '   + self.charge_name
        pdcp_command = self.pdcp_exe + pdcp_args1 + pdcp_args2 + pdcp_args3 + pdcp_args4 + pdcp_log
        try:
            if self.verbose:
                print(pdcp_command)
            os.system(pdcp_command)
        except:
            print(" !!! ERROR: Program pdcp failed!")



    def run_respac(self):

        self.init()
        self.show_basic_settings()
        self.run_pdb2pqr()
        self.generate_apbs_inputs()
        self.run_apbs()
        self.run_surface()
        self.generate_pdc_input()
        self.run_pdc()

        print("")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print(" Calculation finished! Please see results in {}".format(self.charge_name))
        print(" Additional results are also provided in tools.")



if __name__ == '__main__':
    
    if len(sys.argv) == 2:
        x = Respac(sys.argv[1])
        x.run_respac()
    elif len(sys.argv) == 3:
        x = Respac(sys.argv[1], sys.argv[2])
        x.run_respac()
    elif len(sys.argv) == 4:
        x = Respac(sys.argv[1], sys.argv[2], sys.argv[3])
        x.run_respac()        
    else:
        print("Usage: python3 respac_module.py pdb_name pdb_dir out_dir")
