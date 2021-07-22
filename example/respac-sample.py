import sys
import requests

sys.path.append("../")
sys.dont_write_bytecode = True

from respac import Respac


def fetch_pdb(pdb_name, out_dir):
    target_url = 'https://files.rcsb.org/download/' + pdb_name + '.pdb'
    print(target_url)
    req_pdb = requests.get(target_url)

    pdb = open(out_dir + '/' + pdb_name + '.pdb', 'w')
    pdb.write(req_pdb.text)
    pdb.close()
    


def main():

    fetch_pdb('2igd', '../pdb_protein')

    # pdb_dir, out_dir, template_dir is optional arguments
    x = Respac('2igd', pdb_dir = '../pdb_protein', out_dir = '.', template_dir = '../lib/template')

    # respac.py assumes $PATH to programs has already set in default.
    # If not, you can specify path as follows:
    
    # x.pqr_exe        = "path/to/pdb2pqr30 "                                         # default "pdb2pqr30 "
    # x.apbs_exe       = "path/to/APBS-3.0.0.Linux/bin/apbs "                         # default "apbs "
    # x.dxmath_exe     = "path/to/APBS-3.0.0.Linux/share/apbs/tools/dxmath "          # default "dxmath "
    # x.surface_exe    = "path/to/cafemol_3.2.1/utility/RESPAC/surface/bin/surface "  # default "surface "
    # x.pdcp_exe       = "path/to/cafemol_3.2.1/utility/RESPAC/pdcp/bin/pdcp "        # default "pdcp "
    # x.env_ldlib_path = "export LD_LIBRARY_PATH=/path/to/APBS-3.0.0.Linux/lib/:/usr/local/lib/:/usr/lib/:$LD_LIBRARY_PATH; "

    
    # You can change several conditions
    x.ionic_strenght  = 0.15   # default = 0.15 [M]
    x.apbs_box_margin = 20.0   # default = 20.0 [Angstrom]
    x.grid_size       = 0.45   # default = 0.45 [Angstrom]
    x.apbs_radius_A   = 3.0    # default = 3.0  [Angstrom]
    x.apbs_radius_B   = 12.0   # default = 12.0 [Angstrom]


    # Runnning all procedure
    x.run_respac()


    # You can run individual step
    # 0. Make direcotry
    x.init()
    
    # 1. PDB2PQR step
    x.run_pdb2pqr()

    # 2. APBS step
    x.generate_apbs_inputs()
    x.run_apbs()

    # 3. Surface step
    x.run_surface()

    # 4. PDC step
    x.generate_pdc_input()
    x.run_pdc()
    


if __name__ == '__main__':
    main()
