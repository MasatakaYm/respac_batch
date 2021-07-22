# respac_batch

## Description

Simple script to perform [RESPAC](https://pubs.acs.org/doi/abs/10.1021/ct4007162) calculations (including APBS) for coarse-grained
protein charge distributions.

## Dependency

This package contains only a batch script to easily call `pdb2pqr`, `apbs`, `surface`, and `pdcp`. The last two are contained in [CafeMol](https://www.cafemol.org/), and are originally written by Dr. Tsuyoshi Terakawa.
Basically, this script is also derived from a *Perl* version created by Terakawa-san.

Recently Niina-san rewrite the `pdcp` part, which is faster than the original version.  Please ask him for details and download from his git repository.

## Run

### Using respac.py from command line

Simply add the PDB structure file to the `pdb_protein` directory (`aaa.pdb` for example), and run

```sh
$ python3 respac.py aaa
```
where `aaa` is the file name of the protein.
The output will be a file containing partial charges for surface residues in protein. 
The file will be stored in the `results` directory.

Beside, you can specify directories of input PDB, results, and templates of inputs as 
```sh
$ python3 respac.py aaa /path/to/dir_inppdb /path/to/out_results  ./lib/template
```



### Using respac.py from Python

Because the `respac.py` is designed as class, you can invoke it from python scripts. 
```python
# pdb_dir, out_dir, template_dir is optional arguments
x = Respac('2igd', pdb_dir = './pdb_protein', out_dir = '.', template_dir = './lib/template')

# You can change several conditions
x.ionic_strenght  = 0.15   # default = 0.15 [M]
x.apbs_box_margin = 20.0   # default = 20.0 [Angstrom]
x.grid_size       = 0.45   # default = 0.45 [Angstrom]
x.apbs_radius_A   = 3.0    # default = 3.0  [Angstrom]
x.apbs_radius_B   = 12.0   # default = 12.0 [Angstrom]

# Runnning all procedures
x.run_respac()

# You can invoke each step individually.
# 0. Make direcotry for output
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
```


## Other Tools

There are also some tools to make it easy to convert the output into *CafeMol*
input format, and to plot the distribution of charges.  Please see the =tools=
directory.
