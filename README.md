DODO: re<ins>D</ins>esign AlphaF<ins>O</ins>ld2 <ins>D</ins>isorered regi<ins>O</ins>ns  
==============================

## What is DODO?  

DODO is a Python package and command-line utility for taking an AF2 structure and redesigning the disordered regions to make them look more like IDRs. To be clear, the work done by DeepMind to make AlphaFold2 is **AMAZING**, and I do not mean to take away from that in *ANY WAY*. However, for visualizing proteins for presentations, etc. it would be nifty to be able to make the IDRs look more 'IDR-like'. DODO does just that! What it does is identify the IDRs in the structure, predict the end-to-end distance for each IDR from its sequence using ALBATROSS (see https://github.com/idptools/sparrow) and rebuild the structure such that the IDRs are the approximate correct overall dimensions (see example below). If that visualization doesn't work for you, there are other options to make the IDRs more compact or expanded than is predicted from sequence. In addition, you can make a PDB with multiple IDRs in a single 'structure' and keep the folded domains fixed, which when opened in VMD makes something that looks *like* a simulation trajectroy (to be very clear, it is **NOT** the equivalent to an actual simulation trajectory but is really nice for visualizations).  

![DODO_EXAMPLE](https://github.com/ryanemenecker/dodo/blob/main/images/DODO_example.png)

### Current Limitations  

1. Rebuilding employs a simple random walk approach, so the actual IDR conformation is completely random *and not scientifically useful*. However, it is still quite useful for visualization of the protein overall.  

2. At this moment, the rebuilt IDRs only contain alpha carbons whereas folded regions contain all atoms. This is mainly to do with the trickiness of placing all the atoms back after dramatically changing the IDR. I'm working on all-atom IDR generation, which I'm hoping to bring in the future (but this is not a gaurantee because it's tricky business and for visualization is not really necessary to be honest).  

3. Because some IDRs only have alpha carbons, the order of bonds result in an 'unusual bond' warning in VMD.

4. Some visualization modes don't work in VMD. Licorice seems to work fine, tube and trace do not. I'm working on this.  

### How can you use DODO?

DODO is currently usable in Python and from the command-line.  

### Installation

**Note** - to install DODO, you first need to have cython and numpy installed. To install cython and numpy, simpy run:  

```console
pip install cython numpy
  ```

Once you have cython and numpy installed, you should be able to install DODO. 

To install DODO, run the following command from terminal:  
```console
pip install git+https://github.com/ryanemenecker/dodo.git
```

## DODO Python Functions

First import build from DODO.  

```python
from dodo import build
```

You can build new structures from from an existing PDB or just have DODO download the structure from the AF2 database. You can also generate PDBs of IDRs from sequence alone!  
  
### Generating a structure from the name alone

To have DODO download a structure from a protein name and alter the disordered regions, you can use the ``pdb_from_name()`` function. There are two required arguments unless you set ``graph=True``: 1. the protein name as a string, 2. the out_path for where to save the PDB. If you set ``graph=True``, you don't need to specify the outpath and DODO will just show your PDB in a 3D graph using matplotlib. 

```python
build.pdb_from_name('human p53', out_path='/Users/your_user_name/Desktop/my_cool_proteins/my_protein.pdb')
```
**Additional usage:**  

All arguments for ``build.pdb_from_name()`` are as follows:  
**protein_name** - required. The name of your protein as a string. Specifying the organism increases your chance of success.  
  
**out_path** - optional if you set graph=True. Otherwise raises an exception. Where to save your protein structure file. Specify the file name here.  

**mode** - optional. Default: 'predicted'. The ``predicted`` option predicts the end-to-end distance of your disordered regions from sequence and then makes the IDRs fit within that distance. Additional options are ``super_compact``, ``compact``, ``normal``, ``expanded``, ``super_expanded``, ``max_expansion``. These are pretty self explanatory.  

**num_models** - optional. Default: 1. ``num_models`` lets you choose the number of models of IDRs to make for your protein. The folded domains are left in the same location for all models wherease the IDRs vary.  

**linear_placement** - optional. Default: False. Whether to place the folded domains linearly for visualization.  

**beta_for_FD_IDR** - optional. Default: False. Whether to set beta values such that all IDRs = 0 and FDs=100 for visualization.  
  
**include_FD_atoms** - optional. Default: True. Whether to include all atoms for the FDs. Only CA for IDRs for now.  

**CONECT_lines** - optional. Default: True. Whether to included CONECT lines in the generated PDB. Makes visualization generally better.  

**verbose** - optional. Default: True. Whether to show progress as structure is being made.  

**use_metapredict** - optional. Default: False. This option lets you use metapredict to predicte the IDRs and folded regions. Although *fairly accurate*, it doesn't get the exact cutoffs for some regions and fails to predict small loops within large folded regions. The default is to use the number of atoms neighboring each atom in the AF2 structure. The default behavior is slower but works better.  

**graph** - optional. Default: False. Setting this to True will pull up a really rough looking structure of your protein using the 3D graphing functionality in matplotlib. This is something I made when developing this to quickly look at structures. You shouldn't use this, but you can if you want. It's kind of fun TBH.  
   
**attempts_per_region** - optional. Default: 20. Number of times to try and make each region.  

**attempts_per_coord** - optional. Default: 2000. Number of times to try to generate each coordinate for each alpha carbon in the structure.  
  

### Modifying the IDR from an existing PDB file

You can also have DODO alter a pre-existing AF2 pdb file. The AF2 file should have all atom information (though this isn't required). There are two required arguments if not graphing: ``path_to_pdb`` : the path to your pdb file as a string and ``out_path`` : the path and filename of where to save your file. If you set ``graph=True``, you don't need to specify the out_path.

```python
build.pdb_from_pdb('/Users/your_user_name/Desktop/my_AF2_pdb.pdb', out_path='/Users/your_user_name/Desktop/my_AF2_PDB_DODO.pdb')
```
**Additional usage:**  

All arguments for ``build.pdb_from_pdb()`` are as follows:  
**path_to_pdb** - required. The filepath to your pdb as a string.  
  
**out_path** - optional if you set graph=True. Otherwise raises an exception. Where to save your protein structure file. Specify the file name here.  

**mode** - optional. Default: 'predicted'. The ``predicted`` option predicts the end-to-end distance of your disordered regions from sequence and then makes the IDRs fit within that distance. Additional options are ``super_compact``, ``compact``, ``normal``, ``expanded``, ``super_expanded``, ``max_expansion``. These are pretty self explanatory.  

**num_models** - optional. Default: 1. ``num_models`` lets you choose the number of models of IDRs to make for your protein. The folded domains are left in the same location for all models wherease the IDRs vary.  

**linear_placement** - optional. Default: False. Whether to place the folded domains linearly for visualization.  

**beta_for_FD_IDR** - optional. Default: False. Whether to set beta values such that all IDRs = 0 and FDs=100 for visualization.  
  
**include_FD_atoms** - optional. Default: True. Whether to include all atoms for the FDs. Only CA for IDRs for now.  

**CONECT_lines** - optional. Default: True. Whether to included CONECT lines in the generated PDB. Makes visualization generally better.  

**verbose** - optional. Default: True. Whether to show progress as structure is being made.  

**use_metapredict** - optional. Default: False. This option lets you use metapredict to predicte the IDRs and folded regions. Although *fairly accurate*, it doesn't get the exact cutoffs for some regions and fails to predict small loops within large folded regions. The default is to use the number of atoms neighboring each atom in the AF2 structure. The default behavior is slower but works better.  

**graph** - optional. Default: False. Setting this to True will pull up a really rough looking structure of your protein using the 3D graphing functionality in matplotlib. This is something I made when developing this to quickly look at structures. You shouldn't use this, but you can if you want. It's kind of fun TBH.  
   
**attempts_per_region** - optional. Default: 20. Number of times to try and make each region.  
  
**attempts_per_coord** - optional. Default: 2000. Number of times to try to generate each coordinate for each alpha carbon in the structure.  


### Making a PDB for an IDR just from sequence.

If you just have an IDR sequence, you can also generate coordinates for that. Just like before, you can save it or you can graph it in matplotlib.

```python
build.pdb_from_sequence('VQQQGIYNNGTIAVANQVSCQSPNQ', out_path='/Users/your_user_name/Desktop/my_AF2_PDB_DODO.pdb')
```
**Additional usage:**  

All arguments for ``build.pdb_from_sequence()`` are as follows:  
  
**sequence** - The sequence of the IDR as a string.
  
**out_path** - optional if you set ``graph=True``. Otherwise raises an exception. Where to save your protein structure file. Specify the file name here.  
  
**mode** - optional. Default: 'predicted'. The ``predicted`` option predicts the end-to-end distance of your disordered regions from sequence and then makes the IDRs fit within that distance. For IDRs between folded domains, the folded domains are moved apart from each other to accomodate the IDR length. For terminal IDRs, they will just adopt a configuration within the distance equal to the end to end distance. Additional options are ``super_compact``, ``compact``, ``normal``, ``expanded``, ``super_expanded``, ``max_expansion``. These are pretty self explanatory.  
 
**CONECT_LINES** - optional. Default: True. Writes CONECT lines to the PDB so that all alpha carbons are connected.  

**graph** - optional. Default: False. Setting this to True will pull up a really rough looking structure of your protein using the 3D graphing functionality in matplotlib. This is something I made when developing this to quickly look at structures. You shouldn't use this, but you can if you want. It's kind of fun TBH.  
  
**verbose** optional. Default: False. Set to True to get more info on what is happening as your structure is being made.  
  
**attempts_per_idr** - optional. Default: 50. Number of attempts to IDR. Basically if number of attempts per coordinate reaches max and it still hasn't succeeded, this is the number of times it tries again.  
  
**attempts_per_res** - optional. Number of attempts to make each coordinate fit without clashing.  
  
**end_coord** - optional. Default: (0,0,0). The end coordinate for the IDR as a tuple.

  

## DODO command-line functions

Just like in Python, you can generate AF2 structures with re-designed IDRs from the command-line using an existing AF2 PDB file or the name of a protein, and you can also generate PDBs for coordinates of IDRs from sequence alone. The only thing missing from the command-line is the graphing functionality.  

### Generating a structure from a protein name from the command-line

To have DODO download a structure using a protein name and alter the disordered regions, you can use the ``pdb-from-name`` command. There are two required arguments: 1. the protein name, 2. the out_path for where to save the PDB.

```bash
pdb-from-name human p53 -o /Users/your_user_name/Desktop/my_cool_proteins/my_protein.pdb
```
**Additional usage:**  

All arguments for ``pdb-from-name`` are as follows:  
``-o`` or ``--out_path`` : **required**. Where to save your protein structure file. Specify the file name here.  
  
``-m`` or ``--mode`` : optional. Default: **predicted**. ``--mode`` lets you specify how to build the IDR. The default **predicted** option predicts the end-to-end distance of your disordered regions from sequence and then makes the IDRs fit within that distance. Additional options are ``super_compact``, ``compact``, ``normal``, ``expanded``, ``super_expanded``, ``max_expansion``. These are pretty self explanatory.  
  
``-n`` or ``--num_models`` : optional. Default: 1. ``--num_models`` lets you choose the number of models of IDRs to make for your protein. The folded domains are left in the same location for all models wherease the IDRs vary.   

``-l`` or ``--linear_placement`` : optional. Default: False. The ``--linear_placement`` flag lets you place the IDRs linearly across a vector for visualization.  
  
``-b`` or ``--beta_for_FD_IDR`` : optional. Default: False. Whether to set beta values such that all IDRs = 0 and FDs=100 for visualization.  
  
``-c`` or ``--no_CONECT_lines`` : optional. Default: False. The ``--no_CONECT_lines`` flag lets you make files without CONECT lines.  
  
``-f`` or ``--no_FD_atoms`` : optional. Default: False. The ``--no_FD_atoms`` flag lets you make structures with ONLY alpha carbon atoms. By default, the IDRs are only alpha carbons and the FDs are all atom.  
  
``-u`` or ``--use_metapredict`` : optional. Default: False. The ``--use_metapredict`` flag lets use metapredict V2-ff to predict the disordered regions in your structure. By default, the location of all atoms in the structure are used to infer the location of the IDRs.  
  
``-s`` or ``--silent`` : optional. Default: False. The ``--silent`` flag lets you silent most print output to your terminal.  
  
``-apr`` or ``--attempts_per_region`` : optional. Default: 40. ``--attempts_per_region`` lets you specify the number of attempts to make each region of the structure.  
  
``-apc`` or ``--attempts_per_coord`` : optional. Default: 2000. ``--attempts_per_coord`` lets you specify the number of attempts to make to generate each coordinate in your structure.  
  
  
### Modifying the IDR from an existing AF2 PDB file from the command-line
  
To have DODO redesign the disordered regions using an AF2 PDB you have saved locally, you can use the ``pdb-from-pdb`` command. There are two required arguments: 1. The path to your PDB file including the filename and extension and, 2. the out_path for where to save the final PDB.

```bash
pdb-from-pdb /Users/your_user_name/Desktop/my_AF2_proteins/AF2_protein.pdb -o /Users/your_user_name/Desktop/my_cool_proteins/my_protein_with_a_new_IDR.pdb
```
**Additional usage:**  

All arguments for ``pdb-from-pdb`` are as follows:  
``-o`` or ``--out_path`` : **required**. Where to save your protein structure file. Specify the file name here.  
  
``-m`` or ``--mode`` : optional. Default: **predicted**. ``--mode`` lets you specify how to build the IDR. The default **predicted** option predicts the end-to-end distance of your disordered regions from sequence and then makes the IDRs fit within that distance. Additional options are ``super_compact``, ``compact``, ``normal``, ``expanded``, ``super_expanded``, ``max_expansion``. These are pretty self explanatory.  
  
``-n`` or ``--num_models`` : optional. Default: 1. ``--num_models`` lets you choose the number of models of IDRs to make for your protein. The folded domains are left in the same location for all models wherease the IDRs vary.   

``-l`` or ``--linear_placement`` : optional. Default: False. The ``--linear_placement`` flag lets you place the IDRs linearly across a vector for visualization.  
  
``-b`` or ``--beta_for_FD_IDR`` : optional. Default: False. Whether to set beta values such that all IDRs = 0 and FDs=100 for visualization.  
  
``-c`` or ``--no_CONECT_lines`` : optional. Default: False. The ``--no_CONECT_lines`` flag lets you make files without CONECT lines.  
  
``-f`` or ``--no_FD_atoms`` : optional. Default: False. The ``--no_FD_atoms`` flag lets you make structures with ONLY alpha carbon atoms. By default, the IDRs are only alpha carbons and the FDs are all atom.  
  
``-u`` or ``--use_metapredict`` : optional. Default: False. The ``--use_metapredict`` flag lets use metapredict V2-ff to predict the disordered regions in your structure. By default, the location of all atoms in the structure are used to infer the location of the IDRs.  
  
``-s`` or ``--silent`` : optional. Default: False. The ``--silent`` flag lets you silent most print output to your terminal.  
  
``-apr`` or ``--attempts_per_region`` : optional. Default: 40. ``--attempts_per_region`` lets you specify the number of attempts to make each region of the structure.  
  
``-apc`` or ``--attempts_per_coord`` : optional. Default: 2000. ``--attempts_per_coord`` lets you specify the number of attempts to make to generate each coordinate in your structure.
   

### Making a PDB for an IDR just from sequence from the command-line.

If you just have an IDR sequence, you can also generate coordinates for that.

```bash
pdb-from-sequence GRNQNGGGYQNYNNQGYQGHGGQHQNNYNQYPCNYFGPGYNN -o /Users/your_user_name/Desktop/my_cool_proteins/my_cool_IDR.pdb
```
**Additional usage:**  

All arguments for ``pdb-from-pdb`` are as follows:  
``-o`` or ``--out_path`` : **required**. Where to save your protein structure file. Specify the file name here.  
  
``-m`` or ``--mode`` : optional. Default: **predicted**. ``--mode`` lets you specify how to build the IDR. The default **predicted** option predicts the end-to-end distance of your disordered regions from sequence and then makes the IDRs fit within that distance. Additional options are ``super_compact``, ``compact``, ``normal``, ``expanded``, ``super_expanded``, ``max_expansion``. These are pretty self explanatory.  
  
``-n`` or ``--num_models`` : optional. Default: 1. ``--num_models`` lets you choose the number of models of IDRs to make for your protein. The folded domains are left in the same location for all models wherease the IDRs vary.   
  
``-c`` or ``--no_CONECT_lines`` : optional. Default: False. The ``--no_CONECT_lines`` flag lets you make files without CONECT lines.  
  
``-api`` or ``--attempts_per_IDR`` : optional. Default: 50. ``--attempts_per_IDR`` lets you specify the number of attempts to make to the IDR.  

``-apr`` or ``--attempts_per_residue`` : optional. Default: 1000. ``--attempts_per_residue`` lets you specify the number of attempts to make the coordinates for each residue for your IDR.   
   
  
#### Changes

Logging changes below.  

V0.11 - October 30, 2023. Couple small fixes, updated documentation, fixed some embarassingly bad typos.

V0.10 - October 24, 2023. Big changes! Added functionality to generate multiple IDRs for a single PDB so you can make 'simulation-like' visualizations when viewing in PDB! I also added command-line functionality and fixed some more bugs. Note: the multiple models are for visualization only and not equivalent to actual simulations!

V0.06 - October 17, 2023. Added functionality to place folded domains in an approximate linear arrangement for visualization purposes.  

V0.05 - October 5, 2023. Major overhaul to the backend to make structure generation more robust and efficient. Changed some user facing functionality and added ability to generate a PDB of an IDR from sequence alone. Improved documentation.

V0.04 - September 29, 2023. Made it so the atoms of folded domains can be kept in the structure. Only CA for the IDRs for now. Improved performance a bit. Still need to clean up the code a lot.. Made everything callable by the PDBParserObj made from the class pdb_tools.PDBParser 

V0.03 - September 28, 2023. Made it so you can save the IDRs generated from sequence alone. Removed the ability to make all atom structures because it wasn't great. This hopefully will come in the future. 

V0.02 - September 26, 2023. Added generating IDR coords from sequence alone. Added filling in all atom coordinates from alpha carbon coordinate using fixed bond angles and distances, which isn't great but is better than nothing. 

V0.01 - September 25, 2023. Initial release.  


### Copyright

Copyright (c) 2023, Ryan Emenecker - Holehouse Lab


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
