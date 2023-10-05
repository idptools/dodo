DODO
==============================

**D**e-sillifying AlphaF**O**ld2 **D**isorered regi**O**ns  

### What is DODO?  

DODO is a A Python package for taking an AF2 structure and making the disordered regions look less silly (I believe the scientific term is *de-sillifying*). To be clear, the work done by DeepMind to make AlphaFold is AMAZING, and I do not mean to take away from that in ANY WAY. However, for visualizing the proteins for presentations, etc. it would be nifty to be able to make the IDRs look more 'IDR-like'. DODO does just that! What it does is identify the IDRs in the structure, predict their end-to-end distance from the sequence of the IDR (by default, though there are other options..), and rebuild the structure such that the IDRs are the approximate correct overall dimensions. 

#### Current Limitations  

1. Rebuilding employs a simple random walk approach, so the actual IDR conformation is completely random *and not scientifically useful*. However, it is still quite useful for visualization of the protein overall.  

2. At this moment, the rebuilt IDRs only contain alpha carbons whereas folded regions maintain all atoms. This is mainly to do with the trickiness of placing all the atoms back after dramatically changing the IDR. I'm working on all-atom IDR generation, which I'm hoping to bring in the future (but this is not a gaurantee because it's tricky business and for visualization is quite unnecessary to be honest).  

3. Some bonds aren't *quite* the right distance apart, so you might see an 'unusual bond between residues' warning in VMD. I'm working on this. It isn't so bad that it ruins the visualization, but I'll still try to fix it. This is largely because the IDRs only have alpha carbons, which results in alpha carbon to N or C bonds, which wouldn't occur if the IDR had all the atoms. I'm pretty close to getting this resolved.

### How can you use DODO?

DODO is currently only usable in Python.  

### Installation

To install DODO, run the following command from terminal:  

    $ pip install git+https://github.com/ryanemenecker/dodo.git

### DODO Python Functions

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

**attempts** - optional. Default: 25000. Number of attempts to make each coordinate fit without clashing.
  
**end_coord** - optional. Default: (0,0,0). The end coordinate for the IDR as a tuple.

  
  
#### Changes

Logging changes below.  

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
