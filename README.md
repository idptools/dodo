DODO
==============================

**D**e-sillifying AlphaF**O**ld2 disor**D**ered regi**O**ns  


DODO is a A Python package for taking an AF2 structure and making the disordered regions look less silly (I believe the scientific term is *de-sillifying*).  

### How can you use DODO?

DODO is currently under development and only usable in Python.  

### Installation

To install DODO, run the following command from terminal:  

    $ pip install git+https://github.com/ryanemenecker/dodo.git

### DODO Python Functions

DODO is under development, so there's some work to do finishing up functions. Nonetheless, I figured I'd make basic functionality available.  

First import build from DODO.  

```python
from dodo import build
```

You can build new structures from from an existing PDB or just have DODO download the structure from the AF2 database.

**Using a protein name.**  
To have DODO download a structure from a protein name and alter the disordered regions, you can use the ``pdb_from_name()`` function. There is one required argument: the protein name as a string. If you don't specify the ``out_path``, DODO will return a dictionary that really won't make a lot of sense to you because I haven't finished making it useful. So for the time, specify the out_path.  

```python
build.pdb_from_name('human p53', out_path='/Users/your_user_name/Desktop/my_cool_proteins/my_protein.pdb')
```
**Additional usage:**  

All arguments for ``build.pdb_from_name`` are as follows:  
**protein_name** - required. The name of your protein as a string. Specifying the organism increases your chance of success.  
**out_path** - optional. Where to save your protein structure file. Specify the file name here.  
**mode** - optional. Default: 'predicted'. The ``predicted`` option predicts the end-to-end distance of your disordered regions from sequence and then makes the IDRs fit within that distance. For IDRs between folded domains, the folded domains are moved apart from each other to accomodate the IDR length. For terminal IDRs, they will just adopt a configuration within the distance equal to the end to end distance. Additional options are ``super_compact``, ``compact``, ``normal``, ``expanded``, ``super_expanded``, ``max_expansion``. These are pretty self explanatory.  
**use_metapredict** - optional. Default: False. This option lets you use metapredict to predicte the IDRs and folded regions. Although *fairly accurate*, it doesn't get the exact cutoffs for some regions and fails to predict small loops within large folded regions. The default is to use the number of atoms neighboring each atom in the AF2 structure to figure out what is a folded domain, what is a separate folded domain that is proximal to another folded domain, what is an IDR, and what is a disordered loop. The default behavior is slower but works better.  
**graph** - optional. Default: False. Setting this to True will pull up a really rough looking structure of your protein using the 3D graphing functionality in matplotlib. This is something I made when developing this to quickly look at structures. You shouldn't use this, but you can if you want. It's kind of fun TBH.  
**beta_by_region** optional. Default: True. This takes the beta values and changes them such that if you view the amino acids in something like VMD by beta values, you can see what DODO predicted to be an IDR and what it left as a folded domain. Loops will also be a different color if they were present.  
**silent** optional. Default: False. This silences the message telling you what structure is being downloaded. Turn this off at your own risk. DODO does not gaurantee that it will download the structure you wanted from name alone. If you just put 'p53' instead of 'human p53', you might get a different organism based on how Uniprot ranks it.  
**verbose** optional. Default: False. Set to True to get more info on what is happening as your structure is being made.  
**attempts_per_region** - optional. Default: 50. Number of times to try and make each region.  
**attempts_per_coord** - optional. Default: 5000. Number of times to try to generate each coordinate for each alpha carbon in the structure.  
**bond_length** - optional. Default : 3.8. Distance between each bond in angstroms. Probably leave this alone.
**clash_dist** - optional. Default: 3.4. How far atoms must be apart before being considered to be clashing.  
**min_bond_dist** - optional. Default: 2.8. Minimum distance between each bond in the structure before raising an exception. AF2 has some shorter ones, so probably leave this alone.  
**max_bond_dist** - optional. Default: 4.2. Maximum distance between each bond in the structure before raising an exception. AF2 has some longer ones, so probably leave this alone.  


  
  
**Starting from a PDB file.**  
You can also have DODO alter a pre-existing AF2 pdb file. The AF2 file should have all atom information if you don't want to use metapredict to identify the IDRs. Otherwise it must use metapredict. There is obe required argument: ``path_to_pdb`` : the path to your pdb file as a string and. If you don't specify the ``out_path``, DODO will return a dictionary that really won't make a lot of sense to you because I haven't finished making it useful. So for the time, specify the out_path.

```python
build.pdb_from_pdb('/Users/your_user_name/Desktop/my_AF2_pdb.pdb', out_path='/Users/your_user_name/Desktop/my_AF2_PDB_DODO.pdb')
```
**Additional usage:**  

All arguments for ``build.pdb_from_name`` are as follows:  
**path_to_pdb** - required. The filepath to your pdb as a string.  
**out_path** - optional. Where to save your protein structure file. Specify the file name here.  
**mode** - optional. Default: 'predicted'. The ``predicted`` option predicts the end-to-end distance of your disordered regions from sequence and then makes the IDRs fit within that distance. For IDRs between folded domains, the folded domains are moved apart from each other to accomodate the IDR length. For terminal IDRs, they will just adopt a configuration within the distance equal to the end to end distance. Additional options are ``super_compact``, ``compact``, ``normal``, ``expanded``, ``super_expanded``, ``max_expansion``. These are pretty self explanatory.  
**use_metapredict** - optional. Default: False. This option lets you use metapredict to predicte the IDRs and folded regions. Although *fairly accurate*, it doesn't get the exact cutoffs for some regions and fails to predict small loops within large folded regions. The default is to use the number of atoms neighboring each atom in the AF2 structure to figure out what is a folded domain, what is a separate folded domain that is proximal to another folded domain, what is an IDR, and what is a disordered loop. The default behavior is slower but works better.  
**graph** - optional. Default: False. Setting this to True will pull up a really rough looking structure of your protein using the 3D graphing functionality in matplotlib. This is something I made when developing this to quickly look at structures. You shouldn't use this, but you can if you want. It's kind of fun TBH.  
**beta_by_region** optional. Default: True. This takes the beta values and changes them such that if you view the amino acids in something like VMD by beta values, you can see what DODO predicted to be an IDR and what it left as a folded domain. Loops will also be a different color if they were present.  
**silent** optional. Default: False. This silences the message telling you what structure is being downloaded. Turn this off at your own risk. DODO does not gaurantee that it will download the structure you wanted from name alone. If you just put 'p53' instead of 'human p53', you might get a different organism based on how Uniprot ranks it.  
**verbose** optional. Default: False. Set to True to get more info on what is happening as your structure is being made.  
**attempts_per_region** - optional. Default: 50. Number of times to try and make each region.  
**attempts_per_coord** - optional. Default: 5000. Number of times to try to generate each coordinate for each alpha carbon in the structure.  
**bond_length** - optional. Default : 3.8. Distance between each bond in angstroms. Probably leave this alone.  
**clash_dist** - optional. Default: 3.4. How far atoms must be apart before being considered to be clashing.  
**min_bond_dist** - optional. Default: 2.8. Minimum distance between each bond in the structure before raising an exception. AF2 has some shorter ones, so probably leave this alone.  
**max_bond_dist** - optional. Default: 4.2. Maximum distance between each bond in the structure before raising an exception. AF2 has some longer ones, so probably leave this alone.  
  
  
  
**Starting n IDR sequence.**  
If you just have an IDR sequence, you can also generate coordinates for that. Functionality is coming to automatically generate the PDB file as well, but it's near midnight so that's for later. The only required argument is the sequence. 

```python
build.idr('VQQQGIYNNGTIAVANQVSCQSPNQ')
```
**Additional usage:**  

All arguments for ``build.pdb_from_name`` are as follows:  
**mode** - optional. Default: 'predicted'. The ``predicted`` option predicts the end-to-end distance of your disordered regions from sequence and then makes the IDRs fit within that distance. For IDRs between folded domains, the folded domains are moved apart from each other to accomodate the IDR length. For terminal IDRs, they will just adopt a configuration within the distance equal to the end to end distance. Additional options are ``super_compact``, ``compact``, ``normal``, ``expanded``, ``super_expanded``, ``max_expansion``. These are pretty self explanatory.  
**all_atoms** - optional. Default: False. Setting to True will result in DODO adding in atoms other than the alpha carbon (CA). NOTE: This is VERY MUCH in beta. It is NOT sophisticated AT ALL whatsoever. However, for visualization purposes, it's not bad. Given this is largely for visualization, I probably won't make it too much better to be honest but you never know. I'll add this to the other functions today or tomorrow (today?). It's late. 
**graph** - optional. Default: False. Setting this to True will pull up a really rough looking structure of your protein using the 3D graphing functionality in matplotlib. This is something I made when developing this to quickly look at structures. You shouldn't use this, but you can if you want. It's kind of fun TBH.  
**beta_by_region** optional. Default: True. This takes the beta values and changes them such that if you view the amino acids in something like VMD by beta values, you can see what DODO predicted to be an IDR and what it left as a folded domain. Loops will also be a different color if they were present.  
**silent** optional. Default: False. This silences the message telling you what structure is being downloaded. Turn this off at your own risk. DODO does not gaurantee that it will download the structure you wanted from name alone. If you just put 'p53' instead of 'human p53', you might get a different organism based on how Uniprot ranks it.  
**verbose** optional. Default: False. Set to True to get more info on what is happening as your structure is being made.  
**attempts_per_region** - optional. Default: 50. Number of times to try and make each region.  
**attempts_per_coord** - optional. Default: 5000. Number of times to try to generate each coordinate for each alpha carbon in the structure.  
**bond_length** - optional. Default : 3.8. Distance between each bond in angstroms. Probably leave this alone.  
**clash_dist** - optional. Default: 3.4. How far atoms must be apart before being considered to be clashing.  
**min_bond_dist** - optional. Default: 2.8. Minimum distance between each bond in the structure before raising an exception. AF2 has some shorter ones, so probably leave this alone.  
**max_bond_dist** - optional. Default: 4.2. Maximum distance between each bond in the structure before raising an exception. AF2 has some longer ones, so probably leave this alone.  



#### Changes

Logging changes below.

V0.02 - September 26, 2023. Added generating IDR coords from sequence alone. Added filling in ATOM coordinates from alpha carbon coordinate using fixed bond angles and distances, which isn't great but is better than nothing. 

V0.01 - September 25, 2023. Initial release.  

### Copyright

Copyright (c) 2023, Ryan Emenecker - Holehouse Lab


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
