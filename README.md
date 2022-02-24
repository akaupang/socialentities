
## **Social Entities**


#### **Description**

A PyMOL/python3 script that returns a selection of entities (atoms, residues, etc) that have a certain number of other entities (from other objects) around them within a specified distance. Depending on the input data, there may be several applications for such a number- and proximity based "clustering" function.

#### **Example usage - structural water molecules**

The term "structural" water molecules refers to water molecules that are important for the integrity of a structure. Although, in principle, these water molecules are not part of the covalent network, they can contribute to the globularity of a structure in solution, to the packing of a structure in the crystal phase or to the binding of a ligand in a binding pocket. In the former role, such water molecules can be deeply buried and practically non-exchangeable with those of the bulk solvent. Water molecules in the latter roles may be exchangeable, but can be critical to include when setting up systems for molecular dynamics or docking/virtual screening.

One way to identify such water molecules from structures of proteins/protein complexes determined by x-ray crystallography, is to look for water molecules (more precisely; their oxygen atoms) that are present in several crystals e.g. across crystallization conditions, different multimeric complexes, with different ligands, etc. This approach could consist of searching a set of aligned (superimposed) structures, that have their crystallographically determined water molecules intact, for locations in which water molecules appear across many of the structures. This script enables the selection of those water molecules whose positions are corroborated as structurally relevant by the presence of water molecules from other protein structures in their vicinity. The search algorithm is straightforward (brute force), and the result is exhaustive, yet somewhat redundant. To run as efficiently as possible, the script takes advantage of paralellization throught the *multiprocess* package, and provides a progress bar from the *tqdm* package.

- Load the script in PyMOL with:

`run path_to_script/SocialEntities.py`

- The general syntax is:

`socialEntities(string, selectionA, selectionB, fraction, distance, string='Entities', cpus='all')`
                             
- Example:

`socialEntities('allWaters', 'resn HOH', 'resn HOH', 0.9, 1.5, 'Waters', 6)`

This creates a selection ('allWaters') of all the water molecules ("resn HOH", selectionA), and then searches the vicinity (1.5 Å) of each of these for other water molecules ('resn HOH', selectionB) from all the loaded objects. If a 
water molecule in 'allWaters' is part of a cluster in which it is surrounded by water molecules from a  minimum of 90% (0.9) of the other objects (e.g. the other protein structures), at a maximum distance of 1.5 Å, it is added to the output selection.

A part of the name of the produced output selections can be modified to reflect the type of entity you are investigating (e.g. 'Waters'). If a modified name is not provided, the modifiable part of the output selection name will be set to 'Entities'. The search distinguishes between the following three groups;

- *socialEntities* - those entities that have neighbours from at least the specified percentage of the input objects, within the specified distance.
- *pseudosocialEntities* - those entities that have fewer neighbours than specified, within the specified distance.
- *lonelyEntities* - those entities that have no neighbours within the specified distance.
                        
The number of CPUs used (default, cpus='all') can also be set to a number <= the number of available CPUs/threads to allow for running in the background while keeping resources free for other work.


#### **Authors**

Åsmund Kaupang
 
#### **Acknowledgements**
 
Johannes Karwounopoulos, for suggesting the inclusion of multiprocess and tqdm.
