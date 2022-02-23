
Name:           Social Entities
Author:         Åsmund Kaupang, 2018-2022
Thank you:      Johannes Karwounopoulos, for suggesting the inclusion of 
                multiprocess and tqdm!

Description:    A PyMOL/python3 script that returns a selection of entities
                (atoms, residues, etc) that have a certain number of other 
                entities (from other objects) around them within a specified
                distance. Depending on the input data, there may be several 
                applications for such a function.

Here is an example of a usage in which we look for structural water molecules: 
One way to look for structural water molecules, i.e. water molecules that 
appear in several protein crystal structures (perhaps across crystallization 
conditions, unit cells, space groups or multimeric complexes), would be to 
superimpose the  protein structures and see whether the locations of the  
crystallographically determined water molecules coincide. This script helps you 
select those waters whose positions are corroborated as structurally relevant 
by the presence of water molecules from other protein structures in their 
vicinity.

Load the script with:   run SocialEntities.py

The general syntax is:  socialEntities(string, selectionA, selectionB, 
                             fraction, distance, string='Entities', cpus='all')
Example:
socialEntities('allWaters', 'resn HOH', 'resn HOH', 0.9, 1.5, 'Waters', 6)

This creates a selection ('allWaters') of all the water molecules ("resn HOH", 
selectionA), and then searches the vicinity (1.5 Å) of each of these for other 
water molecules ('resn HOH', selectionB) from all the loaded objects. If a 
water molecule in 'allWaters' is part of a cluster in which it is surrounded by 
water molecules from a  minimum of 90% (0.9) of the other objects (e.g. the 
other protein structures), at a maximum distance of 1.5 Å, it is added to the 
output selection.

A part of the name of the produced output selections can be modified to reflect 
the type of entity you are investigating (e.g. 'Waters'). If a modified name is 
not provided, the modifiable part of the output selection name will be set to 
'Entities'. The search distinguishes between the following three groups;

socialEntities          Those entities that have neighbours within the 
                        specified criteria (number, distance).
                        
pseudosocialEntities    Those entities that have fewer neighbours than 
                        specified within the specified distance.
                        
lonelyEntities          Those entities that have no neighbours within the 
                        specified distance.
                        
The number of CPUs used can also be set to a number <= the number of available CPUs/threads to allow for running in the background while keeping resources free for other work.
