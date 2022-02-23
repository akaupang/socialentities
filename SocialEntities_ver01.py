#!/usr/bin/env python3
#
###############################################################################
#
#         Name: Social Entities
#
#       Author: Åsmund Kaupang, 2018-2022
#    Thank you: Johannes Karwounopoulos, for suggesting the inclusion
#               of multiprocess and tqdm!
#
#  Description: A PyMOL/python3 script that returns a selection of entities
#               (atoms, residues, etc) that have a certain number of other
#               entities (from other objects) around them within a specified
#               distance. Depending on the input data, there may be several
#               applications for such a function.
#
#               Here is an example of a usage in which we look for structural
#               water molecules:
#               One way to look for structural water molecules, i.e. water
#               molecules that appear in several protein crystal structures
#               (perhaps across crystallization conditions, unit cells, space
#               groups or multimeric complexes), would be to superimpose the 
#               protein structures and see whether the locations of the 
#               crystallographically determined water molecules coincide. This
#               script helps you select those waters whose positions are
#               corroborated as structurally relevant by the presence of water
#               molecules from other protein structures in their vicinity.
#
#               Load the script with:
#               run SocialEntities.py
#
#               To select those water molecules ("resn HOH") from a selection
#               of all the water molecules of all the loaded objects 
#               ('allWaters'), that are surrounded by a minimum of 5 water
#               molecules belonging to other objects (e.g. other protein
#               structures) at a maximum distance of 2 Å from their position,
#               issue: 
#
#  UPDATE TEXT  socialEntities('allWaters', 'resn HOH', 'resn HOH', 0.90, 1, 'Waters')
#
#               A part of the name of the produced selections can be modified
#               to reflect the type of entity you are investigating ('Waters')
#               or if not provided, it will appear as 'Entities', in the
#               following three groups;
#
#               socialEntities          Those entities that have neighbours
#                                       within the specified criteria (number,
#                                       distance).
#               pseudosocialEntities    Those entities that have fewer
#                                       neighbours than specified within the
#                                       specified distance.
#               lonelyEntities          Those entities that have no neighbours
#                                       within the specified distance.
#
################################################################################

from pymol import cmd
from pymol import stored
from tqdm import tqdm
from itertools import repeat, chain
import multiprocess as mp


def modIndTupleListToSelctionString(modindtuplelist):
    print(f"start to go through {len(modindtuplelist)} items")
    if len(modindtuplelist) == 0:
        selection = f"none"
    else:
        for idx, modind in enumerate(modindtuplelist):
            # First pass; create the start of the selection string
            if idx == 0:
                selection = f"({modind[0]} and index {modind[1]})"
            # Subsequent passes; grow the selection string
            else:
                selection = f"{selection} or ({modind[0]} and index {modind[1]})"
    print("end loop")
    return selection

def neighbourSearch(start, end, objectindexlist, poolSel, numberOfEntitiesAround):
    #print(mp.current_process())

    # Initiate lists for entities with enough found neighbours, for entities
    # with too few neighbours and for entities with no neighbours. 
    nblist = []
    fewnblist = []
    nonblist = []

    for entitynum in tqdm(range(start, end)):

        # Define a selection for the neighbours of the current entity, that are present in the swarm/pool.
        # Recall that poolSel is:
        # '(allEntSelName and nbEntTypeSel) near_to distanceToEntitiesAround of ', e.g.
        # '(allwaters and resn HOH) near_to 2 of '
        nbSel = poolSel + '(' + str(objectindexlist[entitynum][0]) + ' and index ' + str(objectindexlist[entitynum][1]) + ')'
        # which then results in e.g.:
        # '(allwaters and resn HOH) near_to 2 of (2ABC_A and index 3052)', in which 3052 is the current index dictated by entitynum

        # Get the indexes of these neighbours defined by the selection
        nbs = cmd.index(nbSel)

        # If no neighbours are found, add the entity to a list of lonely entities (no neighbours) 
        # If neighbours are found, evaluate if there are enough of them, and if so, add the 
        # (model, index) tuple to a list of social entities (enough neighbours). If there are
        # some, but to few neighbours, add the entity to a list of less lonely entities 
        # (few neighbours).
        if nbs == []:
            nonblist.append(objectindexlist[entitynum])
        elif len(nbs) >= numberOfEntitiesAround:
            nblist.extend(nbs)
        else:
            fewnblist.append(objectindexlist[entitynum])


    return nonblist, fewnblist, nblist

def socialEntities( allEntSelName, allEntTypeSel, nbEntTypeSel, objFractOfEntitiesAround, distanceToEntitiesAround, newSelectionType='Entities', cpus='all' ):

    # Initiate a list in the PyMOL name space to hold the models and indexes
    stored.modelsindexes = []

    # Select all of the pool entities
    cmd.select(allEntSelName, allEntTypeSel)

    # Fill the list with the (model, index) tuples from the user selection
    cmd.iterate(allEntSelName, 'stored.modelsindexes.append((model,index))')
    
    # Patch together a generic opening to the neighbour selection string
    poolSel = '(' + allEntSelName + ' and ' + nbEntTypeSel +  ') near_to ' + str(distanceToEntitiesAround) + ' of '

    # Get the total number of input entities ('the whole swarm')
    # stored.modelsindexes is a list of tuples; ..., (object/model name, index number), ... 
    totnumentities = len(stored.modelsindexes)

    # Get the number of objects (for later)
    numberOfObjects = len(cmd.get_object_list())

    # Evaluate the requested fraction of objects containing neighbours
    numberOfEntitiesAround = int(objFractOfEntitiesAround * numberOfObjects)
   
    # Get the number of CPUs/threads available to divide the workload
    availcpus = mp.cpu_count()

    # Moderate CPU usage if requested by user
    if cpus == 'all':
        cpus = availcpus
    elif cpus > availcpus:
        cpus = availcpus
    else:
        pass

    # Print some information before the calculations start
    print(#
        f"The requested neighbour fraction ({objFractOfEntitiesAround*100}%) corresponds to neighbours from {numberOfEntitiesAround} of {numberOfObjects} objects.\n" \
        f"The total number of entities of type <{allEntTypeSel}> whose neighbourhoods will be searched for neighbours of type <{nbEntTypeSel}> is {totnumentities}.\n" \
        f"{availcpus} CPUs/threads were found and {cpus} will be used in the calculation. "
        )

    # Define the pool based on the CPU count
    pool = mp.Pool(cpus)

    # The map function applies a function (the first argument) to each of the 
    # items in an interable, which is passed as the second argument.
    # If more than one iterable is passed, the number of iterables must 
    # correspond to the number of arguments taken by the function
    # If the iterables are of different lengths, the iteration through them
    # will stop when the shortes iterable is exhausted.
    #
    # The map function returns a map object, which is an iterator. This object
    # can be transformed to a sequence object (a list, a tuple, etc) by calling
    # the respective constructor; list(), tuple(), etc
    # The map object can also be iterated through using a for loop
    #
    # If starmap is used instead of map, the iterators must be "pre-zipped" so
    # that the input is a zip object. The zip object can be visualized when
    # converted to a list, and is then a list of tuples, where the length of 
    # each tuple corresponds to the number of iterators zipped together.
    # The zip function can work to transpose a list of lists

    # Partition the entities equally (the last part is longer or shorter) for
    # each CPU. 
    wholepartslength = totnumentities // cpus
    numwholeparts = (totnumentities // wholepartslength) - 1

    # In practice, the indices that refer to their placement in the
    # stored.modelsindexes are noted down and added to lists of start-
    # and end indices
    startindexes = []
    endindexes = []

    for idx, i in enumerate(range(0, totnumentities, wholepartslength)):
        if idx  < numwholeparts:
            startindexes.append(i)
            endindexes.append(i + wholepartslength - 1)
        elif idx  == numwholeparts:
            startindexes.append(i)
            endindexes.append(totnumentities - 1)
        else:
            pass

    # Now we throw a bone to each dog.
    # The input data structure is the zip transposition of the following lists:
    # [ start1, start2, start3, ...]
    # [   end1,   end2,   end3, ...]
    # [ MI[()],[MI[()],[MI[()], ...]
    # [poolSel,poolSel,poolSel, ...]
    # [ numOEA, numOEA, numOEA, ...]
    #
    # Where '[MI[()]' is [[(m1,i1),(m1,i2),...,(m2,i1),(m2,i2),...], 
    #                     [(m1,i1),(m1,i2),...,(m2,i1),(m2,i2),...],
    #                                                          ...]
    #
    result = pool.starmap(neighbourSearch, zip(startindexes, endindexes, repeat(stored.modelsindexes), repeat(poolSel), repeat(numberOfEntitiesAround)))

    # End the parallelized part of the script
    pool.close()
    #pool.join()

    # The output data structure is a list of lists, whose length in the top
    # level is defined by the number of processes/workers/forks/cpus. Each
    # process has a number of sublists (here, three) defined by the number of
    # output variables from the function.

    # Gather the outputs of each process to reconstitute the outputs of the
    # neighbourSearch function and make new lists containing only the unique
    # entities in each raw listv(as neighbours will have been found multiple
    # times)

    nonblist = list(set(chain(*[result[i][0] for i in range(0, len(result))])))
    fewnblist = list(set(chain(*[result[i][1] for i in range(0, len(result))])))
    nblist = list(set(chain(*[result[i][2] for i in range(0, len(result))])))

    # Print some rudimentary statistics to help tune the search parameters
    frcnonb = len(nonblist)/totnumentities
    frcfewnb = len(fewnblist)/totnumentities
    frcnb = len(nblist)/totnumentities


    print(#
        f"The requested neighbour fraction ({objFractOfEntitiesAround*100}%) corresponds to neighbours from {numberOfEntitiesAround} of {numberOfObjects} objects.\n" \
        f"Out of {totnumentities} entities of type <{allEntTypeSel}>, across {numberOfObjects} objects:\n" \
        f"{frcnb*100:.1f}% had >= {numberOfEntitiesAround} neighbours of type <{nbEntTypeSel}> within {distanceToEntitiesAround} Å\n" \
        f"{frcfewnb*100:.1f}% had < {numberOfEntitiesAround} neighbours of type <{nbEntTypeSel}> within {distanceToEntitiesAround} Å\n" \
        f"{frcnonb*100:.1f}% had no neighbours of type <{nbEntTypeSel}> within {distanceToEntitiesAround} Å" \
          )
    if numberOfEntitiesAround > numberOfObjects:
        print(f"WARNING: The number of entities in the search is larger than the number of objects. Please check your input!")
    else:
        pass

    # Get a selection string from each list of unique entities
    print('Preparing selections and visualization. Please wait...')
    nbsel = modIndTupleListToSelctionString(nblist)
    #fewnbsel = modIndTupleListToSelctionString(fewnblist)
    #nonbsel = modIndTupleListToSelctionString(nonblist)
    print('selecting in pymol')
    # Export these selections to the PyMOL GUI
    cmd.select(f"social{newSelectionType}", nbsel)
    #cmd.select(f"pseudosocial{newSelectionType}", fewnbsel)
    #cmd.select(f"lonely{newSelectionType}", nonbsel)

    #cmd.show_as(f"nb_spheres", f"({nbsel} or {fewnbsel} or {nonbsel})")
    cmd.show_as(f"nb_spheres", nbsel)
    #cmd.color(f"red", nonbsel)
    #cmd.color(f"orange", fewnbsel)
    cmd.color(f"green", nbsel)
    cmd.deselect()

    print('Done')

if __name__ == "__main__":
    socialEntities()

# Make the function runnable in PyMOL
cmd.extend("socialEntities", socialEntities)

