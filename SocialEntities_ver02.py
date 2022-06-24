#!/usr/bin/env python3
import sys
from pymol import cmd, stored
from tqdm.auto import tqdm
from itertools import repeat, chain
import multiprocess as mp

from time import sleep
import psutil

def modIndTupleListToSelectionString(modindtuplelist):
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
    return selection

def incrementalPyMOLSelection(longtuplelist, partlength, finalselectionname):
    partnumber = 1
    tuplelistlength = len(longtuplelist)
    partselnames = []

    for partstart in range(0, tuplelistlength, partlength):
        partend = partstart + partlength # will exceed the tuple list length in
                                         # the last iteration, but this is
                                         # inconsequential

        # Transform the current part of the tuple list into a PyMOL selection 
        print(f"\nSelection '{finalselectionname}', part {partnumber}:")
        print(f"Go through {len(longtuplelist[partstart:partend])} tuples.")
        partsel = modIndTupleListToSelectionString(longtuplelist[partstart:partend])
        print(f"Made selection strings")
####### Strategy 1 - temporary, persistant selections, then fuse and
      # delete intermediates
        print(f"Selecting part {partnumber} in PyMOL")
        cmd.select(f"ps{partnumber}", partsel)
        partselnames.append(f"ps{partnumber}")
        print(f"Done selecting part {partnumber} in PyMOL")
        partnumber = partnumber + 1
    #
    # "Select the selections" to join them
    cmd.select(finalselectionname, ' or '.join(partselnames))
    #
    # Delete the temporary selections
    for tmpsel in partselnames:
        cmd.delete(tmpsel)
#######

####### Strategy 2 - temporary, non-persistant selection, with gradual
      # selection expansion

 

def neighbourSearch(start, end, objectindexlist, poolSel, numberOfEntitiesAround):
    #print(mp.current_process())

    # Initiate lists for entities with enough found neighbours, for entities
    # with too few neighbours and for entities with no neighbours. 
    ewnblist = []
    ewfewnblist = []
    ewnonblist = []

    nblist = []
    fewnblist = []

    for entitynum in tqdm(range(start, end)):
        
        # Select those entities in the pool/swarm that are near to a given single entity
        # Recall that poolSel is:
        # 'allEntSelName near_to distanceToEntitiesAround of ', e.g.
        # 'allwaters near_to 2 of ' , which becomes
        nbSel = f"{poolSel} ({objectindexlist[entitynum][0]} and index {objectindexlist[entitynum][1]})"
        # which then results in e.g.:
        # 'allwaters near_to 2 of (2ABC_A and index 3052)', in which 2ABC_A is the current object and 3052 is the current index, as dictated by objectindexlist[entitynum]
        
        # Get the indexes of these neighbours defined by the selection
        nbs = cmd.index(nbSel)
        
        # If no neighbours are found, add the entity to a list of lonely entities (no neighbours) 
        # If neighbours are found, evaluate if there are enough of them, and if so, add the 
        # (model, index) tuple to a list of social entities (enough neighbours). If there are
        # some, but to few neighbours, add the entity to a list of less lonely entities 
        # (few neighbours).
        
        # Reference-centered (current entity)
        #
        # The entity Has no neighbours - it is lonely
        if nbs == []:
            ewnonblist.append(objectindexlist[entitynum])
        # The entity has enough neighbours - it is social
        elif len(nbs) >= numberOfEntitiesAround:
            ewnblist.append(objectindexlist[entitynum])
        # The entity has fewer neighbours than specified - it is pseudosocial
        else:
            ewfewnblist.append(objectindexlist[entitynum])
 
        # Neighbour-centered
        #
        # There are no neighbours around
        if nbs == []:
            pass
        # There are enough neighbours around. Document these.
        elif len(nbs) >= numberOfEntitiesAround:
            nblist.extend(nbs)
        # There fewer neighbours around. Document these.
        else:
            pass#fewnblist.extend(nbs)

    return ewnonblist, ewfewnblist, ewnblist, nblist

def socialEntities( allEntSelName, allEntTypeSel, nbEntTypeSel, objFractOfEntitiesAround, distanceToEntitiesAround, newSelectionType='Entities', cpus='all', selpartsize=1000, force=False ):

    # check selections for "all", except for in a forced situation
    if force == True:
        print(f"WARNING: Forced mode. Skipping selection sanity checks")
    else:
        if ((allEntSelName.lower() == 'all') or (allEntTypeSel.lower()  == 'all') or (nbEntTypeSel.lower() == 'all')):
            print(f"A selection was set to 'all'. Exiting. Avoid this using force = True at your own peril.")
            sys.exit(1)
        elif objFractOfEntitiesAround == 0:
            sys.exit(1)
        elif distanceToEntitiesAround == 0:
            sys.exit(1)
        else:
            pass

    # Initiate a list in the PyMOL name space to hold the models and indexes
    stored.modelsindexes = []

    # Select all of the pool entities
    cmd.select(allEntSelName, allEntTypeSel)

    # Fill the list with the (model, index) tuples from the user selection
    cmd.iterate(allEntSelName, 'stored.modelsindexes.append((model,index))')
    
    # Patch together a generic opening to the neighbour selection string
    poolSel = f"{allEntSelName} near_to {distanceToEntitiesAround} of "

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
    pool.join()

    # The output data structure is a list of lists, whose length in the top
    # level is defined by the number of processes/workers/forks/cpus. Each
    # process has a number of sublists (here, three) defined by the number of
    # output variables from the function.

    # Gather the outputs of each process to reconstitute the outputs of the
    # neighbourSearch function and make new lists containing only the unique
    # entities in each raw listv(as neighbours will have been found multiple
    # times)
#



#   with tqdm(total=100, desc='cpu%', position=1) as cpubar, tqdm(total=100, desc='ram%', position=0) as rambar:
#       while True:
#           rambar.n=psutil.virtual_memory().percent
#           cpubar.n=psutil.cpu_percent()
#           rambar.refresh()
#           cpubar.refresh()
#           sleep(0.5)

    # Entity-centered
    ewnonblist = list(set(chain(*[result[i][0] for i in range(0, len(result))])))
    ewfewnblist = list(set(chain(*[result[i][1] for i in range(0, len(result))])))
    ewnblist = list(set(chain(*[result[i][2] for i in range(0, len(result))])))

    # Neighbour-centered
    nblist = list(set(chain(*[result[i][3] for i in range(0, len(result))])))
    #fewnblist = list(set(chain(*[result[i][4] for i in range(0, len(result))])))
    print(f"Recorded unique neighbours (to someone): {len(nblist)}")
    print(f"The number of socEnt that figure among these neighbours is {len([ent for ent in nblist if ent in ewnblist])}")

    # Print some rudimentary statistics to help tune the search parameters
    frcnonb = len(ewnonblist)/totnumentities
    frcfewnb = len(ewfewnblist)/totnumentities
    frcnb = len(ewnblist)/totnumentities

    print(f"\n" \
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



    # Name the selections
    ewnbselname = f"social{newSelectionType}"
    ewfewnbselname = f"pseudosocial{newSelectionType}"
    ewnonbselname = f"lonely{newSelectionType}"
    nbselname = f"neighbourhood{newSelectionType}"
    #fewnbselname = f"less-nbSocEnt"
    
    # OLD SELECTION STRATEGY - MEMORY HEAVY AT 35K ENTITIES
    # Get a selection string from each list of unique entities
    #print('Preparing selections and visualization. Please wait...')
    #nbsel = modIndTupleListToSelctionString(nblist)
    #fewnbsel = modIndTupleListToSelctionString(fewnblist)
    #nonbsel = modIndTupleListToSelctionString(nonblist)
    #print('selecting in pymol')

    # Export these selections to the PyMOL GUI
    #cmd.select(nbselname, nbsel)
    #cmd.select(fewnbselname, fewnbsel)
    #cmd.select(nonbselname, nonbsel)

    # NEW INCREMENTAL SELECTION STRATEGY
    # Entity-centered selections
    incrementalPyMOLSelection(ewnblist, selpartsize, ewnbselname)
    #incrementalPyMOLSelection(ewfewnblist, selpartsize, ewfewnbselname)
    #incrementalPyMOLSelection(ewnonblist, selpartsize, ewnonbselname)
    # Neighbour-centered selections
    incrementalPyMOLSelection([item for item in nblist if item not in ewnblist], selpartsize, nbselname)
    #incrementalPyMOLSelection(fewnblist, selpartsize, fewnbselname)
    #cmd.select(f"remneigh", f"{nbselname} and not {ewnbselname}")
    #cmd.delete(f"{nbselname}")
    #cmd.select(f"remlessneigh", f"{fewnbselname} and not {ewfewnbselname}")
    #cmd.delete(f"{fewnbselname}")

    # Show and colour the nb_spheres
    #cmd.show_as(f"nb_spheres", f"({nbselname} or {fewnbselname} or {nonbselname})")
    cmd.show_as(f"nb_spheres", f"{ewnbselname} or {nbselname}")
    cmd.color(f"green", ewnbselname)
    cmd.color(f"yellow", nbselname)
    #cmd.color(f"orange", ewfewnbselname)
    #cmd.color(f"red", ewnonbselname)

    # Deselect all for good measure
    cmd.deselect()

    print('Done')

if __name__ == "__main__":
    socialEntities()

# Make the function runnable in PyMOL
cmd.extend("socialEntities", socialEntities)

