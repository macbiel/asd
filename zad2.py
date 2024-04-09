#Used in all_proper_subsets()
from itertools import chain
from itertools import combinations

#Returns a set of all proper subsets of `x`, sans the empty set.
#Sets within this set are ordered from largest to smallest.
#The relative ordering of `x`'s elements is preserved within each subset.
#Based on powerset() at https://docs.python.org/2/library/itertools.html#recipes
def all_proper_subsets(x):
    return chain.from_iterable(combinations(x, r) for r in reversed(range(1,len(x))))


#Used to test if a set of enhancers is fully non-overlapping.
#Requires entries of `enhancers` to be sorted by end position.
def fully_non_overlapping(enhancers):

    #Loop over `enhancers` backwards and check that each element's starts
    #"after" its previous element ends.
    i = len(enhancers) - 1
    while i > 0:
        if enhancers[i-1][1] >= enhancers[i][0]:
            #If an overlap was found, stop checking further and exit early
            return False
        i -= 1

    #If this point was reached, no overlap was found,
    #so `enhancers` is fully non-overlapping.
    return True


def find_the_best_nonoverlaping_enhancers(enhancers):

    #Check input validity
    if enhancers == None:
        return ([], 0)

    #Check for trivial input
    l = len(enhancers)
    if l == 0:                                      #Empty list
        return ([], 0)
    elif l == 1:                                    #1-element list
        return (enhancers, enhancers[0][2])
    elif fully_non_overlapping(enhancers):          #`enhancers` has no overlaps
        return (enhancers, sum(i[2] for i in enhancers))


    #Prepare to find optimal solution

    #Set of all non-empty proper subsets of `enhancers`.
    #Guaranteed to contain the optimal solution.
    subsets = all_proper_subsets(enhancers)

    #This list will contain subsets which were found to be fully non-overlapping
    fno_subsets = []

    #Variables to keep track of the best solution found so far...
    best_subset = []
    #...along with its binding site total
    best_subset_bs = 0


    #Check if each subset of `enhancers` is fully non-overlapping
    #and add it to `fno_subsets` if so,
    #unless `fno_subsets` contains a previous subset which is a superset
    #of the current subset.
    for subset in subsets:

        #Look through `fno_subsets` to check if a superset of `subset`
        #   exists therein
        superset_exists = False
        for prev_subset in fno_subsets:

            #Somewhat non-obvious optimization: given that
            #1. if X is Y's proper superset, |X| > |Y|, and
            #2. `fno_subsets` elements are ordered from largest to smallest,
            #once this loop hits a previous subset of size <= current subset's,
            #it is certain that no superset exists.
            if len(prev_subset) <= len(subset):
                break

            if set(prev_subset).issuperset(subset):
                superset_exists = True
                break


        #Call fully_non_overlapping() only if no superset was found
        if not superset_exists:
                if fully_non_overlapping(subset):

                    #Add `subset` to `fno_subsets`
                    fno_subsets.append(subset)

                    #Update `best_subset` and `best_subset_bs` if needed
                    subset_bs = sum(i[2] for i in subset)
                    if subset_bs > best_subset_bs:
                        best_subset_bs = subset_bs
                        best_subset = subset

        #If a superset of `subset` exists in `fno_subsets`,
        #`subset` itself won't be added to `fno_subsets`.
        #
        #This is fine, though,
        #as any later subsets for whom the current subset is a superset
        #will also match the current subset's superset
        #(if X is Y's superset and Y is Z's superset, X is Z's superset).
        #
        #Not calculating its binding site count is also fine,
        #since the superset is guaranteed to have a higher binding site count
        #(this would not hold if the binding site count could be negative).


    #Finally, return optimal solution in requested format.
    #Iterating over all_proper_subsets()'s return value
    #produces a tuple, rather than a list, hence a cast to list is required.
    return (list(best_subset), best_subset_bs)