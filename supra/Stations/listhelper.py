
def rmBelow(lst, cutoff):
    ''' 
    Returns list of indicies of points equal to or above the cutoff in a list

    Usage:
    rmBelow([100, 200, 300], 200)
    >>> [1, 2]

    Arguments:
    lst [list] - list to remove points from
    cutoff [float] - number which if the elements of lst are below, to remove

    Returns:
    a list of indicies where lst is equal to or above the cutoff
    '''
    try:

        return [i for i, element in enumerate(lst) if element >= cutoff]

    # If list is not floats    
    except TypeError:
        
        return None



def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than maxgap

        I found this from here, I had something similar, but this code works much better:
        https://stackoverflow.com/questions/14783947/grouping-clustering-numbers-in-python
        it has been modified to work with pick objects

        Usage:
        (where the numbers below are the first point of Pick Objects)
        >>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]
        
        Arguments:
        data [list] - list of Pick Objects to cluster together
        maxgap [float] - maximum cluster size (from first element to last element)

        Returns:
        groups [list] - a list of clusters, each of which is its own list of Pick Objects

    '''
    if len(data) == 0:
        return None

    # First point will be in the first group
    groups = [[data[0]]]

    for x in data[1:]:

        # Join current group if within maxgap
        if abs(x - groups[-1][0]) <= maxgap:
            groups[-1].append(x)
        
        # Create new group
        else:
            groups.append([x])

    return groups

def degroupPicks(groups):

    ''' Takes the output of the cluster function and combines the picks of each group into one pick

    Usage:
    where Pick Objs are shown as P[first_pt, last_pt]
    degroupPicks([[P[1, 3], P[3, 4]], [P[7, 8]], [P[10, 12], P[11, 13]]])
    >>> [P[1, 4], P[7, 8], P[10, 13]]

    Arguments:
    groups [list] - output of the cluster function, a list of groups of clustered pick objects, 
                    each of which is it's own list

    Returns:
    filtered_list [list] - a list of Pick Objects, each cluster of Pick Objects from the input combined
                            into one Pick Object
    '''
    if groups is None:
        return None

    filtered_list = []
    for cluster in groups:

        #join all elements of a cluster together

        # Automatically add the first pick
        result = cluster[0]

        for pick in cluster[1:]:
            result = result + pick

        filtered_list.append(result)
    
    return filtered_list

if __name__ == "__main__":

    #####################
    #### Pick Tests
    #####################

    pick_list = [Pick(0.5, 1.5), Pick(1.6, 2.0), Pick(2.1, 2.2), Pick(2.2, 2.8)]

    def picktestHelperFunc(pick_list, pick_len):

        a = degroupPicks(cluster(pick_list, pick_len))
        lst = []
        for i in a:
            lst.append([i.first_pt, i.last_pt])
        return lst

    assert picktestHelperFunc(pick_list, 1) == [[0.5, 1.5], [1.6, 2.8]]
    assert picktestHelperFunc(pick_list, 0.1) == [[0.5, 1.5], [1.6, 2.0], [2.1, 2.2], [2.2, 2.8]]

    pick_list = [Pick(0.5, 0.6), Pick(0.5, 0.7)]

    assert picktestHelperFunc(pick_list, 0.1) == [[0.5, 0.7]]

    pick_list = [Pick(0.5, 0.7), Pick(0.5, 0.6)]

    assert picktestHelperFunc(pick_list, 0.1) == [[0.5, 0.7]]

    print("Pick tests complete!")

    ####################
    #### Test rmbelow
    ####################

    testing_list = [1, 2, 3, 4, 5, 10, 24, 36, 0, -1]
    bad_list = ["A", "B", 2]

    # Returns elements >= cutoff
    assert rmBelow(testing_list, 3) == [2, 3, 4, 5, 6, 7]
    assert rmBelow(testing_list, 0) == [0, 1, 2, 3, 4, 5, 6, 7, 8]
    assert rmBelow(testing_list, -2) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    assert rmBelow(testing_list, 40) == []
    assert rmBelow(bad_list, 2) is None

    print("Remove below testing complete!")
