import numpy as np
# Creates a list of centers (placement of the forcings)
# and a list of ranges (the nonzero entries that are associated)
# to this forcing. 
def ranges_and_centers(bandwidth, n, offset=0):
    ranges = []
    centers = []
    i = offset
    # Making sure that at one more set fits with sufficient distance
    # to the periodically wrapped-around last element. 
    while i + (2 * bandwidth) < n + offset: 
        if i + (4 * bandwidth + 1) < n + offset:
            ranges.append(range(i, i + 2 * bandwidth + 1))
            centers.append(i + bandwidth)
            i = i + 2 * bandwidth + 1
        else:
            ranges.append(range(i, n))
            # possible can modify to make sure that the center is 
            # actually in the center of the range of influence
            centers.append(i + bandwidth)
            i = n
    return ranges, centers

def construct_cbar(centers, n, pick=None):
    cbars = np.zeros(n)
    if pick is None:
        for c in centers:
            cbars[c] = 1.0
    else:
        cbars[centers[pick]] = 1.0
    return cbars

def deconstruct_sbar(ranges, sbar):
    sbars = np.zeros((len(sbar), len(ranges)))
    for j in range(len(ranges)):
        for i in ranges[j]:
            sbars[i, j] = sbar[i]
    return sbars
