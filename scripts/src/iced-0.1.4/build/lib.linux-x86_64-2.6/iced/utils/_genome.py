import numpy as np


def get_intra_mask(lengths, resolution=10000):
    """
    Returns a mask for intrachromosomal interactions

    Parameters
    ----------
    lengths : ndarray, (n, )
        lengths of the chromosomes

    resolution : int, optional, default: 10000
        resolution in which the mask should be created

    Returns
    -------
    mask : ndarray (m, m)
        boolean mask
    """
    # Copy the lengths in order not to modify the original matrix
    lengths = _change_lengths_resolution(lengths, resolution=resolution)
    mask = np.zeros((lengths.sum(), lengths.sum()))
    begin = 0
    for end in lengths.cumsum():
        mask[begin:end, begin:end] = 1
        begin = end
    return mask.astype(bool)


def get_inter_mask(lengths, resolution=10000):
    """
    Returns a mask for interchromosomal interactions

    Parameters
    ----------
    lengths : ndarray, (n, )
        lengths of the chromosomes

    resolution : int, optional, default: 10000
        resolution in which the mask should be created

    Returns
    -------
    mask : ndarray of dtype boolean
        boolean mask
    """
    intra_mask = get_intra_mask(lengths, resolution=resolution)
    return np.invert(intra_mask)


def get_genomic_distances(lengths, resolution=10000):
    """
    Returns a matrix of the genomic distances

    Inter chromosomal interactions are set to -1

    Parameters
    ----------
    lengths : ndarray (n, )
        lengths of the chromosomes

    resolution : int, optional, default: 10000
        resolution in which to compute the matrix

    Returns
    -------
    dis: ndarray (n, n), dtype: int
        returns the genomic distance matrix, with -1 for inter chromosomal
        interactions
    """
    inter_mask = get_inter_mask(lengths, resolution=resolution)
    lengths = _change_lengths_resolution(lengths, resolution=resolution)
    n = lengths.sum()

    dis = np.concatenate([np.concatenate([np.arange(i, 0, -1),
                                          np.arange(n - i)])
                          for i in range(n)])
    dis = dis.reshape((n, n))
    dis[inter_mask] = -1

    return dis.astype(int)


def undersample_per_chr(X, lengths, resolution=10000):
    """
    Undersample matrix to chromosomes

    Undersample the matrix ununiformaly per chromosome lengths.

    Parameters
    ----------
    X : ndarray (n, n)
        The matrix to undersample

    lengths : ndarray (L, )
        Lengths of the chromosomes

    resolution : integer, optional, default: 10000
        Resolution in which the matrix `X` is provided

    Returns
    -------
    undersampled_X : ndarray (L, L)
        `X` undersampled per chromosome
    """
    lengths = _change_lengths_resolution(lengths, resolution=resolution)
    lengths_cum = lengths.cumsum()
    chr1_begin = 0
    undersampled_X = np.zeros((len(lengths), len(lengths)))
    for i, chr1_end in enumerate(lengths_cum):
        chr2_begin = 0
        for j, chr2_end in enumerate(lengths_cum):
            surface = X[chr1_begin:chr1_end, chr2_begin:chr2_end]
            undersampled_X[i, j] = surface[np.invert(np.isnan(surface))].mean()
            chr2_begin = chr2_end
        chr1_begin = chr1_end

    return undersampled_X


def _change_lengths_resolution(lengths, resolution=10000, copy=True):
    if copy:
        lengths = lengths.copy()
    lengths = lengths.astype(float)
    return np.ceil(lengths / resolution).astype(int)


def get_chromosome_counts(counts, lengths, chromosome):
    """
    """
    if chromosome > len(lengths) - 1:
        raise ValueError("Chromosome %d doesn't exists. Possible values are "
                         "from 0 to %d" % (chromosome, len(lengths) - 1))

    if len(counts) != lengths.sum():
        raise ValueError("The total lengths and the counts matrix shape "
                         "should be the same. They are respectively %d and %d"
                         % (lengths.sum(), len(counts)))

    lengths_cum = np.concatenate([[0], lengths.cumsum()])
    return counts[lengths_cum[chromosome]:lengths_cum[chromosome + 1],
                  lengths_cum[chromosome]:lengths_cum[chromosome + 1]]
