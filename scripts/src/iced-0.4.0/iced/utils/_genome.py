import numpy as np
from .validation import is_symetric_or_tri


def get_intra_mask(lengths):
    """
    Returns a mask for intrachromosomal interactions

    Parameters
    ----------
    lengths : ndarray, (n, )
        lengths of the chromosomes

    Returns
    -------
    mask : ndarray (m, m)
        boolean mask
    """
    # Copy the lengths in order not to modify the original matrix
    mask = np.zeros((lengths.sum(), lengths.sum()))
    begin = 0
    for end in lengths.cumsum():
        mask[begin:end, begin:end] = 1
        begin = end
    return mask.astype(bool)


def get_inter_mask(lengths):
    """
    Returns a mask for interchromosomal interactions

    Parameters
    ----------
    lengths : ndarray, (n, )
        lengths of the chromosomes

    Returns
    -------
    mask : ndarray of dtype boolean
        boolean mask
    """
    intra_mask = get_intra_mask(lengths)
    return np.invert(intra_mask)


def get_genomic_distances(lengths):
    """
    Returns a matrix of the genomic distances

    Inter chromosomal interactions are set to -1

    Parameters
    ----------
    lengths : ndarray (n, )
        lengths of the chromosomes

    Returns
    -------
    dis: ndarray (n, n), dtype: int
        returns the genomic distance matrix, with -1 for inter chromosomal
        interactions
    """
    inter_mask = get_inter_mask(lengths)
    n = lengths.sum()

    dis = np.concatenate([np.concatenate([np.arange(i, 0, -1),
                                          np.arange(n - i)])
                          for i in range(n)])
    dis = dis.reshape((n, n))
    dis[inter_mask] = -1

    return dis.astype(int)


def extract_sub_contact_map(counts, lengths, chromosomes):
    """
    Extract the contact map associated to a given list of chromosome

    Parameters
    ----------
    counts : ndarray (n, n)

    lengths : ndarray (L, )

    chromosomes : list of ids

    Returns
    -------
    sub_counts, sub_lengths : (ndarray, ndarray)

    Examples
    --------

    >>> from iced import datasets
    >>> from iced.utils import extract_sub_contact_map
    >>> counts, lengths = datasets.load_sample_yeast()
    >>> scounts, slengths = extract_sub_contact_map(counts, lengths, [0, 2])
    >>> print len(counts), len(scounts)
    ... # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    350 56
    """
    chromosomes = np.array(chromosomes)
    if chromosomes.max() >= len(lengths):
        raise ValueError(
            "The chromosomes provided are not compatible with the "
            "lengths array. Possible values are"
            " %s" % " ".join("%s" % i for i in np.arange(len(lengths))))
    if lengths.sum() != counts.shape[0]:
        raise ValueError(
            "The lengths provided is incompatible with the counts matrix"
            "shape. The total lengths is %d while the contact count matrix "
            "is %d" % (lengths.sum(), counts.shape[0]))

    is_symetric_or_tri(counts)
    chromosomes.sort()

    new_lengths = lengths[chromosomes]
    new_counts = np.zeros((new_lengths.sum(), new_lengths.sum()))
    begin1, end1 = 0, 0
    for i, l1 in enumerate(lengths):
        end1 += l1
        if i not in chromosomes:
            begin1 = end1
            continue
        # Find position of this pair of chromosome in the matrix
        new_num_chrom = (chromosomes == i).argmax()
        if new_num_chrom == 0:
            new_begin1 = 0
        else:
            new_begin1 = new_lengths.cumsum()[new_num_chrom - 1]
        new_end1 = new_lengths.cumsum()[new_num_chrom]

        begin2, end2 = 0, 0
        for j, l2 in enumerate(lengths):
            end2 += l2
            if j not in chromosomes:
                begin2 = end2
                continue
            # Find position of this pair of chromosome in the matrix
            new_num_chrom = (chromosomes == j).argmax()
            if new_num_chrom == 0:
                new_begin2 = 0
            else:
                new_begin2 = new_lengths.cumsum()[new_num_chrom - 1]
            new_end2 = new_lengths.cumsum()[new_num_chrom]
            new_counts[new_begin1:new_end1,
                       new_begin2:new_end2] = counts[begin1:end1, begin2:end2]
            begin2 = end2

        begin1 = end1

    return new_counts, new_lengths


def undersample_per_chr(X, lengths):
    """
    Undersample matrix to chromosomes

    Undersample the matrix ununiformaly per chromosome lengths.

    Parameters
    ----------
    X : ndarray (n, n)
        The matrix to undersample

    lengths : ndarray (L, )
        Lengths of the chromosomes

    Returns
    -------
    undersampled_X : ndarray (L, L)
        `X` undersampled per chromosome
    """
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


def downsample_resolution(counts, lengths, factor=2):
    """
    Downsamples the resolution of a matrix

    Parameters
    ----------
    counts : ndarray (N, N)
        contact counts matrix to downsample

    lengths : ndarray (L, )
        chromosomes lengths

    coef : int, optionnal, default: 2
        downsample resolution of the counts matrix by `coef`

    Returns
    -------
    target_counts, target_lengths : ndarray
    """
    if factor == 1:
        return counts, lengths
    # FIXME there is probably a better way to do this
    target_lengths = np.ceil(lengths.astype(float) / factor).astype(int)
    target_counts = np.zeros((target_lengths.sum(),
                              target_lengths.sum()))
    begin_i, end_i = 0, 0
    target_begin_i, target_end_i = 0, 0
    for i, length_i in enumerate(lengths):
        end_i += length_i
        target_end_i += target_lengths[i]
        begin_j, end_j = 0, 0
        target_begin_j, target_end_j = 0, 0
        for j, length_j in enumerate(lengths):
            end_j += length_j
            target_end_j += target_lengths[j]

            sub_counts = counts[begin_i:end_i, begin_j:end_j]
            sub_target_counts = target_counts[target_begin_i:target_end_i,
                                              target_begin_j:target_end_j]
            d = np.zeros(sub_target_counts.shape)
            for start in range(factor):
                s = sub_counts[start::factor, start::factor]
                d[:s.shape[0], :s.shape[1]] += np.invert(np.isnan(s))
                s[np.isnan(s)] = 0
                sub_target_counts[:s.shape[0], :s.shape[1]] += s
            sub_target_counts /= d

            begin_j = end_j
            target_begin_j = target_end_j
        begin_i = end_i
        target_begin_i = target_end_i
    return target_counts, target_lengths


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
