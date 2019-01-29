from .fastio_ import loadtxt, savetxt
try:
    from ._io_pandas import load_counts, load_lengths
    from ._io_pandas import write_counts
except ImportError:
    from ._io_else import load_lengths, load_counts
    from ._io_else import write_counts


def write_lengths(filename, lengths, resolution=1):
    """
    Write lengths as bed file
    """
    chromosomes = ["Chr%02d" % (i + 1) for i in range(len(lengths))]
    j = 0
    with open(filename, "w") as bed_file:
        for chrid, l in enumerate(lengths):
            for i in range(l):
                bed_file.write(
                    "%s\t%d\t%d\t%d\n" % (chromosomes[chrid],
                                          i * resolution + 1,
                                          (i + 1) * resolution,
                                          j))
                j += 1
