from .cytobands import hg19_cytobands, hg38_cytobands
from .fragile_sites import hg19_fragile_sites, hg38_fragile_sites
from .gaps import hg19_gaps, hg38_gaps
from .genomes import hg19_chr_lengths, hg38_chr_lengths

class Assembly:
    """
    A class to represent a genomic assembly.

    Attributes
    ----------
    name : str
        The name of the assembly.
    chr_lens : dict
        The lengths of the chromosomes.
    chr_x : str
        The name of the X chromosome.
    chr_y : str
        The name of the Y chromosome.
    cytobands : list
        The cytobands of the chromosomes.
    gaps : list
        The gaps in the chromosomes.
    fragile_sites : list
        The common fragile sites of the chromosomes.
    """

    def __init__(self, name, chr_lens, x_name = "chrX", y_name = "chrY", cytobands=None, gaps=None, fragile_sites=None):
        self.name = name
        self.chr_names = list(chr_lens.keys())
        self.sex_names = [x_name, y_name]
        self.aut_names =  [item for item in chr_lens if item not in self.sex_names]
        self.chr_x = x_name
        self.chr_y = y_name
        self.chr_lens = chr_lens
        self.chr_starts = {}
        i = 0
        self.chr_starts = {}
        for k, v in chr_lens.items():
            self.chr_starts[k] = i
            i += v
        self.gen_len = sum(chr_lens.values())
        self.aut_len = self.gen_len - chr_lens[x_name] - chr_lens[y_name]
        self.cytobands = cytobands
        self.gaps = gaps
        self.fragile_sites = fragile_sites


hg19 = Assembly(
    name="hg19",
    chr_lens=hg19_chr_lengths,
    cytobands=hg19_cytobands,
    gaps=hg19_gaps,
    fragile_sites=hg19_fragile_sites
)
"""
An instance of the Assembly class representing the hg19 genomic assembly.
"""

hg38 = Assembly(
    name="hg38",
    chr_lens=hg38_chr_lengths,
    cytobands=hg38_cytobands,
    gaps=hg38_gaps,
    fragile_sites=hg38_fragile_sites
)
"""
An instance of the Assembly class representing the hg38 genomic assembly.
"""


def get_assembly(assembly_id):
    """
    Retrieve an Assembly instance by its ID.

    Parameters
    ----------
    assembly_id : str
        The ID of the assembly to retrieve. Valid values are "hg19" and "hg38".

    Returns
    -------
    Assembly
        The Assembly instance corresponding to the given ID.

    Raises
    ------
    ValueError
        If the assembly_id is not "hg19" or "hg38".
    """
    if assembly_id == "hg19":
        return hg19
    elif assembly_id == "hg38":
        return hg38
    else:
        raise ValueError(f"Assembly {assembly_id} not found")