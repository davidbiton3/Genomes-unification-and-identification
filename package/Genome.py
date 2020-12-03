import statistics

"""
A class for defining an Genome.
The c'tor receives the following arguments:
    name: genome identifier name from input file.
    path: genome fasta file path form input file.
    completeness: percent of completeness in the genome protein.
    contamination: percent of contamination in the genome protein.
"""


class Genome:

    # Initializer / Instance Attributes
    def __init__(self, name, path, completeness, contamination):
        self.name = name
        self.path = path
        self.completeness = completeness
        self.contamination = contamination
        self.quality = completeness - 5 * contamination
        self.linked = 0
        self.list_path_genome_linked = []  # is tuple (path,ani_score)

    # this function update 2 values - linked, list_path_genome_linked
    def update_linkage_list_path_genome_linked_is_low_linked(self, linked, list_path_genome_linked):
        self.linked = linked
        self.list_path_genome_linked = list_path_genome_linked

    # this function will calc the median ani score of a genome with all his friends.
    # return a median value all ani score with all his friends.
    def calc_median(self):
        if len(self.list_path_genome_linked) == 0:
            return 0

        temp, list_path_genome_linked = zip(*self.list_path_genome_linked)
        ani_score_median = statistics.median(list(list_path_genome_linked))
        return ani_score_median

    # this function will check if path exist in my list_path_genome_linked if yes delete
    def delete_path_from_list_path(self, path):
        index = -1
        for i, tuple_list in enumerate(self.list_path_genome_linked):
            path_list = tuple_list[0]
            if path_list == path:
                index = i
                break
        if index != -1:
            del self.list_path_genome_linked[index]
            return True
        return False
