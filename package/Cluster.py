"""
A class for defining an Cluster.
The c'tor receives the following arguments:
    genome: genome that create the cluster, temporary representative genome.
"""


class Cluster:
    counter = 0  # static class Attribute

    # Initializer / Instance Attributes
    def __init__(self, genome):
        self.id = Cluster.counter
        Cluster.counter += 1
        self.representative_genome = genome
        self.genome_list = []
        self.genome_list.append(genome)
        self.taxonomy = ""

    # this function will select a representative to the cluster between all the genomes in the cluster.
    # all the genomes with the highest link score will be a candidates to be a representative.
    # from the candidates will chose one genome that  his median is the highest of all candidates.
    # return the genome that selected.
    def select_reperesentive(self):
        list_max_link_genome_representative = []
        max_link = 0
        for genome in self.genome_list:
            if genome.linked > max_link:
                max_link = genome.linked
                list_max_link_genome_representative.clear()
                list_max_link_genome_representative.append(genome)
            elif genome.linked == max_link:
                list_max_link_genome_representative.append(genome)

        ani_score_median_max = 0
        for genome_repr in list_max_link_genome_representative:
            if len(genome_repr.list_path_genome_linked):
                ani_score_median = genome_repr.calc_median()
                if ani_score_median > ani_score_median_max:
                    ani_score_median_max = ani_score_median
                    self.representative_genome = genome_repr
