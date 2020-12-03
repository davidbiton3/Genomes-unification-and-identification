from package import functions


# this function will perform the algorithm stating.
# didnt return a value, just run it.
def alg(data_input, n_threads, completeness_threshold, contamination_threshold,
        genome_group_size, ani_threshold, linkage_threshold):

    # state - 1
    # functions.check_m(data_input, n_threads)
    
    # state - 2
    genome_list = functions.init_genomes(data_input, completeness_threshold, contamination_threshold)

    # state - 3
    genome_list = functions.sorted_genome(genome_list)

    # state - 4
    clusters_list = functions.cluster_the_genomes(genome_list, genome_group_size, ani_threshold, n_threads)

    # state - 5
    clusters_list = functions.calc_and_maping(clusters_list, n_threads, ani_threshold, linkage_threshold)

    # state - 6
    representatives_list = functions.select_reperesentive(clusters_list)
    
    # state - 7

    functions.gtdbtk(clusters_list, representatives_list)
    
    # state - 8 
    functions.create_final_file(clusters_list)
