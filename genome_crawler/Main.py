"""
File:           Main.py (for Genome Compilation Module)

Author:         Alexander S. Adranly

Description:    Main.py is the execution script for creating and running the Genome Compilation Module

"""
from crawler.GCM import *

def run():

    genome_compilation_module = GCM()

    genome_compilation_module.crawl_genomes()
    #["Abiotrophia_defectiva","Acetivibrio_cellulolyticus","Acetivibrio_ethanolgignens", "Acaryochloris_sp._CCMEE_5410", "Acaryochloris_sp._CCMEE_5410","zeta_proteobacterium_SCGC_AB-133-C04"]

if __name__ == "__main__":
    run()
