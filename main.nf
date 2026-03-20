/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { SCRANK } from "./modules/local/scrank/main.nf"
include { DOWNSAMPLE } from "./modules/local/downsample_and_split/main.nf"
include { MERGE } from "./modules/local/merge_and_downstream/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    obj = file(params.obj)
    column = params.column
    species = params.species
    n_cells = params.n_cells
    n_cores = params.n_cores
    target = file(params.target)

    DOWNSAMPLE( obj, target, column, species, n_cells )

    DOWNSAMPLE.out.scrank_obj
    .flatten()
    .set { sc_obj }


    SCRANK( sc_obj, n_cores )

    SCRANK.out.rank_obj
    .collect()
    .set { rank_cells  }

    MERGE( obj, target, species, column, rank_cells ) 

}
