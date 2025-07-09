process addParams {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/params_added", mode: 'copy'

    input:
    //val pipeline_run_dir_name
    tuple val(run_name), val(params_file), val(ref_obs), val(pipeline_results)


    output:
    path "*_predicted_meta_combined.tsv", emit:  predicted_meta_combined

    script:
    ref_keys = params.ref_keys.join(' ')
    """
    python $projectDir/bin/add_params.py --run_name ${run_name}  \\
				--pipeline_results ${pipeline_results} \\
				--params_file ${params_file} \\
				--subsample ${params.subsample} \\
                --ref_keys ${ref_keys} \\
                --mapping_file ${params.mapping_file}
						
    """

}

process aggregateResults {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/aggregated_results", mode: 'copy'

    input:
    path predicted_meta_combined

    output:
   // path "f1_results_all_pipeline_runs.tsv", emit: f1_results_aggregated
	path "**predicted_meta_combined.tsv", emit: aggregated_results
	path "multiqc/", emit: multiqc_dir

    script:
    """
    python $projectDir/bin/aggregate_results.py --pipeline_results ${predicted_meta_combined} \\
                                                --ref_keys ${params.ref_keys.join(' ')}
    """

}

process model_correct {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/logreg_results", mode: 'copy'

    input:
    path predicted_meta_combined

    output:
    path "**png", emit: plots
    //path "**_logreg_fits.tsv", emit: logreg_results

    script:
    """
    python $projectDir/bin/model_correct.py --predicted_meta ${predicted_meta_combined} \\
                                             --mapping_file ${params.mapping_file} \\
                                             --ref_keys ${params.ref_keys.join(' ')}
    """
}

process split_by_label {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/split_by_label", mode: 'copy'

    input:
    path predicted_meta_combined

    output:
    path "**predicted_meta_subset.tsv", emit: predicted_meta_labels
    
    script:
    """
    python $projectDir/bin/split_by_label.py --predicted_meta ${predicted_meta_combined} \\
                                                --ref_keys ${params.ref_keys.join(' ')}
    """
}

process model_per_celltype {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/cell_type_models/${cell_type}", mode: 'copy'

    input:
    tuple val(cell_type), path(predicted_meta_subset)

    output:
    path "**png"
    path "**_coefficients.tsv", emit: coefficients

    script:
    """
    python $projectDir/bin/model_correct.py --predicted_meta ${predicted_meta_subset} \\
                                                  --mapping_file ${params.mapping_file} \\
                                                  --ref_keys ${params.ref_keys.join(' ')} \\
                                                    --cell_type ${cell_type}
    """
}

process run_multiqc {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
        path (qc_dir)

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc ${qc_dir} -d --config ${params.multiqc_config}
    """
}

workflow {

    Channel
    .fromPath("${params.results}/*", type: 'dir') // Get all subdirectories
    .map { pipeline_run_dir ->
        def pipeline_run_dirname = pipeline_run_dir.getName().toString()
        def params_file = "${pipeline_run_dir}/params.yaml"
        def ref_obs = "${pipeline_run_dir}/refs/"
        // Collect 'f1_results' directories
        def pipeline_results = []
        pipeline_run_dir.eachDirRecurse { dir ->
                if (dir.getName() == 'scvi' || dir.getName() == 'seurat') {
                        def dir_path = dir.toString()
                        pipeline_results << dir_path
            }
        }

        [pipeline_run_dirname, params_file, ref_obs, pipeline_results.flatten().join(' ')] // Return collected results for this pipeline_run_dir
    }

    .set { all_pipeline_results }

    // add parameters to files  addParams(all_pipeline_results) 
    addParams(all_pipeline_results)
	aggregateResults(addParams.out.predicted_meta_combined.flatten().toList())


	predicted_meta_combined = aggregateResults.out.aggregated_results

    model_correct(predicted_meta_combined)

    split_by_label(predicted_meta_combined)

    split_by_label.out.predicted_meta_labels.flatMap().set { predicted_meta_labels }

    predicted_meta_labels.map { predicted_meta_subset ->
        def cell_type = predicted_meta_subset.getName().replace("_predicted_meta_subset.tsv", "")
        [cell_type, predicted_meta_subset]
    }.set { cell_type_subset }

    
   model_per_celltype(cell_type_subset) 


    aggregateResults.out.multiqc_dir
    .set { qc_dir }

    run_multiqc(qc_dir)
}