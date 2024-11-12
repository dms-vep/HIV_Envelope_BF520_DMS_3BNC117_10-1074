"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""
with open(config["antibody_escape_config"]) as f:
    antibody_escape_config = yaml.safe_load(f)

antibody_escape_pdbs = antibody_escape_config["antibody_escape_PDBs"]

Env_chains_by_pdb = antibody_escape_config["env_chains_by_PDB"]

chains_to_exclude = antibody_escape_config["chains_to_exclude"]

antibody_gv = [x for x in antibody_escape_config['avg_antibody_escape']]

rule spatial_distances:
    """Get spatial distances from PDB."""
    input: 
        pdb="data/env_trimer_6UDJ.pdb",
    output:
        csv="results/spatial_distances/6udj.csv",
    params:
        target_chains=["A", "B", "C"],
    log:
        log="results/logs/spatial_distances.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    script:
        "scripts/spatial_distances.py"

rule color_PDB_structures:
    """Assign b factor values to PDB structures based on escape"""
    input: 
        average_escape_model = rules.avg_escape.output.effect_csv,
        input_pdb_file = lambda wc: os.path.join(
            antibody_escape_config["PDB_structures"], 
            antibody_escape_pdbs[wc.antibody]
        ),
    output:
        output_pdb_file_name = os.path.join("results/escape_PDBs/", "{assay}", "{antibody}" + ".pdb"),
    params: 
        env_chains = lambda wc: Env_chains_by_pdb[antibody_escape_pdbs[wc.antibody]],
        output_file_sub_name = os.path.join("results/escape_PDBs/", "{assay}", "{antibody}"),
    log:
        os.path.join("results/logs/", "antibody_escape_pdbs_{assay}_{antibody}.txt"),
    script:
        "scripts/color_pdb_structures.py"

rule func_effects_reference_sites: 
    """Add reference sites to functional effects dataframe"""
    input: 
        functional_data="results/func_effects/averages/TZM-bl_entry_func_effects.csv",
    output:
        output_file = os.path.join("results/func_effects/averages/TZM-bl_entry_func_effects_references_sites.csv"),
    log:
        os.path.join("results/logs/functional_effects_reference_sites.txt"),
    script:
        "scripts/func_effects_ref_sites.py"

rule format_dms_viz:
    """Format the data for input into dms-viz"""
    input:
        average_escape_model = rules.avg_escape.output.effect_csv,
        input_pdb_file = lambda wc: os.path.join(
            antibody_escape_config["PDB_structures"], 
            antibody_escape_pdbs[wc.antibody]
        ),
        site_map="data/site_numbering_map.csv",
        functional_data="results/func_effects/averages/TZM-bl_entry_func_effects_references_sites.csv",
    output: 
        output_json_file_name = os.path.join("results/dms-viz/", "{assay}", "{antibody}" + ".json"),
    params:
        env_chains = lambda wc: Env_chains_by_pdb[antibody_escape_pdbs[wc.antibody]],
        exclude_chains = lambda wc: chains_to_exclude[wc.antibody],
        name="{antibody}",
    log:
        os.path.join("results/logs/", "dms-viz_file_{assay}_{antibody}.txt"),
    conda:
        "dms-viz.yml"
    shell:
        """
        configure-dms-viz format \
            --name {params.name} \
            --input {input.average_escape_model} \
            --metric  "escape_mean" \
            --metric-name "Escape" \
            --condition "epitope" \
            --condition-name "Epitope" \
            --exclude-amino-acids "*, -" \
            --join-data {input.functional_data} \
            --structure {input.input_pdb_file} \
            --sitemap {input.site_map} \
            --output {output.output_json_file_name} \
            --included-chains "{params.env_chains}" \
            --excluded-chains "{params.exclude_chains}" \
            --tooltip-cols "{{'times_seen': '# Obsv', 'effect': 'Func. Eff.'}}" \
            --filter-cols "{{'times_seen': 'Times Seen', 'effect': 'Functional Effect'}}" \
            --filter-limits "{{'times_seen': [0, 3, 10], 'effect': [-5, -4, 0.1]}}" \
            &> {log}
        """

rule dms_viz_join:
    """Join the dms-viz files into one file"""
    input:
        expand("results/dms-viz/{file}", file=['antibody_escape/10-1074.json', 'antibody_escape/3BNC117.json'])
    output: 
        output_json_file_name = os.path.join("results/dms-viz/TRO11_V3Abs.json"),
    params: 
        file_list = ', '.join(expand("results/dms-viz/{file}", file=['antibody_escape/10-1074.json', 'antibody_escape/3BNC117.json']))
    log:
        os.path.join("results/logs/", "dms-viz_file_join.txt"),
    conda:
        "dms-viz.yml"
    shell:
        """
        configure-dms-viz join \
            --input "{params.file_list}" \
            --output {output.output_json_file_name} \
            &> {log}
        """

rule plot_average_func_scores:
    """Plot average functional scores for variants"""
    input: 
        "results/func_effects/averages/TZM-bl_entry_func_effects.csv",
        nb="notebooks/average_func_score_distributions.ipynb",
    output:
        nb="results/notebooks/average_func_score_distributions.ipynb"
    log:
        "results/logs/average_func_score_distributions.txt"
    shell:
        """
        papermill {input.nb} {output.nb} \
            &> {log}
        """

rule logo_plots:
    """Make logo plots for antibody escape."""
    input:
        "results/antibody_escape/averages/10-1074_mut_effect.csv",
        "results/antibody_escape/averages/3BNC117_mut_effect.csv",
        "data/site_numbering_map.csv",
        nb="notebooks/logoplots.ipynb",
    output:
        "results/logoplots/10-1074_logoplot.svg",
        "results/logoplots/3BNC117_logoplot.svg",
        "results/logoplots/scalebar_figure.svg",
        nb="results/notebooks/logoplots.ipynb",
    log:
        log="results/logs/logoplots.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        """
        papermill {input.nb} {output.nb} \
            &> {log}
        """

# Files (Jupyter notebooks, HTML plots, or CSVs) that you want included in
# the HTML docs should be added to the nested dict `docs`:
docs["Additional analysis-specific files"] = {
    "Reference to sequential site-numbering map": config["site_numbering_map"],
    "Average functional score plots": "results/notebooks/average_func_score_distributions.ipynb", 
}
for file in expand(rules.format_dms_viz.output.output_json_file_name, antibody=antibody_gv, assay='antibody_escape'):
    other_target_files.append(file)
for file in expand(rules.color_PDB_structures.output.output_pdb_file_name, antibody=antibody_gv, assay='antibody_escape'):
    other_target_files.append(file)
other_target_files.append("results/notebooks/logoplots.ipynb") 
