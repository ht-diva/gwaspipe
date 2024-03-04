from pathlib import Path

import numpy as np
import click
import gwaslab as gl
import cistrans_tagger
from cistrans_tagger import *

from configuring import ConfigurationManager
from utils import __appname__, logger


class SumstatsManager:
    def __init__(self, input_path, input_format, input_separator, formatbook_path):
        if formatbook_path.exists():
            gl.options.set_option("formatbook", str(formatbook_path))
        if input_format == 'pickle':
            self.mysumstats = gl.load_pickle(input_path)
        else:
            self.mysumstats = gl.Sumstats(input_path, fmt=input_format,
                                          sep=input_separator)


@click.command()
@click.option("-c", "--config_file", required=True, help="Configuration file path")
@click.option("-i", "--input_file", required=True, help="Input file path")
@click.option(
    "-f",
    "--input_file_format",
    required=True,
    type=click.Choice(["regenie", "fastgwa", "ldsc", "fuma", "pickle"], case_sensitive=False),
    help="Input file format",
)
@click.option("-o", "--output", help="Path where results should be saved")
@click.option("-s", "--input_file_separator", default="\t", help="Input file separator")
@click.option("-q", "--quiet", default=False, is_flag=True, help="Set log verbosity")
def main(config_file, input_file, input_file_format, input_file_separator, output, quiet):
    cm = ConfigurationManager(config_file=config_file, root_path=output)
    log_file = cm.log_file_path

    if quiet:
        logger.add(log_file, level="INFO", retention="30 days")
    else:
        logger.add(log_file, level="DEBUG", retention="30 days")
    logger.info("{} started".format(__appname__.capitalize()))

    input_file_path = Path(input_file)
    input_file_name = input_file_path.name

    mask, sep = cm.filename_settings
    if mask:
        string_list = np.array(input_file_name.split(sep))
        input_file_stem = sep.join(string_list[mask].tolist())
    else:
        input_file_stem = input_file_path.stem

    formatbook_file_path = Path(cm.formatbook_path)

    if input_file_path.exists():
        sm = SumstatsManager(
            input_file_path.as_posix(), input_file_format,
            input_file_separator,
            formatbook_file_path
        )
    else:
        msg = f"{input_file_path} input file not found"
        exit(msg)

    for step in cm.run_sequence:
        params, gl_params = cm.step(step)
        run = params.get('run', False)
        ws = params.get('workspace', 'default')
        ws_subfolder = params.get('workspace_subfolder', False)
        if ws_subfolder:
            workspace_path = Path(cm.root_path, ws, input_file_stem)
        else:
            workspace_path = Path(cm.root_path, ws)
        workspace_path.mkdir(parents=True, exist_ok=True)

        if run:
            logger.info(f"Started {step} step")
            if step == "basic_check":
                sm.mysumstats.basic_check(**gl_params)
            elif step == "infer_build":
                sm.mysumstats.infer_build()
                genome_build = sm.mysumstats.meta["gwaslab"]["genome_build"]
                text = f"\nInferred genome build: {genome_build}\n"
                # with open(report_if_file_path, "a") as fp:
                #     fp.write(text)
            elif step == "fill_data":
                sm.mysumstats.fill_data(**gl_params)
            elif step == "harmonize":
                sm.mysumstats.harmonize(**gl_params)
                sm.mysumstats.flip_allele_stats()
            elif step =="liftover":
                sm.mysumstats.liftover(**gl_params)
            elif step == "report_summary":
                header = f"Summary:"
                summary = sm.mysumstats.summary().to_string()
                print(summary)
                # with open(report_if_file_path, "a") as fp:
                #     fp.write(header)
                #     fp.write(summary)
            elif step == "report_min_pvalues":
                nrows = params['nrows']
                df = sm.mysumstats.data.nsmallest(nrows, 'P', keep="all")
                df.drop(columns=['STATUS'], inplace=True)
                output_path = str(
                    Path(workspace_path, '.'.join([input_file_stem, 'nsmallest.tsv'])))
                df.to_csv(output_path, sep='\t')
            elif step == "report_inflation_factors":
                df = sm.mysumstats.data
                CHISQ = df.Z**2
                max_chisq = str(round(CHISQ.max(), 3))
                mean_chisq = str(round(CHISQ.mean(), 3))
                lambda_GC = str(round(CHISQ.median() / 0.4549, 3))

                report_if_file_path = Path(
                    workspace_path,
                    cm.report_if_filename.replace("placeholder", input_file_stem)
                )
                with open(report_if_file_path, "w") as fp:
                    fp.write("input_file_path\tlambda_GC\tmean_chisq\tmax_chisq\n")
                    fp.write(f"{input_file_path}\t{lambda_GC}\t{mean_chisq}\t{max_chisq}\n")
            elif step == 'write_pickle':
                output_path = str(
                    Path(workspace_path, '.'.join([input_file_stem, 'pkl'])))
                gl.dump_pickle(sm.mysumstats, output_path, overwrite=params['overwrite'])
            elif step in ["write_regenie", "write_ldsc", "write_metal", "write_tsv", "write_fastgwa"]:
                output_path = str(Path(workspace_path, input_file_stem))
                sm.mysumstats.to_format(output_path, **gl_params)
            elif step == 'write_vcf':
                study_name = input_file_stem
                sm.mysumstats.meta["gwaslab"][
                    "study_name"] = study_name
                sm.mysumstats.infer_build()
                output_path = str(Path(workspace_path, input_file_stem))
                sm.mysumstats.to_format(output_path, **gl_params)
            elif step == 'write_same_input_format':
                output_path = str(Path(workspace_path, input_file_stem))
                sm.mysumstats.to_format(output_path, fmt=input_file_format, **gl_params)
            elif step == 'qq_manhattan_plots':
                output_path = str(
                    Path(workspace_path, '.'.join([input_file_stem, 'png'])))
                cut = round(-np.log10(gl_params['sig_level'])) + params['dist']
                sm.mysumstats.plot_mqq(cut=cut, save=output_path, **gl_params)
            elif step == 'cistrans_annotation':
                cistrans_gene_tagger(params, workspace_path)
            logger.info(f"Finished {step} step")
        else:
            logger.info(f"Skipping {step} step")


if __name__ == "__main__":
    main()
