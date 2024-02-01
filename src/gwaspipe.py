from pathlib import Path

import numpy as np
import click
import gwaslab as gl

from configuring import ConfigurationManager
from utils import __appname__, logger


class SumstatsManager:
    def __init__(self, input_path, input_format, input_separator):
        self.mysumstats = gl.Sumstats(input_path, fmt=input_format, sep=input_separator)


@click.command()
@click.option("-c", "--config_file", required=True, help="Configuration file path")
@click.option("-i", "--input_file", required=True, help="Input file path")
@click.option(
    "-f",
    "--input_file_format",
    required=True,
    type=click.Choice(["regenie", "fastgwa", "ldsc", "fuma"], case_sensitive=False),
    help="Input file format",
)
@click.option("-s", "--input_file_separator", default="\t", help="Input file separator")
@click.option("-q", "--quiet", default=False, is_flag=True, help="Set log verbosity")
def main(config_file, input_file, input_file_format, input_file_separator, quiet):
    cm = ConfigurationManager(config_file=config_file)
    log_file = cm.log_file_path

    if quiet:
        logger.add(log_file, level="INFO", retention="30 days")
    else:
        logger.add(log_file, level="DEBUG", retention="30 days")
    logger.info("{} started".format(__appname__.capitalize()))

    input_file_path = Path(input_file)
    input_file_name = input_file_path.name
    input_file_stem = input_file_path.stem

    if input_file_path.exists():
        sm = SumstatsManager(
            input_file_path.as_posix(), input_file_format, input_file_separator
        )
    else:
        msg = f"{input_file_path} input file not found"
        exit(msg)

    for step in cm.run_sequence:
        params, gl_params = cm.step(step)
        run  = params['run']
        workspace = params['workspace']
        workspace_path = Path(cm.root_path,workspace)
        workspace_path.mkdir(parents=True, exist_ok=True)

        if run and step == "basic_check":
            logger.info(f"Started {step} step")
            sm.mysumstats.basic_check(**gl_params)
            logger.info(f"Finished {step} step")
        elif run and step == "infer_build":
            logger.info(f"Started {step} step")
            sm.mysumstats.infer_build()
            genome_build = sm.mysumstats.meta["gwaslab"]["genome_build"]
            text = f"\nInferred genome build: {genome_build}\n"
            # with open(report_if_file_path, "a") as fp:
            #     fp.write(text)
            logger.info(f"Finished {step} step")
        elif run and step == "fill_data":
            logger.info(f"Started {step} step")
            sm.mysumstats.fill_data(**gl_params)
            logger.info(f"Finished {step} step")
        elif run and step == "harmonize":
            logger.info(f"Started {step} step")
            sm.mysumstats.harmonize(**gl_params)
            sm.mysumstats.flip_allele_status()
            logger.info(f"Finished {step} step")
        elif run and step =="liftover":
            logger.info(f"Started {step} step")
            sm.mysumstats.liftover(**gl_params)
            logger.info(f"Finished {step} step")
        elif run and step == "report_summary":
            logger.info(f"Started {step} step")
            header = f"Summary:"
            summary = sm.mysumstats.summary().to_string()
            # with open(report_if_file_path, "a") as fp:
            #     fp.write(header)
            #     fp.write(summary)
            logger.info(f"Finished {step} step")
        elif run and step == "report_inflation_factors":
            logger.info(f"Started {step} step")
            df = sm.mysumstats.data
            CHISQ = df.Z**2
            max_chisq = str(round(CHISQ.max(), 3))
            mean_chisq = str(round(CHISQ.mean(), 3))
            lambda_GC = str(round(CHISQ.median() / 0.4549, 3))

            report_if_file_path = Path(
                workspace_path,
                cm.report_if_filename.replace("placeholder", input_file_name)
            )
            with open(report_if_file_path, "w") as fp:
                fp.write("input_file_path\tlambda_GC\tmean_chisq\tmax_chisq\n")
                fp.write(f"{input_file_path}\t{lambda_GC}\t{mean_chisq}\t{max_chisq}")
            logger.info(f"Finished {step} step")
        elif run and step in ["write_ldsc", "write_metal", "write_vcf"]:
            logger.info(f"Started {step} step")
            output_path = str(Path(workspace_path, input_file_stem))
            sm.mysumstats.to_format(output_path, **gl_params)
            logger.info(f"Finished {step} step")
        elif run and step == 'write_same_input_format':
            logger.info(f"Started {step} step")
            output_path = str(Path(workspace_path, input_file_stem))
            sm.mysumstats.to_format(output_path, fmt=input_file_format, **gl_params)
            logger.info(f"Finished {step} step")
        elif run and step == 'qq_manhattan_plots':
            logger.info(f"Started {step} step")
            output_path = str(Path(cm.save_path, input_file_path.stem))
            cut = round(-np.log10(gl_params['sig_level'])) + params['dist']
            sm.mysumstats.plot_mqq(cut = cut, save = output_path, **gl_params)
            logger.info(f"Finished {step} step")
        else:
            logger.info(f"Skipping {step} step")


if __name__ == "__main__":
    main()
