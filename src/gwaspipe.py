from pathlib import Path

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
    type=click.Choice(["regenie", "vcf"], case_sensitive=False),
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

    report_file_path = Path(
        cm.save_path, cm.report_filename.replace("placeholder", input_file_name)
    )
    with open(report_file_path, "w") as fp:
        fp.write(f"input file path: {input_file_path}")

    if input_file_path.exists():
        sm = SumstatsManager(
            input_file_path.as_posix(), input_file_format, input_file_separator
        )
    else:
        msg = f"{input_file_path} input file not found"
        exit(msg)

    for step in cm.run_sequence:
        run, params = cm.step(step)
        if run and step == "basic_check":
            logger.info(f"Started {step} step")
            sm.mysumstats.basic_check(**params)
            logger.info(f"Ended {step} step")
        elif run and step == "infer_build":
            logger.info(f"Started {step} step")
            sm.mysumstats.infer_build()
            genome_build = sm.mysumstats.meta["gwaslab"]["genome_build"]
            text = f"\nInferred genome build: {genome_build}\n"
            with open(report_file_path, "a") as fp:
                fp.write(text)
            logger.info(f"Ended {step} step")
        elif run and step == "fill_data":
            logger.info(f"Started {step} step")
            sm.mysumstats.fill_data(**params)
            logger.info(f"Ended {step} step")
        elif run and step == "harmonize":
            logger.info(f"Started {step} step")
            sm.mysumstats.harmonize(**params)
            sm.mysumstats.flip_allele_status()
            logger.info(f"Ended {step} step")
        elif run and step == "report_summary":
            logger.info(f"Started {step} step")
            header = f"Summary:"
            summary = sm.mysumstats.summary().to_string()
            with open(report_file_path, "a") as fp:
                fp.write(header)
                fp.write(summary)
            logger.info(f"Ended {step} step")
        elif run and step == "report_inflation_factors":
            logger.info(f"Started {step} step")
            header = f"\nInflation factors:"
            df = sm.mysumstats.data
            CHISQ = df.Z**2
            mean_chisq = CHISQ.mean()
            with open(report_file_path, "a") as fp:
                fp.write(header)
                fp.write("\nLambda GC = " + str(round(CHISQ.median() / 0.4549, 3)))
                fp.write("\nMean chi^2 = " + str(round(mean_chisq, 3)))
                if mean_chisq < 1.02:
                    fp.write(" WARNING: mean chi^2 may be too small.")
                fp.write("\nMax chi^2 = " + str(round(CHISQ.max(), 3)))
            logger.info(f"Ended {step} step")
        elif run and step in ["write_ldsc", "write_metal", "write_vcf"]:
            logger.info(f"Started {step} step")
            output_path = str(Path(cm.save_path, input_file_path.stem))
            sm.mysumstats.to_format(output_path, **params)
            logger.info(f"Ended {step} step")
        else:
            logger.info(f"Skipping {step} step")


if __name__ == "__main__":
    main()
