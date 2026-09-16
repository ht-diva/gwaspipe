import gzip
import json
import os
import re
from datetime import UTC, datetime
from pathlib import Path

import click
import gwaslab as gl
import numpy as np
import polars as pl

from gwaspipe import __appname__, __version__, logger
from gwaspipe.configuring import ConfigurationManager
from gwaspipe.order_alleles import order_alleles as order_alleles_func


class AssemblyValidationError(ValueError):
    """Raised when an assembly declaration cannot be validated safely."""


_ASSEMBLY_ALIASES = {
    "19": "19",
    "37": "19",
    "b37": "19",
    "grch37": "19",
    "hg19": "19",
    "38": "38",
    "b38": "38",
    "grch38": "38",
    "hg38": "38",
}
_ASSEMBLY_LABELS = {"19": "GRCh37", "38": "GRCh38"}
_QC_STEPS = {
    "basic_check",
    "infer_build",
    "harmonize",
    "liftover",
    "canonicalize_effect_alleles",
    "sort_alphabetically",
    "filter_conflicting_snpids",
    "check_ambiguous_snps",
}
_GWASLAB_REMOVAL_PATTERN = re.compile(r"-Removed variants (.+): (\d+)$")
_GWASLAB_COUNT_FIRST_REMOVAL_PATTERN = re.compile(r"-Removed (\d+) variants (.+?)[.:]$")
_ORDER_ALLELES_FLIPPED_PATTERN = re.compile(r"-For Flipped match \((\d+) matches\)")
_EXTREME_P_THRESHOLD = 1e-300


def _basic_check_exclusion_counts(log_delta):
    """Translate GWASLab basic_check removal messages into stable reason codes."""
    exclusion_counts = {}
    position_outlier_count = 0
    for line in log_delta.splitlines():
        match = _GWASLAB_REMOVAL_PATTERN.search(line)
        if match is not None:
            reason, count_text = match.groups()
            count = int(count_text)
        else:
            count_first_match = _GWASLAB_COUNT_FIRST_REMOVAL_PATTERN.search(line)
            if count_first_match is None:
                continue
            count_text, reason = count_first_match.groups()
            count = int(count_text)
        if reason == "in total":
            continue
        if reason == "outliers":
            reason_code = "basic_check_position_out_of_bounds"
            position_outlier_count += count
        elif reason == "with bad positions":
            reason_code = "basic_check_missing_or_invalid_position"
            count = max(count - position_outlier_count, 0)
        elif reason == "with chromosome notations not in CHR list" or reason.startswith("with CHR < "):
            reason_code = "basic_check_invalid_chromosome"
        elif reason == "with NA alleles or alleles that contain bases other than A/C/T/G":
            reason_code = "basic_check_missing_or_invalid_allele"
        elif reason == "with same allele for EA and NEA":
            reason_code = "basic_check_identical_effect_and_non_effect_alleles"
        elif reason == "based on SNPID":
            reason_code = "basic_check_duplicate_snpid"
        elif reason == "based on rsID":
            reason_code = "basic_check_duplicate_rsid"
        elif reason == "based on CHR,POS,EA and NEA":
            reason_code = "basic_check_duplicate_variant_identity"
        elif reason == "multiallelic variants":
            reason_code = "basic_check_multiallelic_variant"
        elif reason.lower() == "with bad/na p":
            reason_code = "basic_check_invalid_or_missing_p_value"
        elif reason.startswith("with NA values in "):
            reason_code = "basic_check_missing_required_value"
        else:
            reason_code = "basic_check_gwaslab_" + re.sub(r"[^a-z0-9]+", "_", reason.lower()).strip("_")
        if count:
            exclusion_counts[reason_code] = exclusion_counts.get(reason_code, 0) + count
    return exclusion_counts


def _order_alleles_flipped_count(log_delta):
    """Return the number of allele pairs swapped by canonicalization."""
    return sum(int(match.group(1)) for match in _ORDER_ALLELES_FLIPPED_PATTERN.finditer(log_delta))


def _normalise_assembly(assembly):
    if assembly is None:
        return None
    return _ASSEMBLY_ALIASES.get(str(assembly).strip().lower())


def _hapmap3_match_counts(sumstats_data):
    """Count input rows matching the HapMap3 coordinates used by GWASLab."""
    from gwaslab.util.util_in_filter_value import _get_hapmap_df_polars

    missing_columns = {"CHR", "POS"}.difference(sumstats_data.columns)
    if missing_columns:
        raise AssemblyValidationError(f"Cannot infer genome assembly; missing columns: {sorted(missing_columns)}")

    positions = (
        pl.from_pandas(sumstats_data[["CHR", "POS"]])
        .filter(pl.col("CHR").is_not_null() & pl.col("POS").is_not_null())
        .with_columns(
            pl.col("CHR").cast(pl.Int64),
            pl.col("POS").cast(pl.Int64),
        )
    )
    return {
        build: positions.join(_get_hapmap_df_polars(build), on=["CHR", "POS"], how="semi").height
        for build in ("19", "38")
    }


def _assembly_audit_config(config):
    validation_config = config.get("assembly_validation", {})
    if validation_config is None:
        validation_config = {}
    if not isinstance(validation_config, dict):
        raise AssemblyValidationError("assembly_validation must be a mapping")

    min_matches = validation_config.get("min_hapmap3_matches", 10000)
    if not isinstance(min_matches, int) or min_matches < 1:
        raise AssemblyValidationError("assembly_validation.min_hapmap3_matches must be a positive integer")

    override_reason = validation_config.get("override_reason")
    allow_override = bool(validation_config.get("allow_override", False))
    if allow_override and not isinstance(override_reason, str):
        raise AssemblyValidationError("assembly_validation.override_reason is required when allow_override is true")

    return {
        "declared_input_assembly": config.get("genome_assembly"),
        "min_hapmap3_matches": min_matches,
        "allow_override": allow_override,
        "override_reason": override_reason,
        "reference_resources": config.get("reference_resources", {}),
    }


def _store_assembly_audit(mysumstats, audit):
    mysumstats.meta.setdefault("gwaspipe", {})["assembly_validation"] = audit


def validate_declared_assembly(mysumstats, config, infer_build_params=None, declared_assembly=None, scope="input"):
    """Run build inference and verify it against a declared assembly."""
    policy = _assembly_audit_config(config)
    declared_assembly = policy["declared_input_assembly"] if declared_assembly is None else declared_assembly
    declared_build = _normalise_assembly(declared_assembly)
    counts = _hapmap3_match_counts(mysumstats.data)
    if counts["19"] > counts["38"]:
        inferred_build = "19"
    elif counts["38"] > counts["19"]:
        inferred_build = "38"
    else:
        inferred_build = "Unknown"

    # Keep GWASLab's STATUS and metadata behaviour, while retaining the evidence
    # needed to audit the decision at the pipeline layer.
    mysumstats.infer_build(**(infer_build_params or {}))

    reasons = []
    if declared_build is None:
        reasons.append("missing or unsupported declared genome_assembly")
    if inferred_build == "Unknown":
        reasons.append("ambiguous inferred genome assembly")
    if max(counts.values()) < policy["min_hapmap3_matches"]:
        reasons.append(f"fewer than {policy['min_hapmap3_matches']} HapMap3 coordinate matches")
    if declared_build and inferred_build != "Unknown" and declared_build != inferred_build:
        reasons.append("declared and inferred genome assemblies disagree")

    overridden = bool(reasons) and policy["allow_override"]
    passed = not reasons or overridden
    audit = {
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "validation_scope": scope,
        "declared_assembly": declared_assembly,
        "declared_build": _ASSEMBLY_LABELS.get(declared_build),
        "inferred_assembly": _ASSEMBLY_LABELS.get(inferred_build, "Unknown"),
        "hapmap3_match_counts": {
            "GRCh37": counts["19"],
            "GRCh38": counts["38"],
        },
        "minimum_hapmap3_matches": policy["min_hapmap3_matches"],
        "inference_method": "HapMap3 CHR/POS coordinate overlap via GWASLab infer_build",
        "reference_resource_versions": policy["reference_resources"],
        "decision": "passed" if not reasons else "overridden" if overridden else "failed",
        "override_reason": policy["override_reason"] if overridden else None,
        "failure_reasons": reasons,
    }
    _store_assembly_audit(mysumstats, audit)
    if not passed:
        raise AssemblyValidationError("Assembly validation failed: " + "; ".join(reasons))
    return audit


def _require_validated_assembly(mysumstats):
    audit = mysumstats.meta.get("gwaspipe", {}).get("assembly_validation")
    if not audit or audit["decision"] not in {"passed", "overridden"}:
        raise AssemblyValidationError(
            "A successful infer_build validation is required before this reference-aware step."
        )
    return audit


def _record_qc_step(mysumstats, step, rows_before, exclusion_counts=None, metrics=None):
    """Record row counts, exclusion reasons, and non-exclusion QC metrics."""
    rows_after = len(mysumstats.data)
    exclusion_counts = exclusion_counts or {}
    exclusions = [
        {"reason_code": reason_code, "row_count": int(row_count)}
        for reason_code, row_count in exclusion_counts.items()
        if row_count
    ]
    accounted_exclusions = sum(exclusion["row_count"] for exclusion in exclusions)
    unaccounted_exclusions = max(rows_before - rows_after - accounted_exclusions, 0)
    if unaccounted_exclusions:
        exclusions.append(
            {
                "reason_code": f"unattributed_removal_by_{step}",
                "row_count": unaccounted_exclusions,
            }
        )

    audit = {
        "step": step,
        "rows_before": rows_before,
        "rows_after": rows_after,
        "rows_excluded": rows_before - rows_after,
        "exclusions": exclusions,
    }
    if metrics:
        audit["metrics"] = metrics
    gwaspipe_meta = mysumstats.meta.setdefault("gwaspipe", {})
    gwaspipe_meta.setdefault("qc_steps", []).append(audit)


def _write_run_provenance(output_path, mysumstats, source_path=None, steps=None):
    """Write a machine-readable sidecar for an output artifact."""
    audit = mysumstats.meta.get("gwaspipe", {}).get("assembly_validation")
    provenance_path = Path(f"{output_path}.provenance.json")
    provenance = {
        "gwaspipe_version": __version__,
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "source_path": str(Path(source_path).resolve()) if source_path is not None else None,
        "output_path": str(Path(output_path).resolve()),
        "assembly_validation": audit,
    }
    if steps is not None:
        provenance["steps"] = steps
    qc_steps = mysumstats.meta.get("gwaspipe", {}).get("qc_steps")
    if qc_steps:
        provenance["qc"] = {"steps": qc_steps}
    with provenance_path.open("w") as output_file:
        json.dump(provenance, output_file, indent=2, sort_keys=True)


class SumstatsManager:
    def _make_gwaslab_snpid(self):
        """Return SNPID in GWASLab format (CHR:POS:EA:NEA)"""
        df = self.mysumstats.data
        return df[["CHR", "POS", "EA", "NEA"]].astype(str).agg(":".join, axis=1)

    def __init__(self, input_path, input_format, input_separator, input_study, formatbook_path, pid, bcfliftover):
        if formatbook_path.exists():
            gl.options.set_option("formatbook", str(formatbook_path))
        if input_format == "pickle":
            self.mysumstats = gl.load_pickle(input_path)
        elif input_format == "vcf":
            self.mysumstats = gl.Sumstats(input_path, fmt=input_format, sep=input_separator, study=input_study)
        else:
            self.mysumstats = gl.Sumstats(input_path, fmt=input_format, sep=input_separator)
        if input_format == "gtex":
            if "CHR" not in self.mysumstats.data.columns:
                self.mysumstats.data["CHR"] = self.mysumstats.data["SNPID"].str.split("_", expand=True)[0]
            if "POS" not in self.mysumstats.data.columns:
                self.mysumstats.data["POS"] = self.mysumstats.data["SNPID"].str.split("_", expand=True)[1].astype(int)
            if "NEA" not in self.mysumstats.data.columns:
                self.mysumstats.data["NEA"] = self.mysumstats.data["SNPID"].str.split("_", expand=True)[2]
            if "EA" not in self.mysumstats.data.columns:
                self.mysumstats.data["EA"] = self.mysumstats.data["SNPID"].str.split("_", expand=True)[3]
        if pid:
            self.mysumstats.data["PREVIOUS_ID_GWASLAB"] = self._make_gwaslab_snpid()
            if "rsID" in self.mysumstats.data.columns:
                self.mysumstats.data["PREVIOUS_rsID"] = self.mysumstats.data["rsID"].astype("string")
            if "SNPID" in self.mysumstats.data.columns:
                self.mysumstats.data["PREVIOUS_ID"] = self.mysumstats.data["SNPID"].astype("string")
            else:
                self.mysumstats.data["PREVIOUS_ID"] = self._make_gwaslab_snpid()
            if bcfliftover:
                self.mysumstats.data["PREVIOUS_ID"] = (
                    self.mysumstats.data["rsID"].astype("string").str.replace("_", ":", regex=False)
                )
        if bcfliftover:
            self.mysumstats.data.drop(columns=["rsID"], inplace=True)

    def fill_mlog10p(self, gl_params) -> bool:
        """
        Fill the MLOG10P column using standard and extreme-value calculations.

        When the argument extreme is False, rows with P values below
        _EXTREME_P_THRESHOLD are processed using extreme-value methods.
        Remaining rows are processed using the P-to-MLOG10P conversion.
        This function replicates the pre-v4 GWASLab behavior for non-extreme P
        values, while handling extremely small P values.
        Note: The current default in v4 is to always use the extreme methods.

        Returns:
            True if MLOG10P was handled by this method; otherwise False.
        """
        if (
            not gl_params.get("extreme")
            and "MLOG10P" in gl_params.get("to_fill", [])
            and "P" in self.mysumstats.data.columns
        ):
            from gwaslab.util.util_in_fill_data import (
                fill_extreme_mlog10p as _fill_extreme_mlog10p,
                fill_mlog10p as _fill_mlog10p,
            )

            data = self.mysumstats.data
            p = data["P"]
            extreme_mask = p.lt(_EXTREME_P_THRESHOLD)
            nonextreme_mask = ~extreme_mask

            if extreme_mask.any():
                extreme_sumstats = data.loc[extreme_mask].copy()
                self.mysumstats.log.write(f"Extremely low P detected: {extreme_mask.sum()}")
                self.mysumstats.log.write("Start filling MLOG10P for extremely low P using extreme-value methods...")
                _fill_extreme_mlog10p(extreme_sumstats, df=None, log=self.mysumstats.log)
                data.loc[extreme_mask, "MLOG10P"] = extreme_sumstats["MLOG10P"].to_numpy()
                self.mysumstats.log.write("Finished filling MLOG10P with extreme-value methods.")

            if nonextreme_mask.any():
                nonextreme_sumstats = data.loc[nonextreme_mask].copy()
                self.mysumstats.log.write("Start filling MLOG10P from P...")
                _fill_mlog10p(nonextreme_sumstats, log=self.mysumstats.log)
                data.loc[nonextreme_mask, "MLOG10P"] = nonextreme_sumstats["MLOG10P"].to_numpy()
                self.mysumstats.log.write("Finished filling MLOG10P from P.")

            self.mysumstats.data = data
            return True
        else:
            return False

    def float_dict_custom(self, gp):
        """Preserve the number of decimals from the input data (statistics)"""
        float_dict = {}
        for col in self.mysumstats.data.columns:
            if str(self.mysumstats.data[col].dtype) in ["Float32", "Float64", "float64", "float32", "float16", "float"]:
                fn = self.mysumstats.data[col].apply(lambda x: len(str(x).split(".")[-1]) if "." in str(x) else 0).max()
                float_dict[col] = "{:." + str(fn) + "f}"
        if "float_formats" in gp:
            float_dict.update({k: v for k, v in gp["float_formats"].items() if k in float_dict})
        return float_dict

    def order_alleles(
        self,
        ea="EA",
        nea="NEA",
        status="STATUS",
        chrom="CHR",
        pos="POS",
        snpid="SNPID",
        format_snpid=True,
        n_cores=1,
        mode="v",
        verbose=True,
    ):
        self.mysumstats.data = order_alleles_func(
            sumstats_data=self.mysumstats.data,
            log=self.mysumstats.log,
            ea=ea,
            nea=nea,
            status=status,
            chrom=chrom,
            pos=pos,
            snpid=snpid,
            format_snpid=format_snpid,
            n_cores=n_cores,
            mode=mode,
            verbose=verbose,
        )


@click.version_option(version=__version__)
@click.command()
@click.option("-c", "--config_file", required=True, help="Configuration file path")
@click.option("-i", "--input_file", required=True, help="Input file path")
@click.option("-b", "--formatbook_file", default=None, help="Formatbook file path")
@click.option(
    "-f",
    "--input_file_format",
    required=True,
    type=click.Choice(
        [
            "plink_pvar",
            "literature_rev",
            "gtex",
            "gwascatalog_hm_custom",
            "ssf_custom",
            "finngen",
            "vcf",
            "decode",
            "gwaslab",
            "regenie",
            "regenie_gene",
            "fastgwa",
            "ldsc",
            "fuma",
            "pickle",
            "metal_het",
            "auto",
        ],
        case_sensitive=False,
    ),
    help="Input file format",
)
@click.option("-o", "--output", help="Path where results should be saved", default="results")
@click.option("-s", "--input_file_separator", default="\t", help="Input file separator")
@click.option("--study_label", default="Study", help="Input study label, valid only for VCF files")
@click.option("-q", "--quiet", default=False, is_flag=True, help="Set log verbosity")
@click.option("--pid", default=False, is_flag=True, help="Preserve ID")
@click.option("--bcfliftover", default=False, is_flag=True, help="Input from BCFtools liftover")
def main(
    config_file,
    input_file,
    formatbook_file,
    input_file_format,
    input_file_separator,
    study_label,
    output,
    quiet,
    pid,
    bcfliftover,
):
    cm = ConfigurationManager(config_file=config_file, formatbook_file=formatbook_file, root_path=output)
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

    if input_file_format == "vcf":
        with gzip.open(input_file_path, "rt") as vcf_f:
            vcf_header = next(line.strip() for line in vcf_f if line.startswith("#CHROM"))
        vcf_cols = vcf_header.split("\t")
        vcf_format_idx = vcf_cols.index("FORMAT")
        try:
            study_label = vcf_cols[vcf_format_idx + 1]
        except IndexError:
            print("No study label found after FORMAT.")

    if "write_snp_mapping" in cm.run_sequence:
        pid = True
    if input_file_path.exists():
        sm = SumstatsManager(
            input_file_path.as_posix(),
            input_file_format,
            input_file_separator,
            study_label,
            formatbook_file_path,
            pid,
            bcfliftover,
        )
    else:
        msg = f"{input_file_path} input file not found"
        exit(msg)

    # Setup cache if needed
    if "harmonize" in cm.run_sequence:
        params, gl_params = cm.step("harmonize")
        run = params.get("run", False)
        preload_cache = params.get("preload_cache", False)
        if run and preload_cache:
            if "ref_infer" in gl_params:
                NUM_WORKERS = cm.config.get("n_cores", None) or int(
                    os.environ.get("SLURM_CPUS_PER_TASK", 1)
                )  # default to 1 if not set. It is used only if cache has to be built
                ref_alt_freq = gl_params.get("ref_alt_freq", None)
                base_path = gl_params["ref_infer"]
                cache_process = gl.cache_manager.CacheProcess(
                    base_path,
                    ref_alt_freq=ref_alt_freq,
                    category=gl.cache_manager.PALINDROMIC_INDEL,
                    n_cores=NUM_WORKERS,
                    log=sm.mysumstats.log,
                    verbose=True,
                )
                cache_process.start()

                # Add cache options to inferstrand_args
                inferstrand_args = gl_params.get("inferstrand_args", {})
                cache_options = inferstrand_args.get("cache_options", {})
                cache_options.update({"cache_process": cache_process})
                inferstrand_args.update({"cache_options": cache_options})
                gl_params["inferstrand_args"] = inferstrand_args

    # EAF floating format in config
    if_eaf_float_format = any(
        y.get("run", False) and "float_formats" in x and "EAF" in x.get("float_formats", {})
        for step in cm.run_sequence
        if step != "write_parquet"
        for y, x in [cm.step(step)]
    )

    for step in cm.run_sequence:
        params, gl_params = cm.step(step)
        run = params.get("run", False)
        ws = params.get("workspace", "default")
        ws_subfolder = params.get("workspace_subfolder", False)
        if ws_subfolder:
            workspace_path = Path(cm.root_path, ws, input_file_stem)
        else:
            workspace_path = Path(cm.root_path, ws)
        workspace_path.mkdir(parents=True, exist_ok=True)

        if run:
            logger.info(f"Started {step} step")
            qc_rows_before = len(sm.mysumstats.data) if step in _QC_STEPS else None
            qc_exclusion_counts = {}
            qc_metrics = {}
            qc_log_text = getattr(sm.mysumstats.log, "log_text", None)
            qc_log_start = (
                len(qc_log_text)
                if step in {"basic_check", "canonicalize_effect_alleles", "sort_alphabetically"}
                and isinstance(qc_log_text, str)
                else None
            )
            if step == "write_snp_mapping":
                output_path = str(Path(workspace_path, "table"))
                sm.mysumstats.data["EQUALS"] = sm.mysumstats.data["SNPID"] == sm.mysumstats.data["PREVIOUS_ID"]
                snp_split = sm.mysumstats.data["SNPID"].str.split(":", expand=True)
                prev_split = sm.mysumstats.data["PREVIOUS_ID_GWASLAB"].str.split(":", expand=True)
                sm.mysumstats.data["FLIPPED"] = (snp_split.iloc[:, -2] != prev_split.iloc[:, -2]) & (
                    snp_split.iloc[:, -1] != prev_split.iloc[:, -1]
                )
                gl_params["float_formats"] = sm.float_dict_custom(gl_params)
                sm.mysumstats.to_format(output_path, **gl_params)
                _write_run_provenance(output_path, sm.mysumstats, input_file_path)
            elif step == "basic_check":
                sm.mysumstats.basic_check(**gl_params)
                if qc_log_start is not None:
                    qc_exclusion_counts = _basic_check_exclusion_counts(sm.mysumstats.log.log_text[qc_log_start:])
                if not if_eaf_float_format and "EAF" in sm.mysumstats.data.columns:
                    sm.mysumstats.data["EAF"] = round(sm.mysumstats.data["EAF"].astype("float64"), 7)
            elif step == "infer_build":
                validate_declared_assembly(sm.mysumstats, cm.config, gl_params)
            elif step == "fill_data":
                mlog10p_handled = sm.fill_mlog10p(gl_params)
                if mlog10p_handled:
                    gl_params["to_fill"] = [x for x in gl_params.get("to_fill", []) if x != "MLOG10P"]
                if gl_params["to_fill"]:
                    sm.mysumstats.fill_data(**gl_params)
            elif step == "harmonize":
                _require_validated_assembly(sm.mysumstats)
                sm.mysumstats.harmonize(**gl_params)
            elif step == "liftover":
                input_audit = _require_validated_assembly(sm.mysumstats)
                sm.mysumstats.liftover(**gl_params)

                sm.mysumstats.log.write("Start to process unmapped variants...")
                unmapped_mask = sm.mysumstats.data["STATUS"].astype(str).str.startswith("97")
                unmapped = sm.mysumstats.data.loc[unmapped_mask].copy()
                mapped = sm.mysumstats.data.loc[~unmapped_mask].copy()
                sm.mysumstats.data = mapped
                output_path = str(Path(workspace_path, ".".join([input_file_stem, "unmapped_variants.tsv.gz"])))
                sm.mysumstats.log.write(f" -Saving {len(unmapped)} unmapped variants to: {output_path}")
                unmapped.to_csv(output_path, sep="\t", index=False, compression="gzip")
                sm.mysumstats.log.write("Unmapped variants saved successfully!")

                if "to_build" not in gl_params:
                    raise AssemblyValidationError("liftover requires gl_params.to_build for post-liftover validation")
                output_audit = validate_declared_assembly(
                    sm.mysumstats,
                    cm.config,
                    declared_assembly=gl_params["to_build"],
                    scope="post_liftover",
                )
                output_audit["input_validation"] = input_audit
                _store_assembly_audit(sm.mysumstats, output_audit)
            elif step == "report_harmonization_summary":
                summary = sm.mysumstats.lookup_status().to_string()
                output_path = str(Path(workspace_path, ".".join([input_file_stem, "harmonization_summary.tsv"])))
                with open(output_path, "w") as fp:
                    fp.write(summary)
            elif step == "report_min_pvalue":
                nrows = params.get("nrows", 1)
                df = sm.mysumstats.data.nlargest(nrows, "MLOG10P", keep="first").reset_index(drop=True)
                snpid = df.at[0, "SNPID"]
                mlog10p = df.at[0, "MLOG10P"]
                output_path = str(Path(workspace_path, ".".join([input_file_stem, "nlargest.txt"])))
                with open(output_path, "w") as fp:
                    fp.write("input_file\tSNPID\tMLOG10P\n")
                    fp.write(f"{input_file_name}\t{snpid}\t{mlog10p}\n")
            elif step == "report_inflation_factors":
                df = sm.mysumstats.data
                CHISQ = df.Z**2
                max_chisq = str(round(CHISQ.max(), 3))
                mean_chisq = str(round(CHISQ.mean(), 3))
                lambda_GC = str(round(CHISQ.median() / 0.4549, 3))

                output_path = str(Path(workspace_path, ".".join([input_file_stem, "if.txt"])))
                with open(output_path, "w") as fp:
                    fp.write("input_file\tlambda_GC\tmean_chisq\tmax_chisq\n")
                    fp.write(f"{input_file_name}\t{lambda_GC}\t{mean_chisq}\t{max_chisq}\n")
            elif step in {"canonicalize_effect_alleles", "sort_alphabetically"}:
                if step == "sort_alphabetically":
                    logger.warning("sort_alphabetically is deprecated; use canonicalize_effect_alleles instead.")
                sm.order_alleles(**gl_params)
                if qc_log_start is not None:
                    qc_metrics["flipped_alleles"] = _order_alleles_flipped_count(
                        sm.mysumstats.log.log_text[qc_log_start:]
                    )
                if not if_eaf_float_format and "EAF" in sm.mysumstats.data.columns:
                    sm.mysumstats.data["EAF"] = round(sm.mysumstats.data["EAF"].astype("float64"), 7)
            elif step == "write_pickle":
                output_path = str(Path(workspace_path, ".".join([input_file_stem, "pkl"])))
                gl.dump_pickle(sm.mysumstats, output_path, overwrite=params["overwrite"])
                _write_run_provenance(output_path, sm.mysumstats, input_file_path)
            elif step in ["write_regenie", "write_ldsc", "write_metal", "write_tsv", "write_fastgwa", "write_parquet"]:
                output_path = str(Path(workspace_path, input_file_stem))
                gl_params["float_formats"] = sm.float_dict_custom(gl_params)
                sm.mysumstats.to_format(output_path, **gl_params)
                _write_run_provenance(output_path, sm.mysumstats, input_file_path)
            elif step == "write_vcf":
                _require_validated_assembly(sm.mysumstats)
                study_name = input_file_stem
                sm.mysumstats.meta["gwaslab"]["study_name"] = study_name
                output_path = str(Path(workspace_path, input_file_stem))
                gl_params["float_formats"] = sm.float_dict_custom(gl_params)
                sm.mysumstats.to_format(output_path, **gl_params)
                _write_run_provenance(output_path, sm.mysumstats, input_file_path)
            elif step == "write_same_input_format":
                output_path = str(Path(workspace_path, input_file_stem))
                gl_params["float_formats"] = sm.float_dict_custom(gl_params)
                sm.mysumstats.to_format(output_path, fmt=input_file_format, **gl_params)
                _write_run_provenance(output_path, sm.mysumstats, input_file_path)
            elif step in {"filter_conflicting_snpids", "check_ambiguous_snps"}:
                if step == "check_ambiguous_snps":
                    logger.warning("check_ambiguous_snps is deprecated; use filter_conflicting_snpids instead.")
                df = sm.mysumstats.data

                # True duplicated SNPs
                dup_mask = df.duplicated(subset=["SNPID", "EAF", "BETA", "SE"], keep="first")
                nr_dup_snps = dup_mask.sum()
                if nr_dup_snps > 0:
                    df = df.loc[~dup_mask].reset_index(drop=True)

                # Ambiguous SNPs
                snp_groups = df.groupby("SNPID")
                ambiguous_mask = snp_groups[["EAF", "BETA", "SE"]].transform("nunique").gt(1).any(axis=1)
                nr_ambiguous_snps = ambiguous_mask.sum()
                if nr_ambiguous_snps > 0:
                    df = df.loc[~ambiguous_mask].reset_index(drop=True)

                # Multi-allelic SNPs
                nr_multiallelic_snps = df.groupby(["CHR", "POS"])["SNPID"].transform("nunique").gt(1).sum()
                nr_multiallelic_loci = df.groupby(["CHR", "POS"])["SNPID"].nunique().gt(1).sum()

                sm.mysumstats.data = df
                qc_exclusion_counts = {
                    "duplicate_snpid_eaf_beta_se": int(nr_dup_snps),
                    "conflicting_snpid_eaf_beta_se": int(nr_ambiguous_snps),
                }
                qc_metrics = {
                    "multiallelic_variant_rows": int(nr_multiallelic_snps),
                    "multiallelic_loci": int(nr_multiallelic_loci),
                }

                sm.mysumstats.log.write("Start to check ambiguous variants...")
                sm.mysumstats.log.write(f" -Dropped duplicated SNPs: {nr_dup_snps}")
                sm.mysumstats.log.write(f" -Dropped ambiguous SNPs: {nr_ambiguous_snps}")
                sm.mysumstats.log.write(f" -Multi-allelic SNPs: {nr_multiallelic_snps}")
                sm.mysumstats.log.write(f" -Multi-allelic positions: {nr_multiallelic_loci}")
                sm.mysumstats.log.write(
                    f" -Current Dataframe shape : {len(sm.mysumstats.data)} x {len(sm.mysumstats.data.columns)}"
                )
            elif step == "qq_manhattan_plots":
                output_path = str(Path(workspace_path, ".".join([input_file_stem, "png"])))
                cut = round(-np.log10(gl_params["sig_level"])) + params["dist"]
                sm.mysumstats.plot_mqq(cut=cut, save=output_path, **gl_params)
            if qc_rows_before is not None:
                _record_qc_step(sm.mysumstats, step, qc_rows_before, qc_exclusion_counts, qc_metrics)
            logger.info(f"Finished {step} step")
        else:
            logger.info(f"Skipping {step} step")


if __name__ == "__main__":
    main()
