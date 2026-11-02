from unittest.mock import patch

import gwaslab as gl
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from gwaspipe.configuring import SingletonConfigurationManager
from gwaspipe.gwaspipe import SumstatsManager, main


def make_manager(chromosomes):
    data = pd.DataFrame(
        {
            "CHR": pd.Categorical(chromosomes),
            "POS": range(1, len(chromosomes) + 1),
            "EA": "A",
            "NEA": "G",
            "BETA": 0.123456789012345,
            "SE": 0.01,
            "P": 1e-8,
        }
    )
    manager = SumstatsManager.__new__(SumstatsManager)
    manager.mysumstats = gl.Sumstats(data, build="38", verbose=False)
    # Exercise the categorical string chromosomes found in input/export data.
    manager.mysumstats.data["CHR"] = data["CHR"]
    return manager


def test_default_parquet_schema_and_row_groups(tmp_path):
    manager = make_manager(["1"] * 65_535 + ["2", "23", "24", None])
    original = manager.mysumstats.data.copy(deep=True)
    output_path = str(tmp_path / "default")

    manager.write_parquet(output_path, fmt="gwaslab", verbose=False)

    parquet = pq.ParquetFile(f"{output_path}.gwaslab.parquet")
    assert parquet.schema_arrow.field("CHR").type == pa.int16()
    assert [parquet.metadata.row_group(i).num_rows for i in range(parquet.num_row_groups)] == [65_536, 3]
    exported = parquet.read()
    assert exported.column("CHR").to_pylist() == [1] * 65_535 + [2, 23, 24, None]
    assert exported.column("POS").to_pylist() == original["POS"].tolist()
    assert exported.column("BETA").to_pylist() == original["BETA"].tolist()
    pd.testing.assert_frame_equal(manager.mysumstats.data, original)


def test_custom_row_groups_and_compression_keep_chromosomes_numeric(tmp_path):
    manager = make_manager(["1", "2", "23", "24", "25", "1", "2"])
    output_path = str(tmp_path / "custom")

    manager.write_parquet(
        output_path,
        fmt="gwaslab",
        verbose=False,
        xymt_number=False,
        chr_prefix="chr",
        to_tabular_kwargs={"row_group_size": 3, "compression": "zstd"},
    )

    parquet = pq.ParquetFile(f"{output_path}.gwaslab.parquet")
    assert parquet.schema_arrow.field("CHR").type == pa.int16()
    assert [parquet.metadata.row_group(i).num_rows for i in range(parquet.num_row_groups)] == [3, 3, 1]
    assert parquet.metadata.row_group(0).column(0).compression == "ZSTD"


@pytest.mark.parametrize("chromosomes", [["1", "unknown"], [1.5, 2], [1, 32768]])
def test_invalid_chromosomes_fail_before_export(tmp_path, chromosomes):
    manager = make_manager(chromosomes)
    original = manager.mysumstats.data.copy(deep=True)

    with pytest.raises((ValueError, TypeError)):
        manager.write_parquet(str(tmp_path / "invalid"), fmt="gwaslab", verbose=False)

    assert not list(tmp_path.glob("*.parquet"))
    pd.testing.assert_frame_equal(manager.mysumstats.data, original)


def test_writer_failure_restores_chromosome_dtype(tmp_path):
    manager = make_manager(["1", "23"])
    original = manager.mysumstats.data.copy(deep=True)

    with patch.object(manager.mysumstats, "to_format", side_effect=OSError("write failed")):
        with pytest.raises(OSError, match="write failed"):
            manager.write_parquet(str(tmp_path / "failed"), fmt="gwaslab")

    pd.testing.assert_frame_equal(manager.mysumstats.data, original)


def test_cli_parquet_export(tmp_path):
    input_path = tmp_path / "test_sumstats.tsv"
    make_manager(["1", "2", "23"]).mysumstats.data.to_csv(input_path, sep="\t", index=False)
    config_path = tmp_path / "config.yml"
    config_path.write_text(
        "log_filename: gwaspipe.log\n"
        "run_sequence:\n"
        "  1: write_parquet\n"
        "steps:\n"
        "  write_parquet:\n"
        "    params:\n"
        "      run: true\n"
        "      workspace: outputs\n"
        "    gl_params:\n"
        "      fmt: gwaslab\n"
        "      to_tabular_kwargs:\n"
        "        row_group_size: 2\n"
    )
    with patch.dict(SingletonConfigurationManager._instances, {}, clear=True):
        result = CliRunner().invoke(
            main,
            [
                "-c",
                str(config_path),
                "-i",
                str(input_path),
                "-f",
                "gwaslab",
                "-o",
                str(tmp_path),
            ],
        )

    assert result.exit_code == 0, result.output
    parquet = pq.ParquetFile(tmp_path / "outputs/test_sumstats.gwaslab.parquet")
    assert parquet.schema_arrow.field("CHR").type == pa.int16()
    assert max(parquet.metadata.row_group(i).num_rows for i in range(parquet.num_row_groups)) <= 2
    assert (tmp_path / "outputs/test_sumstats.provenance.json").exists()
