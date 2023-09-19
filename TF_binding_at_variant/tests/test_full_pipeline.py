import subprocess


def test_run_pipeline():
    """
    Test that the pipeline runs successfully and produces the expected output files
    :return:
    """
    # Run the Snakemake pipeline with subprocess
    result = subprocess.run(["snakemake", "--use-conda"], capture_output=True, text=True)

    # Check if the pipeline completed successfully (exit code 0)
    assert result.returncode == 0

    # Check for expected output files
    assert "output.txt" in result.stdout
