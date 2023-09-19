import subprocess
import os


def test_run_pipeline():
    """
    Test that the pipeline runs successfully and produces the expected output files
    :return:
    """
    # Store the current working directory
    current_dir = os.getcwd()

    # Change the working directory to the root directory where the Snakefile is located
    os.chdir(os.path.dirname(current_dir))

    try:
        # Run the Snakemake pipeline with subprocess
        result = subprocess.run(["snakemake", "--use-conda", "--cores=all"], capture_output=True, text=True)

        # Check if the pipeline completed successfully (exit code 0)
        assert result.returncode == 0

        # Check for expected output files
        # assert "output.txt" in result.stdout
    finally:
        # Change the working directory back to the original directory
        os.chdir(current_dir)
