## How to run the benchmarks for the package comparison

1. Create a virtual environment for Python: `python3 -m venv pyenv`
2. Activate the virtual environment: `source pyenv/bin/activate`
3. Install the requirements: `pip install -r requirements.txt`
4. Load the environment variables: `source .env`
5. Run the benchmarks: `julia --project package_comparison.jl`
6. (Optional) Run Literate.jl to generate the Markdown file: `julia --project -e 'using Literate; Literate.markdown("package_comparison.jl", "../../docs/src/resources/", execute=true)'`
