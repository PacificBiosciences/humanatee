# Installing Humanatee

## From GitHub

It's recommended to install Humanatee into its own virtual environment using the `uv` package manager. Installation instructions for `uv` can be found [here](https://github.com/astral-sh/uv?tab=readme-ov-file#installation).

```bash
# clone the repo
git clone git@github.com:PacificBiosciences/humanatee.git
# move into the humanatee directory
cd humanatee
# create a virtual environment
uv venv --python 3.11
# activate the virtual environment and install humanatee
source .venv/bin/activate; uv pip install .
```

## From Docker

Humanatee can also be run from a Docker container. The Docker image is available on [Quay.io](https://quay.io/repository/pacbio/humanatee) or can be built locally.

```bash
# clone the repo
git clone git@github.com:PacificBiosciences/humanatee.git
# move into the humanatee directory
cd humanatee
# build the docker image
docker build -t humanatee:latest .
# run the docker image
docker run humanatee:latest humanatee version
```
