# Use an official Python runtime based on Debian 10 "buster" as a parent image
FROM continuumio/miniconda3:4.10.3

# The working directory in the docker image
WORKDIR /app

# Add the current directory (on your machine) to the path `/app` in the image
ADD . /app

# Set up Python and pip
RUN conda install -y python=3.8 \
    && conda clean -a -y \
    && pip install --upgrade pip

RUN conda install tqdm -y

RUN pip install DBRetina

# The command that will be run in the container
ENTRYPOINT ["DBRetina"]
