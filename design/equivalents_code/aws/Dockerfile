FROM continuumio/miniconda3:4.7.12

COPY spec-file.txt .
RUN conda install --yes --file spec-file.txt
COPY run_experiment.py .
ENTRYPOINT python run_experiment.py
