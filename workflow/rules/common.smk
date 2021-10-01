from snakemake.utils import validate
import pandas as pd
import os

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

# configfile: "config/config.yaml"
# validate(config, schema="../schemas/config.schema.yaml")

# barcodes = pd.read_csv(config["barcode"], sep="\t").set_index("barcode", drop=False)
# barcodes.index.names = ["barcode"]
# validate(barcodes, schema="../schemas/samples.schema.yaml")

#barcodes = [str(fq_file) for fq_file in os.scandir("./input_data") if os.path.isfile(fq_file)]
