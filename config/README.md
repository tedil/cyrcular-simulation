
# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Region sheet

Add regions to ``config/example_regions.tsv`` (path to the region sheet can be changed in ``config/config.yaml``). For each region, the columns `name`, `chrom`, `start`, `stop` and `length` have to be defined.

Lines can be commented out with `#`.

For each region in the sheet, the respective sequence is extracted from the reference.

# Error model

For simulation of nanopore reads with nanosim, an error model is required. This will be built from the sample specified in the respective config section.

# Additional noise

Additional non-circular WGS nanopore reads can be simulated at the specified coverage. These will be merged with the circular reads before calling.
