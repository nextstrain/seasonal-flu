rule overwrite_auspice_config:
  input:
    auspice_config = "config/{lineage}/{segment}/auspice_config.json",
  output:
    auspice_config = "data/{lineage}/{segment}/auspice_config.json",
  shell:
    """
    jq '.data_provenance = [{{"name": "GenSpectrum", "url": "https://loculus.genspectrum.org"}}]' \
      {input.auspice_config:q} \
      > {output.auspice_config:q}
    """
