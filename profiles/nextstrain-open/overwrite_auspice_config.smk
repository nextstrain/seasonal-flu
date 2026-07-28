rule overwrite_auspice_config:
  """
  Overwrite data_provenance since the built-in merging of Auspice configs
  in Augur does not support removing elements from the base list.
  <https://github.com/nextstrain/augur/blob/c8cc27cade9268fbbeec806f0e7ba5b7ba45c649/augur/util_support/auspice_config.py#L44-L57>
  """
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
