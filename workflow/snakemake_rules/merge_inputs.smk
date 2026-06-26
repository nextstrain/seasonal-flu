"""
This part of the workflow merges inputs based on what is defined in the config.

OUTPUTS:

    metadata  = results/metadata.tsv
    sequences = results/sequences_{segment}.fasta

The config dict is expected to have a top-level `inputs` list that defines the
separate inputs' name, metadata, and sequences. Optionally, the config can have
a top-level `additional-inputs` list that is used to define additional data that
are combined with the default inputs:

```yaml
inputs:
    - name: default
      metadata: <path-or-url>
      sequences: <path-or-url>

additional_inputs:
    - name: private
      metadata: <path-or-url>
      sequences: <path-or-url>
```

Sequences can also be a defined a dict with keys for specific segments:

```yaml
inputs:
    - name: default
      metadata: <path-or-url>
      sequences:
        ha: <path-or-url>
        na: <path-or-url>

additional_inputs:
    - name: private
      metadata: <path-or-url>
      sequences:
        ha: <path-or-url>
        na: <path-or-url>
```

Supports any of the compression formats that are supported by `augur read-file`,
see <https://docs.nextstrain.org/projects/augur/page/usage/cli/read-file.html>
"""

def _parse_config_input(input):
    """
    Parses information from an individual config-defined input, i.e. an element within `config.inputs` or `config.additional_inputs`
    and returns information snakemake rules can use to obtain the underlying data.

    The structure of `input` is a dictionary with keys:
    - name:string (required)
    - metadata:string (optional) - a s3 URI or a local file path
    - sequences:string|dict[string,string] (optional) - either a s3 URI or a local file path, in which case
      it must include a '{segment}' wildcard substring, or a dict of segment → s3 URI or local file path,
      in which case it must not include the wildcard substring.

    Returns a dictionary with optional keys:
    - metadata:string - path or url to the metadata file.
    - sequences:function. Takes in wildcards and returns path or url to the sequences FASTA for the provided
      segment wildcard, or returns `None` if this input doesn't define sequences for the provided segment.

    Raises InvalidConfigError
    """
    info = {
        "name": input["name"],
        "metadata": path_or_url(input["metadata"]) if input.get("metadata") else None,
        "sequences": None,
    }

    if location:=input.get('sequences', False):
        if isinstance(location, dict):
            info['sequences'] = lambda w: path_or_url(location[w.segment]) \
                if w.segment in location \
                else None
        elif isinstance(location, str):
            info['sequences'] = lambda w: path_or_url(location)
        else:
            raise InvalidConfigError(f"Config input for {info['name']} specifies sequences in an unknown format; must be dict or string")

    return info

def _gather_inputs():
    all_inputs = [*config['inputs'], *config.get('additional_inputs', [])]

    if len(all_inputs)==0:
        raise InvalidConfigError("Config must define at least one element in config.inputs or config.additional_inputs lists")
    if not all([isinstance(i, dict) for i in all_inputs]):
        raise InvalidConfigError("All of the elements in config.inputs and config.additional_inputs lists must be dictionaries. "
            "If you've used a command line '--config' double check your quoting.")
    if len({i['name'] for i in all_inputs})!=len(all_inputs):
        raise InvalidConfigError("Names of inputs (config.inputs and config.additional_inputs) must be unique")
    if not all(['name' in i and ('sequences' in i or 'metadata' in i) for i in all_inputs]):
        raise InvalidConfigError("Each input (config.inputs and config.additional_inputs) must have a 'name' and 'metadata' and/or 'sequences'")
    if not any(['metadata' in i for i in all_inputs]):
        raise InvalidConfigError("At least one input must have 'metadata'")
    if not any (['sequences' in i for i in all_inputs]):
        raise InvalidConfigError("At least one input must have 'sequences'")

    available_keys = set(['name', 'metadata', 'sequences'])
    if any([len(set(el.keys())-available_keys)>0 for el in all_inputs]):
        raise InvalidConfigError(f"Each input (config.inputs and config.additional_inputs) can only include keys of {', '.join(available_keys)}")

    return {i['name']: _parse_config_input(i) for i in all_inputs}

input_sources = _gather_inputs()


rule merge_metadata:
    """
    This rule will run different commands depending on the number of inputs:
    - one input = augur read-file
    - otherwise = augur merge
    """
    input:
        **{name: info['metadata'] for name,info in input_sources.items() if info.get('metadata', None)}
    params:
        num_of_inputs = lambda w, input: len(input),
        metadata = lambda w, input: list(map("=".join, input.items()))
    output:
        metadata = "results/metadata.tsv"
    log:
        "logs/merge_metadata.txt",
    benchmark:
        "benchmarks/merge_metadata.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        if [[ {params.num_of_inputs:q} == 1 ]]; then
            echo "[INFO] Reading single metadata input"
            augur read-file {input:q} > {output.metadata:q}
        else
            echo "[INFO] Merging multiple metadata inputs"
            augur merge \
                --metadata {params.metadata:q} \
                --source-columns 'input_{{NAME}}' \
                --output-metadata {output.metadata}
        fi
        """

rule merge_sequences:
    """
    This rule will run different commands depending on the number of inputs.
    - one input = augur read-file
    - otherwise = augur merge
    """
    input:
        lambda w: list(filter(None, [info['sequences'](w) for info in input_sources.values() if info.get('sequences', None)]))
    params:
        num_of_inputs = lambda w, input: len(input),
    output:
        sequences = "results/sequences_{segment}.fasta"
    log:
        "logs/{segment}/merge_sequences.txt",
    benchmark:
        "benchmarks/{segment}/merge_sequences.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        if [[ {params.num_of_inputs:q} == 1 ]]; then
            echo "[INFO] Reading single sequences input"
            augur read-file {input:q} > {output.sequences:q}
        else
            echo "[INFO] Merging multiple sequences inputs"
            augur merge \
                --sequences {input:q} \
                --output-sequences {output.sequences:q}
        fi
        """
