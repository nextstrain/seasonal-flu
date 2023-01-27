
rule all_public:
    input:
        [
            "auspice_renamed/" + build.get("auspice_name", f"{build_name}_{{segment}}").format(segment=segment) + suffix + ".json"
            for build_name, build in config["builds"].items()
            for segment in config["segments"]
            for suffix in ["", "_root-sequence", "_tip-frequencies"]
        ],

def _get_file_by_auspice_name(wildcards):
    for build_name, build_params in config["builds"].items():
        for segment in config["segments"]:
            if build_params.get("auspice_name", f"{build_name}_{{segment}}").format(segment=segment) == wildcards.auspice_name:
                return f"auspice/{build_name}_{segment}.json"

    return ""

rule rename_auspice_main:
    input:
        _get_file_by_auspice_name,
    output:
        "auspice_renamed/{auspice_name}.json",
    shell:
        """
        cp {input} {output}
        """

rule rename_auspice_root_sequence:
    input:
        lambda wildcards: _get_file_by_auspice_name(wildcards).replace(".json", "_root-sequence.json"),
    output:
        "auspice_renamed/{auspice_name}_root-sequence.json",
    shell:
        """
        cp {input} {output}
        """

rule rename_auspice_tip_frequencies:
    input:
        lambda wildcards: _get_file_by_auspice_name(wildcards).replace(".json", "_tip-frequencies.json"),
    output:
        "auspice_renamed/{auspice_name}_tip-frequencies.json",
    shell:
        """
        cp {input} {output}
        """
