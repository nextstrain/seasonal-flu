BUILD_DATE = datetime.date.today().strftime("%Y-%m-%d")

all_private_builds = [
    "auspice_renamed/" + build.get("auspice_name", f"{build_name}_{segment}").format(build_date=BUILD_DATE, segment=segment) + suffix + ".json"
    for build_name, build in config["builds"].items()
    for segment in config["segments"]
    for suffix in ["", "_root-sequence", "_tip-frequencies"]
]

rule all_private:
    input:
        all_private_builds

def _get_file_by_auspice_name(wildcards):
    """Find the original Auspice JSON file that needs to be renamed.
    """
    for build_name, build_params in config["builds"].items():
        for segment in config["segments"]:
            if build_params.get("auspice_name", f"{build_name}_{{segment}}").format(build_date=BUILD_DATE, segment=segment) == wildcards.auspice_name:
                return f"auspice/{build_name}_{segment}.json"

    return ""

rule rename_auspice_main:
    input:
        _get_file_by_auspice_name,
    output:
        "auspice_renamed/{auspice_name}.json",
    shell:
        """
        ln {input} {output}
        """

rule rename_auspice_root_sequence:
    input:
        lambda wildcards: _get_file_by_auspice_name(wildcards).replace(".json", "_root-sequence.json"),
    output:
        "auspice_renamed/{auspice_name}_root-sequence.json",
    shell:
        """
        ln {input} {output}
        """

rule rename_auspice_tip_frequencies:
    input:
        lambda wildcards: _get_file_by_auspice_name(wildcards).replace(".json", "_tip-frequencies.json"),
    output:
        "auspice_renamed/{auspice_name}_tip-frequencies.json",
    shell:
        """
        ln {input} {output}
        """

# auspice_renamed/flu_seasonal_2024-01-05_h3n2_2y_titers_ha.json
# {auspice_name}_{forecast_model}_forecast-tip-frequencies.json
# rule rename_auspice_tip_frequencies:
#     input:
#         lambda wildcards: _get_file_by_auspice_name(wildcards).replace(".json", "_tip-frequencies.json"),
#     output:
#         "auspice_renamed/{auspice_name}_{forecast_model}_forecast-tip-frequencies.json",
#     shell:
#         """
#         ln {input} {output}
#         """
