"""
This part of the workflow handles automatic deployments of nextflu-private builds.
"""

rule all_private:
    input:
        jsons=_get_build_outputs(),
    output:
        json_dir=directory("auspice_renamed"),
    params:
        build_date=config.get("build_date", datetime.date.today().strftime("%Y-%m-%d")),
    shell:
        """
        mkdir -p {output.json_dir};
        for file in {input.jsons}
        do
            ln ${{file}} {output.json_dir}/"flu_seasonal_{params.build_date}_`basename ${{file}}`"
        done
        """

rule deploy_all:
    input:
        json_dir="auspice_renamed",
    params:
        deploy_url = config["deploy_url"]
    shell:
        """
        nextstrain login --no-prompt;
        nextstrain remote upload {params.deploy_url} {input.json_dir}/*.json
        """
