"""
This part of the workflow handles automatic deployments of nextflu-private builds.
Depends on the `all_private` rule from rename.smk
"""

rule deploy_all:
    input: rules.all_private.input
    params:
        deploy_url = config["deploy_url"]
    shell:
        """
        nextstrain login;
        nextstrain remote upload {params.deploy_url} {input}
        """
