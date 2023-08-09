"""
This part of the workflow handles automatic deployments of public builds.
Depends on the `all_public` rule from rename.smk
"""

rule deploy_all:
    input: rules.all_public.input
    params:
        s3_dst = config["deploy_url"]
    shell:
        """
        nextstrain remote upload {params.s3_dst} {input}
        """
