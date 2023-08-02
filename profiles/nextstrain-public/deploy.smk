"""
This part of the workflow handles automatic deployments of public builds.
Depends on the `all_public` rule from rename.smk
"""

rule all_deploy:
    input: rules.all_public.input
    log:
        "logs/all_deploy.txt"
    params:
        s3_dst = config["deploy_url"]
    shell:
        """
        nextstrain remote upload {params.s3_dst} {input} 2>&1 | tee {log}
        """
