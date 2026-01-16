"""
This part of the workflow handles various Slack notifications.
Designed to be used internally by the Nextstrain team with hard-coded paths
to files on AWS S3.

All rules here require two environment variables:
    * SLACK_TOKEN
    * SLACK_CHANNELS
"""


rule notify_on_record_change:
    input:
        ndjson="data/gisaid.ndjson",
    params:
        ndjson_on_s3 = f"{config['s3_src']}/gisaid.ndjson.zst"
    output:
        touch("data/notify-on-record-change.done")
    shell:
        """
        ./vendored/notify-on-record-change {input.ndjson} {params.ndjson_on_s3} "GISAID cache"
        """
