# Curation docs

> WORK IN PROGRESS!

_This guides you through how to change the rules & overrides used by the curation pipeline in order to fix metadata for records._


### High level sketch of the curation experience

* How do you know what needs fixing? Just read `logs/curate.txt` (which parallels the fauna path) or ???

I think the experience goes something like:

1. Download newest cache
2. Run with the old code to get the (full) curated NDJSON & store it somewhere
3. Make code/TSV changes, following instructions here
4. Rerun the curate pipeline & diff against (2). We want a user-friendly diffing approach which identifies what's changed without being too verbose


## Geography

Geography fields "region", "country", "division" and "location" are set via a multi-layered approach. Here's a high-level overview, as this will help identify where a fix can/should be made:

1. The `parse-gisaid-location` script parses the (capital-L) "Location" field from the GISAID record, which is expected to have slash-separated geographic information.
  1a. `defaults/locations.tsv` sets custom fields for epi-isls in the TSV. If the record has an annotation then no more modifications are done by the script (i.e. 1b, 1d and 1d are skipped).
  1b. `defaults/gisaid_location_rules.tsv` remaps the slash-separated Location fields. This is mostly used to correct location values with more than 4 fields so that the script can parse them.
  1c. Sets the slash-separated Location value to "region", "country", "division" and "location"
  1d. If four fields aren't populated then we attempt to extract a further piece of geographic information from the strain name.
2. `augur curate titlecase` changes casing
3. `augur curate apply-geolocation-rules`
  3a. First runs with Augur default fixes
  3b. Then runs with our custom rules from `defaults/geolocation_rules.tsv`. This is to ensure that we use specific regions such as "Japan Korea", as well as to apply other specific fixes.

Metadata which is _wrong_ is best corrected via a manual override for the specific EPI ISL in step (1a).

GISAID Location fields with too much information are best corrected in step (1a) if they're a one-off or in step (1b) if you wish the rule to be applied to multiple records, or to future records.

Spelling mistakes or synonyms are best corrected by adding rules to step (3b).
For instance, we currently have a rule which changes the country "Timor-Leste" to "East Timor", encoded via the line in the TSV `*/Timor-Leste/*/*	*/East Timor/*/*` 

