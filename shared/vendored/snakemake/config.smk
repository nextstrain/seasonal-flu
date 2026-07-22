"""
Shared functions to be used within a Snakemake workflow for handling
workflow configs.
"""
from dependencies import set_min_augur_version
set_min_augur_version("34.1.0")

import os
import sys
import yaml
from augur.config import resolve_filepath
from collections.abc import Callable
from typing import Optional
from textwrap import dedent, indent


# Set search paths
if "AUGUR_SEARCH_PATHS" in os.environ:
    print(dedent(f"""\
        Using existing search paths in AUGUR_SEARCH_PATHS:

            {os.environ["AUGUR_SEARCH_PATHS"]!r}
        """), file=sys.stderr)
else:
    search_paths = [
        # User analysis directory
        Path.cwd(),

        # Workflow defaults folder
        Path(workflow.basedir) / "defaults",

        # Workflow root (contains Snakefile)
        Path(workflow.basedir),
    ]

    # This should work for majority of workflows, but we could consider doing a
    # more thorough search for the nextstrain-pathogen.yaml. This would likely
    # replicate how CLI searches for the root.¹
    # ¹ <https://github.com/nextstrain/cli/blob/d5e184c5/nextstrain/cli/command/build.py#L413-L420>
    repo_root = Path(workflow.basedir) / ".."
    if (repo_root / "nextstrain-pathogen.yaml").is_file():
        search_paths.extend([
            # Pathogen repo root
            repo_root,
        ])

    seen = set()
    normalized_search_paths = []
    for path in search_paths:
        # Skip paths that are not directories
        if not path.is_dir():
            continue

        # Resolve to absolute paths
        resolved = path.resolve()

        # Skip duplicate paths (e.g. often the CWD == workflow.basedir)
        if resolved in seen:
            continue

        seen.add(resolved)
        normalized_search_paths.append(resolved)

    os.environ["AUGUR_SEARCH_PATHS"] = ":".join(map(str, normalized_search_paths))


class InvalidConfigError(Exception):
    pass


def resolve_config_path(path: str) -> Callable:
    """
    Resolve a relative *path* given in a configuration value. Will always try to
    resolve *path* after expanding wildcards with Snakemake's `expand` functionality.

    Returns the path for the first existing file found relative to directories
    defined by the AUGUR_SEARCH_PATHS environment variable.
    """
    global workflow

    def _resolve_config_path(wildcards):
        try:
            expanded_path = expand(path, **wildcards)[0]
        except snakemake.exceptions.WildcardError as e:
            available_wildcards = "\n".join(f"  - {wildcard}" for wildcard in wildcards)
            raise snakemake.exceptions.WildcardError(indent(dedent(f"""\
                {str(e)}

                However, resolve_config_path({{path}}) requires the wildcard.

                Wildcards available for this path are:

                {{available_wildcards}}

                Hint: Check that the config path value does not misspell the wildcard name
                and that the rule actually uses the wildcard name.
                """.lstrip("\n").rstrip()).format(path=repr(path), available_wildcards=available_wildcards), " " * 4))

        search_paths = [Path(p) for p in os.environ["AUGUR_SEARCH_PATHS"].split(":")]
        try:
            return str(resolve_filepath(Path(expanded_path), search_paths))
        except Exception as error:
            raise InvalidConfigError(
                indent(
                    dedent(f"""
                        Unable to resolve the config-provided path {path!r},
                        expanded to {expanded_path!r} after filling in wildcards.

                        """)
                    + str(error)
                    + dedent(f"""

                        Hint: Check that the file {expanded_path!r} exists in your analysis
                        directory or remove the config param to use the workflow defaults.
                        """),
                    " " * 4,
                )
            )

    return _resolve_config_path


def write_config(path, section=None):
    """
    Write Snakemake's 'config' variable, or a section of it, to a file.

    *section* is an optional list of keys to navigate to a specific section of
    config. If provided, only that section will be written.
    """
    global config

    os.makedirs(os.path.dirname(path), exist_ok=True)

    data = config
    section_str = "config"

    if section:
        # Navigate to the specified section
        for key in section:
            # Error if key doesn't exist
            if key not in data:
                raise Exception(f"ERROR: Key {key!r} not found in {section_str!r}.")

            data = data[key]
            section_str += f".{key}"

            # Error if value is not a mapping
            if not isinstance(data, dict):
                raise Exception(f"ERROR: {section_str!r} is not a mapping of key/value pairs.")

    with open(path, 'w') as f:
        yaml.dump(data, f, sort_keys=False, Dumper=NoAliasDumper)

    print(f"Saved {section_str!r} to {path!r}.", file=sys.stderr)


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data):
        return True
