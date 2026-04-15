"""
Shared functions to be used within a Snakemake workflow for handling
workflow configs.
"""
import os
import sys
import yaml
from collections.abc import Callable
from typing import Optional
from textwrap import dedent, indent


# Set search paths for Augur
if "AUGUR_SEARCH_PATHS" in os.environ:
    print(dedent(f"""\
        Using existing search paths in AUGUR_SEARCH_PATHS:

            {os.environ["AUGUR_SEARCH_PATHS"]!r}
        """), file=sys.stderr)
else:
    # Note that this differs from the search paths used in
    # resolve_config_path().
    # This is the preferred default moving forwards, and the plan is to
    # eventually update resolve_config_path() to use AUGUR_SEARCH_PATHS.
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

    search_paths = [path.resolve() for path in search_paths if path.is_dir()]

    os.environ["AUGUR_SEARCH_PATHS"] = ":".join(map(str, search_paths))


class InvalidConfigError(Exception):
    pass


def resolve_config_path(path: str, defaults_dir: Optional[str] = None) -> Callable:
    """
    Resolve a relative *path* given in a configuration value. Will always try to
    resolve *path* after expanding wildcards with Snakemake's `expand` functionality.

    Returns the path for the first existing file, checked in the following order:
    1. relative to the analysis directory or workdir, usually given by ``--directory`` (``-d``)
    2. relative to *defaults_dir* if it's provided
    3. relative to the workflow's ``defaults/`` directory if *defaults_dir* is _not_ provided

    This behaviour allows a default configuration value to point to a default
    auxiliary file while also letting the file used be overridden either by
    setting an alternate file path in the configuration or by creating a file
    with the conventional name in the workflow's analysis directory.
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

        if os.path.exists(expanded_path):
            return expanded_path

        if defaults_dir:
            defaults_path = os.path.join(defaults_dir, expanded_path)
        else:
            # Special-case defaults/… for backwards compatibility with older
            # configs.  We could achieve the same behaviour with a symlink
            # (defaults/defaults → .) but that seems less clear.
            if path.startswith("defaults/"):
                defaults_path = os.path.join(workflow.basedir, expanded_path)
            else:
                defaults_path = os.path.join(workflow.basedir, "defaults", expanded_path)

        if os.path.exists(defaults_path):
            return defaults_path

        raise InvalidConfigError(indent(dedent(f"""\
            Unable to resolve the config-provided path {path!r},
            expanded to {expanded_path!r} after filling in wildcards.
            The workflow does not include the default file {defaults_path!r}.

            Hint: Check that the file {expanded_path!r} exists in your analysis
            directory or remove the config param to use the workflow defaults.
            """), " " * 4))

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
