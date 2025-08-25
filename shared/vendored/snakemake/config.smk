"""
Shared functions to be used within a Snakemake workflow for parsing
workflow configs.
"""
import os.path
from collections.abc import Callable
from snakemake.io import Wildcards
from typing import Optional
from textwrap import dedent, indent


class InvalidConfigError(Exception):
    pass


def resolve_config_path(path: str, defaults_dir: Optional[str] = None) -> Callable[[Wildcards], str]:
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
