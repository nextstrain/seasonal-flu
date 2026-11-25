from packaging import version
from augur.__version__ import __version__ as augur_version
import sys


def set_min_augur_version(min_augur_version: str):
    if version.parse(augur_version) < version.parse(min_augur_version):
      print("This pipeline needs a newer version of augur than you currently have...")
      print(f"Current augur version: {augur_version}. Minimum required: {min_augur_version}")
      sys.exit(1)
