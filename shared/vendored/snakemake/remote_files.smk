"""
Helper functions to set-up storage plugins for remote inputs/outputs. See the
docstring of `path_or_url` for usage instructions.

The errors raised by storage plugins are often confusing. For instance, a HTTP
404 error will result in a `MissingInputException` with little hint as to the
underlying issue. S3 credentials errors are similarly confusing and we attempt
to check these ourselves to improve UX here.
"""
import os
import socket
from urllib.parse import urlparse
from snakemake.io import get_flag_value, AnnotatedString
from snakemake.logging import logger


def is_online(host="8.8.8.8", port=53, timeout=3) -> bool:
    """
    Check if workflow has network connection.
    Based on <https://stackoverflow.com/a/33117579>
    Host: 8.8.8.8 (google-public-dns-a.google.com)
    OpenPort: 53/tcp
    Service: domain (DNS/TCP)
    """
    try:
        socket.setdefaulttimeout(timeout)
        socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect((host, port))
        return True
    except socket.error as err:
        logger.debug(f"Assuming no internet connection with error: {err}")
        return False


# Determine whether workflow should try to fetch or revalidate remote input
# files based on prioritization order:
#     1. Environment variable `NEXTSTRAIN_REMOTE_FILES`
#     2. Per workflow config parameter `use_remote_files`
#     3. Whether the workflow has network connection
USE_REMOTE_FILES = True
if "NEXTSTRAIN_USE_REMOTE_FILES" in os.environ:
    USE_REMOTE_FILES = bool(os.environ["NEXTSTRAIN_USE_REMOTE_FILES"])
elif "use_remote_files" in config:
    USE_REMOTE_FILES = bool(config["use_remote_files"])
else:
    USE_REMOTE_FILES = is_online()


# Keep a list of known public buckets, which we'll allow uncredentialled (unsigned) access to
# We could make this config-definable in the future
PUBLIC_BUCKETS = set(['nextstrain-data'])

# Keep track of registered storage plugins to enable reuse
_storage_registry = {}

class RemoteFilesMissingCredentials(Exception):
    pass

def _storage_s3(*, bucket, keep_local, retries) -> snakemake.storage.StorageProviderProxy:
    """
    Registers and returns an instance of snakemake-storage-plugin-s3. Typically AWS
    credentials are required for _any_ request however we allow requests to known
    public buckets (see `PUBLIC_BUCKETS`) to be unsigned which allows for a nice user
    experience in the common case of downloading inputs from s3://nextstrain-data.

    The intended behaviour for various (S3) URIs supplied to `path_or_url` is:

    |          | S3 buckets                 | credentials present | credentials missing |
    |----------|----------------------------|---------------------|---------------------|
    | download | private / private + public | signed              | Credentials Error   |
    |          | public                     | signed              | unsigned            |
    | upload   | private / private + public | signed              | Credentials Error   |
    |          | public                     | signed              | AccessDenied Error  |
    """
    # If the bucket is public then we may use an unsigned request which has the nice UX
    # of not needing credentials to be present. If we've made other signed requests _or_
    # credentials are present then we just sign everything. This has implications for upload:
    # if you attempt to upload to a public bucket without credentials then we allow that here
    # and you'll get a subsequent `AccessDenied` error when the upload is attempted.
    if bucket in PUBLIC_BUCKETS and \
        "s3_signed" not in _storage_registry and \
        ("s3_unsigned" in _storage_registry or not _aws_credentials_present()):

        if provider:=_storage_registry.get('s3_unsigned', None):
            return provider

        from botocore import UNSIGNED # dependency of snakemake-storage-plugin-s3
        storage s3_unsigned:
            provider="s3",
            signature_version=UNSIGNED,
            retries=retries,
            keep_local=keep_local,

        _storage_registry['s3_unsigned'] = storage.s3_unsigned
        return _storage_registry['s3_unsigned']

    # Resource fetched/uploaded via a signed request, which will require AWS credentials
    if provider:=_storage_registry.get('s3_signed', None):
        return provider

    # Enforce the presence of credentials to paper over <https://github.com/snakemake/snakemake/issues/3663>
    if not _aws_credentials_present():
        raise RemoteFilesMissingCredentials()

    # the tag appears in the local file path, so reference 'signed' to give a hint about credential errors
    storage s3_signed:
        provider="s3",
        retries=retries,
        keep_local=keep_local,

    _storage_registry['s3_signed'] = storage.s3_signed
    return _storage_registry['s3_signed']

def _aws_credentials_present() -> bool:
    import boto3 # dependency of snakemake-storage-plugin-s3
    session = boto3.Session()
    creds = session.get_credentials()
    return creds is not None

def _storage_http(*, keep_local, retries) -> snakemake.storage.StorageProviderProxy:
    """
    Registers and returns an instance of snakemake-storage-plugin-http
    """
    if provider:=_storage_registry.get('http', None):
        return provider

    storage:
        provider="http",
        allow_redirects=True,
        supports_head=True,
        keep_local=keep_local,
        retries=retries,

    _storage_registry['http'] = storage.http
    return _storage_registry['http']


def path_or_url(uri, *, keep_local=True, retries=2) -> str | AnnotatedString:
    """
    Intended for use in Snakemake inputs / outputs to transparently use remote
    resources. Returns the URI wrapped by an applicable storage plugin. Local
    filepaths will be returned unchanged.

    For example, the following rule will download inputs from HTTPs and upload
    the output to S3:

        rule filter:
            input:
                sequences = path_or_url("https://data.nextstrain.org/..."),
                metadata = path_or_url("https://data.nextstrain.org/..."),
            output:
                sequences = path_or_url("s3://...")
            shell:
                r'''
                augur filter \
                    --sequences {input.sequences:q} \
                    --metadata {input.metadata:q} \
                    --metadata-id-columns accession \
                    --output-sequences {output.sequences:q}
                '''

    If *keep_local* is True (the default) then downloaded/uploaded files will
    remain in `.snakemake/storage/`. The presence of a previously downloaded
    file (via `keep_local=True`) does not guarantee that the file will not be
    re-downloaded if the storage plugin decides the local file is out of date.

    Depending on the *uri* authentication may be required. See the specific
    helper functions (such as `_storage_s3`) for more details.

    See <https://snakemake.readthedocs.io/en/stable/snakefiles/storage.html> for
    more information on Snakemake storage plugins. Note: various snakemake
    plugins will be required depending on the URIs provided.
    """
    info = urlparse(uri)

    if info.scheme=='': # local
        return uri      # no storage wrapper

    if info.scheme=='s3':
        try:
            return _local_or_storage(_storage_s3(bucket=info.netloc, keep_local=keep_local, retries=retries)(uri))
        except RemoteFilesMissingCredentials as e:
            raise Exception(f"AWS credentials are required to access {uri!r}") from e

    if info.scheme=='https':
        return _local_or_storage(_storage_http(keep_local=keep_local, retries=retries)(uri))
    elif info.scheme=='http':
        raise Exception(f"HTTP remote file support is not implemented in nextstrain workflows (attempting to access {uri!r}).\n"
            "Please use an HTTPS address instead.")

    if info.scheme in ['gs', 'gcs']:
        raise Exception(f"Google Storage is not yet implemented for nextstrain workflows (attempting to access {uri!r}).\n"
            "Please get in touch if you require this functionality and we can add it to our workflows")

    raise Exception(f"Input address {uri!r} (scheme={info.scheme!r}) is from a non-supported remote")


def _local_or_storage(wrapped: AnnotatedString) -> str | AnnotatedString:
    """
    Given a storage-wrapped input (an AnnotatedString carrying a `storage_object`
    flag), transparently fall back to the local cache copy as determined by
    USE_REMOTE_FILES:
      * If USE_REMOTE_FILES == True, then return object unchanged and let the
        allow normal Snakemake retrieval, revalidation, or
        MissingInputException behaviour.
      * If USE_REMOTE_FILES == False and cache does not exist, raise an error.
      * If USE_REMOTE_FILES == False and cache does exist, then return the
        plain local path. Snakemake then treats it as an ordinary local input
        (existence via `os.path.exists`, no network), allowing the workflow to
        run offline from cache.
    """
    so = get_flag_value(wrapped, "storage_object")
    if so is None:
        raise Exception("_local_or_storage must be called with a storage object AnnotatedString. Provided argument: ", wrapped)

    # Return StorageObject as-is since remote files are expected and allow it
    # to fetch or revalidate cached files as needed.
    # Any connection errors will be raised downstream
    if USE_REMOTE_FILES:
        return wrapped

    # The keep_local flag is managed by the workflow author and this cannot be
    # changed by the user. So output warning for user and return the StorageObject as-is.
    if not so.keep_local:
        logger.warning(
            f"Unable to use cached copy of {so.query!r} because the workflow "
            "author has opted out of keeping a local copy for this input. "
            "Please contact the workflow author to enable use of local cached files."
        )
        return wrapped

    local = so.local_path()
    if not local.exists():
        raise Exception(
            f"Unable to find local cached file for {so.query!r}. "
            "Please make sure you have run the workflow at least once with "
            "internet connection and config.use_remote_files set to True to "
            "to fetch and save the remote file locally.")

    logger.warning(
        f"Using the cached copy of {so.query} at {local} without revalidation. "
        "Set config.use_remote_files to True if you want to force revalidation."
    )

    return str(local)
