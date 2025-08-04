from shlex import (
    quote as shquote,       # shquote() is used in this file and also other workflow files
    split as shsplitwords,
)

def shquotewords(s: str) -> str:
    """
    Split string *s* into (POSIX) shell words, quote each word, and join them
    back into a string.

    This is suitable for properly quoting multi-word, user-defined values which
    should follow shell quoting and escaping semantics (e.g. to allow spaces in
    single words) but not allow shell features like variable interpolation,
    command substition, redirection, piping, etc.

    For example, quote a query string used as input to augur filter like this:
    f"--query {shquote(query)}".

    See usage in https://github.com/nextstrain/ncov for more examples.
    """
    return " ".join(shquote(word) for word in shsplitwords(s))
