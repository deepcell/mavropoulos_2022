"""Module containing string processing functionality."""

from __future__ import annotations

from itertools import count, cycle
from typing import Iterator, Sequence


def prefixes(s: str) -> Iterator[str]:
    """Yield all prefixes of the given string.

    Results include trivial prefixes like the empty string and s, because the empty string is a prefix of all strings
    and s is a prefix of itself.

    Args:
        s: any string

    Yields:
        prefixes from shortest to longest

    Examples:
        >>> list(prefixes(''))
        ['']
        >>> list(prefixes('A'))
        ['', 'A']
        >>> list(prefixes('AB'))
        ['', 'A', 'AB']
        >>> list(prefixes('ABCDEFG'))
        ['', 'A', 'AB', 'ABC', 'ABCD', 'ABCDE', 'ABCDEF', 'ABCDEFG']

    """
    for i in range(len(s) + 1):
        yield s[:i]


def non_trivial_prefixes(s: str) -> Iterator[str]:
    """Yield all non-trivial prefixes of the given string.

    Results exclude trivial prefixes like the empty string and s.  These are trivial because the empty string is a
    prefix of all strings and s is a prefix of itself.  So this command yields the equivalent of::

        list(prefixes(s))[1:-1]

    Args:
        s: any string

    Yields:
        prefixes from shortest to longest

    Examples:
        >>> list(non_trivial_prefixes(''))
        []
        >>> list(non_trivial_prefixes('A'))
        []
        >>> list(non_trivial_prefixes('AB'))
        ['A']
        >>> list(non_trivial_prefixes('ABCDEFG'))
        ['A', 'AB', 'ABC', 'ABCD', 'ABCDE', 'ABCDEF']

    """
    for i in range(1, len(s)):
        yield s[:i]


def suffixes(s: str) -> Iterator[str]:
    """Yield all suffixes of the given string.

    Results include trivial suffixes like the empty string and s, because the empty string is a suffix of all strings
    and s is a suffix of itself.

    Args:
        s: any string

    Yields:
        suffixes from shortest to longest

    Examples:
        >>> list(suffixes(''))
        ['']
        >>> list(suffixes('A'))
        ['', 'A']
        >>> list(suffixes('AB'))
        ['', 'B', 'AB']
        >>> list(suffixes('ABCDEFG'))
        ['', 'G', 'FG', 'EFG', 'DEFG', 'CDEFG', 'BCDEFG', 'ABCDEFG']

    """
    for i in range(len(s), -1, -1):
        yield s[i:]


def non_trivial_suffixes(s: str) -> Iterator[str]:
    """Yield all non-trivial suffixes of the given string.

    Results exclude trivial suffixes like the empty string and s.  These are trivial because the empty string is a
    suffix of all strings and s is a suffix of itself.  So this command yields the equivalent of::

        list(suffixes(s))[1:-1]

    Args:
        s: any string

    Yields:
        suffixes from shortest to longest

    Examples:
        >>> list(non_trivial_suffixes(''))
        []
        >>> list(non_trivial_suffixes('A'))
        []
        >>> list(non_trivial_suffixes('AB'))
        ['B']
        >>> list(non_trivial_suffixes('ABCDEFG'))
        ['G', 'FG', 'EFG', 'DEFG', 'CDEFG', 'BCDEFG']

    """
    for i in range(len(s) - 1, 0, -1):
        yield s[i:]


def rotations(s: str) -> Iterator[str]:
    """Yield rotations of the given string.

    A rotation of a string s is any other string s[i:]+s[:i] for any i in range(len(s)).  Rotations may not always be
    distinct.  If distinct rotations are required, then use::

        set(rotations(s))

    Args:
        s: any string

    Yields:
        rotations of s

    Examples:
        >>> list(rotations(''))
        ['']
        >>> list(rotations('A'))
        ['A']
        >>> list(rotations('AB'))
        ['AB', 'BA']
        >>> list(rotations('AA'))
        ['AA', 'AA']
        >>> list(rotations('ABC'))
        ['ABC', 'BCA', 'CAB']

    """
    if not s:
        yield ""

    for i in range(len(s)):
        yield s[i:] + s[:i]


def is_rotation(s1: str, s2: str) -> bool:
    """Test whether s1 and s2 are rotations of each other.

    A rotation of a string s is any other string s[i:]+s[:i] for any i in range(len(s)).

    Args:
        s1: any string
        s2: any string

    Returns:
        True if s1 and s2 are rotations of each other, else False

    Examples:
        >>> is_rotation('A', 'B')
        False
        >>> is_rotation('A', 'A')
        True
        >>> is_rotation('AB', 'BA')
        True
        >>> is_rotation('CAB', 'CBA')
        False
        >>> is_rotation('ABCDEFG', 'GABCDEF')
        True
        >>> all(is_rotation('ABCDEFGH', r) for r in rotations('ABCDEFGH'))
        True

    """
    if len(s1) != len(s2):
        return False

    return s1 in s2 + s2


def minimum_rotation(s: str) -> str:
    """Return the minimum rotation of string s using an inefficient algorithm.

    Args:
        s: input string

    Returns:
        minimum rotation of s

    Examples:
       >>> minimum_rotation('ABC')
       'ABC'
       >>> minimum_rotation('BCA')
       'ABC'
       >>> minimum_rotation('CAB')
       'ABC'

    """
    return min(rotations(s))


def rotate(s: str, i: int) -> str:
    """Rotate string by i characters.

    args:
        s: input string to rotate
        i: characters to rotate; i > 0 to rotate right, i < 0 to rotate left

    Examples:
        >>> rotate('ABC', 0)
        'ABC'
        >>> rotate('ABC', 1)
        'BCA'
        >>> rotate('ABC', 2)
        'CAB'
        >>> rotate('ABC', 3)
        'ABC'
        >>> rotate('ABC', 4)
        'BCA'
        >>> rotate('ABC', 5)
        'CAB'
        >>> rotate('ABC', 6)
        'ABC'
        >>> rotate('ABC', -1)
        'CAB'
        >>> rotate('ABC', -2)
        'BCA'
        >>> rotate('ABC', -3)
        'ABC'
        >>> rotate('ABC', -4)
        'CAB'
        >>> rotate('ABC', -5)
        'BCA'
        >>> rotate('ABC', -6)
        'ABC'

    """
    n = len(s)

    if not i or n <= 1:
        return s

    i = i % n

    if i > 0:
        s = s[i:] + s[:i]
    elif i < 0:
        i = abs(i)
        s = s[:i] + s[i:]

    return s


def cyclic_match_right(s: str, pattern: str) -> int:
    """Match pattern cyclically within s from the right and return number of matching characters.

    Args:
        s: input string to search
        pattern: cyclic pattern to find in s

    Returns:
        count of matching characters (0..len(s))

    Examples:
        >>> cyclic_match_right('AAAA', '')
        0
        >>> cyclic_match_right('AAAA', 'A')
        4
        >>> cyclic_match_right('AAAA', 'B')
        0
        >>> cyclic_match_right('AAABA', 'A')
        3
        >>> cyclic_match_right('ABAB', 'A')
        1
        >>> cyclic_match_right('ABAB', 'AB')
        4
        >>> cyclic_match_right('ABABA', 'AB')
        5

    """
    if not s or not pattern:
        return 0

    for a, b, match in zip(s, cycle(pattern), count()):
        if a != b:
            return match

    return match + 1


def cyclic_match_left(s: str, pattern: str) -> int:
    """Match pattern cyclically within s from the left and return number of matching characters.

    Args:
        s: input string to search
        pattern: cyclic pattern to find in s

    Returns:
        count of matching characters (0..len(s))

    Examples:
        >>> cyclic_match_left('AAAA', '')
        0
        >>> cyclic_match_left('AAAA', 'A')
        4
        >>> cyclic_match_left('AAAA', 'B')
        0
        >>> cyclic_match_left('ABAAA', 'A')
        3
        >>> cyclic_match_left('BABA', 'A')
        1
        >>> cyclic_match_left('ABAB', 'AB')
        4
        >>> cyclic_match_left('ABAB', 'BA')
        0
        >>> cyclic_match_left('ABABA', 'BA')
        5

    """
    if not s or not pattern:
        return 0

    for a, b, match in zip(reversed(s), cycle(reversed(pattern)), count()):
        if a != b:
            return match

    return match + 1


def trim_common_prefixes(strs: Sequence[str], max_trim: int | None = None) -> tuple[str, list[str]]:
    """Trim common prefixes from a list of strings.

    Args:
        strs: list of strings
        max_trim: maximum length common prefix to trim; default None is unlimited length

    Returns:
        tuple of common prefix and prefix-trimmed list of strings

    Examples:
        >>> trim_common_prefixes([])
        ('', [])
        >>> trim_common_prefixes(['AA'])
        ('', ['AA'])
        >>> trim_common_prefixes(['A','AB'])
        ('A', ['', 'B'])
        >>> trim_common_prefixes(['A','BA'])
        ('', ['A', 'BA'])
        >>> trim_common_prefixes(['A','AB','AAAB'])
        ('A', ['', 'B', 'AAB'])
        >>> trim_common_prefixes(['AB','ABA','ABAA'])
        ('AB', ['', 'A', 'AA'])
        >>> trim_common_prefixes(['AB','ABA','ABAA','C'])
        ('', ['AB', 'ABA', 'ABAA', 'C'])
        >>> trim_common_prefixes(['AAAAAAAAA','AAAAAAAAA','AAAAAAAAA','AAAAAAAAA'])
        ('AAAAAAAAA', ['', '', '', ''])
        >>> trim_common_prefixes(['AAAAAAAAAB','AAAAAAAAAC','AAAAAAAAAD','AAAAAAAAA'])
        ('AAAAAAAAA', ['B', 'C', 'D', ''])
        >>> trim_common_prefixes(['AAAAAB', 'AAAAB'], max_trim=2)
        ('AA', ['AAAB', 'AAB'])

    """
    if not isinstance(strs, list):
        strs = list(strs)

    if len(strs) < 2:
        return "", strs

    istrs = iter(strs)
    first = next(istrs)
    trim = len(first)

    for s in istrs:
        trim = min(trim, len(s))
        i = 0
        while i < trim and first[i] == s[i]:
            i += 1
        trim = i
        if not trim:
            break

    if not trim:
        return "", strs

    if max_trim is not None and trim > max_trim:
        trim = max_trim

    return first[:trim], [s[trim:] for s in strs]


def trim_common_suffixes(strs: Sequence[str], max_trim: int | None = None) -> tuple[str, list[str]]:
    """Trim common suffixes from a list of strings.

    Args:
        strs: list of strings
        max_trim: maximum length common suffix to trim; default None is unlimited length

    Returns:
        tuple of common suffix and suffix-trimmed list of strings

    Examples:
        >>> trim_common_suffixes([])
        ('', [])
        >>> trim_common_suffixes(['AA'])
        ('', ['AA'])
        >>> trim_common_suffixes(['A','BA'])
        ('A', ['', 'B'])
        >>> trim_common_suffixes(['A','AB'])
        ('', ['A', 'AB'])
        >>> trim_common_suffixes(['A','BA','BAAA'])
        ('A', ['', 'B', 'BAA'])
        >>> trim_common_suffixes(['BA','ABA','AABA'])
        ('BA', ['', 'A', 'AA'])
        >>> trim_common_suffixes(['AB','ABA','ABAA','C'])
        ('', ['AB', 'ABA', 'ABAA', 'C'])
        >>> trim_common_suffixes(['AAAAAAAAA','AAAAAAAAA','AAAAAAAAA','AAAAAAAAA'])
        ('AAAAAAAAA', ['', '', '', ''])
        >>> trim_common_suffixes(['BAAAAAAAAA','CAAAAAAAAA','DAAAAAAAAA','AAAAAAAAA'])
        ('AAAAAAAAA', ['B', 'C', 'D', ''])

    """
    if not isinstance(strs, list):
        strs = list(strs)

    if len(strs) < 2:
        return "", strs

    istrs = iter(strs)
    first = next(istrs)
    trim = len(first)

    for s in istrs:
        trim = min(trim, len(s))
        i = -1
        while i >= -trim and first[i] == s[i]:
            i -= 1
        trim = -i - 1
        if not trim:
            break

    if not trim:
        return "", strs

    if max_trim is not None and trim > max_trim:
        trim = max_trim

    return first[-trim:], [s[:-trim] for s in strs]


def strip_prefix(s: str, prefix: str) -> str:
    """Strip a given prefix from a string, if it exists, and then return the result.

    Args:
         s: string with or without prefix
         prefix: prefix to strip from string

    Returns:
        string with prefix removed or original string if prefix was not present

    Examples:
        >>> strip_prefix('aaaaaaaaaa', 'a')
        'aaaaaaaaa'
        >>> strip_prefix('Sample_Sample_ample', 'Sample_')
        'Sample_ample'

    """
    return s[len(prefix) :] if s.startswith(prefix) else s


def strip_suffix(s: str, suffix: str) -> str:
    """Strip a given suffix from a string, if it exists, and then return the result.

    Args:
         s: string with or without suffix
         suffix: suffix to strip from string

    Returns:
        string with suffix removed or original string if suffix was not present

    Examples:
        >>> strip_suffix('aaaaaaaaaa', 'a')
        'aaaaaaaaa'
        >>> strip_suffix('Sample_Sample_ample', '_ample')
        'Sample_Sample'
        >>> strip_suffix('Sample_8009.vcf', '.vcf')
        'Sample_8009'

    """
    return s[: -len(suffix)] if s.endswith(suffix) else s


def to_optional_bool(value: str | None) -> bool | None:
    """Convert a string to a bool or None.

    'true'    (case insensitive) -> True
    'false'   (case insensitive) -> False
    'unknown' (case insensitive) -> None
    None or ''                   -> None
    otherwise raises ValueError

    It should be noted that the above conversion rules give 'unknown' a Falsy value. This may need to be revisited if a
    a clearer distinction is needed between 'false' and 'unknown' values.

    Args:
        value: string to cast

    Returns:
        True, False, or None

    Raises:
        ValueError: unknown boolean value

    Examples:
        >>> to_optional_bool('true')
        True
        >>> to_optional_bool('TRUE')
        True
        >>> to_optional_bool('False')
        False
        >>> to_optional_bool('unknown') is None
        True
        >>> to_optional_bool(None) is None
        True
        >>> to_optional_bool('') is None
        True
        >>> to_optional_bool('error')
        Traceback (most recent call last):
            ...
        ValueError: unknown boolean value: error

    """
    if not value:
        return None

    lvalue = value.lower()

    if lvalue == "true":
        return True
    elif lvalue == "false":
        return False
    elif lvalue == "unknown":
        return None

    raise ValueError(f"unknown boolean value: {value}")
