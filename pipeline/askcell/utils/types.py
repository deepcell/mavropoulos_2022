"""Utilities to facilitate managing string-keyed nested-mapping data structures and schemas."""

from __future__ import annotations

import logging
import math

from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import date, datetime
from types import MappingProxyType
from typing import (
    IO,
    Any,
    Generic,
    Iterable,
    Iterator,
    Mapping,
    MutableMapping,
    Sequence,
    TypeVar,
    Union,
)

import numpy as np

from .. import Path, PathLike, parse_path
from .expression import _join_keys


Value = Any
ImmutableStrMapping = Mapping[str, Any]
MutableStrMapping = MutableMapping[str, Any]
StrMapping = Union[MutableStrMapping, ImmutableStrMapping]
MemoCache = dict[int, Any]

SerializationType = TypeVar("SerializationType")
NativeType = TypeVar("NativeType")


class SchemaType(ABC, Generic[SerializationType, NativeType]):
    """Base class for StrMapping schema types.

    Attributes:
        default: default value, defaults to None
        required: value can be None if False, defaults to False

    ``default`` provides a default value when applying a schema to a data fragment or building default data from a
    schema. ``required`` is used for validation and prohibits None values if set.

    """

    def __init__(self, **kwargs: Any) -> None:
        """Parse attributes from a schema definition.

        Args:
            kwargs: schema attributes

        Supported key-worded arguments are ``default`` and ``required``. They are used to set corresponding attributes.

        """
        self.required = kwargs.get("required", False)
        self.default = self.parse(kwargs.get("default"))

    @abstractmethod
    def parse(self, value: SerializationType | None) -> NativeType | None:
        """Parse value of serialization type to native type."""

    @abstractmethod
    def unparse(self, value: NativeType | None) -> SerializationType | None:
        """Unparse value of native type to serialization type."""

    @abstractmethod
    def validate(self, value: NativeType | None) -> bool:
        """Validate value is a valid native type."""

    @abstractmethod
    def parse_string(self, value: str) -> NativeType | str | None:
        """Parse string to native type."""

    @abstractmethod
    def normalize(self, value: Any) -> Any:
        """Normalize a compatible value to native type."""

    def __repr__(self) -> str:
        """Return string representation of this type."""
        return type(self).__name__

    def __eq__(self, other: object) -> bool:
        """Return True if self and other are equal."""
        # make it clear to the type checker that other must be a SchemaType to be considered below
        if not isinstance(other, SchemaType):
            return False
        return type(self) is type(other) and self.required == other.required and self.default == other.default


AtomTypeTuple = (str, int, bool, float, date, datetime, Path, SchemaType, type(None))
AtomType = Union[str, int, bool, float, date, datetime, Path, type[SchemaType[Any, Any]], None]

# Type groups for normalization
IntTypes = Union[int, np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64, np.bool_]
FloatTypes = Union[
    int, np.int8, np.int16, np.int32, np.uint8, np.uint16, np.uint32, np.bool_, float, np.float32, np.float64
]
BoolTypes = Union[bool, np.bool_]
DateTypes = Union[str, date]
DateTimeTypes = Union[str, date, datetime]

# Type tuples for normalization
IntTypesTuple = (int, np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64, np.bool_)
FloatTypesTuple = (
    int,
    np.int8,
    np.int16,
    np.int32,
    np.uint8,
    np.uint16,
    np.uint32,
    np.bool_,
    float,
    np.float32,
    np.float64,
)
BoolTypesTuple = (bool, np.bool_)
DateTypesTuple = (str, date)
DateTimeTypesTuple = (str, date, datetime)


class StrSchemaType(SchemaType[str, str]):
    """Schema type for converting serialized strings to native strings and vice-versa.

    Attributes:
        default: default value, defaults to None
        required: value can be None if False, defaults to False

    """

    def parse(self, value: str | None) -> str | None:
        """Parse value of serialization type to native type."""
        return value

    def unparse(self, value: str | None) -> str | None:
        """Unparse value of native type to serialization type."""
        return value

    def validate(self, value: str | None) -> bool:
        """Validate value is a valid native type."""
        return (not self.required and value is None) or type(value) is str

    def parse_string(self, value: str) -> str | None:
        """Parse string to native type."""
        return value or None

    def normalize(self, value: str | None) -> str | None:
        """Normalize a compatible value to native type."""
        return value


class IntSchemaType(SchemaType[int, int]):
    """Schema type for converting serialized integers to native integers and vice-versa.

    This class also supports converting strings and numpy integers to int.

    Attributes:
        default: default value, defaults to None
        required: value can be None if False, defaults to False

    """

    def parse(self, value: int | None) -> int | None:
        """Parse value of serialization type to native type."""
        return value

    def unparse(self, value: int | None) -> int | None:
        """Unparse value of native type to serialization type."""
        return value

    def validate(self, value: int | None) -> bool:
        """Validate value is a valid native type."""
        return (not self.required and value is None) or type(value) is int

    def parse_string(self, value: str) -> int | str | None:
        """Parse string to native type."""
        if not value:
            return None

        try:
            return int(value)
        except (TypeError, ValueError):
            return value

    def normalize(self, value: IntTypes | None) -> int | None:
        """Normalize a compatible value to native type."""
        if isinstance(value, IntTypesTuple):
            return int(value)
        return value


class FloatSchemaType(SchemaType[float, float]):
    """Schema type for converting serialized floats to native floats and vice-versa.

    This class also supports converting strings, numpy floats, numpy integers, and int to float.

    Attributes:
        default: default value, defaults to None
        required: value can be None if False, defaults to False

    """

    def parse(self, value: float | None) -> float | None:
        """Parse value of serialization type to native type."""
        return value

    def unparse(self, value: float | None) -> float | None:
        """Unparse value of native type to serialization type."""
        return value

    def validate(self, value: float | None) -> bool:
        """Validate value is a valid native type."""
        return (not self.required and value is None) or type(value) is float

    def parse_string(self, value: str) -> float | str | None:
        """Parse string to native type."""
        if not value:
            return None

        try:
            return float(value)
        except (TypeError, ValueError):
            return value

    def normalize(self, value: FloatTypes | None) -> float | None:
        """Normalize a compatible value to native type."""
        if isinstance(value, FloatTypesTuple):
            return float(value)
        return value


class BoolSchemaType(SchemaType[bool, bool]):
    """Schema type for converting serialized booleans to native booleans and vice-versa.

    This class also supports converting strings and numpy's bool to bool.

    Attributes:
        default: default value, defaults to None
        required: value can be None if False, defaults to False

    """

    def parse(self, value: bool | None) -> bool | None:
        """Parse value of serialization type to native type."""
        return value

    def unparse(self, value: bool | None) -> bool | None:
        """Unparse value of native type to serialization type."""
        return value

    def validate(self, value: bool | None) -> bool:
        """Validate value is a valid native type."""
        return (not self.required and value is None) or type(value) is bool

    def parse_string(self, value: str) -> bool | str | None:
        """Parse string to native type."""
        if not value:
            return None
        elif not isinstance(value, str):
            return value

        lower = value.lower()
        if lower == "false":
            return False
        elif lower == "true":
            return True
        return value

    def normalize(self, value: BoolTypes | None) -> bool | None:
        """Normalize a compatible value to native type."""
        # FIXME: Should bool normalize also evaluate all non-None instances as truth values?
        # e.g. set_value(schema, data, 'a.b', a or b) looks like a boolean expression
        if isinstance(value, BoolTypesTuple):
            return bool(value)
        return value


class DateSchemaType(SchemaType[date, date]):
    """Schema type for converting serialized dates to native dates and vice-versa.

    This class also supports converting strings to datetime.date.

    Attributes:
        default: default value, defaults to None
        required: value can be None if False, defaults to False

    """

    def parse(self, value: date | None) -> date | None:
        """Parse value of serialization type to native type."""
        return value

    def unparse(self, value: date | None) -> date | None:
        """Unparse value of native type to serialization type."""
        return value

    def validate(self, value: date | None) -> bool:
        """Validate value is a valid native type."""
        return (not self.required and value is None) or type(value) is date

    def parse_string(self, value: str) -> date | str | None:
        """Parse string to native type."""
        if not value:
            return None

        try:
            return date.fromisoformat(value)
        except (TypeError, ValueError):
            return value

    def normalize(self, value: DateTypes) -> date | str | None:
        """Normalize a compatible value to native type."""
        if value and isinstance(value, str):
            return self.parse_string(value)
        return value


class DateTimeSchemaType(SchemaType[datetime, datetime]):
    """Schema type for converting serialized datetimes to native datetimes and vice-versa.

    This class also supports converting strings and datetime.date to datetime.datetime.

    Attributes:
        default: default value, defaults to None
        required: value can be None if False, defaults to False

    """

    def parse(self, value: datetime | None) -> datetime | None:
        """Parse value of serialization type to native type."""
        return value

    def unparse(self, value: datetime | None) -> datetime | None:
        """Unparse value of native type to serialization type."""
        return value

    def validate(self, value: datetime | None) -> bool:
        """Validate value is a valid native type."""
        return (not self.required and value is None) or type(value) is datetime

    def parse_string(self, value: str) -> datetime | str | None:
        """Parse string to native type."""
        if not value:
            return None

        try:
            return datetime.fromisoformat(value)
        except (TypeError, ValueError):
            return value

    def normalize(self, value: DateTimeTypes | None) -> datetime | str | None:
        """Normalize a compatible value to native type."""
        if isinstance(value, datetime):
            return value
        elif isinstance(value, date):
            return datetime(value.year, value.month, value.day)
        elif value and isinstance(value, str):
            return self.parse_string(value)
        return value


class PathSchemaType(SchemaType[str, Path]):
    """Schema type for converting serialized strings to native Paths and vice-versa.

    Empty strings are treated as missing values (:code:`None`).  Otherwise, :code:`Path('')` becomes a
    synonym for :code:`Path('.')`

    Attributes:
        default: default value, defaults to None
        required: value can be None if False, defaults to False

    """

    def parse(self, value: str | None) -> Path | None:
        """Parse value of serialization type to native type."""
        if not isinstance(value, str):
            return value
        return parse_path(value) if value else None

    def unparse(self, value: Path | None) -> str | None:
        """Unparse value of native type to serialization type."""
        if not isinstance(value, Path):
            return value
        return str(value)

    def validate(self, value: Path | None) -> bool:
        """Validate value is a valid native type."""
        # may be None or any Path subtype (e.g. PosixPath)
        return (not self.required and value is None) or isinstance(value, Path)

    def parse_string(self, value: str) -> Path | None:
        """Parse string to native type."""
        if not isinstance(value, str):
            return value
        return parse_path(value) if value else None

    def normalize(self, value: PathLike | None) -> Path | None:
        """Normalize a compatible value to native type."""
        if not isinstance(value, str):
            return value
        return parse_path(value) if value else None


ATOM_TYPES: dict[str, type[SchemaType[Any, Any]]] = {
    "str": StrSchemaType,
    "int": IntSchemaType,
    "float": FloatSchemaType,
    "bool": BoolSchemaType,
    "date": DateSchemaType,
    "datetime": DateTimeSchemaType,
    "path": PathSchemaType,
}


class ListSchemaType(SchemaType[Any, Any]):
    """Schema type for list of another atomic element type.

    Individual list elements are always required, and no list element is allowed to be None. Nested lists are not
    supported.

    Attributes:
        default: default sequence, defaults to None
        element_type: atomic SchemaType such as StrSchemaType, IntSchemaType, etc.
        min_length: minimum length, defaults to 0
        max_length: maximum length, defaults to None (unbounded)
        required: value can be None if False (ignores length constraints), defaults to False

    ``required``, ``min_length``, and ``max_length`` are used to validate a sequence. ``element_type`` is used to parse,
    unparse, normalize, and validate elements of a sequence.

    """

    def __init__(self, *, type: str, **kwargs: Any) -> None:
        """Parse attributes from a schema definition.

        Args:
            type: element schema type. One of {str, int, float, bool, date, datetime, Path}
            kwargs: schema attributes

        Supported key-worded arguments are ``default``, ``required``, ``min_length``, and ``max_length``. They are used
        to set corresponding attributes. ``length`` if provided sets the attributes ``min_length`` and ``max_length``.
        It must be consistent with the arguments ``min_length`` and ``max_length`` or an error is raised.

        """
        length = kwargs.pop("length", None)
        default = kwargs.pop("default", None)
        self.min_length = kwargs.pop("min_length", length or 0)
        self.max_length = kwargs.pop("max_length", length)
        self.required = kwargs.pop("required", False)

        if length is not None and not (self.min_length == length == self.max_length):
            raise ValueError("schema value definition has contradictory lengths")

        self.element_type = ATOM_TYPES[type](required=True, **kwargs)
        self.default = self.parse(default)

    def parse(self, values: Sequence[Any] | None) -> Sequence[Any] | None:
        """Parse values of serialization type to native type."""
        if values is None:
            return None
        elif not values:
            return []
        elif isinstance(values, Iterable) and not isinstance(values, str):
            return [self.element_type.parse(value) for value in values]
        return values

    def unparse(self, values: Sequence[Any] | None) -> Sequence[Any] | None:
        """Unparse values of native type to serialization type."""
        if values is None:
            return None
        elif not values:
            return []
        return [self.element_type.unparse(value) for value in values]

    def validate(self, values: Sequence[Any] | None) -> bool:
        """Validate length constraints and values are valid native type."""
        if values is None:
            return not self.required

        max_length = self.max_length
        if max_length is None:
            max_length = math.inf

        if not (self.min_length <= len(values) <= max_length):
            return False

        return all(self.element_type.validate(value) for value in values)

    def parse_string(self, value: str) -> Sequence[Any] | None:
        """Parse string to native type."""
        if not value:
            return None
        elif not isinstance(value, str):
            return value
        return [self.element_type.parse_string(v.strip()) for v in value.split(",")]

    def normalize(self, values: Sequence[Any] | None) -> Sequence[Any] | None:
        """Normalize compatible values to native type."""
        # prevent negating a numpy array from raising an error
        if isinstance(values, np.ndarray):
            values = values.tolist()

        if values is None:
            return None
        elif not values:
            return []
        elif isinstance(values, Iterable) and not isinstance(values, str):
            return [self.element_type.normalize(value) for value in values]
        return values

    def __repr__(self) -> str:
        """Return string representation of this type."""
        return f"ListSchemaType[{self.element_type!r}]"

    def __eq__(self, other: object) -> bool:
        """Return True if self and other are equal."""
        # make it clear to the type checker that other must be a ListSchemaType to be considered below
        if not isinstance(other, ListSchemaType):
            return False

        return (
            self.required == other.required
            and self.default == other.default
            and self.element_type == other.element_type
            and self.min_length == other.min_length
            and self.max_length == other.max_length
        )


def parse_schema(data: StrMapping) -> MutableStrMapping:
    """Parse mappings with schema definitions.

    This function takes string-keyed nested mapping data and converts schema definition values
    into in-memory representations using SchemaType classes.

    Args:
        data: string-keyed nested mappings

    Examples:
        >>> d = {'a': {'type': 'str'}}
        >>> parse_schema(d)
        {'a': StrSchemaType}

        >>> d = {'a': {'b': {'type': 'str'}}}
        >>> parse_schema(d)
        {'a': {'b': StrSchemaType}}

    """
    if not isinstance(data, (dict, MappingProxyType)):
        raise TypeError(f"invalid schema: found {type(data)!r}")

    result: MutableStrMapping = {}
    for key, value in data.items():
        if not isinstance(value, (dict, MappingProxyType)):
            raise TypeError(f"invalid schema: found {type(data)!r}")
        elif all(isinstance(v, (dict, MappingProxyType)) for v in value.values()):
            result[key] = parse_schema(value)
        else:
            result[key] = parse_schema_type(**value)

    return result


def parse_schema_type(*, type: str, **kwargs: Any) -> SchemaType[Any, Any]:
    """Parse schema type.

    Args:
        type: schema type
        kwargs: schema attributes

    Returns:
        Python type.

    Examples:
        >>> parse_schema_type(type='str')
        StrSchemaType
        >>> parse_schema_type(type='list[str]')
        ListSchemaType[StrSchemaType]

    """
    if type.startswith("list[") and type.endswith("]"):
        # drop list[ and ] at the ends
        return ListSchemaType(type=type[5:-1], **kwargs)
    else:
        return ATOM_TYPES[type](**kwargs)


def unparse_data(schema: StrMapping, data: StrMapping) -> MutableStrMapping:
    """Replace values with ones that are TOML-valid.

    Before mappings can be written as TOML invalid values must be cast to valid values.
    None values are outright removed as they are not supported.

    Args:
        schema: atomized string-keyed mappings
        data:   string-keyed schema mappings

    Returns:
        Updated atomized string-keyed value mappings.

    Examples:
        >>> schema = {
        ...    'a': FloatSchemaType(),
        ...    'b': IntSchemaType(),
        ... }
        >>> data = {'a': 9.0, 'b': 1}
        >>> unparse_data(schema, data)
        {'a': 9.0, 'b': 1}

        >>> schema = {'a': PathSchemaType()}
        >>> data = {'a': parse_path('/foo/bar')}
        >>> unparse_data(schema, data)
        {'a': '/foo/bar'}

    """
    result = {}
    for key, value in data.items():
        if value is None:
            continue

        if isinstance(value, (dict, MappingProxyType)):
            result[key] = unparse_data(schema[key], data[key])
        else:
            result[key] = schema[key].unparse(value)

    return result


def build_data_from_schema(schema: StrMapping) -> MutableStrMapping:
    """Build a data-only representation of a schema based on types and default values.

    Args:
        schema: string-keyed schema mappings

    Examples:
        >>> d = {'a': StrSchemaType()}
        >>> build_data_from_schema(d)
        {'a': None}

        >>> d = {'a': StrSchemaType(default='foo')}
        >>> build_data_from_schema(d)
        {'a': 'foo'}

    """
    return {
        key: value.default if isinstance(value, SchemaType) else build_data_from_schema(value)
        for key, value in schema.items()
    }


def apply_schema(schema: StrMapping, data: StrMapping) -> MutableStrMapping:
    """Apply schema definitions.

    Enforces schema definitions on separate atomized mapping object.
    Values that do not exist in atomized mappings are created with defaults (or None).
    PathLike elements defined as Path type are converted to Path objects.

    Args:
        schema: string-keyed schema mappings
        data:   atomized string-keyed mappings

    Returns:
        Updated atomized string-keyed mappings.

    Examples:
        >>> schema = {'a': StrSchemaType()}
        >>> data   = {}
        >>> apply_schema(schema, data)
        {'a': None}

        Schema default values are used to fill in missing data elements:

        >>> schema = {'a': StrSchemaType(required=True, default='foo')}
        >>> data   = {}
        >>> apply_schema(schema, data)
        {'a': 'foo'}

        Schema Path type declarations are enforced:

        >>> schema = {'a': PathSchemaType(required=True, default=parse_path('/foo/bar'))}
        >>> data = apply_schema(schema, data)
        >>> data
        {'a': PosixPath('/foo/bar')}

    """
    result = {}
    for key, stype in schema.items():
        if isinstance(stype, SchemaType):
            result[key] = stype.parse(data[key]) if key in data else stype.default
        else:
            result[key] = apply_schema(stype, data.get(key, {}))
    return result


@dataclass(init=False, frozen=True)
class ValidationError:
    """Schema validation failure error.

    Attributes:
        msg:   explanation of why schema validation failed
        path:  keys traversed in data and schema where validation failed
        value: data value that failed schema validation
        level: python logging level. One of {logging.CRITICAL, logging.ERROR, etc.}

    """

    msg: str
    path: str
    value: Any
    level: int

    def __init__(self, msg: str, path: str | list[str], value: Any = None, level: int = logging.ERROR) -> None:
        """Set frozen attributes.

        Args:
            msg:   explanation of why schema validation failed
            path:  keys traversed in data and schema where validation failed
            value: data value that failed schema validation
            level: python logging level. One of {logging.CRITICAL, logging.ERROR, etc.}

        """
        set_frozen_attr = object.__setattr__
        set_frozen_attr(self, "msg", msg)
        set_frozen_attr(self, "path", path if isinstance(path, str) else _join_keys(path))
        set_frozen_attr(self, "value", value)
        set_frozen_attr(self, "level", level)

    def __str__(self) -> str:
        """Return summary message."""
        return f"{logging.getLevelName(self.level)}: {self.msg}"


@dataclass(init=False)
class ValidationException(ValueError):
    """:class:`ValidtationError` container and ValueError subclass.

    Attributes:
        errors: container of :class:`ValidationError`

    """

    def __init__(self, errors: Iterable[ValidationError]) -> None:
        """Set exception message and save validation errors.

        Args:
            errors:  validation errors

        """
        super().__init__("Schema error")
        self.errors = list(errors)

    def log_to(self, logger: logging.Logger) -> None:
        """Log each :class:`ValidationError` in errors.

        Args:
            logger: logging.Logger instance

        """
        for error in self.errors:
            logger.log(error.level, error.msg)

    def write_to(self, file: IO[str]) -> None:
        """Write each :class:`ValidationError` in errors to ``file``.

        Args:
            file: output I/O stream

        """
        for error in self.errors:
            file.write(f"{error}\n")

    def __bool__(self) -> bool:
        """Return True if errors is not empty, False otherwise."""
        return bool(self.errors)


def validate_data(
    schema: StrMapping,
    data: StrMapping,
    *,
    errors: str = "raise",
    normalize: bool = True,
) -> ValidationException | None:
    """Handle :class:`ValidationError` yielded by :func:`validate_data_helper`.

    This function invokes :func:`validate_data_helper` to normalize and validate ``data`` using ``schema``. The list of
    :class:`ValidationError` are stored in the container :class:`ValidationException`. If :class:`ValidationException`
    is empty (i.e. no :class:`ValidationError` were yielded by :func:`validate_data_helper`), this function returns
    None. Otherwise, ``errors`` determines if the :class:`ValidationException` is ignored, returned, or raised.

    Args:
        schema:    string-keyed schema mappings
        data:      atomized string-keyed mappings
        errors:    one of {ignore, return, raise}
        normalize: normalize data using schema if True

    Returns:
        :class:`ValidationException` if not empty and ``errors`` is return

    Raises:
        :class:`ValidationException` if not empty and ``errors`` is raise

    Examples:
        >>> schema = {'a': StrSchemaType(required=True)}
        >>> data   = {'a': None}
        >>> validate_data(schema, data, errors='ignore')

        >>> validate_data(schema, data, errors='return')
        ValidationException()
        >>> validate_data(schema, data, errors='raise')
        Traceback (most recent call last):
        ...
        askcell.utils.types.ValidationException: Schema error

    """
    if errors not in ["ignore", "raise", "return"]:
        raise ValueError(f"Invalid validate_data errors: {errors!r}")

    exception = ValidationException(validate_data_helper(schema, data, normalize=normalize))

    if not exception:
        return None
    elif errors == "raise":
        raise exception
    elif errors == "return":
        return exception
    return None  # errors == 'ignore'


def validate_data_helper(
    schema: StrMapping,
    data: StrMapping,
    *,
    normalize: bool = True,
    path: list[str] | None = None,
) -> Iterator[ValidationError]:
    r"""Verify that data values are valid based on schema after (optionally) data value normalization.

    Utilizes attribute schema definitions to validate present values.
    Ensures that every schema key mapping exists in data.
    If normalize is True, the values in data are normalized using schema. An error will be raised if normalize is True
    and data is an immutable string mapping.

    Args:
        schema:    string-keyed schema mappings
        data:      atomized string-keyed mappings
        normalize: normalize data using schema if True
        path:      iterable of nested keys traversed to get schema and data

    Yields:
        ValidationError for each failed validation condition

    Examples:
        >>> schema = {'a': StrSchemaType(), 'b': IntSchemaType()}
        >>> data   = {'a': 'foo', 'b': 1}
        >>> list(validate_data_helper(schema, data))
        []

        >>> schema = {'a': IntSchemaType()}
        >>> data   = {'a': '{a.b.c}'}
        >>> list(validate_data_helper(schema, data))    # doctest: +NORMALIZE_WHITESPACE
        [ValidationError(msg="cannot validate {a.b.c} of type <class 'str'> using IntSchemaType",
        path='a', value='{a.b.c}', level=40)]

        >>> schema = {'a': IntSchemaType(required=True)}
        >>> data   = {'a': 5}
        >>> list(validate_data_helper(schema, data))
        []

        >>> schema = {'a': IntSchemaType(required=True)}
        >>> data   = {'a': None}
        >>> list(validate_data_helper(schema, data))  # doctest: +NORMALIZE_WHITESPACE
        [ValidationError(msg="cannot validate None of type <class 'NoneType'> using IntSchemaType",
        path='a', value=None, level=40)]

        >>> schema = {'a': IntSchemaType()}
        >>> data   = {'b': 1}
        >>> list(validate_data_helper(schema, data))  # doctest: +NORMALIZE_WHITESPACE
        [ValidationError(msg="key 'b' was not in schema", path='b', value=None, level=40),
        ValidationError(msg="key 'a' was not in data", path='a', value=None, level=40)]

    """
    if path is None:
        path = []

    fail = False
    if not isinstance(schema, (dict, MappingProxyType)):
        yield ValidationError("schema is not a mapping", path)
        fail = True
    elif not isinstance(data, (dict, MappingProxyType)):
        yield ValidationError("data is not a mapping", path)
        fail = True

    if fail:
        return

    for key in data.keys() - schema.keys():
        yield ValidationError(f"key {key!r} was not in schema", path + [key])

    for key, stype in schema.items():
        if key not in data:
            yield ValidationError(f"key {key!r} was not in data", path + [key])
            continue

        value = data[key]

        # Validate mappings recursively
        if not isinstance(stype, SchemaType):
            yield from validate_data_helper(stype, value, normalize=normalize, path=path + [key])
            continue

        # Optionally normalize and validate values
        if normalize:
            value = stype.normalize(value)

        if not stype.validate(value):
            yield ValidationError(f"cannot validate {value} of type {type(value)} using {stype}", path + [key], value)
        elif normalize:
            data[key] = value  # type: ignore
