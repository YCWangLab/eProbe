"""
Result type for functional error handling.

This module implements a Result pattern (similar to Rust's Result<T, E>) 
for explicit error handling without exceptions in pure functions.

Usage:
    >>> result = do_something()
    >>> if result.is_ok():
    ...     value = result.unwrap()
    >>> else:
    ...     error = result.unwrap_err()

    >>> # Or use pattern matching style
    >>> match result:
    ...     case Ok(value):
    ...         print(f"Success: {value}")
    ...     case Err(error):
    ...         print(f"Error: {error}")
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import TypeVar, Generic, Callable, Union, Any

T = TypeVar("T")  # Success type
E = TypeVar("E")  # Error type
U = TypeVar("U")  # Transformed type


@dataclass(frozen=True, slots=True)
class Ok(Generic[T]):
    """Represents a successful result containing a value."""
    
    value: T
    
    def is_ok(self) -> bool:
        return True
    
    def is_err(self) -> bool:
        return False
    
    def unwrap(self) -> T:
        """Returns the contained value."""
        return self.value
    
    def unwrap_or(self, default: T) -> T:
        """Returns the contained value (ignores default)."""
        return self.value
    
    def unwrap_err(self) -> Any:
        """Raises ValueError since this is Ok."""
        raise ValueError("Called unwrap_err on Ok value")
    
    def map(self, fn: Callable[[T], U]) -> Ok[U]:
        """Applies function to contained value."""
        return Ok(fn(self.value))
    
    def map_err(self, fn: Callable[[Any], Any]) -> Ok[T]:
        """Returns self (no error to map)."""
        return self
    
    def and_then(self, fn: Callable[[T], Result[U, Any]]) -> Result[U, Any]:
        """Chains another Result-returning function."""
        return fn(self.value)
    
    def __repr__(self) -> str:
        return f"Ok({self.value!r})"


@dataclass(frozen=True, slots=True)
class Err(Generic[E]):
    """Represents a failed result containing an error."""
    
    error: E
    
    def is_ok(self) -> bool:
        return False
    
    def is_err(self) -> bool:
        return True
    
    def unwrap(self) -> Any:
        """Raises ValueError with the error message."""
        raise ValueError(f"Called unwrap on Err: {self.error}")
    
    def unwrap_or(self, default: T) -> T:
        """Returns the default value."""
        return default
    
    def unwrap_err(self) -> E:
        """Returns the contained error."""
        return self.error
    
    def map(self, fn: Callable[[Any], Any]) -> Err[E]:
        """Returns self (no value to map)."""
        return self
    
    def map_err(self, fn: Callable[[E], U]) -> Err[U]:
        """Applies function to contained error."""
        return Err(fn(self.error))
    
    def and_then(self, fn: Callable[[Any], Any]) -> Err[E]:
        """Returns self (short-circuits on error)."""
        return self
    
    def __repr__(self) -> str:
        return f"Err({self.error!r})"


# Type alias for Result
Result = Union[Ok[T], Err[E]]


def try_except(fn: Callable[..., T], *args, **kwargs) -> Result[T, str]:
    """
    Wraps a function call in try/except and returns a Result.
    
    Args:
        fn: Function to call
        *args: Positional arguments for fn
        **kwargs: Keyword arguments for fn
        
    Returns:
        Ok(result) on success, Err(error_message) on exception
    """
    try:
        return Ok(fn(*args, **kwargs))
    except Exception as e:
        return Err(f"{type(e).__name__}: {str(e)}")


def collect_results(results: list[Result[T, E]]) -> Result[list[T], E]:
    """
    Collects a list of Results into a Result of list.
    
    Returns Err on first error encountered, otherwise Ok with all values.
    
    Args:
        results: List of Result objects
        
    Returns:
        Ok(list of values) if all Ok, else first Err encountered
    """
    values = []
    for result in results:
        if result.is_err():
            return result  # type: ignore
        values.append(result.unwrap())
    return Ok(values)
