"""
Tests for eprobe.core.result module.
"""

import pytest
from eprobe.core.result import Ok, Err, try_except, collect_results


class TestResult:
    """Tests for Result type."""
    
    def test_ok_is_ok(self):
        """Test that Ok result returns True for is_ok()."""
        result = Ok(42)
        assert result.is_ok() is True
        assert result.is_err() is False
    
    def test_err_is_err(self):
        """Test that Err result returns True for is_err()."""
        result = Err("error message")
        assert result.is_err() is True
        assert result.is_ok() is False
    
    def test_ok_unwrap(self):
        """Test unwrapping Ok value."""
        result = Ok("success")
        assert result.unwrap() == "success"
    
    def test_err_unwrap_raises(self):
        """Test that unwrapping Err raises ValueError."""
        result = Err("error")
        with pytest.raises(ValueError, match="error"):
            result.unwrap()
    
    def test_err_unwrap_err(self):
        """Test unwrapping Err error."""
        result = Err("error message")
        assert result.unwrap_err() == "error message"
    
    def test_ok_unwrap_err_raises(self):
        """Test that unwrap_err on Ok raises ValueError."""
        result = Ok(42)
        with pytest.raises(ValueError):
            result.unwrap_err()
    
    def test_ok_unwrap_or(self):
        """Test unwrap_or returns value for Ok."""
        result = Ok(42)
        assert result.unwrap_or(0) == 42
    
    def test_err_unwrap_or(self):
        """Test unwrap_or returns default for Err."""
        result = Err("error")
        assert result.unwrap_or(0) == 0
    
    def test_ok_map(self):
        """Test map transforms Ok value."""
        result = Ok(5)
        mapped = result.map(lambda x: x * 2)
        assert mapped.is_ok()
        assert mapped.unwrap() == 10
    
    def test_err_map(self):
        """Test map does not transform Err."""
        result = Err("error")
        mapped = result.map(lambda x: x * 2)
        assert mapped.is_err()
        assert mapped.unwrap_err() == "error"
    
    def test_ok_and_then(self):
        """Test and_then chains Ok operations."""
        def double_if_positive(x):
            if x > 0:
                return Ok(x * 2)
            return Err("not positive")
        
        result = Ok(5)
        chained = result.and_then(double_if_positive)
        assert chained.is_ok()
        assert chained.unwrap() == 10
    
    def test_err_and_then(self):
        """Test and_then passes through Err."""
        def double_if_positive(x):
            return Ok(x * 2)
        
        result = Err("original error")
        chained = result.and_then(double_if_positive)
        assert chained.is_err()
        assert chained.unwrap_err() == "original error"


class TestTryExcept:
    """Tests for try_except function."""
    
    def test_success_returns_ok(self):
        """Test that successful function returns Ok."""
        def divide(a, b):
            return a / b
        
        result = try_except(divide, 10, 2)
        assert result.is_ok()
        assert result.unwrap() == 5.0
    
    def test_exception_returns_err(self):
        """Test that exception returns Err."""
        def divide(a, b):
            return a / b
        
        result = try_except(divide, 10, 0)
        assert result.is_err()
        assert "ZeroDivision" in result.unwrap_err()


class TestCollectResults:
    """Tests for collect_results function."""
    
    def test_all_ok(self):
        """Test collecting all Ok results."""
        results = [Ok(1), Ok(2), Ok(3)]
        collected = collect_results(results)
        assert collected.is_ok()
        assert collected.unwrap() == [1, 2, 3]
    
    def test_one_err(self):
        """Test that one Err fails the collection."""
        results = [Ok(1), Err("failed"), Ok(3)]
        collected = collect_results(results)
        assert collected.is_err()
        assert collected.unwrap_err() == "failed"
    
    def test_empty_list(self):
        """Test collecting empty list."""
        results = []
        collected = collect_results(results)
        assert collected.is_ok()
        assert collected.unwrap() == []
