# Contributing to PICOTA

Thank you for your interest in contributing to PICOTA! This document provides guidelines and instructions for contributing.

## Code of Conduct

Be respectful, inclusive, and professional in all interactions with other contributors and maintainers.

## Getting Started

### 1. Fork & Clone
```bash
git clone https://github.com/YOUR_USERNAME/picota.git
cd picota
```

### 2. Set Up Development Environment
```bash
conda create -n picota-dev python=3.9
conda activate picota-dev
pip install -r requirements.txt
pip install pytest pytest-cov
```

### 3. Create a Feature Branch
```bash
git checkout -b feature/your-feature-name
```

## Development Workflow

### Code Style
- Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) guidelines
- Use descriptive variable and function names
- Add docstrings to all public functions
- Add type hints for better code clarity

### Example:
```python
def calculate_jaccard_similarity(set_a: set, set_b: set) -> float:
    """Calculate Jaccard similarity between two sets.
    
    Args:
        set_a: First set
        set_b: Second set
        
    Returns:
        float: Jaccard similarity (0.0 to 1.0)
    """
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    if union == 0:
        return 0.0
    return intersection / union
```

### Testing
Before submitting PR, ensure all tests pass:

```bash
# Run all tests
pytest tests/ -v

# Run with coverage report
pytest tests/ --cov=picota --cov-report=html

# Run specific test file
pytest tests/test_cycle_finder.py -v
```

### Adding Tests
When adding new features, include corresponding tests:

```python
def test_my_new_feature():
    """Test description"""
    result = my_new_function(test_input)
    assert result == expected_output
```

## Commit Message Guidelines

Write clear, descriptive commit messages in English:

```
feat: Add k-mer similarity filtering with MinHash

- Implement MinHash-based locality-sensitive hashing
- Reduce O(n²) similarity comparisons to O(n log n)
- Add unit tests for hash collision detection
- Benchmark: ~10x speedup on 100+ cycles

Fixes #42
```

### Format:
```
<type>: <subject>

<body (optional)>

<footer (optional)>
```

**Types:** `feat`, `fix`, `docs`, `style`, `refactor`, `perf`, `test`, `chore`

## Submitting a Pull Request

1. **Update your branch** with latest main:
   ```bash
   git fetch origin
   git rebase origin/main
   ```

2. **Run full test suite**:
   ```bash
   pytest tests/ -v
   pytest tests/ --cov=picota
   ```

3. **Push your branch**:
   ```bash
   git push origin feature/your-feature-name
   ```

4. **Create PR** at GitHub with:
   - Clear title and description
   - Reference related issues (#42)
   - Screenshots/benchmarks if applicable
   - Checklist:
     - [ ] Tests pass locally
     - [ ] Code follows PEP 8
     - [ ] Added type hints
     - [ ] Updated documentation if needed

## Documentation

- Update README.md for new features
- Add docstrings following Google style guide
- Update CHANGELOG.md with significant changes
- Include examples for new command-line options

## Performance Optimization Guidelines

When optimizing code:
1. Measure before & after (use `timeit` or `cProfile`)
2. Document performance improvements in commit message
3. Ensure optimization doesn't break existing tests
4. Consider memory usage and trade-offs

Example:
```python
import timeit

# Before
time_before = timeit.timeit(lambda: old_function(data), number=1000)

# After optimization
time_after = timeit.timeit(lambda: new_function(data), number=1000)
speedup = time_before / time_after
print(f"Speedup: {speedup:.1f}x")
```

## Reporting Issues

When reporting bugs:
1. Include Python version and OS
2. Provide minimal reproducible example
3. Show error message and traceback
4. List external tools and versions (SPAdes, BLAST+, etc.)

## Questions?

- Open a GitHub Issue for bugs/features
- Discussions section for general questions
- Email maintainers for sensitive security issues

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

**Thank you for contributing to PICOTA!** 🎉
