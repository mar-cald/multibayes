# Print method for `pd_adjust` objects

Summarises the constant quantities (`pi0`, `m`) in a header and prints
the per-test table. The full set of columns remains available for
programmatic access via the usual `$` and `[`.

## Usage

``` r
# S3 method for class 'pd_adjust'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A `pd_adjust` object returned by
  [`pd.adjust`](https://mar-cald.github.io/multibayes/reference/pd.adjust.md).

- digits:

  Number of digits used when printing. Defaults to 4.

- ...:

  Further arguments passed to `print.data.frame`.

## Value

`x`, invisibly.
