# Print method for `pd_adjust` objects

Prints a header with the number of tests \\m\\, \\\pi_0\\, the
unadjusted threshold \\c\\ and the adjusted threshold \\c^\*\\ (a
direction is claimed when \\pd \> c^\*\\), followed by the per-test
table and a short note on `pd.adj`. Columns that take the same value in
every row (`null.value`, `direction`) are reported once in the header.
The full set of columns remains available for programmatic access via
the usual `$` and `[`.

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
