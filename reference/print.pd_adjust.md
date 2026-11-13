# Print method for `pd_adjust` objects

Prints the number of tests \\m\\, \\\bar{\pi}\\, the threshold \\c\\ and
the value \\c^\*\\ that *pd* must reach, then one row per test with `pd`
and `pd.adj`. Columns with the same value in every row (`null.value`,
`direction`) are shown once, above the table.

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
