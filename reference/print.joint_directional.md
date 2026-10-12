# Print method for joint directional statements

Prints the table returned by
[`joint`](https://mar-cald.github.io/multibayes/reference/joint.md) with
`interval = FALSE`. Columns that take the same value in every row
(`null.value`, `direction`) are reported once in a header. The full set
of columns remains available for programmatic access via the usual `$`
and `[`.

## Usage

``` r
# S3 method for class 'joint_directional'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A `joint_directional` object returned by
  [`joint`](https://mar-cald.github.io/multibayes/reference/joint.md).

- digits:

  Number of digits used when printing. Defaults to 4.

- ...:

  Further arguments passed to `print.data.frame`.

## Value

`x`, invisibly.
