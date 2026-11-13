# Print method for joint directional statements

Prints the table returned by
[`joint`](https://mar-cald.github.io/multibayes/reference/joint.md) with
`interval = FALSE`. Columns with the same value in every row
(`null.value`, `direction`) are shown once, above the table.

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
