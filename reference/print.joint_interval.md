# Print method for simultaneous credible intervals

Prints the table returned by
[`joint`](https://mar-cald.github.io/multibayes/reference/joint.md) with
`interval = TRUE`. The joint probability is shown once, above the table.

## Usage

``` r
# S3 method for class 'joint_interval'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A `joint_interval` object returned by
  [`joint`](https://mar-cald.github.io/multibayes/reference/joint.md).

- digits:

  Number of decimal places used when printing. Defaults to 4.

- ...:

  Further arguments passed to `print.data.frame`.

## Value

`x`, invisibly.
