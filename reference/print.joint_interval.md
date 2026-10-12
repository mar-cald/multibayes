# Print method for simultaneous credible intervals

Prints the table returned by
[`joint`](https://mar-cald.github.io/multibayes/reference/joint.md) with
`interval = TRUE`, reporting the joint coverage probability once in a
header. The full set of columns remains available for programmatic
access via the usual `$` and `[`.

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
