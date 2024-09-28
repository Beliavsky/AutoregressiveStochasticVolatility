# AutoregressiveStochasticVolatility
Simulate from and fit a discrete-time autoregressive log stochastic volatility model. Run with `python xarsv.py`.

```
Sample output:

#obs: 1000000

            mu       phi     sigma
true -0.500000  0.900000  0.400000
 fit -0.497048  0.910509  0.377125

                 mean       std      skew  kurtosis        min        max
vol          0.865088  0.417543  1.540771  4.381076   0.068404   6.841719
log(vol)    -0.249821  0.458181 -0.004259  0.003998  -2.682318   1.923039
returns      0.000676  0.961069  0.016876  4.004919 -13.000362  17.599214
returns/vol  0.000224  1.001080  0.002882 -0.004051  -4.770352   4.496121

            mu       phi     sigma
true -0.500000  0.900000  0.400000
 fit -0.498815  0.916295  0.365594

                 mean       std      skew  kurtosis        min        max
vol          0.865183  0.420412  1.564852  4.576003   0.088911   5.967455
log(vol)    -0.250893  0.460690 -0.002258 -0.001021  -2.420119   1.786320
returns     -0.000413  0.962733  0.015700  3.903745 -10.188885  11.650475
returns/vol -0.001248  1.000968 -0.002871  0.001489  -4.787957   4.829210

            mu       phi     sigma
true -0.500000  0.900000  0.400000
 fit -0.503347  0.889827  0.422405

                 mean       std      skew  kurtosis        min        max
vol          0.865523  0.419989  1.544392  4.325622   0.075386   5.731296
log(vol)    -0.250369  0.460482 -0.003971 -0.003401  -2.585135   1.745942
returns      0.000182  0.960359  0.004223  3.915706 -12.880456  11.390660
returns/vol -0.000037  0.998975 -0.000440 -0.000098  -4.861752   5.047978
```
