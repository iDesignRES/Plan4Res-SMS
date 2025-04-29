# SMSppThroughJulia.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.com//SMSppThroughJulia.jl.svg?branch=master)](https://travis-ci.com//SMSppThroughJulia.jl)
[![codecov.io](http://codecov.io/github//SMSppThroughJulia.jl/coverage.svg?branch=master)](http://codecov.io/github//SMSppThroughJulia.jl?branch=master)

# How to use
Set up a directory cpp - where the dynamically compiled library emx_smspp_lib must go and some dependencies must go
Then follow up with 
using SMSppThroughJulia
SMSppThroughJulia.test()

Expected output:
[855.0, 829.9999999999993, 804.9999999999987, 779.9999999999983, 754.9999999999978, 729.9999999999974, 704.9999999999969, 679.9999999999965, 654.999999999996, 629.9999999999956, 604.9999999999951, 579.9999999999947, 554.9999999999942, 529.9999999999937, 504.9999999999933, 479.9999999999929, 454.9999999999925, 429.9999999999921, 404.9999999999917, 379.99999999999136, 354.999999999991, 329.9999999999907, 304.99999999999034, 279.99999999999]
