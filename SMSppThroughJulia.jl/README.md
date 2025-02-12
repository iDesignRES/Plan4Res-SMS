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
[855.0, 830.0, 805.0, 780.0, 755.0, 730.0, 705.0, 680.0, 655.0, 630.0, 605.0, 580.0, 555.0, 530.0, 505.0, 480.0, 455.0, 430.0, 405.0, 330.99999999999994, 312.99999999999994, 294.99999999999994, 276.99999999999994, 258.99999999999994]
