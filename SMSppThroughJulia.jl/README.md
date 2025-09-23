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

SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
, "test/data/NBlock_ramp.nc4")

Expected output:

[855.0, 829.9999999999993, 804.9999999999987, 779.9999999999983, 754.9999999999978, 729.9999999999974, 704.9999999999969, 679.9999999999965, 654.999999999996, 629.9999999999956, 604.9999999999951, 579.9999999999947, 554.9999999999942, 529.9999999999937, 504.9999999999933, 479.9999999999929, 454.9999999999925, 429.9999999999921, 404.9999999999917, 379.99999999999136, 354.999999999991, 329.9999999999907, 304.99999999999034, 279.99999999999]

SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
, "test/data/NBlock_ramp2.nc4")
 
(In fact the same data but only the nature of the unit has changed which should change nothing to the output)

SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
, "test/data/NBlock_mindown.nc4")

[854.9999999999994, 829.9999999999994, 804.9999999999994, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
 
SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
, "test/data/NBlock_mindown2.nc4")

[855.0, 829.9999999999994, 804.999999999999, 800.0, 0.0, 0.0, 0.0, 0.0, 800.0, 830.0, 860.0, 890.0, 920.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0]
 
 SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
, "test/data/Nblock_mindown.nc4")

[855.0, 829.9999999999994, 804.999999999999, 800.0, 800.0, 800.0, 830.0, 860.0, 890.0, 920.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0]
 
 SMSppThroughJulia.value_nuclear_on_price( [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
, "test/data/NBlock_mod.nc4")

[887.0, 886.0, 885.0, 884.0, 883.0, 882.0, 881.0, 856.0, 855.0, 854.0, 853.0, 828.0, 827.0, 826.0, 825.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

 SMSppThroughJulia.value_nuclear_on_price( [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
, "test/data/NBlock_mod2.nc4")

[905.0, 910.0, 909.0, 884.0, 883.0, 882.0, 881.0, 856.0, 855.0, 854.0, 853.0, 828.0, 827.0, 826.0, 825.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
 
 SMSppThroughJulia.value_nuclear_on_price( [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
, "test/data/Nblock_mod3.nc4")
 
[910.0000000000005, 911.0000000000005, 912.000000000001, 913.000000000001, 925.0, 925.0, 920.000000000033, 910.0, 900.0, 890.0, 880.0, 855.0, 845.0, 835.0, 825.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Alternative - automated test

Alternatively one can run

using SMSppThroughJulia

followed up by 

SMSppThroughJulia.test()

then one should see displayed

Testing : 8 / 8 - passed => indicating succes


