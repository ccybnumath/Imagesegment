#!/bin/bash
## do stuff
R --slave <<EOF
  ## R code
  Rcpp::sourceCpp('~/Rstudio_ccy/Imagesegment/ImageGibbs.cpp')
  set.seed(101)
  x<-rnorm(1)
  save.image("~/Rstudio_ccy/script.RData")
EOF
python << EOF
print 'Hello $1'
EOF
python mail.py
