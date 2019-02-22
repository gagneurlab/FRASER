## ----style, eval=TRUE, echo=FALSE, results="asis"--------------------------
BiocStyle::latex2()

## ----knitr, echo=FALSE, cache=FALSE, results="hide"------------------------
library("knitr")
opts_chunk$set(
  tidy=FALSE,
  dev="png",
  #fig.show="hide",
  fig.width=4, fig.height=4.5,
  dpi = 300,
  cache=TRUE
  #message=FALSE,
  #warning=FALSE
)

