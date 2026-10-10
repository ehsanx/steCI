# Column names used in non-standard evaluation (model formulas passed to
# glm()/aalen(), and a few dplyr expressions) that R CMD check would otherwise
# flag as "no visible binding for global variable". Declared here to keep the
# check clean.
globalVariables(c(
  "A", "A.baseline", "A.lead1", "Alag1", "Alag1.baseline",
  "Anext.equal.to.baseline", "L", "L.baseline", "L.lead1", "Tstop",
  "event", "id", "ipw.stab", "ipw.stab.cum", "ipw.stab.lag", "ones",
  "rownum", "time", "time.new", "time.stop", "time.stop.new", "trial",
  "visit", "wt.denom", "wt.num"
))
