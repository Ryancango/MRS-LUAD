# Shared configuration for the MRS-LUAD code release.
# Run scripts from the repository root or set environment variables explicitly.

mrs_project_root <- function() {
  Sys.getenv("MRS_PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
}

scfea_project_root <- function() {
  Sys.getenv("SCFEA_ROOT", unset = file.path(dirname(mrs_project_root()), "scFEA"))
}

mrs_final_genes <- c("RRAS", "GPM6B", "FMO3", "CHRDL1", "ADRB2", "LPL", "ANGPT1", "SPOCK2", "P2RY14", "ADH1B")

mrs_coefficients <- c(
  RRAS = 0.342398338,
  GPM6B = -0.249540066,
  FMO3 = 0.034609870,
  CHRDL1 = -0.205775214,
  ADRB2 = 0.033408756,
  LPL = -0.092544047,
  ANGPT1 = 0.007075607,
  SPOCK2 = 0.022607685,
  P2RY14 = -0.063043353,
  ADH1B = 0.091332717
)
